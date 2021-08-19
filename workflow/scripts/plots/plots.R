library(dplyr)
library(ggplot2)
library(pbmcapply)
library(kableExtra)
library(latex2exp)
library(ggpubr)

load('../results/Tag_SNP_Selection/setting2.rda')
remove(list = setdiff(ls(),'region_info'))

LoadDosage <- function(chr=NULL,pos=NULL,REF=NULL,ALT=NULL,vcf_path,software_dir,incl_samples = NULL,chr_prefix='',names=F){
  if(is.null(chr) & is.null(pos)){
    dosage_data <- system(paste0('bcftools plugin dosage -- -t "GT" > ',vcf_path,'.chr',chr,'.dosage'),
                          intern = F)
    if(names){
      rsid <-system(paste0('bcftools query -f "%ID\n" ',vcf_path),
                    intern = T)
    }
  }else if(is.null(pos)){
    dosage_data <- system(paste0('bcftools view -r ',paste(paste0(chr_prefix,chr),collapse = ','),
                                 ' ',vcf_path,' | ','bcftools +dosage -- -t "GT" > ',vcf_path,'.chr',chr,'.dosage'),
                          intern = F)
    if(names){
      rsid <-system(paste0('bcftools view -r ',paste(paste0(chr_prefix,chr),collapse = ','),
                           ' ',vcf_path,' | ','bcftools query -f "%ID\n" '),
                    intern = T)
    }
  }else{
    pos_vect <- paste(paste0(chr_prefix,chr),pos,sep = ':')
    dosage_data <- system(paste0('bcftools view -r ',paste(pos_vect,collapse = ','),
                                 ' ',vcf_path,' | ','bcftools +dosage -- -t "GT" '),
                          intern = T)
    if(names){
      rsid <-system(paste0('bcftools view -r ',paste(pos_vect,collapse = ','),
                           ' ',vcf_path,' | ','bcftools query -f "%ID\n" '),
                    intern = T)
    }
  }
  #Parse dosage matrix
  if(!is.null(pos) & !is.null(chr)){
    dosage_matrix <- data.table::fread(text = dosage_data)
  }else{
    dosage_matrix <- data.table::fread(file = paste0(vcf_path,'.chr',chr,'.dosage'))
  }
  #Parse column names
  colnames(dosage_matrix) <- sapply(colnames(dosage_matrix),function(x) gsub(x=x,pattern = "\\[.*\\]", replacement = "",))
  #Remove uneeded columns
  if(!is.null(incl_samples)){
    ind = match(incl_samples,colnames(dosage_matrix))
  }else{
    ind = sort(which(!colnames(dosage_matrix) %in% c('#CHROM','POS','REF','ALT')))
  }
  if(!is.null(REF) & !is.null(ALT)){
    if(dosage_matrix$REF != REF | dosage_matrix$ALT != ALT){
      return(rep(NA,length(ind)))
    }
  }
  if(names==T){
    row_names <- rsid
    row_names[rsid=='.'] <- paste0(dosage_matrix$`#CHROM`,':',dosage_matrix$POS,'[b37]',dosage_matrix$REF,',',dosage_matrix$ALT)[rsid=='.']
    col_names <- colnames(dosage_matrix)[ind]
    dosage_matrix <- matrix(as.integer(as.matrix(dosage_matrix)[,ind,drop=F]),nrow = length(rsid),ncol = length(ind))
    rownames(dosage_matrix) <- row_names
    colnames(dosage_matrix) <- col_names
  }else{
    dosage_matrix <- matrix(as.integer(as.matrix(dosage_matrix)[,ind,drop=F]),nrow = length(dosage_matrix$POS),ncol = length(ind))
  }
  dosage_matrix[dosage_matrix==-1] <- NA
  return(dosage_matrix)
}
BinDf <- function(df,type,cut_off_points,metric_list){
  start_cut_off <- max(c(min(which(cut_off_points > 0.05)) - 1),1)
  cut_off_points <- cut_off_points[start_cut_off:length(cut_off_points)]
  bins <- pbmcmapply(function(x,y) {
    curBin <- dplyr::filter(df,MAF > x & MAF <= y)
    if(type=='info'){
      curBin$INFO
    }else if(type=='mismatch'){
      curIndex <- match(curBin$ID,metric_list$snp)
      metric_list$mismatch[curIndex]
    }else if(type =='cor'){
      curIndex <- match(curBin$ID,metric_list$snp)
      metric_list$cor[curIndex]
    }
  },cut_off_points[1:(length(cut_off_points) - 1)],cut_off_points[2:(length(cut_off_points))],SIMPLIFY = F,mc.cores = 5)
  return(list(bins=bins,cut_off_points=cut_off_points))
}
BinInfo <- function(df,non_h3a_snps,cut_off_points = NULL,non_h3a_cor_vector = NULL){
  non_h3a_snps <- unique(non_h3a_snps)
  non_h3a_matching_index <- match(non_h3a_snps,df$ID)
  df_non_h3a <- df[non_h3a_matching_index,]
  if(is.null(cut_off_points)){
    cut_off_points <- quantile(df_non_h3a$MAF,probs = c(seq(0.1,0.4,0.1),seq(0.51,1,0.07)))
  }
  non_h3a_cor_bin <- BinDf(df_non_h3a,'info',cut_off_points,non_h3a_cor_vector)
  return(non_h3a_cor_bin)
}
BinCor <- function(df,non_h3a_snps,cut_off_points = NULL,metric_list){
  non_h3a_snps <- unique(non_h3a_snps)
  non_h3a_matching_index <- match(non_h3a_snps,df$ID)
  df_non_h3a <- df[non_h3a_matching_index,]
  if(is.null(cut_off_points)){
    cut_off_points <- quantile(df_non_h3a$MAF,probs = c(seq(0.1,0.4,0.1),seq(0.51,1,0.07)))
  }
  non_h3a_cor_bin <- BinDf(df_non_h3a,'cor',cut_off_points,metric_list)
  return(non_h3a_cor_bin)
}
GetImputationDiff <- function(baseline_path,add_on_path,random_path,regions_improved,illumina_info,setting,addtl_tags,addtl_tags_random){
  if(setting == 'setting1'){
    tb_regions_tagged_snps <- names(unlist(sapply(regions_improved$tag_info,function(x) sapply(x,function(y)y$tagged_snps))))
    
    baseline_tb_regions <- readRDS(baseline_path)
    matched_ind_baseline <- unlist(pbmclapply(regions_improved$region_info,function(x) which(sapply(baseline_tb_regions$snps_in_gene_regions_parsed,function(y) any(x$candidate_snps$ID %in% y$wgs_snps$ID))),mc.cores = 20)) #match(sapply(regions_improved$region_info,function(x) x$region_id),sapply(baseline_tb_regions$snps_in_gene_regions_parsed,function(x) x$region_id))
    if(length(matched_ind_baseline) != length(regions_improved$region_info)){
      stop('Length Mismatch')
    }
    baseline_tb_regions <- baseline_tb_regions$snps_in_gene_regions_parsed[matched_ind_baseline]
    baseline_tb_regions <- baseline_tb_regions[sapply(baseline_tb_regions,function(x) x$region_id) != ""]
    existing_tags <- do.call(rbind,lapply(baseline_tb_regions, function(x) x$h3a_tags)) %>% distinct(ID,.keep_all = T)
    
    addon_tb_regions <- readRDS(add_on_path)
    matched_ind_addon <- unlist(pbmclapply(regions_improved$region_info,function(x) which(sapply(addon_tb_regions$snps_in_gene_regions_parsed,function(y) any(x$candidate_snps$ID %in% y$wgs_snps$ID))),mc.cores = 20)) #match(sapply(regions_improved$region_info,function(x) x$region_id),sapply(addon_tb_regions$snps_in_gene_regions_parsed,function(x) x$region_id))
    if(length(matched_ind_addon) != length(regions_improved$region_info)){
      stop('Length Mismatch')
    }
    addon_tb_regions <- addon_tb_regions$snps_in_gene_regions_parsed[matched_ind_addon]
    addon_tb_regions <- addon_tb_regions[sapply(addon_tb_regions,function(x) x$region_id) != ""]
    
    random_tb_regions <- readRDS(random_path)
    matched_ind_random <- unlist(pbmclapply(regions_improved$region_info,function(x) which(sapply(random_tb_regions$snps_in_gene_regions_parsed,function(y) any(x$candidate_snps$ID %in% y$wgs_snps$ID))),mc.cores = 20)) #match(sapply(regions_improved$region_info,function(x) x$region_id),sapply(random_tb_regions$snps_in_gene_regions_parsed,function(x) x$region_id))
    if(length(matched_ind_random) != length(regions_improved$region_info)){
      stop('Length Mismatch')
    }
    random_tb_regions <- random_tb_regions$snps_in_gene_regions_parsed[matched_ind_random]
    random_tb_regions <- random_tb_regions[sapply(random_tb_regions,function(x) x$region_id) != ""]
    
  }else if(setting == 'setting2'){
    tb_regions_tagged_snps <- names(unlist(sapply(regions_improved$tag_info,function(x) x$tagged_snps)))
    
    baseline_tb_regions <- readRDS(baseline_path)
    baseline_tb_regions <- baseline_tb_regions$snps_in_gene_regions_parsed[unique(sapply(regions_improved$tag_info,function(x) x$region_id))]
    existing_tags <- do.call(rbind,lapply(baseline_tb_regions, function(x) rbind(x$addtl_snps,x$h3a_tags))) %>% distinct(ID,.keep_all = T)
      
    addon_tb_regions <- readRDS(add_on_path)
    addon_tb_regions <- addon_tb_regions$snps_in_gene_regions_parsed[unique(sapply(regions_improved$tag_info,function(x) x$region_id))]
    
    random_tb_regions <- readRDS(random_path)
    random_tb_regions <- random_tb_regions$snps_in_gene_regions_parsed[unique(sapply(regions_improved$tag_info,function(x) x$region_id))]
    
  }else{
    stop('Invalid Setting')
  }
  baseline_tb_regions_imp_snps_raw <- do.call(rbind,lapply(baseline_tb_regions,function(x) x$imp_snps)) 
  baseline_tb_regions_imp_snps <- baseline_tb_regions_imp_snps_raw %>% dplyr::filter(!is.na(MAF))  %>% dplyr::filter(ID %in% tb_regions_tagged_snps) %>% dplyr::distinct(CHROM,POS,REF,ALT,.keep_all = T) #%>% dplyr::filter(MAF > 0.05)
  
  addon_tb_regions_imp_snps_raw <- do.call(rbind,lapply(addon_tb_regions,function(x) x$imp_snps)) 
  addon_tb_regions_imp_snps <- dplyr::inner_join(baseline_tb_regions_imp_snps %>% dplyr::select('CHROM'='CHROM','POS'='POS','REF'='REF','ALT'='ALT'),addon_tb_regions_imp_snps_raw,by=c('CHROM'='CHROM','POS'='POS','REF'='REF','ALT'='ALT')) %>% dplyr::distinct(CHROM,POS,REF,ALT,.keep_all = T)
  
  random_tb_regions_imp_snps_raw <- do.call(rbind,lapply(random_tb_regions,function(x) x$imp_snps)) 
  random_tb_regions_imp_snps <- dplyr::inner_join(baseline_tb_regions_imp_snps %>% dplyr::select('CHROM'='CHROM','POS'='POS','REF'='REF','ALT'='ALT'),random_tb_regions_imp_snps_raw,by=c('CHROM'='CHROM','POS'='POS','REF'='REF','ALT'='ALT')) %>% dplyr::distinct(CHROM,POS,REF,ALT,.keep_all = T)
  
  
  baseline_bins_info <- BinInfo(baseline_tb_regions_imp_snps,baseline_tb_regions_imp_snps$ID) 
  addon_bins_info <- BinInfo(addon_tb_regions_imp_snps,addon_tb_regions_imp_snps$ID,cut_off_points = baseline_bins_info$cut_off_points) 
  random_bins_info <- BinInfo(random_tb_regions_imp_snps,random_tb_regions_imp_snps$ID,cut_off_points = baseline_bins_info$cut_off_points) 
  
  mean_cut_off_points_baseline_info <- (baseline_bins_info$cut_off_points[1:(length(baseline_bins_info$cut_off_points) - 1)] + baseline_bins_info$cut_off_points[2:(length(baseline_bins_info$cut_off_points))]) / 2
  mean_cut_off_points_addon_info <- (addon_bins_info$cut_off_points[1:(length(addon_bins_info$cut_off_points) - 1)] + addon_bins_info$cut_off_points[2:(length(addon_bins_info$cut_off_points))]) / 2
  mean_cut_off_points_random_info <- (random_bins_info$cut_off_points[1:(length(random_bins_info$cut_off_points) - 1)] + random_bins_info$cut_off_points[2:(length(random_bins_info$cut_off_points))]) / 2
  
  
  GetImprovInfo <- function(baseline_tb_regions_imp_snps_non_h3a,addon_tb_regions_imp_snps_non_h3a,random_tb_regions_imp_snps_non_h3a,addon_tb_regions,existing_tags,illumina_info,addtl_tags,addtl_tags_random,type){
    baseline_tb_regions_imp_snps_non_h3a <- baseline_tb_regions_imp_snps_non_h3a %>% dplyr::distinct(CHROM,POS,REF,ALT,.keep_all = T)
    if(type == 'ALL'){
      addon_tb_regions_imp_snps_non_h3a_filt <- addon_tb_regions_imp_snps_non_h3a %>% dplyr::distinct(CHROM,POS,REF,ALT,.keep_all = T)
      random_tb_regions_imp_snps_non_h3a_filt <- random_tb_regions_imp_snps_non_h3a %>% dplyr::distinct(CHROM,POS,REF,ALT,.keep_all = T)
    }else{
      addon_tb_regions_imp_snps_non_h3a_filt <- addon_tb_regions_imp_snps_non_h3a[addon_tb_regions_imp_snps_non_h3a$Source == type,] %>% dplyr::distinct(CHROM,POS,REF,ALT,.keep_all = T)
      random_tb_regions_imp_snps_non_h3a_filt <- random_tb_regions_imp_snps_non_h3a[random_tb_regions_imp_snps_non_h3a$Source == type,] %>% dplyr::distinct(CHROM,POS,REF,ALT,.keep_all = T)
    }
    addtl_snps<- do.call(rbind,lapply(addon_tb_regions, function(x) x$wgs_snps)) %>% dplyr::filter(ID %in% addtl_tags)
    addtl_snps_random <- do.call(rbind,lapply(addon_tb_regions, function(x) x$wgs_snps)) %>% dplyr::filter(ID %in% addtl_tags_random)
    
    
    snps_added_df <- dplyr::left_join(addtl_snps,addon_tb_regions_imp_snps_non_h3a,by=c('CHROM'='CHROM','POS'='POS','REF'='REF','ALT'='ALT')) %>% 
      dplyr::left_join(illumina_info %>% dplyr::select(Locus_Name,Final_Score,Assay_Type) %>% dplyr::mutate(N_Beads = ifelse(Assay_Type=='InfiniumII',1,2)),by=c('ID.x'='Locus_Name')) %>% distinct(ID.x,.keep_all = T)
    snps_added_df$Source[is.na(snps_added_df$Source)] <- 'WGS'
    
    snps_added_df_random <- dplyr::left_join(addtl_snps_random,addon_tb_regions_imp_snps_non_h3a,by=c('CHROM'='CHROM','POS'='POS','REF'='REF','ALT'='ALT')) %>% 
      dplyr::left_join(illumina_info %>% dplyr::select(Locus_Name,Final_Score,Assay_Type) %>% dplyr::mutate(N_Beads = ifelse(Assay_Type=='InfiniumII',1,2)),by=c('ID.x'='Locus_Name')) %>% distinct(ID.x,.keep_all = T)
    snps_added_df_random$Source[is.na(snps_added_df$Source)] <- 'WGS'
    
    if(type == 'ALL'){
      snps_added <- nrow(snps_added_df)
    }else if (type == 'WGS'){
      snps_added <- sum(snps_added_df$Source == 'WGS')
      snps_added_df <- dplyr::filter(snps_added_df,Source == 'WGS')
      snps_added_df_random <- dplyr::filter(snps_added_df_random,Source == 'WGS')
      
    }else if (type =='AFGR'){
      snps_added <- sum(snps_added_df$Source == 'AFGR')
      snps_added_df <- dplyr::filter(snps_added_df,Source == 'AFGR')
      snps_added_df_random <- dplyr::filter(snps_added_df_random,Source == 'AFGR')
    }

    snps_improved_df <- dplyr::left_join(addon_tb_regions_imp_snps_non_h3a_filt,baseline_tb_regions_imp_snps_non_h3a,by=c('CHROM'='CHROM','POS'='POS','REF'='REF','ALT'='ALT')) %>% dplyr::mutate(improv = ifelse(INFO.x > INFO.y,T,F)) 
    n_improved <- length(which(is.na(snps_improved_df$improv) | snps_improved_df$improv == T))
    
    snps_improved_0p8_df <- dplyr::left_join(addon_tb_regions_imp_snps_non_h3a_filt,baseline_tb_regions_imp_snps_non_h3a,by=c('CHROM'='CHROM','POS'='POS','REF'='REF','ALT'='ALT')) %>% dplyr::mutate(improv = ifelse(INFO.x > INFO.y & INFO.x > 0.8 & INFO.y < 0.8,T,F)) 
    n_improved_0p8 <- length(which((is.na(snps_improved_0p8_df$improv) & snps_improved_0p8_df$INFO.x > 0.8) | snps_improved_0p8_df$improv == T))
    
    
    snps_improved_random_df <- dplyr::left_join(random_tb_regions_imp_snps_non_h3a_filt,baseline_tb_regions_imp_snps_non_h3a,by=c('CHROM'='CHROM','POS'='POS','REF'='REF','ALT'='ALT')) %>% dplyr::mutate(improv = ifelse(INFO.x > INFO.y,T,F)) 
    n_improved_random <- length(which(is.na(snps_improved_random_df$improv) | snps_improved_random_df$improv == T))
    
    snps_improved_0p8_random_df <- dplyr::left_join(random_tb_regions_imp_snps_non_h3a_filt,baseline_tb_regions_imp_snps_non_h3a,by=c('CHROM'='CHROM','POS'='POS','REF'='REF','ALT'='ALT')) %>% dplyr::mutate(improv = ifelse(INFO.x > INFO.y & INFO.x > 0.8 & INFO.y < 0.8,T,F)) 
    n_improved_0p8_random <- length(which((is.na(snps_improved_0p8_random_df$improv) & snps_improved_0p8_random_df$INFO.x > 0.8) | snps_improved_0p8_random_df$improv == T))
    
    return(data.frame(n_snps_added=snps_added,n_improved=n_improved,n_improved_0p8 = n_improved_0p8,
                      n_improved_random=n_improved_random,n_improved_0p8_random=n_improved_0p8_random,mean_design_score = mean(snps_added_df$Final_Score),n_beads_added = sum(snps_added_df$N_Beads)))
  }
  
  all_snps_improv <- GetImprovInfo(baseline_tb_regions_imp_snps_raw,addon_tb_regions_imp_snps_raw,random_tb_regions_imp_snps_raw,addon_tb_regions,existing_tags,illumina_info,addtl_tags,addtl_tags_random,'ALL')
  
  wgs_snps_improv <- GetImprovInfo(baseline_tb_regions_imp_snps_raw,addon_tb_regions_imp_snps_raw,random_tb_regions_imp_snps_raw,addon_tb_regions,existing_tags,illumina_info,addtl_tags,addtl_tags_random,'WGS')
  
  afgr_snps_improv <- GetImprovInfo(baseline_tb_regions_imp_snps_raw,addon_tb_regions_imp_snps_raw,random_tb_regions_imp_snps_raw,addon_tb_regions,existing_tags,illumina_info,addtl_tags,addtl_tags_random,'AFGR')
  
  #Info Score Plot
  INFO_df <- data.frame(MAF = c(mean_cut_off_points_baseline_info,mean_cut_off_points_addon_info,mean_cut_off_points_random_info),
                        INFO = c(sapply(baseline_bins_info$bins,mean),sapply(addon_bins_info$bins,mean),sapply(random_bins_info$bins,mean)),
                        Chip = c(rep('H3Africa\n',length(baseline_bins_info$bins)),rep('H3Africa and \nAdd-On\n',length(addon_bins_info$bins)),rep('H3Africa and \nRandom\n',length(random_bins_info$bins))))
  INFO_plot <- ggplot(data = INFO_df) +
    aes(x=MAF,y=INFO,color = Chip) + geom_line() + geom_point() + guides(color=guide_legend(title="Array Content")) + ylim(0,1)
  
  #R2 Plot
  baseline_cor_afgr <- readRDS('../results/Imputation_Eval/afgr_non_h3a_cor.rds') 
  baseline_cor_wgs <- readRDS('../results/Imputation_Eval/wgs_imp_non_h3a_cor.rds')
  if(grepl(x=tolower(add_on_path),pattern = 'setting1')){
    addon_cor_afgr <- readRDS('../results/Imputation_Eval/afgr_setting1_non_h3a_cor.rds') 
    addon_cor_wgs <- readRDS('../results/Imputation_Eval/wgs_imp_setting1_non_h3a_cor.rds') 
    random_cor_afgr <- readRDS('../results/Imputation_Eval/afgr_setting1_random_non_h3a_cor.rds') 
    random_cor_wgs <- readRDS('../results/Imputation_Eval/wgs_imp_setting1_random_non_h3a_cor.rds') 
    
  }else if(grepl(x=tolower(add_on_path),pattern = 'setting2')){
    addon_cor_afgr <- readRDS('../results/Imputation_Eval/afgr_setting2_non_h3a_cor.rds') 
    addon_cor_wgs <- readRDS('../results/Imputation_Eval/wgs_imp_setting2_non_h3a_cor.rds') 
    random_cor_afgr <- readRDS('../results/Imputation_Eval/afgr_setting2_random_non_h3a_cor.rds') 
    random_cor_wgs <- readRDS('../results/Imputation_Eval/wgs_imp_setting2_random_non_h3a_cor.rds') 
  }else{
    stop('Cannot find cor file')
  }
  baseline_tb_regions_imp_snps_afgr <- dplyr::filter(baseline_tb_regions_imp_snps,Source =='AFGR')
  baseline_tb_regions_imp_snps_wgs <- dplyr::filter(baseline_tb_regions_imp_snps,Source =='WGS')
  addon_tb_regions_imp_snps_afgr <- dplyr::filter(addon_tb_regions_imp_snps,Source =='AFGR')
  addon_tb_regions_imp_snps_wgs <- dplyr::filter(addon_tb_regions_imp_snps,Source =='WGS')
  random_tb_regions_imp_snps_afgr <- dplyr::filter(random_tb_regions_imp_snps,Source =='AFGR')
  random_tb_regions_imp_snps_wgs <- dplyr::filter(random_tb_regions_imp_snps,Source =='WGS')
  
  
  MatchSNPs <- function(imp_df,cor_list){
    index_match <- match(imp_df$ID,cor_list$snps)
    index_match <- index_match[!is.na(index_match)]
    matched <- list(snps = cor_list$snps[index_match],cor = cor_list$cor[index_match])
    return(matched)
  }
  wgs_snps <- data.table::rbindlist(lapply(addon_tb_regions,function(x)x$wgs_snps)) %>% dplyr::filter(MAF > 0.05)
  
  baseline_tb_regions_imp_snps_to_incl <- dplyr::inner_join(wgs_snps,baseline_tb_regions_imp_snps,
                                                            by=c('CHROM'='CHROM','POS'='POS','REF'='REF','ALT'='ALT')) %>% dplyr::select(CHROM=CHROM,POS=POS,ID=ID.y,REF=REF,ALT=ALT,INFO=INFO,AF=AF.x,Source=Source,MAF=MAF.x)
  addon_tb_regions_imp_snps_to_incl <- dplyr::inner_join(wgs_snps,addon_tb_regions_imp_snps,
                                                            by=c('CHROM'='CHROM','POS'='POS','REF'='REF','ALT'='ALT')) %>% dplyr::select(CHROM=CHROM,POS=POS,ID=ID.y,REF=REF,ALT=ALT,INFO=INFO,AF=AF.x,Source=Source,MAF=MAF.x)
  random_tb_regions_imp_snps_to_incl <- dplyr::inner_join(wgs_snps,random_tb_regions_imp_snps,
                                                            by=c('CHROM'='CHROM','POS'='POS','REF'='REF','ALT'='ALT')) %>% dplyr::select(CHROM=CHROM,POS=POS,ID=ID.y,REF=REF,ALT=ALT,INFO=INFO,AF=AF.x,Source=Source,MAF=MAF.x)
  
  baseline_cor_afgr_matched <- MatchSNPs(dplyr::filter(baseline_tb_regions_imp_snps_to_incl,Source=='AFGR'),baseline_cor_afgr)
  baseline_cor_wgs_matched <- MatchSNPs(dplyr::filter(baseline_tb_regions_imp_snps_to_incl,Source=='WGS'),baseline_cor_wgs)
  
  addon_cor_afgr_matched <- MatchSNPs(dplyr::filter(addon_tb_regions_imp_snps_to_incl,Source=='AFGR'),addon_cor_afgr)
  addon_cor_wgs_matched <- MatchSNPs(dplyr::filter(addon_tb_regions_imp_snps_to_incl,Source=='WGS'),addon_cor_wgs)
  
  random_cor_afgr_matched <- MatchSNPs(dplyr::filter(random_tb_regions_imp_snps_to_incl,Source=='AFGR'),random_cor_afgr)
  random_cor_wgs_matched <- MatchSNPs(dplyr::filter(random_tb_regions_imp_snps_to_incl,Source=='WGS'),random_cor_wgs)
  
  
  
  baseline_bins_cor <- BinCor(baseline_tb_regions_imp_snps_to_incl,baseline_tb_regions_imp_snps_to_incl$ID,cut_off_points = NULL,metric_list = list(snps = c(baseline_cor_afgr_matched$snps,baseline_cor_wgs_matched$snps),cor = c(baseline_cor_afgr_matched$cor^2,baseline_cor_wgs_matched$cor^2))) #BinCor(baseline_tb_regions_imp_snps_to_incl,setdiff(baseline_tb_regions_imp_snps_to_incl$ID,h3a_snps_addon$V1),cut_off_points = NULL,metric_list = list(snps = c(baseline_cor_afgr_matched$snps,baseline_cor_wgs_matched$snps),cor = c(baseline_cor_afgr_matched$cor^2,baseline_cor_wgs_matched$cor^2)))
  
  addon_bins_cor <- BinCor(addon_tb_regions_imp_snps_to_incl,addon_tb_regions_imp_snps_to_incl$ID,cut_off_points = baseline_bins_cor$cut_off_points,metric_list = list(snps = c(addon_cor_afgr_matched$snps,addon_cor_wgs_matched$snps),cor = c(addon_cor_afgr_matched$cor^2,addon_cor_wgs_matched$cor^2)))
  
  random_bins_cor <- BinCor(random_tb_regions_imp_snps_to_incl,random_tb_regions_imp_snps_to_incl$ID,cut_off_points = baseline_bins_cor$cut_off_points,metric_list = list(snps = c(random_cor_afgr_matched$snps,random_cor_wgs_matched$snps),cor = c(random_cor_afgr_matched$cor^2,random_cor_wgs_matched$cor^2)))
  
  mean_cut_off_points_baseline_cor <- (baseline_bins_cor$cut_off_points[1:(length(baseline_bins_cor$cut_off_points) - 1)] + baseline_bins_cor$cut_off_points[2:(length(baseline_bins_cor$cut_off_points))]) / 2
  mean_cut_off_points_addon_cor <- (addon_bins_cor$cut_off_points[1:(length(addon_bins_cor$cut_off_points) - 1)] + addon_bins_cor$cut_off_points[2:(length(addon_bins_cor$cut_off_points))]) / 2
  mean_cut_off_points_random_cor <- (random_bins_cor$cut_off_points[1:(length(random_bins_cor$cut_off_points) - 1)] + random_bins_cor$cut_off_points[2:(length(random_bins_cor$cut_off_points))]) / 2
  
  R2_df <- data.frame(MAF = c(mean_cut_off_points_baseline_cor,mean_cut_off_points_addon_cor,mean_cut_off_points_random_cor),
                      r2 = c(sapply(baseline_bins_cor$bins,function(x) mean(x,na.rm = T)),
                             sapply(addon_bins_cor$bins,function(x) mean(x,na.rm = T)),
                             sapply(random_bins_cor$bins,function(x) mean(x,na.rm = T))),
                      Chip = c(rep('H3Africa\n',length(baseline_bins_cor$bins)),rep('H3Africa and\nAdd-On\n',length(addon_bins_cor$bins)),rep('H3Africa and\nRandom\n',length(random_bins_cor$bins))))
  R2_plot <- ggplot(data = R2_df) +
    aes(x=MAF,y=r2,color = Chip) + geom_line() + geom_point() + guides(color=guide_legend(title="Array Content")) + ylim(0,1)
  

  #Improvement Table
  type = c('All SNPs','WGS SNPs','AFGR SNPs')
  df_tb <- rbind(all_snps_improv,wgs_snps_improv,afgr_snps_improv)
  df_tb$type <- type
  df_tb <- df_tb %>% dplyr::select(type,n_snps_added=n_snps_added,mean_design_score=mean_design_score,
                                   n_beads_added=n_beads_added,n_improved=n_improved,n_improved_0p8 = n_improved_0p8,
                                   n_improved_random=n_improved_random,n_improved_0p8_random=n_improved_0p8_random)
  tbl <- kable(df_tb, format = "latex", longtable = TRUE,booktabs = T,caption = "Tag SNPs added and Improvement Statistics.") %>%
    kable_styling(latex_options = "HOLD_position")
  
  return(list(INFO_plot=INFO_plot,R2_plot=R2_plot,tbl=df_tb,R2_df=R2_df,INFO_df=INFO_df,
              addon_tb_regions=addon_tb_regions,random_tb_regions=random_tb_regions))
  
}

GetImputationDiffTagger <- function(baseline_path,tagger_path,regions_improved,illumina_info,setting,addtl_tags){
  if(setting == 'setting1'){
    baseline_tb_regions <- readRDS(baseline_path)
    matched_ind_baseline <- unlist(pbmclapply(regions_improved$region_info,function(x) which(sapply(baseline_tb_regions$snps_in_gene_regions_parsed,function(y) any(x$candidate_snps$ID %in% y$wgs_snps$ID))),mc.cores = 20)) #match(sapply(regions_improved$region_info,function(x) x$region_id),sapply(baseline_tb_regions$snps_in_gene_regions_parsed,function(x) x$region_id))
    if(length(matched_ind_baseline) != length(regions_improved$region_info)){
      stop('Length Mismatch')
    }
    baseline_tb_regions <- baseline_tb_regions$snps_in_gene_regions_parsed[matched_ind_baseline]
    baseline_tb_regions <- baseline_tb_regions[sapply(baseline_tb_regions,function(x) x$region_id) != ""]
    existing_tags <- do.call(rbind,lapply(baseline_tb_regions, function(x) x$h3a_tags)) %>% distinct(ID,.keep_all = T)
    
    tagger_tb_regions <- readRDS(tagger_path)
    matched_ind_tagger <- unlist(pbmclapply(regions_improved$region_info,function(x) which(sapply(tagger_tb_regions$snps_in_gene_regions_parsed,function(y) any(x$candidate_snps$ID %in% y$wgs_snps$ID))),mc.cores = 20)) #match(sapply(regions_improved$region_info,function(x) x$region_id),sapply(tagger_tb_regions$snps_in_gene_regions_parsed,function(x) x$region_id))
    if(length(matched_ind_tagger) != length(regions_improved$region_info)){
      stop('Length Mismatch')
    }
    tagger_tb_regions <- tagger_tb_regions$snps_in_gene_regions_parsed[matched_ind_tagger]
    tagger_tb_regions <- tagger_tb_regions[sapply(tagger_tb_regions,function(x) x$region_id) != ""]
    
  }else if(setting == 'setting2'){
    baseline_tb_regions <- readRDS(baseline_path)
    baseline_tb_regions <- baseline_tb_regions$snps_in_gene_regions_parsed[unique(sapply(regions_improved$tag_info,function(x) x$region_id))]
    existing_tags <- do.call(rbind,lapply(baseline_tb_regions, function(x) rbind(x$addtl_snps,x$h3a_tags))) %>% distinct(ID,.keep_all = T)
    
    tagger_tb_regions <- readRDS(tagger_path)
    tagger_tb_regions <- tagger_tb_regions$snps_in_gene_regions_parsed[unique(sapply(regions_improved$tag_info,function(x) x$region_id))]
    
  }else{
    stop('Invalid Setting')
  }
  baseline_tb_regions_imp_snps_raw <- do.call(rbind,lapply(baseline_tb_regions,function(x) x$imp_snps)) 
  tagger_tb_regions_imp_snps_raw <- do.call(rbind,lapply(tagger_tb_regions,function(x) x$imp_snps)) 

  
  GetImprovInfo <- function(baseline_tb_regions_imp_snps_non_h3a,tagger_tb_regions_imp_snps_non_h3a,existing_tags,illumina_info,tagger_tb_regions,addtl_tags,type){
    baseline_tb_regions_imp_snps_non_h3a <- baseline_tb_regions_imp_snps_non_h3a %>% dplyr::distinct(CHROM,POS,REF,ALT,.keep_all = T)
    if(type == 'ALL'){ 
      tagger_tb_regions_imp_snps_non_h3a_filt <- tagger_tb_regions_imp_snps_non_h3a %>% dplyr::distinct(CHROM,POS,REF,ALT,.keep_all = T)
    }else{
      tagger_tb_regions_imp_snps_non_h3a_filt <- tagger_tb_regions_imp_snps_non_h3a[tagger_tb_regions_imp_snps_non_h3a$Source == type,] %>% dplyr::distinct(CHROM,POS,REF,ALT,.keep_all = T)
      
    }
    addtl_snps<- do.call(rbind,lapply(tagger_tb_regions, function(x) x$wgs_snps)) %>% dplyr::filter(ID %in% addtl_tags)
    
    snps_added_df <- dplyr::left_join(addtl_snps,tagger_tb_regions_imp_snps_non_h3a,by=c('CHROM'='CHROM','POS'='POS','REF'='REF','ALT'='ALT')) %>% 
      dplyr::left_join(illumina_info %>% dplyr::select(Locus_Name,Final_Score,Assay_Type) %>% dplyr::mutate(N_Beads = ifelse(Assay_Type=='InfiniumII',1,2)),by=c('ID.x'='Locus_Name')) %>% distinct(ID.x,.keep_all = T)
    snps_added_df$Source[is.na(snps_added_df$Source)] <- 'WGS'
    
    if(type == 'ALL'){
      snps_added <- nrow(snps_added_df)
    }else if (type == 'WGS'){
      snps_added <- sum(snps_added_df$Source == 'WGS')
      snps_added_df <- dplyr::filter(snps_added_df,Source == 'WGS')
    }else if (type =='AFGR'){
      snps_added <- sum(snps_added_df$Source == 'AFGR')
      snps_added_df <- dplyr::filter(snps_added_df,Source == 'AFGR')
    }
    
    snps_improved_tagger_df <- dplyr::left_join(tagger_tb_regions_imp_snps_non_h3a_filt %>% dplyr::distinct(CHROM,POS,REF,ALT,.keep_all = T),baseline_tb_regions_imp_snps_non_h3a,by=c('CHROM'='CHROM','POS'='POS','REF'='REF','ALT'='ALT')) %>% dplyr::mutate(improv = ifelse(INFO.x > INFO.y,T,F)) 
    n_improved_tagger <- length(which(is.na(snps_improved_tagger_df$improv) | snps_improved_tagger_df$improv == T))
    
    snps_improved_0p8_tagger_df <- dplyr::left_join(tagger_tb_regions_imp_snps_non_h3a_filt,baseline_tb_regions_imp_snps_non_h3a,by=c('CHROM'='CHROM','POS'='POS','REF'='REF','ALT'='ALT')) %>% dplyr::mutate(improv = ifelse(INFO.x > INFO.y & INFO.x > 0.8 & INFO.y < 0.8,T,F)) 
    n_improved_0p8_tagger <- length(which((is.na(snps_improved_0p8_tagger_df$improv) & snps_improved_0p8_tagger_df$INFO.x > 0.8) | snps_improved_0p8_tagger_df$improv == T))
    
    return(data.frame(n_snps_added=snps_added,
                      n_improved_tagger,n_improved_0p8_tagger,mean_design_score = mean(snps_added_df$Final_Score),n_beads_added = sum(snps_added_df$N_Beads)))
  }
  
  all_snps_improv <- GetImprovInfo(baseline_tb_regions_imp_snps_raw,tagger_tb_regions_imp_snps_raw,existing_tags,illumina_info,tagger_tb_regions,addtl_tags,'ALL')
  
  wgs_snps_improv <- GetImprovInfo(baseline_tb_regions_imp_snps_raw,tagger_tb_regions_imp_snps_raw,existing_tags,illumina_info,tagger_tb_regions,addtl_tags,'WGS')
  
  afgr_snps_improv <- GetImprovInfo(baseline_tb_regions_imp_snps_raw,tagger_tb_regions_imp_snps_raw,existing_tags,illumina_info,tagger_tb_regions,addtl_tags,'AFGR')
  
  #Improvement Table
  type = c('All SNPs','WGS SNPs','AFGR SNPs')
  df_tb <- rbind(all_snps_improv,wgs_snps_improv,afgr_snps_improv)
  df_tb$type <- type
  df_tb <- df_tb %>% dplyr::select(type,n_snps_added=n_snps_added,mean_design_score=mean_design_score,n_beads_added=n_beads_added,n_improved_tagger,n_improved_0p8_tagger)
  tbl <- kable(df_tb, format = "latex", longtable = TRUE,booktabs = T,caption = "Tag SNPs added and Improvement Statistics.") %>%
    kable_styling(latex_options = "HOLD_position")
  
  return(list(tbl=df_tb,tagger_tb_regions=tagger_tb_regions))
  
}

PlotCovExample <- function(snps_in_gene_regions_parsed_baseline,snps_in_gene_regions_parsed_setting1,snps_in_gene_regions_parsed_setting1_random){
  #Plot SNPs in example region where imputation improved
  PlotRegion <- function(snps_in_gene_regions_parsed,index,MAF_Plot=T,title=F,r2=F,df_only=T){
    cur_imp_snps <- snps_in_gene_regions_parsed[[index]]$imp_snps
    cur_imp_snps <- dplyr::filter(cur_imp_snps,MAF > 0.05)
    df <- data.frame(Position=cur_imp_snps$POS,r2=cur_imp_snps$r2,INFO=cur_imp_snps$INFO,SNP_Type=ifelse(cur_imp_snps$POS %in% snps_in_gene_regions_parsed[[index]]$h3a_tags$POS,'H3Africa',
                                                                                                         ifelse(cur_imp_snps$POS %in% snps_in_gene_regions_parsed[[index]]$addtl_snps$POS,'Add-On','Imputed'))) %>% dplyr::arrange(desc(SNP_Type))
    if(df_only){
      return(df)
    }
    cbp1 <- c('Imputed' = "#999999", 'H3Africa' = "#E69F00", 'Add-On' = "#D55E00")
    if(title){
      if(r2){
        p1 <- ggplot(data=df) + aes(x=Position,y=r2,color = SNP_Type) + geom_point(alpha = 0.8) + 
          xlim(c(snps_in_gene_regions_parsed[[index]]$start_region,snps_in_gene_regions_parsed[[index]]$end_region)) + 
          ggtitle(paste0('Imputed SNPs: Region',index)) + ylim(0,1) + scale_colour_manual(values = cbp1)
        
      }else{
        p1 <- ggplot(data=df) + aes(x=Position,y=INFO,color = SNP_Type) + geom_point(alpha = 0.8) + 
          xlim(c(snps_in_gene_regions_parsed[[index]]$start_region,snps_in_gene_regions_parsed[[index]]$end_region)) +
          ggtitle(paste0('Imputed SNPs: Region',index)) + ylim(0,1) + scale_colour_manual(values = cbp1)
        
      }
    }else{
      if(r2){
        p1 <- ggplot(data=df) + aes(x=Position/1e6,y=r2,color = SNP_Type) + geom_point(alpha = 0.8) + 
          xlim(c(snps_in_gene_regions_parsed[[index]]$start_region/1e6,snps_in_gene_regions_parsed[[index]]$end_region/1e6)) + 
          ggtitle('') + ylim(0,1) + scale_colour_manual(values = cbp1) + guides(color=guide_legend(title="SNP Type")) + xlab(paste0('Chr ',unique(cur_imp_snps$CHROM),' Position (Mb)'))
        
      }else{
        p1 <- ggplot(data=df) + aes(x=Position/1e6,y=INFO,color = SNP_Type) + geom_point(alpha = 0.8) + 
          xlim(c(snps_in_gene_regions_parsed[[index]]$start_region/1e6,snps_in_gene_regions_parsed[[index]]$end_region/1e6)) + 
          ggtitle('') + ylim(0,1) + scale_colour_manual(values = cbp1) + guides(color=guide_legend(title="SNP Type")) + xlab(paste0('Chr ',unique(cur_imp_snps$CHROM),' Position (Mb)'))
        
      }
      
    }
    if(MAF_Plot){
      wgs_snps <- snps_in_gene_regions_parsed[[index]]$wgs_snps
      wgs_snps <- dplyr::filter(wgs_snps,MAF > 0.05)
      all_snps <- rbind(data.frame(POS = c(wgs_snps$POS,cur_imp_snps$POS),MAF = c(wgs_snps$MAF,cur_imp_snps$MAF),stringsAsFactors = F))
      p2_df <- data.frame(Position=wgs_snps$POS,MAF=wgs_snps$MAF,SNP_Type=ifelse(wgs_snps$POS %in% snps_in_gene_regions_parsed[[index]]$h3a_tags$POS,'H3Africa',
                                                                                 ifelse(wgs_snps$POS %in% snps_in_gene_regions_parsed[[index]]$addtl_snps$POS,'Add-On',
                                                                                        ifelse(wgs_snps$POS %in% snps_in_gene_regions_parsed[[index]]$afgr_sites$POS,'AFGR SNPs','WGS SNPs'))))
      cbp1 <- c('AFGR SNPs' = "#999999", 'H3Africa' = "#E69F00", 'WGS SNPs' = "#56B4E9", 'Add-On' = "#D55E00")
      p2_df$SNP_Type <- factor(p2_df$SNP_Type,levels = c('AFGR SNPs','H3A','WGS SNPs','Add-On'))
      p2 <- ggplot(p2_df) + 
        aes(x=Position,y=MAF,color=SNP_Type) + geom_point(alpha = 0.8) + 
        xlim(c(snps_in_gene_regions_parsed[[index]]$start_region,snps_in_gene_regions_parsed[[index]]$end_region)) +  ggtitle('Candidate Tags') + ylim(0,0.5) + scale_colour_manual(values = cbp1)
      return(ggpubr::ggarrange(p1,p2,nrow = 2))
    }else{
      return(p1)
    }
  }
  
  
  #Get r2
  baseline_cor_afgr <- readRDS('../results/Imputation_Eval/afgr_non_h3a_cor.rds') 
  baseline_cor_wgs <- readRDS('../results/Imputation_Eval/wgs_imp_non_h3a_cor.rds') 
  setting1_cor_afgr <- readRDS('../results/Imputation_Eval/afgr_setting1_non_h3a_cor.rds') 
  setting1_cor_wgs <- readRDS('../results/Imputation_Eval/wgs_imp_setting1_non_h3a_cor.rds') 
  setting1_random_cor_afgr <- readRDS('../results/Imputation_Eval/afgr_setting1_random_non_h3a_cor.rds') 
  setting1_random_cor_wgs <- readRDS('../results/Imputation_Eval/wgs_imp_setting1_random_non_h3a_cor.rds') 
  
  
  baseline_imp_snps <- snps_in_gene_regions_parsed_baseline[[1024]]$imp_snps
  baseline_afgr_ids <- dplyr::filter(baseline_imp_snps,Source == 'AFGR' & MAF > 0.05)$ID
  baseline_wgs_ids <- dplyr::filter(baseline_imp_snps,Source == 'WGS' & MAF > 0.05)$ID
  baseline_r2_afgr <- baseline_cor_afgr$cor[match(baseline_afgr_ids,baseline_cor_afgr$snps)]
  baseline_r2_afgr[is.na(baseline_r2_afgr)] <- 1 #Tag SNPs, r2 were not calculated, but should be 1. 
  baseline_r2_wgs <- baseline_cor_wgs$cor[match(baseline_wgs_ids,baseline_cor_wgs$snps)]
  baseline_r2_wgs[is.na(baseline_r2_wgs)] <- 1 #Tag SNPs, r2 were not calculated, but should be 1. 
  baseline_imp_snps_r2 <- data.frame(ID = c(baseline_afgr_ids,baseline_wgs_ids),r2 = c(baseline_r2_afgr,baseline_r2_wgs))
  
  
  setting1_imp_snps <- snps_in_gene_regions_parsed_setting1[[1024]]$imp_snps
  setting1_afgr_ids <- dplyr::filter(setting1_imp_snps,Source == 'AFGR' & MAF > 0.05)$ID
  setting1_wgs_ids <- dplyr::filter(setting1_imp_snps,Source == 'WGS' & MAF > 0.05)$ID
  setting1_r2_afgr <- setting1_cor_afgr$cor[match(setting1_afgr_ids,setting1_cor_afgr$snps)]
  setting1_r2_afgr[is.na(setting1_r2_afgr)] <- 1 #Tag SNPs, r2 were not calculated, but should be 1. 
  setting1_r2_wgs <- setting1_cor_wgs$cor[match(setting1_wgs_ids,setting1_cor_wgs$snps)]
  setting1_r2_wgs[is.na(setting1_r2_wgs)] <- 1 #Tag SNPs, r2 were not calculated, but should be 1. 
  setting1_imp_snps_r2 <- data.frame(ID = c(setting1_afgr_ids,setting1_wgs_ids),r2 = c(setting1_r2_afgr,setting1_r2_wgs))
  
  setting1_random_imp_snps <- snps_in_gene_regions_parsed_setting1_random[[1024]]$imp_snps
  setting1_random_afgr_ids <- dplyr::filter(setting1_random_imp_snps,Source == 'AFGR' & MAF > 0.05)$ID
  setting1_random_wgs_ids <- dplyr::filter(setting1_random_imp_snps,Source == 'WGS' & MAF > 0.05)$ID
  setting1_random_r2_afgr <- setting1_random_cor_afgr$cor[match(setting1_random_afgr_ids,setting1_random_cor_afgr$snps)]
  setting1_random_r2_afgr[is.na(setting1_random_r2_afgr)] <- 1 #Tag SNPs, r2 were not calculated, but should be 1. 
  setting1_random_r2_wgs <- setting1_random_cor_wgs$cor[match(setting1_random_wgs_ids,setting1_random_cor_wgs$snps)]
  setting1_random_r2_wgs[is.na(setting1_random_r2_wgs)] <- 1 #Tag SNPs, r2 were not calculated, but should be 1. 
  setting1_random_imp_snps_r2 <- data.frame(ID = c(setting1_random_afgr_ids,setting1_random_wgs_ids),r2 = c(setting1_random_r2_afgr,setting1_random_r2_wgs))
  
  
  snps_in_gene_regions_parsed_baseline[[1024]]$imp_snps <- dplyr::left_join(snps_in_gene_regions_parsed_baseline[[1024]]$imp_snps,baseline_imp_snps_r2)
  baseline_plot_info <- PlotRegion(snps_in_gene_regions_parsed_baseline,1024,MAF_Plot=F)
  baseline_plot_r2 <- PlotRegion(snps_in_gene_regions_parsed_baseline,1024,MAF_Plot=F,r2=T)
  
  snps_in_gene_regions_parsed_setting1[[1024]]$imp_snps <- dplyr::left_join(snps_in_gene_regions_parsed_setting1[[1024]]$imp_snps,setting1_imp_snps_r2)
  setting1_plot_info <- PlotRegion(snps_in_gene_regions_parsed_setting1,1024,MAF_Plot=F)
  setting1_plot_r2 <- PlotRegion(snps_in_gene_regions_parsed_setting1,1024,MAF_Plot=F,r2=T)
  
  snps_in_gene_regions_parsed_setting1_random[[1024]]$imp_snps <- dplyr::left_join(snps_in_gene_regions_parsed_setting1_random[[1024]]$imp_snps,setting1_random_imp_snps_r2)
  setting1_random_plot_info <- PlotRegion(snps_in_gene_regions_parsed_setting1_random,1024,MAF_Plot=F)
  setting1_random_plot_r2 <- PlotRegion(snps_in_gene_regions_parsed_setting1_random,1024,MAF_Plot=F,r2=T)
  
  
  
  baseline_plot_info$Pass <- 'H3Africa'
  setting1_plot_info$Pass <- 'H3Africa + Add-On'
  setting1_random_plot_info$Pass <- 'H3Africa + Random'
  
  cbp1 <- c('Imputed' = "#999999", 'H3Africa' = "#E69F00", 'Add-On' = "#D55E00")
  
  merged_info <- rbind(baseline_plot_info,setting1_random_plot_info,setting1_plot_info)
  merged_info$Pass <- factor(merged_info$Pass,levels = c('H3Africa','H3Africa + Random','H3Africa + Add-On'))
  
  p5 <- ggplot(data=merged_info) + aes(x=Position/100000,y=INFO,color = SNP_Type) + 
    geom_point(alpha = 0.8) + ylim(0,1) + scale_colour_manual(values = cbp1) + 
    facet_grid(~Pass) + xlab('Chr10 Position (Mb)') + ylab('INFO Score')  + 
    scale_x_continuous(breaks = seq(519.4,520.2,0.4), labels = seq(519.4,520.2,0.4),limits = c(519.3,520.2))  + 
    labs(colour="SNP Type")

  
  baseline_plot_r2$Pass <- 'H3Africa'
  setting1_plot_r2$Pass <- 'H3Africa + Add-On'
  setting1_random_plot_r2$Pass <- 'H3Africa + Random'
  
  merged_r2 <- rbind(baseline_plot_r2,setting1_random_plot_r2,setting1_plot_r2)
  merged_r2$Pass <- factor(merged_r2$Pass,levels = c('H3Africa','H3Africa + Random','H3Africa + Add-On'))
  
  p6 <- ggplot(data=merged_r2) + aes(x=Position/100000,y=r2,color = SNP_Type) + 
    geom_point(alpha = 0.8) + ylim(0,1) + scale_colour_manual(values = cbp1) + facet_grid(~Pass) + 
    xlab('Chr10 Position (Mb)') + ylab(TeX('$r^{2}$')) + 
    scale_x_continuous(breaks = seq(519.4,520.2,0.4), labels = seq(519.4,520.2,0.4),limits = c(519.3,520.2)) + 
    labs(colour="SNP Type")
  
  return(ggpubr::ggarrange(p5,p6,nrow = 2,common.legend = T))
}

#### Fig 3

regions_improved_setting1 <- readRDS('../results/WGS_Addtl_Tags/addlt.tags.setting1.rds')
addtl_tags_setting1 <- setdiff(data.table::fread('../results/WGS_Addtl_Tags/h3a.with.addlt.tags.setting1.txt',header = F)$V1,
                      data.table::fread('../data/WGS_Host_Data/joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz.pos',header=F)$V3)
addtl_tags_random_setting1 <- setdiff(data.table::fread('../results/WGS_Addtl_Tags/h3a.with.addlt.tags.setting1.random.txt',header = F)$V1,
                      data.table::fread('../data/WGS_Host_Data/joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz.pos',header=F)$V3)

Setting1_Plot = GetImputationDiff('../results/Baseline/Tag_SNP_Selection/TB_gene_regions_baseline.rds',
                            '../results/WGS_Addtl_Tags/Setting1/TB_gene_regions_Setting1.rds',
                            '../results/WGS_Addtl_Tags/Setting1Random/TB_gene_regions_Setting1Random.rds',regions_improved_setting1,
                            rbind(data.table::fread('../results/Tag_SNP_Selection/WGS.Baseline.Illumina.Scores.csv'),
                                  data.table::fread('../results/Tag_SNP_Selection/WGS.KIR.HLA.Illumina.Scores.csv'),
                                  data.table::fread('../results/Tag_SNP_Selection/WGS.TBSNPs.Illumina.Scores.csv')),'setting1',addtl_tags_setting1,addtl_tags_random_setting1)

tags_incl <- data.table::fread('../results/WGS_Addtl_Tags/h3a.with.addlt.tags.setting2.txt',header = F)
regions_improved_setting2 <- readRDS('../results/WGS_Addtl_Tags/addlt.tags.setting2.rds')
regions_incl <- which(regions_improved_setting2$tags_to_add %in% tags_incl$V1)

regions_improved_setting2$tags_to_add <- regions_improved_setting2$tags_to_add[regions_incl]
regions_improved_setting2$tag_info <- regions_improved_setting2$max_tag_info[regions_incl]

addtl_tags_setting2 <- setdiff(data.table::fread('../results/WGS_Addtl_Tags/h3a.with.addlt.tags.setting2.txt',header = F)$V1,
                      data.table::fread('../results/WGS_Addtl_Tags/h3a.with.addlt.tags.setting1.txt',header=F)$V1)

addtl_tags_random_setting2 <- setdiff(data.table::fread('../results/WGS_Addtl_Tags/h3a.with.addlt.tags.setting2.random.txt',header = F)$V1,
                             data.table::fread('../results/WGS_Addtl_Tags/h3a.with.addlt.tags.setting1.random.txt',header=F)$V1)


Setting2_Plot = GetImputationDiff('../results/WGS_Addtl_Tags/Setting2/TB_gene_regions_Setting2_baseline.rds',
                                  '../results/WGS_Addtl_Tags/Setting2/TB_gene_regions_Setting2.rds',
                                  '../results/WGS_Addtl_Tags/Setting2Random/TB_gene_regions_Setting2Random.rds',regions_improved_setting2,
                                  rbind(data.table::fread('../results/Tag_SNP_Selection/WGS.MT.Illumina.Scores.csv'),
                                        data.table::fread('../results/Tag_SNP_Selection/WGS.Setting2.Illumina.Scores.csv'),
                                        data.table::fread('../results/Tag_SNP_Selection/WGS.ExtraSNPs.Illumina.Scores.csv'),
                                        data.table::fread('../results/Tag_SNP_Selection/WGS.KIR.HLA.Illumina.Scores.csv'),
                                        data.table::fread('../results/Tag_SNP_Selection/WGS.TBSNPs.Illumina.Scores.csv'),
                                        data.table::fread('../results/Tag_SNP_Selection/WGS.Baseline.Illumina.Scores.csv')),'setting2',addtl_tags_setting2,addtl_tags_random_setting2)

library(latex2exp)
Setting1_Plot$INFO_df$Setting <- 'Setting1'
Setting2_Plot$INFO_df$Setting <-  'Setting2'
Setting1_Plot$R2_df$Setting <- 'Setting1'
Setting2_Plot$R2_df$Setting <-  'Setting2'

INFO_df <- rbind(Setting1_Plot$INFO_df,Setting2_Plot$INFO_df)
R2_df <- rbind(Setting1_Plot$R2_df,Setting2_Plot$R2_df)


r2 = TeX('Mean $r^{2}$')
info_plot_maf <- ggplot(data = INFO_df) +
  aes(x=MAF,y=INFO,color = Chip) + geom_line() + geom_point() + guides(color=guide_legend(title="Array Content")) + ylim(0,1) + xlab('Minor Allele Frequency (MAF)') + ylab('Mean INFO Score') + facet_grid(~Setting)
r2_plot_maf <- ggplot(data = R2_df) +
  aes(x=MAF,y=r2,color = Chip) + geom_line() + geom_point() + guides(color=guide_legend(title="Array Content")) + ylim(0,1) + xlab('Minor Allele Frequency (MAF)') + ylab(r2) + facet_grid(~Setting)

fig3 <- ggpubr::ggarrange(info_plot_maf,r2_plot_maf,
                         common.legend = T,nrow = 2)

## Fig 4
fig4 <- PlotCovExample(readRDS('../results/Baseline/Tag_SNP_Selection/TB_gene_regions_baseline.rds')$snps_in_gene_regions_parsed,
                         readRDS('../results/WGS_Addtl_Tags/Setting1/TB_gene_regions_Setting1.rds')$snps_in_gene_regions_parsed,
                         readRDS('../results/WGS_Addtl_Tags/Setting1Random/TB_gene_regions_Setting1Random.rds')$snps_in_gene_regions_parsed)


#### Tbl 1 and 2
regions_improved_setting1 <- readRDS('../results/WGS_Addtl_Tags/addlt.tags.setting1.rds')
addtl_tags_tagger_setting1 <- setdiff(data.table::fread('../results/Tagger/Addtl_Tags/h3a.with.addlt.tags.setting1.tagger.txt',header = F)$V1,
                      data.table::fread('../data/WGS_Host_Data/joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz.pos',header=F)$V3)

Setting1_Tagger_tbl = GetImputationDiffTagger('../results/Baseline/Tag_SNP_Selection/TB_gene_regions_baseline.rds',
                            '../results/Tagger/Addtl_Tags/Setting1/TB_gene_regions_Setting1_Tagger.rds',regions_improved_setting1,
                            rbind(data.table::fread('../results/Tag_SNP_Selection/WGS.Baseline.Illumina.Scores.csv'),
                                  data.table::fread('../results/Tag_SNP_Selection/WGS.KIR.HLA.Illumina.Scores.csv'),
                                  data.table::fread('../results/Tag_SNP_Selection/WGS.TBSNPs.Illumina.Scores.csv')),'setting1',addtl_tags_tagger_setting1)

addtl_tags_tagger_setting2 <- setdiff(data.table::fread('../results/Tagger/Addtl_Tags/h3a.with.addlt.tags.setting2.tagger.txt',header = F)$V1,
                             data.table::fread('../results/Tagger/Addtl_Tags/h3a.with.addlt.tags.setting1.tagger.txt',header=F)$V1)

Setting2_Tagger_tbl = GetImputationDiffTagger('../results/Tagger/Addtl_Tags/Setting2/TB_gene_regions_Setting2_Tagger_baseline.rds',
                                              '../results/Tagger/Addtl_Tags/Setting2/TB_gene_regions_Setting2_Tagger.rds',regions_improved_setting2,
                                              rbind(data.table::fread('../results/Tag_SNP_Selection/WGS.MT.Illumina.Scores.csv'),
                                                    data.table::fread('../results/Tag_SNP_Selection/WGS.Setting2.Illumina.Scores.csv'),
                                                    data.table::fread('../results/Tag_SNP_Selection/WGS.ExtraSNPs.Illumina.Scores.csv'),
                                                    data.table::fread('../results/Tag_SNP_Selection/WGS.KIR.HLA.Illumina.Scores.csv'),
                                                    data.table::fread('../results/Tag_SNP_Selection/WGS.TBSNPs.Illumina.Scores.csv'),
                                                    data.table::fread('../results/Tag_SNP_Selection/WGS.Baseline.Illumina.Scores.csv')),'setting2',addtl_tags_tagger_setting2)

setting1_all_snps <- Setting1_Plot$tbl[Setting1_Plot$tbl$type == 'All SNPs',]
setting1_tagger_all_snps <- Setting1_Tagger_tbl$tbl[Setting1_Tagger_tbl$tbl$type == 'All SNPs',]
setting1_AFGR <- Setting1_Plot$tbl[Setting1_Plot$tbl$type == 'AFGR SNPs',]
setting1_tagger_AFGR <- Setting1_Tagger_tbl$tbl[Setting1_Tagger_tbl$tbl$type == 'AFGR SNPs',]
setting1_WGS <- Setting1_Plot$tbl[Setting1_Plot$tbl$type == 'WGS SNPs',]
setting1_tagger_WGS <- Setting1_Tagger_tbl$tbl[Setting1_Tagger_tbl$tbl$type == 'WGS SNPs',]

setting2_all_snps <- Setting2_Plot$tbl[Setting2_Plot$tbl$type == 'All SNPs',]
setting2_tagger_all_snps <- Setting2_Tagger_tbl$tbl[Setting2_Tagger_tbl$tbl$type == 'All SNPs',]
setting2_AFGR <- Setting2_Plot$tbl[Setting2_Plot$tbl$type == 'AFGR SNPs',]
setting2_tagger_AFGR <- Setting2_Tagger_tbl$tbl[Setting2_Tagger_tbl$tbl$type == 'AFGR SNPs',]
setting2_WGS <- Setting2_Plot$tbl[Setting2_Plot$tbl$type == 'WGS SNPs',]
setting2_tagger_WGS <- Setting2_Tagger_tbl$tbl[Setting2_Tagger_tbl$tbl$type == 'WGS SNPs',]


tbl1 <- rbind(data.frame(n_beads_added = c(setting1_all_snps$n_beads_added,setting1_tagger_all_snps$n_beads_added,NA),
                            mean_design_score = c(setting1_all_snps$mean_design_score,setting1_tagger_all_snps$mean_design_score,NA),
                            n_improved_per_bead = c(setting1_all_snps$n_improved / setting1_all_snps$n_beads_added,setting1_tagger_all_snps$n_improved_tagger/setting1_tagger_all_snps$n_beads_added,NA),
                            n_improved_per_tag = c(setting1_all_snps$n_improved / setting1_all_snps$n_snps_added,setting1_tagger_all_snps$n_improved_tagger/setting1_tagger_all_snps$n_snps_added, setting1_all_snps$n_improved_random / setting1_all_snps$n_snps_added),
                            n_improved_AFGR = c(setting1_AFGR$n_improved / setting1_all_snps$n_improved,setting1_tagger_AFGR$n_improved_tagger/setting1_tagger_all_snps$n_improved_tagger, setting1_AFGR$n_improved_random / setting1_all_snps$n_improved_random),
                            n_improved_WGS = c(setting1_WGS$n_improved / setting1_all_snps$n_improved,setting1_tagger_WGS$n_improved_tagger/setting1_tagger_all_snps$n_improved_tagger, setting1_WGS$n_improved_random / setting1_all_snps$n_improved_random),
                            n_improved_0p8_per_bead =  c(setting1_all_snps$n_improved_0p8 / setting1_all_snps$n_beads_added,setting1_tagger_all_snps$n_improved_0p8_tagger/setting1_tagger_all_snps$n_beads_added,NA),
                            n_improved_0p8_per_tag = c(setting1_all_snps$n_improved_0p8 / setting1_all_snps$n_snps_added,setting1_tagger_all_snps$n_improved_0p8_tagger/setting1_tagger_all_snps$n_snps_added, setting1_all_snps$n_improved_0p8_random / setting1_all_snps$n_snps_added),
                            n_improved_0p8_AFGR = c(setting1_AFGR$n_improved_0p8 / setting1_all_snps$n_improved_0p8,setting1_tagger_AFGR$n_improved_0p8_tagger/setting1_tagger_all_snps$n_improved_0p8_tagger, setting1_AFGR$n_improved_0p8_random / setting1_all_snps$n_improved_0p8_random),
                            n_improved_0p8_WGS = c(setting1_WGS$n_improved_0p8 / setting1_all_snps$n_improved_0p8,setting1_tagger_WGS$n_improved_0p8_tagger/setting1_tagger_all_snps$n_improved_0p8_tagger, setting1_WGS$n_improved_0p8_random / setting1_all_snps$n_improved_0p8_random)),
              data.frame(n_beads_added = c(setting2_all_snps$n_beads_added,setting2_tagger_all_snps$n_beads_added,NA),
                         mean_design_score = c(setting2_all_snps$mean_design_score,setting2_tagger_all_snps$mean_design_score,NA),
                         n_improved_per_bead = c(setting2_all_snps$n_improved / setting2_all_snps$n_beads_added,setting2_tagger_all_snps$n_improved_tagger/setting2_tagger_all_snps$n_beads_added,NA),
                         n_improved_per_tag = c(setting2_all_snps$n_improved / setting2_all_snps$n_snps_added,setting2_tagger_all_snps$n_improved_tagger/setting2_tagger_all_snps$n_snps_added, setting2_all_snps$n_improved_random / setting2_all_snps$n_snps_added),
                         n_improved_AFGR = c(setting2_AFGR$n_improved / setting2_all_snps$n_improved,setting2_tagger_AFGR$n_improved_tagger/setting2_tagger_all_snps$n_improved_tagger, setting2_AFGR$n_improved_random / setting2_all_snps$n_improved_random),
                         n_improved_WGS = c(setting2_WGS$n_improved / setting2_all_snps$n_improved,setting2_tagger_WGS$n_improved_tagger/setting2_tagger_all_snps$n_improved_tagger, setting2_WGS$n_improved_random / setting2_all_snps$n_improved_random),
                         n_improved_0p8_per_bead =  c(setting2_all_snps$n_improved_0p8 / setting2_all_snps$n_beads_added,setting2_tagger_all_snps$n_improved_0p8_tagger/setting2_tagger_all_snps$n_beads_added,NA),
                         n_improved_0p8_per_tag = c(setting2_all_snps$n_improved_0p8 / setting2_all_snps$n_snps_added,setting2_tagger_all_snps$n_improved_0p8_tagger/setting2_tagger_all_snps$n_snps_added, setting2_all_snps$n_improved_0p8_random / setting2_all_snps$n_snps_added),
                         n_improved_0p8_AFGR = c(setting2_AFGR$n_improved_0p8 / setting2_all_snps$n_improved_0p8,setting2_tagger_AFGR$n_improved_0p8_tagger/setting2_tagger_all_snps$n_improved_0p8_tagger, setting2_AFGR$n_improved_0p8_random / setting2_all_snps$n_improved_0p8_random),
                         n_improved_0p8_WGS = c(setting2_WGS$n_improved_0p8 / setting2_all_snps$n_improved_0p8,setting2_tagger_WGS$n_improved_0p8_tagger/setting2_tagger_all_snps$n_improved_0p8_tagger, setting2_WGS$n_improved_0p8_random / setting2_all_snps$n_improved_0p8_random)
              ))
rownames(tbl1) <- c('Setting 1 - Proposed Approach','Setting 1 - Tagger','Setting 1 - Random',
                             'Setting 2 - Proposed Approach','Setting 2 - Tagger','Setting 2 - Random')


#### Fig S3
# PlotInfoH3A <- function(info_rds){
#   info_file <- readRDS(info_rds)
#   info_file <- lapply(info_file,function(x) {x$MAF = pmin(x$AFGR_AF,1-x$AFGR_AF); return(x)})
#   MAF_bins <- quantile(info_file$TBDAR$MAF,probs = c(seq(0,0.4,0.1),seq(0.51,1,0.07)))
# 
#   Binned_INFO <- lapply(info_file,function(df) pbmcmapply(function(x,y) {
#     curBin <- dplyr::filter(df,MAF > x & MAF <= y)
#     return(curBin$INFO)
#   },MAF_bins[1:(length(MAF_bins) - 1)],MAF_bins[2:(length(MAF_bins))],SIMPLIFY = F,mc.cores = 5))
# 
#   Binned_INFO_df <- lapply(Binned_INFO,function(x) data.frame(MAF = (MAF_bins[1:(length(MAF_bins) - 1)] + MAF_bins[2:(length(MAF_bins))]) / 2,
#                                                               INFO = sapply(x,function(y) sum(y > 0.8)/length(y))))
#   Binned_INFO_df_bind <- dplyr::bind_rows(Binned_INFO_df,.id = "Population")
#   Binned_INFO_df_bind$Population[Binned_INFO_df_bind$Population == 'TBDAR'] <- 'TB-DAR'
#   Binned_INFO_df_bind$Population <- factor(Binned_INFO_df_bind$Population,levels = c('TB-DAR',setdiff(Binned_INFO_df_bind$Population,'TB-DAR')))
#   cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
#             "#0072B2", "#D55E00", "#CC79A7")
# 
#   return(ggplot(Binned_INFO_df_bind,aes(x=MAF,y=INFO,color = Population)) + geom_point() + geom_line() + xlab('Minor Allele Frequency (MAF)') + ylab('Fraction of Variant Sites Imputed\n (INFO > 0.8)') + scale_color_manual(values = cbp2))
# 
# }
# autosomes_info <- PlotInfoH3A('../results/1000_Genomes/info_file.rds')
# x_info <- PlotInfoH3A('../results/1000_Genomes/info_file_X.rds')
# ggarrange(autosomes_info + ggtitle('Autosomes'),x_info + ggtitle('X Chromosome'),labels = c('A)','B)'))
