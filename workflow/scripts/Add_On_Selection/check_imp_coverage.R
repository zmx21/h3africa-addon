#Check imputation accuracy within defined regions. 
library(dplyr)
library(biomaRt)
library(pbmcapply)
library(ggplot2)
library(ggrepel)
library(IRanges)
load('../results/Tag_SNP_Selection/setting2.rda')
remove(list = setdiff(ls(),'snps_in_region'))

#Get prioritized regions (Associated with TB)
GetTBRegion <- function(tb_gwas_dir = '../data/TB_Targets/'){
  #Open Targest association threshold
  assoc_score_filt <- 0.1
  tb_open_targets <- data.table::fread(paste0(tb_gwas_dir,'targets_associated_with_Tuberculosis.csv')) %>%
    dplyr::filter(association_score.overall > assoc_score_filt)
  tb_open_targets_broad <- data.table::fread(paste0(tb_gwas_dir,'targets_associated_with_Tuberculosis.csv')) %>%
    dplyr::filter(association_score.overall < assoc_score_filt & association_score.overall > 0.01)
  
  tb_gwas <- data.table::fread(paste0(tb_gwas_dir,'gwas-association-downloaded_2019-08-06-Orphanet_3389.tsv')) %>% dplyr::filter(`DISEASE/TRAIT`!="Tuberculosis (SNP x SNP interaction)")
  hla_genes <- data.table::fread(paste0(tb_gwas_dir,'hla_genes.csv'),header = T)
  kir_genes <- c('KIR2DL1', 'KIR2DL2', 'KIR2DL3', 'KIR2DL4', 'KIR2DL5A', 
                 'KIR2DL5B','KIR2DS1', 'KIR2DS2', 'KIR2DS3', 'KIR2DS4', 'KIR2DS5',
                 'KIR3DL1', 'KIR3DL2', 'KIR3DL3','KIR3DS1')
  tb_gwas_snps <- unlist(lapply(tb_gwas$SNPS,function(x) strsplit(x,split = ' x ')[[1]]),recursive = T)
  write(tb_gwas_snps,file = '../data/tb_gwas_snps.txt')
  #Get gene positions using Biomart (GRCh37)
  grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
  
  #Defind regions of interst, 5KB upstream and downstream of gene of interest
  #Get opentargets region
  offset <- 5000
  tb_open_target_genes <- tb_open_targets$target.id
  tb_open_target_genes_broad <- tb_open_targets_broad$target.id
  
  hla_genes <- hla_genes$`Approved symbol`
  tb_open_targets_regions <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol","chromosome_name", "start_position","end_position","strand"),
                                   filters = "ensembl_gene_id", values = tb_open_target_genes,mart = grch37)%>% dplyr::mutate(Focused = T)
  tb_open_targets_regions_broad <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol","chromosome_name", "start_position","end_position","strand"),
                                   filters = "ensembl_gene_id", values = tb_open_target_genes_broad,mart = grch37) %>% dplyr::mutate(Focused = F)
  
  #Get HLA regions, those which are not in open targets
  hla_regions <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol","chromosome_name", "start_position","end_position","strand"),
                                   filters = "hgnc_symbol", values = hla_genes,mart = grch37)
  hla_regions <- dplyr::filter(hla_regions,chromosome_name=='6') %>% dplyr::filter(!hgnc_symbol %in% tb_open_targets_regions$hgnc_symbol) %>% dplyr::mutate(Focused = T)
  
  #Get KIR region, those which are not in open targets
  kir_regions <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol","chromosome_name", "start_position","end_position","strand"),
                       filters = "hgnc_symbol", values = kir_genes,mart = grch37)
  kir_regions <- dplyr::filter(kir_regions,chromosome_name=='19') %>% dplyr::filter(!hgnc_symbol %in% tb_open_targets_regions$hgnc_symbol) %>% dplyr::mutate(Focused = T)

  tb_open_targets_regions <- rbind(tb_open_targets_regions,hla_regions,kir_regions,tb_open_targets_regions_broad)
  tb_open_targets_regions <- dplyr::mutate(tb_open_targets_regions,start_region = start_position - offset,end_region = end_position + offset)
  
  #Parse gwas genes into correct format
  gwas_reported_genes <- dplyr::filter(tb_gwas,`REPORTED GENE(S)` != 'NR' & `REPORTED GENE(S)` != 'intergenic')$`REPORTED GENE(S)`
  gwas_reported_genes <- unlist(lapply(gwas_reported_genes,function(x) strsplit(x,split = ' x ')[[1]]),recursive = T)
  gwas_reported_genes <- unlist(lapply(gwas_reported_genes,function(x) strsplit(x,split = ', ')[[1]]),recursive = T)
  gwas_reported_genes <- sapply(gwas_reported_genes,function(x) gsub(x=x,pattern = ' ',replacement = ''),USE.NAMES = F)
  gwas_reported_genes <- unique(gwas_reported_genes)

  gwas_mapped_genes <- unlist(lapply(tb_gwas$MAPPED_GENE,function(x) strsplit(x,split = ' - ')[[1]]),recursive = T)
  gwas_mapped_genes <- unlist(lapply(gwas_mapped_genes,function(x) strsplit(x,split = ' x ')[[1]]),recursive = T)
  gwas_mapped_genes <- unlist(lapply(gwas_mapped_genes,function(x) strsplit(x,split = ', ')[[1]]),recursive = T)
  
  #Add GWAS genes and genes from recent study (https://www.nature.com/articles/s41467-019-11664-1)
  #Add TB susceptibility and interaction genes from review paper (https://doi.org/10.3389/fgene.2019.00865)
  gwas_genes <- union(gwas_reported_genes,
                      c(gwas_mapped_genes,'ATP1B3','UBLCP1','CCDC6','IL12A','TOX',
                        'IKBKG','GC','ATG9B','ATG2B','LAMP3','WIPI1','ATG4C','LRIF1','IL6','TYK2'))
  gwas_regions <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol","chromosome_name", "start_position","end_position","strand"),
                        filters = "hgnc_symbol", values = gwas_genes,mart = grch37)
  gwas_regions <- dplyr::filter(gwas_regions,chromosome_name %in% c(as.character(1:22),'X','Y'))
  gwas_regions <- dplyr::mutate(gwas_regions,start_region = start_position - offset,end_region = end_position + offset) %>% dplyr::mutate(Focused = T)
  
  #Use both open targets and GWAS to extract relevant imputed SNPs
  tb_regions <- rbind(tb_open_targets_regions,dplyr::filter(gwas_regions,!hgnc_symbol %in% tb_open_targets_regions$hgnc_symbol))
  
  #Add intergenic SNPs in TB GWAS 
  tb_gwas_SNPs <- unlist(lapply(tb_gwas$SNPS,function(x) strsplit(x,split = ' x ')[[1]]),recursive = T)
  tb_gwas_SNPs <- c(tb_gwas_SNPs,'rs73226617','rs72993272')
  grch37_snps = useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_snp")
  tb_gwas_SNPs_df <- getBM(attributes = c('refsnp_id','chr_name','chrom_strand','chrom_start','chrom_strand'), 
        filters = c('snp_filter'), 
        values = tb_gwas_SNPs, 
        mart = grch37_snps)
  tb_gwas_SNPs_df <- dplyr::filter(tb_gwas_SNPs_df,chr_name %in% c('X','Y',as.character(seq(1,22,1))))
  tb_gwas_SNPs_df <- data.frame(ensembl_gene_id=tb_gwas_SNPs_df$refsnp_id,
                                           hgnc_symbol=tb_gwas_SNPs_df$refsnp_id,
                                           chromosome_name=tb_gwas_SNPs_df$chr_name,
                                           start_position=tb_gwas_SNPs_df$chrom_start,
                                           end_position=tb_gwas_SNPs_df$chrom_start,
                                           strand = tb_gwas_SNPs_df$chrom_strand.1,
                                           start_region = tb_gwas_SNPs_df$chrom_start - offset,
                                           end_region = tb_gwas_SNPs_df$chrom_start + offset) %>% dplyr::mutate(Focused = T)
  tb_regions <- rbind(tb_regions,tb_gwas_SNPs_df)
  
  #Make sure regions are non_overlapping
  tb_regions_unique <- tb_regions[0,]
  for(i in 1:length(unique(tb_regions$chromosome_name))){
    cur_chr <- dplyr::filter(tb_regions,chromosome_name == unique(tb_regions$chromosome_name)[i])
    cur_chr_non_overlapping <- cur_chr[0,]
    overlaps <- IRanges(cur_chr$start_region,cur_chr$end_region)
    group_assignment <- subjectHits(findOverlaps(overlaps, reduce(overlaps)))
    for(j in group_assignment){
      cur_group <- cur_chr[which(group_assignment==j),]
      cur_group_unique <- cur_group[1,]
      cur_group_unique$start_region <- min(cur_group$start_region)
      cur_group_unique$end_region <- max(cur_group$end_region)
      cur_group_unique$ensembl_gene_id <- paste0(cur_group$ensembl_gene_id,collapse = ',')
      cur_group_unique$hgnc_symbol <- paste0(cur_group$hgnc_symbol,collapse = ',')
      cur_group_unique$chromosome_name <- unique(tb_regions$chromosome_name)[i]
      cur_group_unique$start_position <- paste0(cur_group$start_position,collapse = ',')
      cur_group_unique$end_position <- paste0(cur_group$end_position,collapse = ',')
      cur_group_unique$strand <- paste0(cur_group$strand,collapse = ',')
      cur_group_unique$Focused <- any(cur_group_unique$Focused)
      cur_chr_non_overlapping <- rbind(cur_chr_non_overlapping,cur_group_unique)
    }
    tb_regions_unique <- rbind(tb_regions_unique,cur_chr_non_overlapping)
  }
  tb_regions_unique_distinct <- dplyr::arrange(tb_regions_unique, desc(Focused)) %>% dplyr::distinct(chromosome_name,start_region,end_region,.keep_all=T)
  tb_regions_unique_distinct <- dplyr::filter(tb_regions_unique_distinct,chromosome_name %in% c(as.character(seq(1,22,1)),'X'))
  
  #Parse out centromere and x-pseudoautosomal region
  #data(pseudoautosomal.hg19)
  tb_regions_unique_distinct_filt <- tb_regions_unique_distinct #%>% dplyr::filter(chromosome_name!= 'X' | 
                                                                #               (chromosome_name == 'X' & !((start_region >= pseudoautosomal.hg19$start.base[3] & end_region <= pseudoautosomal.hg19$end.base[3]) | 
                                                                #               (start_region >= pseudoautosomal.hg19$start.base[2] & end_region <= pseudoautosomal.hg19$end.base[2]))))
  return(list(tb_regions_unique=tb_regions_unique_distinct_filt,tb_gwas_snps=tb_gwas_snps))
}

#Extract imputed SNPs within defined regions (From both Tanz and AFGR reference panel). 
#Also extract existing tag SNPs and any sequenced SNPs
GetTargetRegionCoverage <- function(region_info,tb_gwas_snps,afgr_imp_path,addtl_snp_path = '',h3a_path,wgs_path,wgs_excl,wgs_dup,merged = F,n_cores=20){
  afgr_sites_path <- '../data/AFGRTesting.vcfs/agr.minAC1.sites.vcf.gz'
  if(!file.exists(afgr_sites_path)){
    system(paste0('wget -p ',afgr_sites_path,' ftp://ngs.sanger.ac.uk/production/hrc/sangerimpute/agr/sites/v1/agr.minAC1.sites.vcf.gz'))
    system(paste0('~/Software/bcftools index --threads 5 ',afgr_sites_path))
    
  }
  
  if(addtl_snp_path != ''){
    addtl_tags_in_region <- pbmclapply(1:nrow(region_info),function(i) 
      system(paste0('~/Software/bcftools view -r ',region_info$chromosome_name[i],':',region_info$start_region[i],'-',region_info$end_region[i],' ',addtl_snp_path,
                    " | ~/Software/bcftools query -f '%CHROM %POS %ID %REF %ALT\n'"),intern = T),mc.cores = n_cores)
    
  }else{
    addtl_tags_in_region <- lapply(1:nrow(region_info),function(x) c())
  }
  
  if(merged == F){
    afgr_imputed_snps_in_gene_regions <- pbmclapply(1:nrow(region_info),function(i) system(paste0('~/Software/bcftools view -r ',region_info$chromosome_name[i],':',region_info$start_region[i],'-',region_info$end_region[i],' ',afgr_imp_path,
                                                                                                        " | ~/Software/bcftools query -f '%CHROM %POS %ID %REF %ALT %INFO/INFO %INFO/RefPanelAF\n'"),intern = T),mc.cores = n_cores)
  }else{
    afgr_imputed_snps_in_gene_regions <- pbmclapply(1:nrow(region_info),
                                                    function(i){ 
                                                      system(paste0("awk '{ if(($2 >= ",as.character(region_info$start_region[i]),') && ($2 <= ',as.character(region_info$end_region[i]),")) { print } }' ",
                                                                    gsub(x=afgr_imp_path,pattern = '.txt',replacement = paste0('_chr',region_info$chromosome_name[i],'.txt'))),intern=T)
                                                    }
  ,mc.cores = n_cores)
  }
  h3a_tags_in_gene_regions <- pbmclapply(1:nrow(region_info),function(i) system(paste0('~/Software/bcftools view -r ',region_info$chromosome_name[i],':',region_info$start_region[i],'-',region_info$end_region[i],' ',h3a_path,
                                                                                      " | ~/Software/bcftools query -f '%CHROM %POS %ID %REF %ALT\n'"),intern = T),mc.cores = n_cores)
  
  afgr_sites_in_gene_regions <- pbmclapply(1:nrow(region_info),function(i) system(paste0('~/Software/bcftools view -r ',region_info$chromosome_name[i],':',region_info$start_region[i],'-',region_info$end_region[i],' ',afgr_sites_path,
                                                                                        " | ~/Software/bcftools query -f '%CHROM %POS %ID %REF %ALT %INFO/AF\n'"),intern = T),mc.cores = n_cores)
  
  wgs_snps_in_gene_regions <- pbmclapply(1:nrow(region_info),function(i) system(paste0('~/Software/bcftools view -r ',region_info$chromosome_name[i],':',region_info$start_region[i],'-',region_info$end_region[i],' ',wgs_path,
                                                                                      " | ~/Software/bcftools +fill-tags -- -t AF | ~/Software/bcftools query -f '%CHROM %POS %ID %REF %ALT %INFO/AF\n'"),intern = T),mc.cores = n_cores)
  
  # grch38_snps <- data.table::fread('~/G2G_TB/WGS_Host_Data/grch38_wgs_snps.txt',header=F)
  # colnames(grch38_snps) <- c('ID','GRCH38.REF','GRCH38.ALT')
  #Parse ~/Software/bcftools output to correct format
  snps_in_gene_regions_parsed <- pbmclapply(1:nrow(region_info),function(i){
    start_region <- region_info$start_region[i]
    end_region <- region_info$end_region[i]
    region_length <- end_region - start_region
    region_id <- region_info$hgnc_symbol[i]
    focused <- region_info$Focused[i]
    if(length(afgr_imputed_snps_in_gene_regions[[i]] > 0)){
      imp_snps <- data.table::fread(input = paste0(afgr_imputed_snps_in_gene_regions[[i]],'\n',collapse = '\n'),fill = T,stringsAsFactors = F) %>% dplyr::filter(!is.na(V2))
      if(merged == F){
        colnames(imp_snps) <- c('CHROM','POS','ID','REF','ALT','INFO','AF')
      }else{
        colnames(imp_snps) <- c('CHROM','POS','ID','REF','ALT','INFO','AF','Source')
      }
      imp_snps$CHROM <- as.character(imp_snps$CHROM)
      imp_snps$MAF <- sapply(imp_snps$AF,function(x) min(x,1-x))
      #imp_snps <-  dplyr::filter(imp_snps,MAF > 0.05)
      missing_id <- which(imp_snps$ID == '.')
      if(length(missing_id) > 0){
        imp_snps$ID[missing_id] <- sapply(missing_id,function(i) paste0(imp_snps$CHROM[i],':',imp_snps$POS[i],'[b37]',imp_snps$REF[i],',',imp_snps$ALT[i]))
      }
    }else{
      if(merged == F){
        imp_snps <- data.frame('CHROM' = as.character(NA),'POS' = as.integer(NA),'ID' = as.character(NA),'REF'=as.character(NA),'ALT'=as.character(NA),'INFO'=as.integer(NA),'AF'=as.integer(NA),'MAF'=as.integer(NA),stringsAsFactors = F)
      }else{
        imp_snps <- data.frame('CHROM' = as.character(NA),'POS' = as.integer(NA),'ID' = as.character(NA),'REF'=as.character(NA),'ALT'=as.character(NA),'INFO'=as.integer(NA),'AF'=as.integer(NA),'MAF'=as.integer(NA),'Source' = as.character(NA),stringsAsFactors = F)
      }
    }
    
    if(length(afgr_sites_in_gene_regions[[i]] > 0)){
      afgr_sites <- data.table::fread(input = paste0(afgr_sites_in_gene_regions[[i]],'\n',collapse = '\n'),fill = T,stringsAsFactors = F) %>% dplyr::filter(!is.na(V2))
      colnames(afgr_sites) <- c('CHROM','POS','ID','REF','ALT','AF')
      afgr_sites$CHROM <- as.character(afgr_sites$CHROM)
      afgr_sites$MAF <- sapply(afgr_sites$AF,function(x) min(x,1-x))
      afgr_sites <-  dplyr::filter(afgr_sites,MAF > 0.05)
    }else{
      afgr_sites <- data.frame('CHROM' = as.character(NA),'POS' = as.integer(NA),'ID' = as.character(NA),'REF'=as.character(NA),'ALT'=as.character(NA),'INFO'=as.integer(NA),'AF'=as.integer(NA),'MAF'=as.integer(NA),stringsAsFactors = F)
    }
    
    if(length(h3a_tags_in_gene_regions[[i]] > 0)){
      h3a_tags <- data.table::fread(input = paste0(h3a_tags_in_gene_regions[[i]],'\n',collapse = '\n'),fill = T,stringsAsFactors = F) %>% dplyr::filter(!is.na(V2))
      colnames(h3a_tags) <- c('CHROM','POS','ID','REF','ALT')
      h3a_tags$CHROM <- as.character(h3a_tags$CHROM)
      missing_id <- which(h3a_tags$ID == '.')
      if(length(missing_id) > 0){
        h3a_tags$ID[missing_id] <- sapply(missing_id,function(i) paste0(h3a_tags$CHROM[i],':',h3a_tags$POS[i],'[b37]',h3a_tags$REF[i],',',h3a_tags$ALT[i]))
      }
    }else{
      h3a_tags <- data.frame('CHROM' = as.character(NA),'POS' = as.integer(NA),'ID' = as.character(NA),'REF'=as.character(NA),'ALT'=as.character(NA),stringsAsFactors = F)
    }
    
    if(length(wgs_snps_in_gene_regions[[i]] > 0)){
      wgs_snps <- data.table::fread(input = paste0(wgs_snps_in_gene_regions[[i]],'\n',collapse = '\n'),fill = T,stringsAsFactors = F) %>% dplyr::filter(!is.na(V2))
      colnames(wgs_snps) <- c('CHROM','POS','ID','REF','ALT','AF')
      wgs_snps$CHROM <- as.character(wgs_snps$CHROM)
      missing_id <- which(wgs_snps$ID == '.')
      if(length(missing_id) > 0){
        wgs_snps$ID[missing_id] <- sapply(missing_id,function(i) paste0(gsub(x = wgs_snps$CHROM[i],pattern = 'chr',replacement = ''),':',wgs_snps$POS[i],'[b37]',wgs_snps$REF[i],',',wgs_snps$ALT[i]))
      }
      cur_region_tb_gwas <- intersect(tb_gwas_snps,wgs_snps$ID) 
      
      wgs_snps$MAF <- sapply(wgs_snps$AF,function(x) min(x,1-x))
      #wgs_snps <- wgs_snps %>% dplyr::filter(MAF > 0.05)
      
      wgs_mismapped_snps <- data.table::fread(wgs_excl,header = F)$V1
      wgs_dup_snps <- data.table::fread(wgs_dup,header = F)
      wgs_snps <- wgs_snps %>% dplyr::filter(!(ID %in% wgs_mismapped_snps | ID %in% wgs_dup_snps))
      #Only include SNPs with rsid
      wgs_snps <- wgs_snps[which(sapply(wgs_snps$ID,function(x) grepl(x=x,'rs'))),]
      #Get GRCh38 POS,REF,ALT
      # grch38_snps <- useMart(biomart="ENSEMBL_MART_SNP", path="/biomart/martservice",dataset="hsapiens_snp")
      # wgs_grch38 <- getBM(attributes=c('refsnp_id','allele'), filters = 'snp_filter', values = wgs_snps$ID, mart = grch38_snps)
      # wgs_grch38 <- dplyr::select(wgs_grch38,ID='refsnp_id',GRCh38.Allele='allele')
      # wgs_snps <- dplyr::left_join(wgs_snps,grch38_snps,by=c('ID'='ID'))
      
      #Add tags outside of H3A
      if(addtl_snp_path != '' & length(addtl_tags_in_region[[i]] > 0)){
        addtl_snps <- data.table::fread(input = paste0(addtl_tags_in_region[[i]],'\n',collapse = '\n'),fill = T,stringsAsFactors = F) %>% dplyr::filter(!is.na(V2))
        colnames(addtl_snps) <- c('CHROM','POS','ID','REF','ALT')
        addtl_snps$CHROM <- as.character(addtl_snps$CHROM)
        addtl_snps <- dplyr::anti_join(addtl_snps,h3a_tags,by=c('CHROM'='CHROM','POS'='POS','REF'='REF','ALT'='ALT'))
        if(nrow(addtl_snps)==0){
          addtl_snps <- data.frame('CHROM' = as.character(NA),'POS' = as.integer(NA),'ID' = as.character(NA),'REF'=as.character(NA),'ALT'=as.character(NA),stringsAsFactors = F)
        }
      }else{
        addtl_snps <- data.frame('CHROM' = as.character(NA),'POS' = as.integer(NA),'ID' = as.character(NA),'REF'=as.character(NA),'ALT'=as.character(NA),stringsAsFactors = F)
      }
    }else{
      wgs_snps <- data.frame('CHROM' = as.character(NA),'POS' = as.integer(NA),'ID' = as.character(NA),'REF'=as.character(NA),'ALT'=as.character(NA),'INFO'=as.integer(NA),'AF'=as.integer(NA),'MAF'=as.integer(NA),'GRCh38.REF'=as.character(NA),'GRCh38.ALT'=as.character(NA),stringsAsFactors = F)
      addtl_snps <- data.frame('CHROM' = as.character(NA),'POS' = as.integer(NA),'ID' = as.character(NA),'REF'=as.character(NA),'ALT'=as.character(NA),stringsAsFactors = F)
    }
    
    return(list(start_region = start_region,end_region = end_region,region_length = region_length,region_id=region_id,focused=focused,
                imp_snps = imp_snps,afgr_sites=afgr_sites,h3a_tags=h3a_tags,wgs_snps=wgs_snps,addtl_snps=addtl_snps,
                cur_region_tb_gwas=cur_region_tb_gwas))
  }
  ,mc.cores = 15)
  return(snps_in_gene_regions_parsed)
}

#Percentage of Imputed SNPs with INFO > 0.8 (and MAF > 0.05)
CalcPercAbove <- function(snps_in_gene_regions_parsed){
  sapply(snps_in_gene_regions_parsed,function(x){
    x$imp_snps <- dplyr::filter(x$imp_snps,MAF > 0.05)
    if(all(is.na(x$imp_snps$ID))){
      return(0)
    }else{
      sum(x$imp_snps$INFO > 0.8)/nrow(x$imp_snps)
    }
  })
}

#Spread of well imputed and common imputed SNPs (INFO > 0.8 and MAF > 0.05)
CalcSpread <- function(snps_in_gene_regions_parsed){
  spread <- sapply(snps_in_gene_regions_parsed,function(x){
    x$imp_snps <- dplyr::filter(x$imp_snps,MAF > 0.05)
    x$wgs_snps <- dplyr::filter(x$wgs_snps,MAF > 0.05)
    
    if(all(is.na(x$imp_snps$ID))){
      return(0)
    }else{
      sd(x$imp_snps$POS[x$imp_snps$INFO > 0.8])/max(c(sd(x$wgs_snps$POS),sd(x$imp_snps$POS)),na.rm = T)
    }
  })
}

#Plot imputed, existing tag, and sequenced SNPs in a defined region. 
PlotRegion <- function(snps_in_gene_regions_parsed,index,MAF_Plot=T,title=F){
  cur_imp_snps <- snps_in_gene_regions_parsed[[index]]$imp_snps
  cur_imp_snps <- dplyr::filter(cur_imp_snps,MAF > 0.05)
  df <- data.frame(Position=cur_imp_snps$POS,INFO=cur_imp_snps$INFO,SNP_Type=ifelse(cur_imp_snps$POS %in% snps_in_gene_regions_parsed[[index]]$h3a_tags$POS,'H3A Tags',
                                                                                    ifelse(cur_imp_snps$POS %in% snps_in_gene_regions_parsed[[index]]$addtl_snps$POS,'Add-On Tags','Imputed')))
  cbp1 <- c('Imputed' = "#999999", 'H3A Tags' = "#E69F00", 'Add-On Tags' = "#D55E00")
  if(title){
    p1 <- ggplot(data=df) + aes(x=Position,y=INFO,color = SNP_Type) + geom_point(alpha = 0.8) + 
      xlim(c(snps_in_gene_regions_parsed[[index]]$start_region,snps_in_gene_regions_parsed[[index]]$end_region)) + ggtitle(paste0('Imputed SNPs: Region',index)) + ylim(0,1) + scale_colour_manual(values = cbp1)
  }else{
    p1 <- ggplot(data=df) + aes(x=Position/1e6,y=INFO,color = SNP_Type) + geom_point(alpha = 0.8) + 
      xlim(c(snps_in_gene_regions_parsed[[index]]$start_region/1e6,snps_in_gene_regions_parsed[[index]]$end_region/1e6)) + ggtitle('') + ylim(0,1) + scale_colour_manual(values = cbp1) + guides(color=guide_legend(title="SNP Type")) + xlab(paste0('Chr ',unique(cur_imp_snps$CHROM),' Position (Mb)'))
    
  }
  if(MAF_Plot){
    wgs_snps <- snps_in_gene_regions_parsed[[index]]$wgs_snps
    wgs_snps <- dplyr::filter(wgs_snps,MAF > 0.05)
    all_snps <- rbind(data.frame(POS = c(wgs_snps$POS,cur_imp_snps$POS),MAF = c(wgs_snps$MAF,cur_imp_snps$MAF),stringsAsFactors = F))
    p2_df <- data.frame(Position=wgs_snps$POS,MAF=wgs_snps$MAF,SNP_Type=ifelse(wgs_snps$POS %in% snps_in_gene_regions_parsed[[index]]$h3a_tags$POS,'H3A Tags',
                                                                               ifelse(wgs_snps$POS %in% snps_in_gene_regions_parsed[[index]]$addtl_snps$POS,'Add-On Tags',
                                                                                      ifelse(wgs_snps$POS %in% snps_in_gene_regions_parsed[[index]]$afgr_sites$POS,'AFGR SNPs','WGS SNPs'))))
    cbp1 <- c('AFGR SNPs' = "#999999", 'H3A Tags' = "#E69F00", 'WGS SNPs' = "#56B4E9", 'Add-On Tags' = "#D55E00")
    p2_df$SNP_Type <- factor(p2_df$SNP_Type,levels = c('AFGR SNPs','H3A Tags','WGS SNPs','Add-On Tags'))
    p2 <- ggplot(p2_df) + 
      aes(x=Position,y=MAF,color=SNP_Type) + geom_point(alpha = 0.8) + 
      xlim(c(snps_in_gene_regions_parsed[[index]]$start_region,snps_in_gene_regions_parsed[[index]]$end_region)) +  ggtitle('Candidate Tags') + ylim(0,0.5) + scale_colour_manual(values = cbp1)
    return(ggpubr::ggarrange(p1,p2,nrow = 2))
  }else{
    return(p1)
  }
}
PlotAllPasses <- function(no_pass,first_pass,second_pass,region_id){
  no_pass_plot <- PlotRegion(no_pass,region_id,MAF_Plot = T) 
  first_pass_plot <- PlotRegion(first_pass,region_id,MAF_Plot = T)
  second_pass_plot <- PlotRegion(second_pass,region_id,MAF_Plot = T)
  ggpubr::ggarrange(plotlist = list(no_pass_plot,first_pass_plot,second_pass_plot),ncol = 3)
}

#Get prioritized regions. 
# tb_regions_unique <- GetTBRegion()
# TB_gene_regions_first_pass <- readRDS('/mnt/data1/zmxu/TB_array_design/Tag_SNP_Selection/TB_gene_regions_first_pass.rds')
# 
# data.table::fwrite(data.frame(hgnc_symbol=sapply(TB_gene_regions_first_pass$snps_in_gene_regions_parsed,function(x) x$region_id),
#                               chromosome_name = sapply(TB_gene_regions_first_pass$snps_in_gene_regions_parsed,function(x) unique(x$imp_snps$CHROM)),
#                               start_region = sapply(TB_gene_regions_first_pass$snps_in_gene_regions_parsed,function(x) x$start_region),
#                               end_region = sapply(TB_gene_regions_first_pass$snps_in_gene_regions_parsed,function(x) x$end_region)),'../data/Prioritized_regions.csv')

prioritized_regions <- data.table::fread('../data/Prioritized_regions.csv')
tb_gwas_snps <- data.table::fread('../data/tb_gwas_snps.txt',header = F)$V1

args <- commandArgs(trailingOnly = T)
setting <- args[[1]]  
if(setting == 'baseline'){
  wgs_afgr_imp_path <-  args[[2]] #"../results/Baseline/WGS_AFGR_Imputed/WGS_AFGR_Merged_Info.txt" 
  out_dir <- args[[3]] #"../results/Baseline/Tag_SNP_Selection/" 
  raw_vcf <- args[[4]] #"../data/WGS_Host_Data/WGS_Fellay.hg38.joint.118h-1947437863.genotyped.renamed.consensus.filt.vcf.gz" 
  
  h3a_path <- args[[5]] #'../data/WGS_Host_Data/joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz' 
  wgs_path <- args[[6]] #'../data/WGS_Host_Data/joined.hg19.nodup.nomismap.vcf.gz' 
  wgs_excl <- args[[7]]  #'../data/WGS_Host_Data/mismapped_snps.txt'
  wgs_dup <- args[[8]]  #'../data/WGS_Host_Data/joined.hg19.dupvars.txt'
  n_cores <- args[[9]] #20
  
  snps_in_gene_regions_baseline <- GetTargetRegionCoverage(prioritized_regions,
                                                                    tb_gwas_snps = tb_gwas_snps,
                                                                    afgr_imp_path = wgs_afgr_imp_path,
                                                                    addtl_snp_path = '',
                                                                    h3a_path=h3a_path,
                                                                    wgs_path=wgs_path,
                                                                    wgs_excl=wgs_excl,
                                                                    wgs_dup=wgs_dup,
                                                                    merged = T,
                                                                    n_cores = n_cores)
  region_ids <- sapply(snps_in_gene_regions_baseline,function(x) x$region_id)
  SpreadBaseline <- CalcSpread(snps_in_gene_regions_baseline)
  PercAboveBaseline <- CalcPercAbove(snps_in_gene_regions_baseline)
  spread_vs_perc_above <- ggplot(data.frame(Spread=SpreadBaseline,PercAbove=PercAboveBaseline),aes(x=Spread,y=PercAbove))+
    geom_point() +
    geom_label_repel(data= dplyr::filter(data.frame(Spread=SpreadBaseline,PercAbove=PercAboveBaseline,
                                                    Region=seq(1,length(SpreadBaseline),1)),(Spread < 0.8 | PercAbove < 0.8) & (Spread >0 & PercAbove >0)),aes(label=Region))+
    xlim(0,1.2) + ggtitle('H3A Tags')
  
  saveRDS(list(snps_in_gene_regions_parsed=snps_in_gene_regions_baseline,
               PercAbove=PercAboveBaseline,Spread=SpreadBaseline),glue::glue("{out_dir}TB_gene_regions_baseline.rds"))
  
  first_pass_plot <- PlotRegion(snps_in_gene_regions_baseline,1024,MAF_Plot=F)
  hla_regions <- sapply(snps_in_gene_regions_baseline,function(x) grepl(x=x$region_id,pattern = 'HLA'))
  kir_regions <- sapply(snps_in_gene_regions_baseline,function(x) grepl(x=x$region_id,pattern = 'KIR'))

  all_cand_snps <- do.call(rbind,lapply(snps_in_gene_regions_baseline[!hla_regions & !kir_regions],function(x) dplyr::filter(x$wgs_snps,MAF>0.05)))
  data.table::fwrite(dplyr::select(all_cand_snps,Locus_Name=ID),glue::glue("{out_dir}first_pass_rsid.csv"))

  kir_hla_snps <- do.call(rbind,lapply(snps_in_gene_regions_baseline[hla_regions | kir_regions],function(x) dplyr::filter(x$wgs_snps,MAF>0.01)))
  data.table::fwrite(dplyr::select(kir_hla_snps,Locus_Name=ID),glue::glue("{out_dir}kir_hla_rsid.csv"))


  data.table::fwrite(data.frame(Locus_Name=unique(tb_gwas_snps)),glue::glue("{out_dir}tb_snps.csv"))
  mt_snps <- system(paste0("~/Software/bcftools view -r chrM ",raw_vcf," | ~ query -f '%ID %AF\n' | awk '{if($2 > 0.05 && $2 < 0.95) {print $1}}' | grep -v '\\.' "),intern = T)
  mt_snps_phylogeny <- data.table::fread('../data/MT_Y_Haplogroups/Add_On_MT_Major_rsid.txt',header = F)$V1
  data.table::fwrite(data.frame(Locus_Name=unique(c(mt_snps,mt_snps_phylogeny))),glue::glue("{out_dir}mt_snps.csv"))
  y_phylogeny <- data.table::fread('../data/MT_Y_Haplogroups/Add_On_Y_Major_rsid.txt',header = F)$V1
  data.table::fwrite(data.frame(Locus_Name=unique(y_phylogeny)),glue::glue("{out_dir}y_snps.csv"))
}

if(setting == 'Setting1'){
  wgs_afgr_imp_path <-  args[[2]] #"../results/WGS_Addtl_Tags/Setting1/WGS_AFGR_Imputed/WGS_AFGR_Merged_Info.txt" 
  addlt_tags <-  args[[3]] #"../results/WGS_Addtl_Tags/addlt.tags.setting1.random.vcf.gz" 
  out_dir <- args[[4]] #"../results/WGS_Addtl_Tags/Setting1/" 
  raw_vcf <- args[[5]] #"../data/WGS_Host_Data/WGS_Fellay.hg38.joint.118h-1947437863.genotyped.renamed.consensus.filt.vcf.gz" 
  
  h3a_path <- args[[6]] #'../data/WGS_Host_Data/joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz'
  wgs_path <- args[[7]] # '../data/WGS_Host_Data/joined.hg19.nodup.nomismap.vcf.gz'
  wgs_excl <- args[[8]]  #'../data/WGS_Host_Data/mismapped_snps.txt' 
  wgs_dup <- args[[9]] #'../data/WGS_Host_Data/joined.hg19.dupvars.txt'
  n_cores <- args[[10]]
    
    snps_in_gene_regions_setting1 <- GetTargetRegionCoverage(prioritized_regions,
                                                             tb_gwas_snps = tb_gwas_snps,
                                                             afgr_imp_path = wgs_afgr_imp_path,
                                                             addtl_snp_path=addlt_tags,
                                                             h3a_path=h3a_path,
                                                             wgs_path=wgs_path,
                                                             wgs_excl=wgs_excl,
                                                             wgs_dup=wgs_dup,
                                                             merged = T)
  region_ids <- sapply(snps_in_gene_regions_setting1,function(x) x$region_id)
  SpreadBaseline <- CalcSpread(snps_in_gene_regions_setting1)
  PercAboveBaseline <- CalcPercAbove(snps_in_gene_regions_setting1)
  spread_vs_perc_above <- ggplot(data.frame(Spread=SpreadBaseline,PercAbove=PercAboveBaseline),aes(x=Spread,y=PercAbove))+
    geom_point() +
    geom_label_repel(data= dplyr::filter(data.frame(Spread=SpreadBaseline,PercAbove=PercAboveBaseline,
                                                    Region=seq(1,length(SpreadBaseline),1)),(Spread < 0.8 | PercAbove < 0.8) & (Spread >0 & PercAbove >0)),aes(label=Region))+
    xlim(0,1.2) + ggtitle('H3A Tags')
  
  saveRDS(list(snps_in_gene_regions_parsed=snps_in_gene_regions_setting1,
               PercAbove=PercAboveBaseline,Spread=SpreadBaseline),glue::glue("{out_dir}TB_gene_regions_Setting1.rds"))
  
}

if(setting == 'Setting1Random'){
  wgs_afgr_imp_path <-  args[[2]] #"../results/WGS_Addtl_Tags/Setting1Random/WGS_AFGR_Imputed/WGS_AFGR_Merged_Info.txt" #args[[2]]
  addlt_tags <- args[[3]] #"../results/WGS_Addtl_Tags/Setting1Random/addlt.tags.setting1.random.testing.vcf.gz"  
  out_dir <-  args[[4]] #"../results/WGS_Addtl_Tags/Setting1Random/"  
  raw_vcf <-  args[[5]]  #"../data/WGS_Host_Data/WGS_Fellay.hg38.joint.118h-1947437863.genotyped.renamed.consensus.filt.vcf.gz"  
  h3a_path <- args[[6]] #'../data/WGS_Host_Data/joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz'  
  wgs_path <- args[[7]] #'../data/WGS_Host_Data/joined.hg19.nodup.nomismap.vcf.gz'  
  wgs_excl <- args[[8]] #'../data/WGS_Host_Data/mismapped_snps.txt'  
  wgs_dup <- args[[9]] #'../data/WGS_Host_Data/joined.hg19.dupvars.txt' 
  n_cores <- args[[10]]  #20 
  
  snps_in_gene_regions_setting1_random <- GetTargetRegionCoverage(prioritized_regions,
                                                           tb_gwas_snps = tb_gwas_snps,
                                                           afgr_imp_path = wgs_afgr_imp_path,
                                                           addtl_snp_path=addlt_tags,
                                                           h3a_path=h3a_path,
                                                           wgs_path=wgs_path,
                                                           wgs_excl=wgs_excl,
                                                           wgs_dup=wgs_dup,
                                                           merged = T)
  region_ids <- sapply(snps_in_gene_regions_setting1_random,function(x) x$region_id)
  SpreadBaseline <- CalcSpread(snps_in_gene_regions_setting1_random)
  PercAboveBaseline <- CalcPercAbove(snps_in_gene_regions_setting1_random)
  spread_vs_perc_above <- ggplot(data.frame(Spread=SpreadBaseline,PercAbove=PercAboveBaseline),aes(x=Spread,y=PercAbove))+
    geom_point() +
    geom_label_repel(data= dplyr::filter(data.frame(Spread=SpreadBaseline,PercAbove=PercAboveBaseline,
                                                    Region=seq(1,length(SpreadBaseline),1)),(Spread < 0.8 | PercAbove < 0.8) & (Spread >0 & PercAbove >0)),aes(label=Region))+
    xlim(0,1.2) + ggtitle('H3A Tags')
  
  saveRDS(list(snps_in_gene_regions_parsed=snps_in_gene_regions_setting1_random,
               PercAbove=PercAboveBaseline,Spread=SpreadBaseline),glue::glue("{out_dir}TB_gene_regions_Setting1Random.rds"))
  
}

if(setting == 'Setting1Tagger'){
  wgs_afgr_imp_path <-  args[[2]] 
  addlt_tags <-  args[[3]] 
  out_dir <- args[[4]]
  raw_vcf <- args[[5]]
  
  h3a_path <- args[[6]]
  wgs_path <- args[[7]]
  wgs_excl <- args[[8]] 
  wgs_dup <- args[[9]] 
  n_cores <- args[[10]]
  
  snps_in_gene_regions_setting1 <- GetTargetRegionCoverage(prioritized_regions,
                                                                  tb_gwas_snps = tb_gwas_snps,
                                                                  afgr_imp_path = wgs_afgr_imp_path,
                                                                  addtl_snp_path=addlt_tags,
                                                                  h3a_path=h3a_path,
                                                                  wgs_path=wgs_path,
                                                                  wgs_excl=wgs_excl,
                                                                  wgs_dup=wgs_dup,
                                                                  merged = T)
  region_ids <- sapply(snps_in_gene_regions_setting1,function(x) x$region_id)
  SpreadBaseline <- CalcSpread(snps_in_gene_regions_setting1)
  PercAboveBaseline <- CalcPercAbove(snps_in_gene_regions_setting1)
  spread_vs_perc_above <- ggplot(data.frame(Spread=SpreadBaseline,PercAbove=PercAboveBaseline),aes(x=Spread,y=PercAbove))+
    geom_point() +
    geom_label_repel(data= dplyr::filter(data.frame(Spread=SpreadBaseline,PercAbove=PercAboveBaseline,
                                                    Region=seq(1,length(SpreadBaseline),1)),(Spread < 0.8 | PercAbove < 0.8) & (Spread >0 & PercAbove >0)),aes(label=Region))+
    xlim(0,1.2) + ggtitle('H3A Tags')
  
  saveRDS(list(snps_in_gene_regions_parsed=snps_in_gene_regions_setting1,
               PercAbove=PercAboveBaseline,Spread=SpreadBaseline),glue::glue("{out_dir}TB_gene_regions_Setting1_Tagger.rds"))
  
}


if(setting == 'Setting2Baseline'){

  wgs_afgr_imp_path <-  args[[2]] #"../results/WGS_Addtl_Tags/Setting2Random/WGS_AFGR_Imputed/WGS_AFGR_Merged_Info.txt" 
  addlt_tags <-  args[[3]] # "../results/WGS_Addtl_Tags/addlt.tags.setting2.random.vcf.gz" 
  out_dir <- args[[4]] # "../results/WGS_Addtl_Tags/Setting2Random/" 
  raw_vcf <- args[[5]] # "../data/WGS_Host_Data/WGS_Fellay.hg38.joint.118h-1947437863.genotyped.renamed.consensus.filt.vcf.gz" 
  
  h3a_path <- args[[6]] # '../data/WGS_Host_Data/joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz' 
  wgs_path <- args[[7]] # '../dat/WGS_Host_Data/joined.hg19.nodup.nomismap.vcf.gz' 
  wgs_excl <- args[[8]] # '../data/WGS_Host_Data/mismapped_snps.txt' 
  wgs_dup <- args[[9]] # '../data/WGS_Host_Data/joined.hg19.dupvars.txt'
  n_cores <- args[[10]] #20
  
  setting2_regions <- data.frame(start_region = sapply(snps_in_region,function(x) x$start_region),
                                 end_region = sapply(snps_in_region,function(x) x$end_region),
                                 chromosome_name = sapply(snps_in_region,function(x) unique(x$wgs_snps$CHROM)))
  
  snps_in_gene_regions_setting2_baseline <- GetTargetRegionCoverage(setting2_regions,
                                                                  tb_gwas_snps = tb_gwas_snps,
                                                                  afgr_imp_path = wgs_afgr_imp_path,
                                                                  addtl_snp_path=addlt_tags,
                                                                  h3a_path=h3a_path,
                                                                  wgs_path=wgs_path,
                                                                  wgs_excl=wgs_excl,
                                                                  wgs_dup=wgs_dup,
                                                                  merged = T)
  region_ids <- sapply(snps_in_gene_regions_setting2_baseline,function(x) x$region_id)
  SpreadBaseline <- CalcSpread(snps_in_gene_regions_setting2_baseline)
  PercAboveBaseline <- CalcPercAbove(snps_in_gene_regions_setting2_baseline)
  spread_vs_perc_above <- ggplot(data.frame(Spread=SpreadBaseline,PercAbove=PercAboveBaseline),aes(x=Spread,y=PercAbove))+
    geom_point() +
    geom_label_repel(data= dplyr::filter(data.frame(Spread=SpreadBaseline,PercAbove=PercAboveBaseline,
                                                    Region=seq(1,length(SpreadBaseline),1)),(Spread < 0.8 | PercAbove < 0.8) & (Spread >0 & PercAbove >0)),aes(label=Region))+
    xlim(0,1.2) + ggtitle('H3A Tags')
  
  saveRDS(list(snps_in_gene_regions_parsed=snps_in_gene_regions_setting2_baseline,
               PercAbove=PercAboveBaseline,Spread=SpreadBaseline),glue::glue("{out_dir}TB_gene_regions_Setting2_baseline.rds"))
  
}


if(setting == 'Setting2Random'){

  wgs_afgr_imp_path <-  args[[2]] #"../results/WGS_Addtl_Tags/Setting2Random/WGS_AFGR_Imputed/WGS_AFGR_Merged_Info.txt" 
  addlt_tags <-  args[[3]] # "../results/WGS_Addtl_Tags/addlt.tags.setting2.random.vcf.gz" 
  out_dir <- args[[4]] # "../results/WGS_Addtl_Tags/Setting2Random/" 
  raw_vcf <- args[[5]] # "../data/WGS_Host_Data/WGS_Fellay.hg38.joint.118h-1947437863.genotyped.renamed.consensus.filt.vcf.gz" 
  
  h3a_path <- args[[6]] # '../data/WGS_Host_Data/joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz' 
  wgs_path <- args[[7]] # '../dat/WGS_Host_Data/joined.hg19.nodup.nomismap.vcf.gz' 
  wgs_excl <- args[[8]] # '../data/WGS_Host_Data/mismapped_snps.txt' 
  wgs_dup <- args[[9]] # '../data/WGS_Host_Data/joined.hg19.dupvars.txt'
  n_cores <- args[[10]] #20
  setting2_regions <- data.frame(start_region = sapply(snps_in_region,function(x) x$start_region),
                                 end_region = sapply(snps_in_region,function(x) x$end_region),
                                 chromosome_name = sapply(snps_in_region,function(x) unique(x$wgs_snps$CHROM)))
  
  
  
  snps_in_gene_regions_setting2_random <- GetTargetRegionCoverage(setting2_regions,
                                                           tb_gwas_snps = tb_gwas_snps,
                                                           afgr_imp_path = wgs_afgr_imp_path,
                                                           addtl_snp_path=addlt_tags,
                                                           h3a_path=h3a_path,
                                                           wgs_path=wgs_path,
                                                           wgs_excl=wgs_excl,
                                                           wgs_dup=wgs_dup,
                                                           merged = T)
  region_ids <- sapply(snps_in_gene_regions_setting2_random,function(x) x$region_id)
  SpreadBaseline <- CalcSpread(snps_in_gene_regions_setting2_random)
  PercAboveBaseline <- CalcPercAbove(snps_in_gene_regions_setting2_random)
  spread_vs_perc_above <- ggplot(data.frame(Spread=SpreadBaseline,PercAbove=PercAboveBaseline),aes(x=Spread,y=PercAbove))+
    geom_point() +
    geom_label_repel(data= dplyr::filter(data.frame(Spread=SpreadBaseline,PercAbove=PercAboveBaseline,
                                                    Region=seq(1,length(SpreadBaseline),1)),(Spread < 0.8 | PercAbove < 0.8) & (Spread >0 & PercAbove >0)),aes(label=Region))+
    xlim(0,1.2) + ggtitle('H3A Tags')
  
  saveRDS(list(snps_in_gene_regions_parsed=snps_in_gene_regions_setting2_random,
               PercAbove=PercAboveBaseline,Spread=SpreadBaseline),glue::glue("{out_dir}TB_gene_regions_Setting2Random.rds"))
  
}

if(setting == 'Setting2'){

  wgs_afgr_imp_path <-  args[[2]] #"../results/WGS_Addtl_Tags/Setting2/WGS_AFGR_Imputed/WGS_AFGR_Merged_Info.txt" 
  addlt_tags <-  args[[3]] # "../results/WGS_Addtl_Tags/addlt.tags.setting2.random.vcf.gz" 
  out_dir <- args[[4]] # "../results/WGS_Addtl_Tags/Setting2/" 
  raw_vcf <- args[[5]] # "../data/WGS_Host_Data/WGS_Fellay.hg38.joint.118h-1947437863.genotyped.renamed.consensus.filt.vcf.gz" 
  
  h3a_path <- args[[6]] # '../data/WGS_Host_Data/joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz' 
  wgs_path <- args[[7]] # '../dat/WGS_Host_Data/joined.hg19.nodup.nomismap.vcf.gz' 
  wgs_excl <- args[[8]] # '../data/WGS_Host_Data/mismapped_snps.txt' 
  wgs_dup <- args[[9]] # '../data/WGS_Host_Data/joined.hg19.dupvars.txt'
  n_cores <- args[[10]] #20
  setting2_regions <- data.frame(start_region = sapply(snps_in_region,function(x) x$start_region),
                             end_region = sapply(snps_in_region,function(x) x$end_region),
                             chromosome_name = sapply(snps_in_region,function(x) unique(x$wgs_snps$CHROM)))

  
  
  snps_in_gene_regions_setting2 <- GetTargetRegionCoverage(setting2_regions,
                                                           tb_gwas_snps = tb_gwas_snps,
                                                           afgr_imp_path = wgs_afgr_imp_path,
                                                           addtl_snp_path=addlt_tags,
                                                           h3a_path=h3a_path,
                                                           wgs_path=wgs_path,
                                                           wgs_excl=wgs_excl,
                                                           wgs_dup=wgs_dup,
                                                           merged = T)
  region_ids <- sapply(snps_in_gene_regions_setting2,function(x) x$region_id)
  SpreadBaseline <- CalcSpread(snps_in_gene_regions_setting2)
  PercAboveBaseline <- CalcPercAbove(snps_in_gene_regions_setting2)
  spread_vs_perc_above <- ggplot(data.frame(Spread=SpreadBaseline,PercAbove=PercAboveBaseline),aes(x=Spread,y=PercAbove))+
    geom_point() +
    geom_label_repel(data= dplyr::filter(data.frame(Spread=SpreadBaseline,PercAbove=PercAboveBaseline,
                                                    Region=seq(1,length(SpreadBaseline),1)),(Spread < 0.8 | PercAbove < 0.8) & (Spread >0 & PercAbove >0)),aes(label=Region))+
    xlim(0,1.2) + ggtitle('H3A Tags')
  
  saveRDS(list(snps_in_gene_regions_parsed=snps_in_gene_regions_setting2,
               PercAbove=PercAboveBaseline,Spread=SpreadBaseline),glue::glue("{out_dir}TB_gene_regions_Setting2.rds"))
  
}
if(setting == 'Setting2TaggerBaseline'){
  
  wgs_afgr_imp_path <-  args[[2]] #"../results/WGS_Addtl_Tags/Setting2Random/WGS_AFGR_Imputed/WGS_AFGR_Merged_Info.txt" 
  addlt_tags <-  args[[3]] # "../results/WGS_Addtl_Tags/addlt.tags.setting2.random.vcf.gz" 
  out_dir <- args[[4]] # "../results/WGS_Addtl_Tags/Setting2Random/" 
  raw_vcf <- args[[5]] # "../data/WGS_Host_Data/WGS_Fellay.hg38.joint.118h-1947437863.genotyped.renamed.consensus.filt.vcf.gz" 
  
  h3a_path <- args[[6]] # '../data/WGS_Host_Data/joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz' 
  wgs_path <- args[[7]] # '../dat/WGS_Host_Data/joined.hg19.nodup.nomismap.vcf.gz' 
  wgs_excl <- args[[8]] # '../data/WGS_Host_Data/mismapped_snps.txt' 
  wgs_dup <- args[[9]] # '../data/WGS_Host_Data/joined.hg19.dupvars.txt'
  n_cores <- args[[10]] #20
  
  setting2_regions <- data.frame(start_region = sapply(snps_in_region,function(x) x$start_region),
                                 end_region = sapply(snps_in_region,function(x) x$end_region),
                                 chromosome_name = sapply(snps_in_region,function(x) unique(x$wgs_snps$CHROM)))
  
  snps_in_gene_regions_setting2_baseline <- GetTargetRegionCoverage(setting2_regions,
                                                                    tb_gwas_snps = tb_gwas_snps,
                                                                    afgr_imp_path = wgs_afgr_imp_path,
                                                                    addtl_snp_path=addlt_tags,
                                                                    h3a_path=h3a_path,
                                                                    wgs_path=wgs_path,
                                                                    wgs_excl=wgs_excl,
                                                                    wgs_dup=wgs_dup,
                                                                    merged = T)
  region_ids <- sapply(snps_in_gene_regions_setting2_baseline,function(x) x$region_id)
  SpreadBaseline <- CalcSpread(snps_in_gene_regions_setting2_baseline)
  PercAboveBaseline <- CalcPercAbove(snps_in_gene_regions_setting2_baseline)
  spread_vs_perc_above <- ggplot(data.frame(Spread=SpreadBaseline,PercAbove=PercAboveBaseline),aes(x=Spread,y=PercAbove))+
    geom_point() +
    geom_label_repel(data= dplyr::filter(data.frame(Spread=SpreadBaseline,PercAbove=PercAboveBaseline,
                                                    Region=seq(1,length(SpreadBaseline),1)),(Spread < 0.8 | PercAbove < 0.8) & (Spread >0 & PercAbove >0)),aes(label=Region))+
    xlim(0,1.2) + ggtitle('H3A Tags')
  
  saveRDS(list(snps_in_gene_regions_parsed=snps_in_gene_regions_setting2_baseline,
               PercAbove=PercAboveBaseline,Spread=SpreadBaseline),glue::glue("{out_dir}TB_gene_regions_Setting2_Tagger_baseline.rds"))
}

if(setting == 'Setting2Tagger'){
  wgs_afgr_imp_path <-  args[[2]] #"../results/WGS_Addtl_Tags/Setting2/WGS_AFGR_Imputed/WGS_AFGR_Merged_Info.txt" 
  addlt_tags <-  args[[3]] # "../results/WGS_Addtl_Tags/addlt.tags.setting2.random.vcf.gz" 
  out_dir <- args[[4]] # "../results/WGS_Addtl_Tags/Setting2/" 
  raw_vcf <- args[[5]] # "../data/WGS_Host_Data/WGS_Fellay.hg38.joint.118h-1947437863.genotyped.renamed.consensus.filt.vcf.gz" 
  
  h3a_path <- args[[6]] # '../data/WGS_Host_Data/joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz' 
  wgs_path <- args[[7]] # '../dat/WGS_Host_Data/joined.hg19.nodup.nomismap.vcf.gz' 
  wgs_excl <- args[[8]] # '../data/WGS_Host_Data/mismapped_snps.txt' 
  wgs_dup <- args[[9]] # '../data/WGS_Host_Data/joined.hg19.dupvars.txt'
  n_cores <- args[[10]] #20
  setting2_regions <- data.frame(start_region = sapply(snps_in_region,function(x) x$start_region),
                                 end_region = sapply(snps_in_region,function(x) x$end_region),
                                 chromosome_name = sapply(snps_in_region,function(x) unique(x$wgs_snps$CHROM)))
  
  
  snps_in_gene_regions_setting2 <- GetTargetRegionCoverage(setting2_regions,
                                                           tb_gwas_snps = tb_gwas_snps,
                                                           afgr_imp_path = wgs_afgr_imp_path,
                                                           addtl_snp_path=addlt_tags,
                                                           h3a_path=h3a_path,
                                                           wgs_path=wgs_path,
                                                           wgs_excl=wgs_excl,
                                                           wgs_dup=wgs_dup,
                                                           merged = T)
  region_ids <- sapply(snps_in_gene_regions_setting2,function(x) x$region_id)
  SpreadBaseline <- CalcSpread(snps_in_gene_regions_setting2)
  PercAboveBaseline <- CalcPercAbove(snps_in_gene_regions_setting2)
  spread_vs_perc_above <- ggplot(data.frame(Spread=SpreadBaseline,PercAbove=PercAboveBaseline),aes(x=Spread,y=PercAbove))+
    geom_point() +
    geom_label_repel(data= dplyr::filter(data.frame(Spread=SpreadBaseline,PercAbove=PercAboveBaseline,
                                                    Region=seq(1,length(SpreadBaseline),1)),(Spread < 0.8 | PercAbove < 0.8) & (Spread >0 & PercAbove >0)),aes(label=Region))+
    xlim(0,1.2) + ggtitle('H3A Tags')
  
  saveRDS(list(snps_in_gene_regions_parsed=snps_in_gene_regions_setting2,
               PercAbove=PercAboveBaseline,Spread=SpreadBaseline),glue::glue("{out_dir}TB_gene_regions_Setting2_Tagger.rds"))
  
}
