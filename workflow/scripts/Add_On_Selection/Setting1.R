#Select Add-on SNPs in Prioritized Regions
library(infotheo)
library(filematrix)
library(dplyr)
library(snpStats)
library(pbmcapply)
library(future.apply)
library(glue)
set.seed(888)

#Calculates pariwise MI between two matrices
MI_Calculation <- function(source,target,source_names,target_names,method = 'MI',n_cores=40){
  plan(multiprocess, workers = n_cores)
  if(method == 'MI'){
    MI_mat <- future_apply(source,2,function(s) apply(target,2,function(t) mutinformation(s,t)))
    #MI_mat <- apply(source,2,function(s) apply(target,2,function(t) mutinformation(s,t)))
    
  }else if(method == 'cor'){
    MI_mat <- future_apply(source,2,function(s) apply(target,2,function(t) cor(s,t,use = 'complete.obs')))
    #MI_mat <- apply(source,2,function(s) apply(target,2,function(t) cor(s,t,use = 'complete.obs')))
  }
  MI_mat <- t(MI_mat)
  if(nrow(MI_mat) != length(source_names) | ncol(MI_mat) != length(target_names)){
    MI_mat <- t(MI_mat)
    rownames(MI_mat) <- source_names
    colnames(MI_mat) <- target_names
    return(MI_mat)
  }else{
    rownames(MI_mat) <- source_names
    colnames(MI_mat) <- target_names
    return(MI_mat)
  }
}
#Find snp which maximizes sum of MI, correcting for previous tag snps
FindBestTagToAdd <- function(candidate_snps,candidate_tag_snps,pairwise_candidate_MI,candidate_max_MI,probe_correction=F){
  #For each candidate tag, calculate sum of MI with other candidates, if higher than prev tag's MI. 
  sum_of_MI <- rep(NA,length(candidate_tag_snps$ID))
  n_improv <- rep(NA,length(candidate_tag_snps$ID))
  MI_info <- list()
  names(sum_of_MI) <- candidate_tag_snps$ID
  tagged_snps <- list()
  for(i in 1:length(candidate_tag_snps$ID)){
    curSNP <- candidate_tag_snps$ID[i]
    #Calculate MI of tag SNPs against all candidate SNPs
    MI_against_candidates <- pairwise_candidate_MI[curSNP,candidate_snps]
    #take current MI, if existing tags has worst MI against target of interest
    MI_diff <- MI_against_candidates - candidate_max_MI
    MI_diff <- sapply(MI_diff,function(x) max(0,x))
    sum_of_MI[i] <- sum(MI_diff)
    tagged_snps[[i]] <- MI_diff[which(MI_diff > 0)]
    #Two probes required for A/T and C/G SNPs (InfiniumI)
    if(probe_correction){
      if(candidate_tag_snps$Assay_Type[i] == 'InfiniumI'){
        sum_of_MI[i] <- sum_of_MI[i]/2
      }
    }
    n_improv[i] <- sum(MI_diff > 0)
  }
  if(all(sum_of_MI == 0)){
    maxIndex <- which.max(candidate_tag_snps$Final_Score)
    if(length(maxIndex) > 1){
      maxIndex <- maxIndex[which.max(candidate_tag_snps$MAF[maxIndex])]
    }
  }else{
    illuminaProbeScoresRank <- rank(candidate_tag_snps$Final_Score[which(sum_of_MI != 0)])
    MI_Rank <- rank(sum_of_MI[which(sum_of_MI != 0)])
    avgRank <- rowMeans(cbind(illuminaProbeScoresRank,MI_Rank))
    maxIndex <- which(avgRank == max(avgRank))
    maxIndex <- which(sum_of_MI != 0)[maxIndex]
    if(length(maxIndex) > 1){
      maxIndex <- maxIndex[which.max(candidate_tag_snps$MAF[maxIndex])]
    }
  }
  bestTag <- candidate_tag_snps$ID[maxIndex]
  
  return(list(bestTag = bestTag,bestMIDiff = sum_of_MI[maxIndex],bestProbeScore=candidate_tag_snps$Final_Score[maxIndex],n_improv = n_improv[maxIndex],tagged_snps = tagged_snps[[maxIndex]]))
}
#Find tags to add to a region
SelectTags <- function(candidate_snps,candidate_tag_snps,MI_matrix,pairwise_candidate_MI,candidate_max_MI,verbose=F,improvThresh=-Inf,max_iter=NULL,returnTagInfo = F){
  iter <- 1
  if(is.null(max_iter)){
    max_iter = length(candidate_snps)
  }
  n_improv_thresh = max(1,round(0.005 * length(candidate_snps)))
  prevBestMIDiff <- 0
  tagSearchResult <- list()
  while(iter <= max_iter){
    #Invoke Search for best tag
    tagSearchResult[[iter]] <- FindBestTagToAdd(candidate_snps,candidate_tag_snps,pairwise_candidate_MI,candidate_max_MI,probe_correction=T)
    bestTag <- tagSearchResult[[iter]]$bestTag
    bestMIDiff <- tagSearchResult[[iter]]$bestMIDiff
    if(verbose){
      print(1-(bestMIDiff/prevBestMIDiff))
    }
    if((tagSearchResult[[iter]]$n_improv < n_improv_thresh) | ((iter > 1) & (1-(bestMIDiff/prevBestMIDiff) < improvThresh))){
      if(returnTagInfo){
        return(tagSearchResult)
      }else{
        return(MI_matrix)
      }
    }
    #Add a new column with the best tag
    bestTagAgainstOthersMI <- pairwise_candidate_MI[bestTag,rownames(MI_matrix)]
    bestTagAgainstOthersMI <- as.matrix(bestTagAgainstOthersMI)
    colnames(bestTagAgainstOthersMI) <- bestTag
    MI_matrix <- cbind(MI_matrix,bestTagAgainstOthersMI)
    #Find best tag so far for each candidate SNP
    candidate_max_MI <- apply(MI_matrix,1,function(x) max(x))
    iter = iter + 1
    prevBestMIDiff <- bestMIDiff
  }
  if(returnTagInfo){
    return(tagSearchResult)
  }else{
    return(MI_matrix)
  }
}
ExtractDosage <- function(dosage_matrix,row_indices,col_indicies,col_names,row_names){
  sub_matrix <- dosage_matrix[row_indices,col_indicies]
  rownames(sub_matrix)<- row_names
  colnames(sub_matrix) <- col_names
  return(sub_matrix)
}
#Get Candidate Target SNPs in a region. 
GetCandidateSNPs <- function(candidate_sites = 'afgr',region_id,training_samples,snps_in_gene_regions_parsed,wgs_mat,wgs_colnames,wgs_rownames,illumina_scores,INFO_thresh = 0.8,excl_snps=c(),MAF_thresh = 0.05){
  #Get snps within gene region (Imputed,WGS,tags)
  gene_regions_obj <- snps_in_gene_regions_parsed$snps_in_gene_regions_parsed
  if(is.null(gene_regions_obj)){
    gene_regions_obj <- snps_in_gene_regions_parsed
  }
  cur_region <- gene_regions_obj[[region_id]]
  #Add TB SNPs as tags SNPs, if they are not already tag SNPs
  tb_snps <- c()
  if(length(cur_region$cur_region_tb_gwas) != 0){
    tb_snps <- dplyr::filter(cur_region$wgs_snps,ID %in% cur_region$cur_region_tb_gwas) %>% dplyr::select(-AF,-MAF)
    tb_snps$CHROM <- sapply(tb_snps$CHROM,function(x) gsub(x=x,pattern = 'chr',replacement = ''))
    cur_region$addtl_snps <- rbind(cur_region$addtl_snps,tb_snps) %>% dplyr::filter(!is.na(POS)) %>% dplyr::distinct(.keep_all = T)
    tb_snps <- tb_snps$ID
  }

  #print(region_id)
  #Get Candidate snps which are subject to improvement(INFO < 0.8 or not imputed)
  if(candidate_sites == 'afgr' | candidate_sites == 'afgr_imp'){
    if(all(is.na(cur_region$afgr_sites))){
      afgr_snps <- cur_region$afgr_sites %>% dplyr::select(CHROM=CHROM,POS,ID=ID,REF=REF,ALT=ALT,MAF=MAF)
      candidate_snps <- afgr_snps %>% dplyr::select(CHROM=CHROM,POS=POS,ID=ID,REF=REF,ALT=ALT,MAF=MAF)
    }else{
      if(candidate_sites == 'afgr'){
        afgr_snps <- dplyr::left_join(cur_region$afgr_sites,cur_region$imp_snps,by=c('POS'='POS','REF'='REF','ALT'='ALT')) %>%
          dplyr::filter(INFO < INFO_thresh | is.na(INFO)) %>% dplyr::select(CHROM=CHROM.x,POS,ID=ID.y,REF=REF,ALT=ALT,INFO,MAF=MAF.x)
        candidate_snps <- dplyr::inner_join(afgr_snps,cur_region$wgs_snps,by=c('POS'='POS','REF'='REF','ALT'='ALT')) %>% 
          dplyr::select(CHROM=CHROM.y,POS=POS,ID=ID.y,REF=REF,ALT=ALT,MAF=MAF.x)
      }else if(candidate_sites =='afgr_imp'){
        afgr_snps <- cur_region$imp_snps %>% dplyr::filter(INFO < INFO_thresh) 
        candidate_snps <- dplyr::inner_join(afgr_snps,cur_region$wgs_snps,by=c('POS'='POS','REF'='REF','ALT'='ALT')) %>% 
          dplyr::select(CHROM=CHROM.y,POS=POS,ID=ID.y,REF=REF,ALT=ALT,MAF=MAF.x)
      }
      candidate_tag_snps <- dplyr::left_join(cur_region$afgr_sites,cur_region$imp_snps,by=c('POS'='POS','REF'='REF','ALT'='ALT')) %>%
        dplyr::select(CHROM=CHROM.x,POS,ID=ID.y,REF=REF,ALT=ALT,INFO,MAF=MAF.x) %>% dplyr::inner_join(cur_region$wgs_snps,by=c('POS'='POS','REF'='REF','ALT'='ALT')) %>% 
        dplyr::select(CHROM=CHROM.y,POS=POS,ID=ID.y,REF=REF,ALT=ALT,MAF=MAF.x)
    }
  }else if(candidate_sites == 'wgs'){
    #Candidate snps are: Info Score < Thresh or Not imputed SNP with MAF > MAF_thresh in either AFGR or sequenced and not in excl snps
    candidate_snps <- dplyr::left_join(cur_region$wgs_snps,cur_region$imp_snps,by=c('POS'='POS','REF'='REF','ALT'='ALT')) %>% 
      dplyr::filter((INFO < INFO_thresh & (MAF.y > MAF_thresh | MAF.x > MAF_thresh)) | (is.na(INFO) & MAF.x > MAF_thresh)) %>%
      dplyr::select(CHROM=CHROM.x,POS,ID=ID.x,REF=REF,ALT=ALT,INFO,MAF = MAF.x,MAF.AFGR = MAF.y) %>% dplyr::filter(!ID %in% excl_snps)
    
    #Candidate tags are: any sequences snps with MAF > MAF_thresh in sequeced data and not in excl snps
    candidate_tag_snps <- dplyr::left_join(cur_region$wgs_snps,cur_region$imp_snps,by=c('POS'='POS','REF'='REF','ALT'='ALT')) %>%  
    dplyr::filter(MAF.x > MAF_thresh) %>%
    dplyr::select(CHROM=CHROM.x,POS,ID=ID.x,REF=REF,ALT=ALT,INFO,MAF = MAF.x)
    candidate_tag_snps <- dplyr::left_join(candidate_tag_snps,illumina_scores,by=c('ID'='Locus_Name')) %>% dplyr::filter(!is.na(Assay_Type)) %>% dplyr::filter(!ID %in% excl_snps) %>% dplyr::filter(Final_Score > 0.3)
  }
  #H3A Tags and Addtl Tags (if any) as tags which already exist
  existing_tag_snps <- rbind(cur_region$h3a_tags,cur_region$addtl_snps) %>% dplyr::filter(!is.na(POS)) %>% dplyr::distinct(.keep_all = T)
  #Remove existing h3a snps from candidate snps
  candidate_tag_snps <- dplyr::filter(candidate_tag_snps,!POS%in%existing_tag_snps$POS)

  #Load genotype for existing tag snps, candidate tag snps, and candidate snps
  trainingIndices <- match(training_samples,wgs_rownames)
  if(all(is.na(existing_tag_snps$ID))){
    existing_geno <- as.character(NA)
  }else{
    existing_geno <- ExtractDosage(dosage_matrix = wgs_mat,
                                   row_indices = trainingIndices,
                                   col_indicies = match(existing_tag_snps$ID,wgs_colnames),
                                   row_names = training_samples,
                                   col_names = existing_tag_snps$ID)
    #Filter for missingess in genotype
    existing_geno_missingness <- apply(existing_geno,2,function(x) sum(is.na(x)) / length(x))
    existing_geno <- existing_geno[,which(existing_geno_missingness < 0.5),drop = F]
    existing_tag_snps <- existing_tag_snps[which(existing_geno_missingness < 0.5),]
  }
  if(all(is.na(candidate_snps$ID))){
    candidate_geno <- as.character(NA)
  }else{
    candidate_geno <- ExtractDosage(dosage_matrix = wgs_mat,
                                    row_indices = trainingIndices,
                                    col_indicies = match(candidate_snps$ID,wgs_colnames),
                                    row_names = training_samples,
                                    col_names = candidate_snps$ID)
    #Filter for missingess in genotype
    candidate_geno_missingness <- apply(candidate_geno,2,function(x) sum(is.na(x)) / length(x))
    candidate_geno <- candidate_geno[,which(candidate_geno_missingness < 0.5),drop = F]
    candidate_snps <- candidate_snps[which(candidate_geno_missingness < 0.5),]
    
  }
  if(all(is.na(candidate_tag_snps$ID))){
    candidate_tags_geno <- as.character(NA)
  }else{
    candidate_tags_geno <- ExtractDosage(dosage_matrix = wgs_mat,
                                         row_indices = trainingIndices,
                                         col_indicies = match(candidate_tag_snps$ID,wgs_colnames),
                                         row_names = training_samples,
                                         col_names = candidate_tag_snps$ID)
    #Filter for missingess in genotype
    candidate_tags_geno_missingness <- apply(candidate_tags_geno,2,function(x) sum(is.na(x)) / length(x))
    candidate_tags_geno <- candidate_tags_geno[,which(candidate_tags_geno_missingness < 0.5),drop = F]
    candidate_tag_snps <- candidate_tag_snps[which(candidate_tags_geno_missingness < 0.5),]
  }
  return(list(candidate_geno=candidate_geno,existing_geno=existing_geno,
              candidate_tags_geno=candidate_tags_geno,candidate_tag_snps=candidate_tag_snps,
              candidate_snps=candidate_snps,existing_tag_snps=existing_tag_snps,
              tb_snps = tb_snps,region_id=cur_region$region_id,focused = cur_region$focused))
}
WriteAddOnFile <- function(lift_over_vcf,phased_vcf,out_prefix,out_path,random){
    addtl_tags_dir <- paste0(out_path)
    system(paste0('mkdir -p ',addtl_tags_dir))
  
    #Write unphased vcf for AFGR Imputation
    system(paste0('~/Software/plink2 --vcf ',gsub(lift_over_vcf,pattern = '.gz',replacement = ''),' --extract ',out_path,'h3a.with.',out_prefix,'.txt',' --export vcf bgz --out ',
                  paste0(addtl_tags_dir,out_prefix)))
    system(paste0('~/Software/bcftools index -t --threads 5 ',addtl_tags_dir,out_prefix,'.vcf.gz'))
    testing_samples <- readRDS(paste0('../data/WGS_Host_Data/testing_set/testing_split.rds'))
    system(paste0('~/Software/bcftools view -O v -s \"',
                  paste(testing_samples,collapse = ','),'\" ',paste0(addtl_tags_dir,out_prefix,'.vcf.gz'),' > ',
                  paste0(addtl_tags_dir,out_prefix,'.testing.vcf')))
    system(paste0('bgzip -c ',paste0(addtl_tags_dir,out_prefix,'.testing.vcf'),'>',paste0(addtl_tags_dir,out_prefix,'.testing.vcf.gz')))
    system(paste0('~/Software/bcftools index -t --threads 5 ',paste0(addtl_tags_dir,out_prefix,'.testing.vcf.gz')))
    #Write phased VCF for WGS Imputation
    chr <- c(seq(1,22,1),'X','PAR')
    if(random){
      addtl_tags_dir <- paste0(out_path,'/Setting1Random/')
    }else{
      addtl_tags_dir <- paste0(out_path,'/Setting1/')
    }
    system(paste0('mkdir -p ',addtl_tags_dir))
    #Write files with h3a and addtl tags for each chr (including X)
    for(cur_chr in chr){
      system(paste0('gunzip -c ',phased_vcf,'chr',cur_chr,'.phased.testing.vcf.gz > ',phased_vcf,'chr',cur_chr,'.phased.testing.vcf'))
      system(paste0('~/Software/plink2 --vcf ',phased_vcf,'chr',cur_chr,'.phased.testing.vcf',' --extract ',out_path,
                    '/h3a.with.',out_prefix,'.txt',' --export vcf bgz --out ',
                    paste0(addtl_tags_dir,'chr',cur_chr,'.',out_prefix,'.testing')))
      system(paste0('~/Software/bcftools index -t --threads 5 ',addtl_tags_dir,'chr',cur_chr,'.',out_prefix,'.testing.vcf.gz'))
    }

}
WriteSelectedTags <- function(out_prefix,site_type,snps_in_gene_regions_parsed_path,keep_samples_path,existing_tags,lift_over_vcf,phased_vcf,SpreadThresh,PercAboveThresh,random = F,n_cores = 40){
  #Run Tag Selection
  data_dir_path <- '../results/'
  
  training_samples <- readRDS('../data/WGS_Host_Data/training_set/training_split.rds')
  
  if(!random){
    wgs_mat <- fm.open('../results/Imputation_Eval/wgs_mat')
    wgs_colnames <- data.table::fread('../results/Imputation_Eval/wgs_mat.nmscol.txt',sep = "",header = F)$V1
    wgs_rownames <- rownames(wgs_mat)
    
    snps_in_gene_regions_parsed <- readRDS(snps_in_gene_regions_parsed_path)
    wgs_regions_to_update <- which((snps_in_gene_regions_parsed$Spread < SpreadThresh | snps_in_gene_regions_parsed$PercAbove < PercAboveThresh) |
                                     (sapply(snps_in_gene_regions_parsed$snps_in_gene_regions_parsed,function(x) x$focused)))
    hla_kir_regions <- wgs_regions_to_update[sapply(wgs_regions_to_update,function(i){
      cur_region_id <- snps_in_gene_regions_parsed$snps_in_gene_regions_parsed[[i]]$region_id
      return(grepl(x=cur_region_id,pattern = 'HLA') | grepl(x=cur_region_id,pattern = 'KIR'))})]
    #Should upload to illumina to get scores for candidate tags ('~/G2G_TB/Tag_SNP_Selection/first_pass_rsid.csv')
    #Should upload to illumina to get scores for candidate tags ('~/G2G_TB/Tag_SNP_Selection/kir_hla_rsid.csv')
    
    illumina_scores <- rbind(data.table::fread('../results/Tag_SNP_Selection/WGS.Baseline.Illumina.Scores.csv'),data.table::fread('../results/Tag_SNP_Selection/WGS.KIR.HLA.Illumina.Scores.csv'))
    illumina_scores <- illumina_scores %>% dplyr::select(Locus_Name,Assay_Type,Final_Score)
    
    region_info_raw <- pbmclapply(wgs_regions_to_update,function(i) GetCandidateSNPs(candidate_sites = site_type,
                                                                                     region_id = i,training_samples = training_samples,
                                                                                     snps_in_gene_regions_parsed=snps_in_gene_regions_parsed,
                                                                                     wgs_mat = wgs_mat,
                                                                                     wgs_colnames = wgs_colnames,
                                                                                     wgs_rownames = wgs_rownames,
                                                                                     illumina_scores=illumina_scores,
                                                                                     INFO_thresh = 0.8,
                                                                                     MAF_thresh = ifelse(i %in% hla_kir_regions,0.01,0.05)),mc.cores = 20)
    
    #Only keep regions where there is atleast 1 candidate snp
    region_info <- region_info_raw[sapply(region_info_raw, function(x) nrow(x$candidate_snps)) > 0]
    print(paste0('Regions To Add Tags:',length(region_info)))
    tags_to_add <- list()
    tag_info <- list()
    for(i in 1:length(region_info)){
      print(i)
      #Calculate MI of candidate snps against existing tags
      cur_region <- region_info[[i]]
      if(all(is.na(cur_region$candidate_tag_snps)) | all(is.na(cur_region$candidate_snps$ID))){
        tags_to_add[[i]] <- NA
        next
      }
      #Add TB snps to tags to add (added to existing tag snps previously, to force inclusion here)
      if(length(cur_region$tb_snps) > 0){
        tags_to_add[[i]] <- cur_region$tb_snps
      }else{
        tags_to_add[[i]] <- character()
      }
      if(!all(is.na(cur_region$existing_tag_snps$ID))){
        existing_MI <- MI_Calculation(cur_region$candidate_geno,cur_region$existing_geno,cur_region$candidate_snps$ID,cur_region$existing_tag_snps$ID,n_cores)
        #Find best tag so far for each candidate SNP
        candidate_max_MI <- apply(existing_MI,1,function(x) max(x))
      }else{
        existing_MI <- matrix(nrow = ncol(cur_region$candidate_geno),ncol = 0)
        rownames(existing_MI) <- colnames(cur_region$candidate_geno)
        #No tag for each candidate SNP
        candidate_max_MI <- rep(0,ncol(cur_region$candidate_geno))
        names(candidate_max_MI) <- colnames(cur_region$candidate_geno)
      }
      #Find best tags to add
      pairwise_candidate_MI <- MI_Calculation(cur_region$candidate_tags_geno,cur_region$candidate_geno,cur_region$candidate_tag_snps$ID,cur_region$candidate_snps$ID,n_cores)
      tag_info[[i]] <- SelectTags(cur_region$candidate_snps$ID,cur_region$candidate_tag_snps,existing_MI,pairwise_candidate_MI,candidate_max_MI,returnTagInfo=T)
      MI_matrix <- SelectTags(cur_region$candidate_snps$ID,cur_region$candidate_tag_snps,existing_MI,pairwise_candidate_MI,candidate_max_MI)
      tags_to_add[[i]] <- c(tags_to_add[[i]],setdiff(colnames(MI_matrix),cur_region$existing_tag_snps$ID))
    }
    saveRDS(list(wgs_regions_to_update=wgs_regions_to_update,region_info=region_info,tags_to_add=tags_to_add,tag_info=tag_info),
            file = paste0(data_dir_path,'WGS_Addtl_Tags/',out_prefix,'.rds'))
    
    tags_to_add_merged <- unlist(tags_to_add)
    tags_to_add_merged <- tags_to_add_merged[!is.na(tags_to_add_merged)]
    #Add TB snps (even if they are not selected tags)
    print(paste0('tags added:',length(tags_to_add_merged)))
    all_tb_snps <- unlist(lapply(snps_in_gene_regions_parsed$snps_in_gene_regions_parsed,function(x) x$cur_region_tb_gwas))
    tags_to_add_merged <- union(tags_to_add_merged,all_tb_snps)
    print(paste0('tags added:',length(tags_to_add_merged)))
    h3a_rsid <- system(paste0('~/Software/bcftools query -f "%ID\n" ',existing_tags),intern = T)
    print(paste0('h3a added:',length(h3a_rsid)))
    print(paste0('h3a and tags added:',length(union(h3a_rsid,tags_to_add_merged))))
    write(union(h3a_rsid,tags_to_add_merged),file = paste0(data_dir_path,'WGS_Addtl_Tags/h3a.with.',out_prefix,'.txt'),sep = '\n',append = F)
  }else{
    tag_summary <- readRDS(paste0(data_dir_path,'WGS_Addtl_Tags/',out_prefix,'.rds'))
    snps_in_gene_regions_parsed <- readRDS(snps_in_gene_regions_parsed_path)
    h3a_rsid <- system(paste0('~/Software/bcftools query -f "%ID\n" ',existing_tags),intern = T)
    
    random_tags <- c()
    match_ind <- sapply(tag_summary$region_info,function(x) which(sapply(snps_in_gene_regions_parsed$snps_in_gene_regions_parsed,function(y) any(x$candidate_snps$ID %in% y$wgs_snps$ID))))
    region_info_matched <- snps_in_gene_regions_parsed$snps_in_gene_regions_parsed[match_ind]
    for(i in 1:length(region_info_matched)){
      cur_num_tags <- length(setdiff(union(tag_summary$tags_to_add[[i]],region_info_matched[[i]]$cur_region_tb_gwas),h3a_rsid))
      candidate_random_tags <- setdiff(region_info_matched[[i]]$wgs_snps$ID,
                                       union(union(tag_summary$tags_to_add[[i]],region_info_matched[[i]]$cur_region_tb_gwas),h3a_rsid))
      if(length(candidate_random_tags) < cur_num_tags){
        random_tags <- c(random_tags,candidate_random_tags,sample(union(tag_summary$tags_to_add[[i]],region_info_matched[[i]]$cur_region_tb_gwas),
                                                                  size = cur_num_tags - length(candidate_random_tags),replace = F))
        print('Selected Tags Added')
      }else{
        random_tags <- c(random_tags,sample(candidate_random_tags,size = cur_num_tags,replace = F))
      }
    }
    
    region_info_matched_no_tags_selected <- snps_in_gene_regions_parsed$snps_in_gene_regions_parsed[-match_ind]
    for(i in 1:length(region_info_matched_no_tags_selected)){
      cur_tb_snps <- setdiff(region_info_matched_no_tags_selected[[i]]$cur_region_tb_gwas,h3a_rsid)
      if(length(cur_tb_snps) > 0){
        random_tags <- c(random_tags,sample(setdiff(region_info_matched_no_tags_selected[[i]]$wgs_snps$ID,h3a_rsid),size = length(cur_tb_snps),replace = F))
      }
    }
      
    write(union(h3a_rsid,random_tags),file = paste0(data_dir_path,'WGS_Addtl_Tags/h3a.with.',out_prefix,'.random.txt'),sep = '\n',append = F)
  }

  #Write add-on snps
  if(!random){
    WriteAddOnFile(lift_over_vcf,phased_vcf,out_prefix,paste0(data_dir_path,'WGS_Addtl_Tags/'),random)
  }else{
    #Write random snps
    WriteAddOnFile(lift_over_vcf,phased_vcf,paste0(out_prefix,'.random'),paste0(data_dir_path,'WGS_Addtl_Tags/'),random)
  }
}
args <- commandArgs(trailingOnly = T)

#Write out selected tags along with H3Africa tags as VCF files. 
WriteSelectedTags(out_prefix = args[[1]],
                  site_type = args[[2]],
                  snps_in_gene_regions_parsed_path = args[[3]],
                  keep_samples_path = args[[4]],
                  existing_tags = args[[5]],
                  lift_over_vcf = args[[6]],
                  phased_vcf = args[[7]],
                  SpreadThresh=0.9,
                  PercAboveThresh=0.9,
                  random = F,
                  n_cores = args[[8]])

WriteSelectedTags(out_prefix = args[[1]],
                  site_type = args[[2]],
                  snps_in_gene_regions_parsed_path = args[[3]],
                  keep_samples_path = args[[4]],
                  existing_tags = args[[5]],
                  lift_over_vcf = args[[6]],
                  phased_vcf = args[[7]],
                  SpreadThresh=0.9,
                  PercAboveThresh=0.9,
                  random = T,
                  n_cores = args[[8]])

#Write out selected tags along with H3Africa tags as VCF files. 
# WriteSelectedTags(out_prefix = 'addlt.tags.setting1',
#                   site_type = 'wgs',
#                   snps_in_gene_regions_parsed_path = '../results/Baseline/Tag_SNP_Selection/TB_gene_regions_baseline.rds',
#                   keep_samples_path = '../data/WGS_Host_Data/testing_set/testing_split.txt',
#                   existing_tags = '../data/WGS_Host_Data/joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz',
#                   lift_over_vcf = '../data/WGS_Host_Data/joined.hg19.nodup.nomismap.vcf.gz',
#                   phased_vcf = '../results/testing_vcfs/',
#                   SpreadThresh=0.9,
#                   PercAboveThresh=0.9,
#                   random = F,
#n_cores = n_cores)
# WriteSelectedTags(out_prefix = 'addlt.tags.setting1',
#                   site_type = 'wgs',
#                   snps_in_gene_regions_parsed_path = '../results/Baseline/Tag_SNP_Selection/TB_gene_regions_baseline.rds',
#                   keep_samples_path = '../data/WGS_Host_Data/testing_set/testing_split.txt',
#                   existing_tags = '../data/WGS_Host_Data/joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz',
#                   lift_over_vcf = '../data/WGS_Host_Data/joined.hg19.nodup.nomismap.vcf.gz',
#                   phased_vcf = '../results/testing_vcfs/',
#                   SpreadThresh=0.9,
#                   PercAboveThresh=0.9,
#                   random = T,
#                   n_cores = n_cores)