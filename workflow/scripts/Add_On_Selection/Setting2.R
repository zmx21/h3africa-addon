#Select Add-on SNPs in all other regions
library(infotheo)
library(filematrix)
library(dplyr)
library(snpStats)
library(here)
library(pbmcapply)
library(future.apply)
library(glue)
library(GenomicRanges)
library(IRanges)
library(GWASTools)
#Extract genome wdie any snp which are poorly imputed
GetPoorlyImputedSNPs <- function(afgr_info_path,wgs_info_path,info_thresh,MAF_thresh){
  afgr_info <- data.table::fread(afgr_info_path,header = T)
  wgs_pos <- data.table::fread(wgs_info_path,header=F)
  colnames(wgs_pos) <- c('CHROM','POS','ID','REF','ALT','AF')
  
  #Filter AFGR first to limit to poorly imputed SNPs, and above MAF threshold
  afgr_info_filt <- dplyr::filter(afgr_info,INFO < info_thresh & AF > MAF_thresh & AF < 1 - MAF_thresh)
  wgs_pos_filt <- dplyr::filter(wgs_pos,AF > MAF_thresh & AF < 1 - MAF_thresh)
  wgs_pos_filt$CHROM <- sapply(wgs_pos_filt$CHROM,function(x) gsub(x=x,pattern = 'chr',replacement = ''))
  #Filter for only AFGR sites that was sequenced
  afgr_info_filt <- dplyr::inner_join(afgr_info_filt,wgs_pos_filt,by=c('CHROM','POS','REF','ALT')) %>% dplyr::select(CHROM,POS,ID=ID.y,REF,ALT,AF=AF.x,INFO)
  afgr_info_filt$MAF <- sapply(afgr_info_filt$AF,function(x) min(c(x,1-x)))
  return(afgr_info_filt)
}
#Get haploblocks which overlap with poorly imputed SNPs
AssignHaploBlock <- function(poorly_imputed_snps,haploblock_df,offset){
  offset_start_region <- lapply(poorly_imputed_snps$POS,function(x) x - offset)
  offset_end_region <- lapply(poorly_imputed_snps$POS,function(x) x + offset)
  
  start_end_regions <- pbmclapply(1:nrow(poorly_imputed_snps),function(i) {
    cur_snp <- poorly_imputed_snps[i,]
    cur_haplo_block <- dplyr::filter(haploblock_df,cur_snp$CHROM == haploblock_df$CHROM & 
                                       cur_snp$POS >= START & cur_snp$POS <= END)
    if(nrow(cur_haplo_block) == 0){
      start_region <- offset_start_region[i]
      end_region <- offset_end_region[i]
    }else if(cur_haplo_block$END - cur_haplo_block$START < offset){
      start_region <- offset_start_region[i]
      end_region <- offset_end_region[i]
    }else{
      start_region <- cur_haplo_block$START
      end_region<- cur_haplo_block$END
    }
    return(list(start_region=start_region,end_region=end_region))
  },mc.cores = 3)
  poorly_imputed_snps$start_region <- unlist(lapply(start_end_regions,function(x) x$start_region),recursive = T)
  poorly_imputed_snps$end_region <- unlist(lapply(start_end_regions,function(x) x$end_region),recursive = T)
  poorly_imputed_snps$chromosome_name <- poorly_imputed_snps$CHROM
  return(poorly_imputed_snps)
}
#For overlapping haploblocks/regions, merge into one region. 
MergeSNPs <- function(poorly_imputed_snps_grouped){
  #Make sure regions are non_overlapping
  poorly_imputed_snps_grouped_unique <- poorly_imputed_snps_grouped[0,]
  for(i in 1:length(unique(poorly_imputed_snps_grouped$CHROM))){
    cur_chr <- dplyr::filter(poorly_imputed_snps_grouped,CHROM == unique(poorly_imputed_snps_grouped$CHROM)[i])
    cur_chr_non_overlapping <- cur_chr[0,]
    overlaps <- IRanges(cur_chr$start_region,cur_chr$end_region)
    group_assignment <- subjectHits(findOverlaps(overlaps, reduce(overlaps)))
    for(j in group_assignment){
      cur_group <- cur_chr[which(group_assignment==j),]
      cur_group_unique <- cur_group[1,]
      cur_group_unique$start_region <- min(cur_group$start_region)
      cur_group_unique$end_region <- max(cur_group$end_region)
      cur_group_unique$ID <- paste0(cur_group$ID,collapse = ',')
      cur_group_unique$REF <- paste0(cur_group$REF,collapse = ',')
      cur_group_unique$ALT <- paste0(cur_group$ALT,collapse = ',')
      cur_group_unique$INFO <- paste0(cur_group$INFO,collapse = ',')
      cur_group_unique$MAF <- paste0(cur_group$MAF,collapse = ',')
      cur_group_unique$CHROM <- unique(poorly_imputed_snps_grouped$CHROM)[i]
      cur_group_unique$POS <- paste0(cur_group$POS,collapse = ',')
      cur_group_unique$n_snps <- length(cur_group$ID)
      cur_chr_non_overlapping <- rbind(cur_chr_non_overlapping,cur_group_unique)
    }
    poorly_imputed_snps_grouped_unique <- rbind(poorly_imputed_snps_grouped_unique,cur_chr_non_overlapping)
  }
  poorly_imputed_snps_grouped_unique <- dplyr::distinct(poorly_imputed_snps_grouped_unique,CHROM,start_region,end_region,.keep_all=T)
  poorly_imputed_snps_grouped_unique$hgnc_symbol <- poorly_imputed_snps_grouped_unique$ID
  return(poorly_imputed_snps_grouped_unique)
}
#Get candidate target/tag SNPs and sequenced SNPs within a region
GetRegionInfo <- function(poorly_imputed_snps_grouped_merged,snps_in_region = NULL,INFO_thresh=0.5,excl_snps=c()){
  tb_regions_unique <- GetTBRegion()
  #Get Avail SNPs in each defined region
  if(is.null(snps_in_region)){
    snps_in_region <- GetTargetRegionCoverage(region_info = poorly_imputed_snps_grouped_merged,
                                              tb_gwas_snps=tb_regions_unique$tb_gwas_snps,
                                              afgr_imp_path = '~/G2G_TB/WGS_Host_Data/WGS_Addtl_Tags/AFGRAndWGSMergedImpSecondPass/WGS_AFGR_Merged_Info.txt',
                                              addtl_snp_path = '~/G2G_TB/WGS_Host_Data/WGS_Addtl_Tags/second_pass/joined.hg19.nodup.nomismap.addlt.WGS.tags.second.pass.0p9.testing.vcf.gz',
                                              merged = T,
                                              pass=2)
    return(snps_in_region)
  }
  #Get Genotypes of Avail SNPs in each region
  software_dir <- paste0(here(),'/software/')
  data_dir_path <- '/mnt/data2/xu/G2G_TB/WGS_Host_Data/'
  
  training_samples <- readRDS('~/G2G_TB/WGS_Host_Data/training_set/training_split.rds')
  wgs_mat <- fm.open('~/G2G_TB/Imputation_Eval/wgs_mat')
  wgs_colnames <- data.table::fread('~/G2G_TB/Imputation_Eval/wgs_mat.nmscol.txt',sep = "",header = F)$V1
  wgs_rownames <- rownames(wgs_mat)
  
  region_info <- pbmclapply(1:length(snps_in_region),function(i){
    if(unique(snps_in_region[[i]]$wgs_snps$CHROM) == 'X'){
      return(GetCandidateSNPs(candidate_sites = 'wgs',
                       region_id = i,
                       training_samples = training_samples,
                       snps_in_gene_regions_parsed=snps_in_region,
                       wgs_mat = wgs_mat,
                       wgs_colnames = wgs_colnames,
                       wgs_rownames = wgs_rownames,
                       illumina_scores = data.table::fread('~/G2G_TB/Tag_SNP_Selection/WGS.Setting2.Illumina.Scores.csv'),
                       INFO_thresh = INFO_thresh,
                       excl_snps = excl_snps))
      
    }else{
      return(GetCandidateSNPs(candidate_sites = 'wgs',
                       region_id = i,
                       training_samples = training_samples,
                       snps_in_gene_regions_parsed=snps_in_region,
                       wgs_mat = wgs_mat,
                       wgs_colnames = wgs_colnames,
                       wgs_rownames = wgs_rownames,
                       illumina_scores = data.table::fread('~/G2G_TB/Tag_SNP_Selection/WGS.Setting2.Illumina.Scores.csv'),
                       INFO_thresh = INFO_thresh,
                       excl_snps = c()))
    }
    },mc.cores = 8)
  
  return(region_info)
}
#Calculate MI between two SNPs
GetMIInfo <- function(region_info){
  MI_info <- list()
  for(i in 1:length(region_info)){
    print(i)
    #Calculate MI of candidate snps against existing tags
    cur_region <- region_info[[i]]
    if(all(is.na(cur_region$candidate_tag_snps)) | all(is.na(cur_region$candidate_snps))){
      MI_info[[i]] <- list(existing_MI=NA,candidate_max_MI=NA,pairwise_candidate_MI=NA)
      next
    }
    if(!all(is.na(cur_region$existing_tag_snps$ID))){
      existing_MI <- MI_Calculation(cur_region$candidate_geno,cur_region$existing_geno,cur_region$candidate_snps$ID,cur_region$existing_tag_snps$ID)
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
    pairwise_candidate_MI <- MI_Calculation(cur_region$candidate_tags_geno,cur_region$candidate_geno,cur_region$candidate_tag_snps$ID,cur_region$candidate_snps$ID)
    MI_info[[i]] <- list(existing_MI=existing_MI,candidate_max_MI=candidate_max_MI,pairwise_candidate_MI=pairwise_candidate_MI)
    # region_best_tag <- SelectTags(candidate_snps=cur_region$candidate_snps$ID,
    #                            candidate_tag_snps=cur_region$candidate_tag_snps$ID,
    #                            MI_matrix = existing_MI,
    #                            pairwise_candidate_MI = pairwise_candidate_MI,
    #                            candidate_max_MI = candidate_max_MI,
    #                            max_iter = 1,
    #                            returnTagInfo = T)
    # MI_info$BestTag[i] <- region_best_tag[[1]]$bestTag
    # MI_info$BestMIDiff[i] <- as.numeric(region_best_tag[[1]]$bestMIDiff)
    # MI_info$N_Improv[i] <- as.numeric(region_best_tag[[1]]$n_improv)
  }
  return(MI_info)
}
#Select best add-on tag SNP, concurrently across all regions. 
SelectBestGenomeWideTags <- function(region_info,MI_info,method){
  max_iter <- 5000
  iter <- 1
  tags_to_add <- c()
  max_tag_info <- list()
  best_tag_info <- pbmclapply(1:length(region_info),function(i){
    if(all(is.na(region_info[[i]]$candidate_snps$ID)) | all(is.na(region_info[[i]]$candidate_tag_snps$ID)) ){
      return(list(bestTag=NA,bestMIDiff=NA,n_improv=NA))
    }else{
      best_tag <- SelectTags(candidate_snps=region_info[[i]]$candidate_snps$ID,
                             candidate_tag_snps=region_info[[i]]$candidate_tag_snps,
                             MI_matrix = MI_info[[i]]$existing_MI,
                             pairwise_candidate_MI = MI_info[[i]]$pairwise_candidate_MI,
                             candidate_max_MI = MI_info[[i]]$candidate_max_MI,
                             max_iter = 1,
                             returnTagInfo = T)
      return(best_tag[[1]])
    }
  },mc.cores = 5,mc.preschedule = F)
  
  while (iter <= max_iter) {
    print(iter)
    if(method == 'Avg'){
      MI_diff <- sapply(best_tag_info,function(x) ifelse(is.null(x$bestMIDiff),NA,ifelse(x$n_improv==1,0,x$bestMIDiff/x$n_improv)))
    }else if(method == 'Sum'){
      MI_diff <- sapply(best_tag_info,function(x) ifelse(is.null(x$bestMIDiff),as.numeric(NA),x$bestMIDiff))
    }
    MI_diff_rank <- rank(MI_diff[which(MI_diff != 0 & !is.na(MI_diff))])
    TagScore <- sapply(best_tag_info,function(x) ifelse(is.null(x$bestProbeScore),NA,x$bestProbeScore))
    TagScoreRank <- rank(TagScore[which(MI_diff != 0 & !is.na(MI_diff))])
    AvgRank <- rowMeans(cbind(MI_diff_rank,TagScoreRank))
    max_index <- which(MI_diff != 0 & !is.na(MI_diff))[which.max(AvgRank)]
    max_info <- best_tag_info[[max_index]]
    max_info$region_id <- max_index
    max_tag_info[[iter]] <- max_info
    region_info[[max_index]]$existing_tag_snps <- rbind(region_info[[max_index]]$existing_tag_snps,
                                                                dplyr::filter(region_info[[max_index]]$candidate_tag_snps,ID == max_info$bestTag) %>% dplyr::select(colnames(region_info[[max_index]]$existing_tag_snps)))
    if(all(is.na(region_info[[max_index]]$existing_geno))){
      region_info[[max_index]]$existing_geno <- region_info[[max_index]]$candidate_tags_geno[,colnames(region_info[[max_index]]$candidate_tags_geno) == max_info$bestTag,drop = F]
    }else{
      region_info[[max_index]]$existing_geno <- cbind(region_info[[max_index]]$existing_geno,region_info[[max_index]]$candidate_tags_geno[,colnames(region_info[[max_index]]$candidate_tags_geno) == max_info$bestTag,drop = F])
    }

    region_info[[max_index]]$candidate_snps <- dplyr::filter(region_info[[max_index]]$candidate_snps,ID != max_info$bestTag)
    region_info[[max_index]]$candidate_tag_snps <- dplyr::filter(region_info[[max_index]]$candidate_tag_snps,ID != max_info$bestTag)
    region_info[[max_index]]$candidate_geno <- region_info[[max_index]]$candidate_geno[,colnames(region_info[[max_index]]$candidate_geno) != max_info$bestTag,drop=F]
    region_info[[max_index]]$candidate_tags_geno <- region_info[[max_index]]$candidate_tags_geno[,colnames(region_info[[max_index]]$candidate_tags_geno) != max_info$bestTag,drop = F]
    
    existing_MI <- MI_info[[max_index]]$existing_MI
    pairwise_candidate_MI <- MI_info[[max_index]]$pairwise_candidate_MI
    candidate_max_MI <- MI_info[[max_index]]$candidate_max_MI
    new_existing_MI <- cbind(existing_MI,t(pairwise_candidate_MI[max_info$bestTag,rownames(existing_MI),drop=F]))
    new_existing_MI <- new_existing_MI[rownames(new_existing_MI) != max_info$bestTag,,drop=F]
    new_pairwise_candidate_MI <- pairwise_candidate_MI[which(rownames(pairwise_candidate_MI) != max_info$bestTag),which(colnames(pairwise_candidate_MI) != max_info$bestTag),drop=F]
    new_candidate_max_MI <- apply(new_existing_MI,1,function(x) max(x))
    #new_MI_info <- unlist(GetMIInfo(region_info[max_index]),recursive = F)
    MI_info[[max_index]] <- list(existing_MI=new_existing_MI,
                                 pairwise_candidate_MI=new_pairwise_candidate_MI,
                                 candidate_max_MI=new_candidate_max_MI)
    
    tags_to_add <- c(tags_to_add,max_info$bestTag)
    best_tag_info[[max_index]] <- unlist(lapply(max_index,function(i){
      if(all(is.na(region_info[[i]]$candidate_snps$ID)) | all(is.na(region_info[[i]]$candidate_tag_snps$ID)) ){
        return(list(bestTag=NA,bestMIDiff=NA,n_improv=NA))
      }else{
        best_tag <- SelectTags(candidate_snps=region_info[[i]]$candidate_snps$ID,
                               candidate_tag_snps=region_info[[i]]$candidate_tag_snps,
                               MI_matrix = MI_info[[i]]$existing_MI,
                               pairwise_candidate_MI = MI_info[[i]]$pairwise_candidate_MI,
                               candidate_max_MI = MI_info[[i]]$candidate_max_MI,
                               max_iter = 1,
                               returnTagInfo = T)
        return(best_tag[[1]])
      }
    }),recursive = F)
    # print(tags_to_add)
    iter <- iter + 1
  }
  return(list(tags_to_add=tags_to_add,max_tag_info=max_tag_info))
}
WriteAddOnFile <- function(lift_over_vcf,phased_vcf,out_prefix,out_path,random){
  addtl_tags_dir <- out_path
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
    addtl_tags_dir <- paste0(out_path,'/Setting2Random/')
  }else{
    addtl_tags_dir <- paste0(out_path,'/Setting2/')
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

#Write out add-on Tags, including Setting 1/2 and MT/Y SNPs. 
WriteGenomeWideTags <- function(lift_over_vcf,phased_vcf,out_prefix,out_path,best_tags,N_Probes_to_add = 5000,extra_snps,random_tags,random=T){
  data_dir_path <- '../data/WGS_Host_Data/'
  #All H3A IDs
  h3a_rsid<- system(paste0('~/Software/bcftools query -f "%ID\n" ',data_dir_path,'joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz'),intern = T)
  
  #H3A IDs which were not sequenced or did not pass WGS filters
  h3a_v1 <- data.table::fread('../data/H3Africa/v1/H3Africa_2017_20021485_A2.csv',skip=5)
  h3a_v2 <- data.table::fread('../data/H3Africa/v2/H3Africa_2019_20037295_A1.csv',skip=5)
  excl_h3a_tags <- setdiff(data.table::fread('../data/Chip_Data/h3achip_dbsnp150.tsv')$rsid,h3a_rsid)
  v2_specific <- dplyr::filter(h3a_v2,IlmnID %in% setdiff(h3a_v2$IlmnID,h3a_v1$IlmnID))
  v2_specific_rs <- sapply(v2_specific$Name,function(x){
    if(grepl(pattern='rs',x = x)){
      if(grepl(pattern='h3a_37_',x = x)){
        split_str <- strsplit(x,split = '_')[[1]]
        return(split_str[length(split_str)])
      }else{
        return(x)
      }
    }else{
      return(NA)
    }
  })
  excl_h3a_tags <- union(excl_h3a_tags,setdiff(v2_specific_rs,h3a_rsid))
  #Tags added by setting1 (TB Regions)
  addtl_tags <- data.table::fread('../results/WGS_Addtl_Tags/h3a.with.addlt.tags.setting1.txt',header = F) %>% dplyr::select(ID=V1)
  addtl_tags_rsid <- setdiff(addtl_tags$ID,h3a_rsid)
  #Calculate number of probes used in addtl tags, and the number of probes avaliable
  addtl_tags_tbl <- data.frame(ID=addtl_tags_rsid,stringsAsFactors = F) %>% dplyr::left_join(rbind(data.table::fread('../results/Tag_SNP_Selection/WGS.Baseline.Illumina.Scores.csv',header = T),
                                                                                                   data.table::fread('../results/Tag_SNP_Selection/WGS.TBSNPs.Illumina.Scores.csv',header = T),
                                                                                                   data.table::fread('../results/Tag_SNP_Selection/WGS.KIR.HLA.Illumina.Scores.csv',header = T)) %>% dplyr::distinct(Locus_Name,.keep_all=T),by=c('ID'='Locus_Name'))
  addtl_tags_tbl$N_Probes <- sapply(addtl_tags_tbl$Assay_Type,function(x) ifelse(x=='InfiniumII',1,2))
  probes_added <- sum(addtl_tags_tbl$N_Probes,na.rm = T)
  probes_remaining <- N_Probes_to_add - probes_added
  #Add Specified extra SNPs
  best_tags$tags_to_add <- c(extra_snps,best_tags$tags_to_add)
  random_tags <- c(extra_snps,random_tags)
  #Exclude tags added in second pass
  #Add all SNPs in Mitochondria with MAF > 0.05
  existing_h3a_mt <- dplyr::filter(h3a_v2,Chr == 'MT')
  system(paste0("~/Software/bcftools view -r chrM ",data_dir_path,"WGS_Fellay.hg38.joint.118h-1947437863.genotyped.renamed.consensus.filt.vcf.gz | ~/Software/bcftools query -f '%ID %AF %POS %REF %ALT\n' > ",data_dir_path,"MT_pos.txt"),intern = T)
  mt_rsid <- data.table::fread(paste0(data_dir_path,"MT_pos.txt"))
  colnames(mt_rsid) <- c('ID','AF','POS','REF','ALT')
  mt_rsid_filt <- dplyr::filter(mt_rsid,AF > 0.05 & AF < 0.95) %>% dplyr::filter(ID != '.') %>% dplyr::filter(!POS %in% existing_h3a_mt$MapInfo)
  mt_rsid_tbl <- dplyr::left_join(data.frame(ID=mt_rsid_filt$ID),
                                  data.table::fread('../results/Tag_SNP_Selection/WGS.MT.Illumina.Scores.csv',header = T),by=c('ID'='Locus_Name')) %>% dplyr::filter(!is.na(Assay_Type))
  mt_rsid_tbl$N_Probes <- sapply(mt_rsid_tbl$Assay_Type,function(x) ifelse(x=='InfiniumII',1,2))
  probes_remaining <- probes_remaining - sum(mt_rsid_tbl$N_Probes)
  #Exclude tags added in MT
  best_tags$tags_to_add <- best_tags$tags_to_add[sort(which(is.na(match(best_tags$tags_to_add,mt_rsid_tbl$ID))))]
  random_tags <- random_tags[sort(which(is.na(match(random_tags,mt_rsid_tbl$ID))))]
  
  #Exclude h3a tags which were not sequenced (wil be added anyways)
  random_tags <- random_tags[sort(which(is.na(match(best_tags$tags_to_add,excl_h3a_tags))))]
  best_tags$tags_to_add <- best_tags$tags_to_add[sort(which(is.na(match(best_tags$tags_to_add,excl_h3a_tags))))]

  
  #Add Tags by Setting 2
  tags_to_add_tbl<- dplyr::left_join(data.frame(ID=best_tags$tags_to_add,stringsAsFactors = F),
                                     rbind(data.table::fread('../results/Tag_SNP_Selection/WGS.Setting2.Illumina.Scores.csv',header = T),
                                           data.table::fread('../results/Tag_SNP_Selection/WGS.ExtraSNPs.Illumina.Scores.csv',header = T),
                                           data.table::fread('../results/Tag_SNP_Selection/WGS.MT.Illumina.Scores.csv',header = T),
                                           data.table::fread('../results/Tag_SNP_Selection/WGS.Y.Illumina.Scores.csv',header = T)) %>% 
                                       dplyr::distinct(Locus_Name,.keep_all=T),by=c('ID'='Locus_Name'))
  tags_to_add_tbl$N_Probes <- sapply(tags_to_add_tbl$Assay_Type,function(x) ifelse(x=='InfiniumII',1,2))
  cum_N_Probes <- cumsum(tags_to_add_tbl$N_Probes)
  last_tag_index <- ifelse(length(which(cum_N_Probes == probes_remaining)) == 0,which(cum_N_Probes - probes_remaining == -1),which(cum_N_Probes == probes_remaining))
  all_tags <- c(mt_rsid_tbl$ID,addtl_tags$ID,best_tags$tags_to_add[1:(last_tag_index)])

  if(random){
    random_addtl_tags <- data.table::fread('../results/WGS_Addtl_Tags/h3a.with.addlt.tags.setting1.random.txt',header = F) %>% dplyr::select(ID=V1)
    random_addtl_tags_rsid <- setdiff(random_addtl_tags$ID,h3a_rsid)
    
    all_random_tags <- c(mt_rsid_tbl$ID,random_addtl_tags$ID,random_tags[1:(last_tag_index)])
    write(all_random_tags,file = paste0('../results/WGS_Addtl_Tags/h3a.with.',out_prefix,'.txt'))
    write(c(mt_rsid_tbl$ID,random_addtl_tags_rsid,random_tags[1:(last_tag_index)]),file = paste0('../results/WGS_Addtl_Tags/',out_prefix,'.txt'))

    WriteAddOnFile(lift_over_vcf,phased_vcf,out_prefix,out_path,random)
      
  }else{
    write(all_tags,file = paste0('../results/WGS_Addtl_Tags/h3a.with.',out_prefix,'.txt'))
    write(c(mt_rsid_tbl$ID,addtl_tags_rsid,best_tags$tags_to_add[1:(last_tag_index)]),file = paste0('../results/WGS_Addtl_Tags/',out_prefix,'.txt'))
    write(mt_rsid_tbl$ID,paste0('../results/WGS_Addtl_Tags/mt.tags.setting2.txt'))
    
    WriteAddOnFile(lift_over_vcf,phased_vcf,out_prefix,out_path,random)
    
  }
  
}

#Get Haploblocks
pbmclapply(1:23,function(i) system(paste0('./scripts/Add_On_Selection/estimate_haploblock.sh ',i),intern = T),mc.cores = 23)
haploblock_df <- do.call(rbind,(lapply(1:23,function(x) data.table::fread(paste0('~/G2G_TB/WGS_Host_Data/joined.hg19.nodup.nomismap.',x,'.blocks.det'),header = T)))) %>% dplyr::select(CHROM=CHR,START=BP1,END=BP2)
haploblock_df$CHROM[which(haploblock_df$CHROM == 23)] <- 'X'

#Get poorly imputed snps genomewide
afgr_info_path <- '../results/WGS_Addtl_Tags/Setting1/WGS_AFGR_Imputed/WGS_AFGR_Merged_Info.txt'
wgs_info_path <- '../data/WGS_Host_Data/joined.hg19.nodup.nomismap.vcf.gz.info'
poorly_imputed_snps <- GetPoorlyImputedSNPs(afgr_info_path,wgs_info_path,info_thresh = 0.8,MAF_thresh = 0.05)
#Merge overlapping regions which contain poorly imputed SNPs
poorly_imputed_snps_grouped <- AssignHaploBlock(poorly_imputed_snps,haploblock_df,offset = 5000)
save.image('../results/Tag_SNP_Selection/setting2_step1.rda')
poorly_imputed_snps_grouped_merged <- MergeSNPs(poorly_imputed_snps_grouped)
save.image('../results/Tag_SNP_Selection/setting2_step2.rda')
#Get sequenced SNPs in regions.
snps_in_region <- GetRegionInfo(dplyr::filter(poorly_imputed_snps_grouped_merged,n_snps > 1),snps_in_region = NULL,INFO_thresh=0.8,excl_snps = data.table::fread('~/G2G_TB/WGS_Host_Data/phased_no_geno_filt_correct_ref/excl_hap_het_id.txt',header = F)$V1)
save.image('../results/Tag_SNP_Selection/setting2_step3.rda')
#Remove snps at the end or start of chr (and pseudoautosomal for X and Y)
snps_in_region_chr <- sapply(snps_in_region,function(x) unique(x$imp_snps$CHROM))
h3a_pos <- data.table::fread('../data/H3Africa/v2/H3Africa_2019_20037295_A1.csv',skip = 5) %>% dplyr::select(CHR = Chr,POS = MapInfo)
h3a_pos_ranges <- lapply(unique(snps_in_region_chr),function(x) range(dplyr::filter(h3a_pos,CHR==x)$POS))
names(h3a_pos_ranges) <- unique(snps_in_region_chr)
index_to_remove <- c()
for(i in 1:length(snps_in_region)){
  print(i/length(snps_in_region))
  cur_chr <- snps_in_region_chr[i]
  cur_range <- h3a_pos_ranges[[cur_chr]]
  cur_coord <- snps_in_region[[i]]$start_region:snps_in_region[[i]]$end_region
  if(cur_chr == 'X'){
    data(pseudoautosomal.hg19)
    if(any(cur_coord < cur_range[1] | cur_coord > cur_range[2])){
      index_to_remove <- c(index_to_remove,i)
    }else if (any(cur_coord >= pseudoautosomal.hg19$start.base[2] & cur_coord <= pseudoautosomal.hg19$end.base[2])) {
      index_to_remove <- c(index_to_remove,i)
    }else if(any(cur_coord >= pseudoautosomal.hg19$start.base[1] & cur_coord <= pseudoautosomal.hg19$end.base[1])){
      index_to_remove <- c(index_to_remove,i)
    }
  }else{
    if(any(cur_coord < cur_range[1] | cur_coord > cur_range[2])){
      index_to_remove <- c(index_to_remove,i)
    }
  }
}
snps_in_region <- snps_in_region[-index_to_remove]
#Get candidate target SNPs which are illumina probes
wgs_candidate_snps <- do.call(rbind,lapply(snps_in_region,function(x) dplyr::filter(x$wgs_snps,MAF > 0.05)))
#Upload this file to illumina, and name the result as '../results/Tag_SNP_Selection/WGS.Setting2.Illumina.Scores.csv'
data.table::fwrite(dplyr::select(wgs_candidate_snps,Locus_Name=ID),file = '../results/Tag_SNP_Selection/setting2_rsid.csv')
save.image('../results/Tag_SNP_Selection/setting2_step4.rda')
#Get candidate target/tag SNPs in regions.
region_info <- GetRegionInfo(dplyr::filter(poorly_imputed_snps_grouped_merged,n_snps > 2),snps_in_region = snps_in_region,
                             excl_snps = data.table::fread('~/G2G_TB/WGS_Host_Data/phased_no_geno_filt_correct_ref/excl_hap_het_id.txt',header = F)$V1)
save.image('../results/Tag_SNP_Selection/setting2_step5.rda')
#Within each region, calculate MI between SNPs
MI_info <- GetMIInfo(region_info)
save.image('../results/Tag_SNP_Selection/setting2.rda')
#Select best add-on tags concurrently across all regions.
best_tags <- SelectBestGenomeWideTags(region_info,MI_info,method = 'Sum')
saveRDS(best_tags,file = '../results/WGS_Addtl_Tags/addlt.tags.setting2.rds')

set.seed(888)
#Select Random Tags
addtl_tags_random_setting1 <- setdiff(data.table::fread('../results/WGS_Addtl_Tags/h3a.with.addlt.tags.setting1.random.txt',header = F)$V1,
                                      data.table::fread('../data/WGS_Host_Data/joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz.pos',header=F)$V3)
random_tags <- c()
for(i in 1:length(best_tags$max_tag_info)){
  cur_tag <- best_tags$max_tag_info[[i]]$bestTag
  excl_tags <- c(cur_tag,snps_in_region[[best_tags$max_tag_info[[i]]$region_id]]$h3a_tags$ID,addtl_tags_random_setting1,random_tags)
  candidate_random_tags <- setdiff(snps_in_region[[best_tags$max_tag_info[[i]]$region_id]]$wgs_snps$ID,excl_tags)
  random_tags <- c(random_tags,sample(candidate_random_tags,size = 1,replace = F))
}


#Add known SNP according to (https://www.frontiersin.org/articles/10.3389/fimmu.2019.02248/full)
extra_snps <- 'rs763780'
#Add known SNPs according to (https://doi.org/10.3389/fgene.2019.00865)
extra_snps <- c('rs34536443','rs2278589','rs9637876','rs148722713','rs73226617',extra_snps)
#Add known snps according to GWAS
extra_snps <- c('rs72993272',extra_snps)
#Y phylogeny
y_phylogeny <- data.table::fread('../results/Tag_SNP_Selection/WGS.Y.Illumina.Scores.csv')$Locus_Name
mt_phylogeny <- data.table::fread('../results/Tag_SNP_Selection/WGS.MT.Illumina.Scores.csv')$Locus_Name
extra_snps <- c(extra_snps,y_phylogeny,mt_phylogeny)

write(c('Locus_Name',extra_snps),file = '../results/Tag_SNP_Selection/extra_snps.csv')

#Write out all add-on tag SNPs (Setting 1/2, Y, and MT)
best_tags_sum <- readRDS('../results/WGS_Addtl_Tags/addlt.tags.setting2.rds')
WriteGenomeWideTags('../data/WGS_Host_Data/joined.hg19.nodup.nomismap.vcf.gz','../results/testing_vcfs/',
                    'addlt.tags.setting2','../results/WGS_Addtl_Tags/',best_tags_sum,5000,extra_snps,random_tags,F)

WriteGenomeWideTags('../data/WGS_Host_Data/joined.hg19.nodup.nomismap.vcf.gz','../results/testing_vcfs/',
                    'addlt.tags.setting2.random','../results/WGS_Addtl_Tags/',best_tags_sum,5000,extra_snps,random_tags,T)
