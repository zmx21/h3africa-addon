library(dplyr)
library(glue)
library(pbmcapply)
library(GWASTools)

RunTaggerSetting1 <- function(Setting1_Out_Dir = '../results/Tagger/Setting1/'){
  WriteAddOnFile <- function(lift_over_vcf,phased_vcf,out_prefix,out_path,random){
    addtl_tags_dir <- out_path
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
  
  set.seed(888)
  Setting1_Region_Info <- readRDS('../results/WGS_Addtl_Tags/addlt.tags.setting1.rds')
  Setting1_TB_Regions_Raw <- readRDS( '../results/Baseline/Tag_SNP_Selection/TB_gene_regions_baseline.rds')
  matched_ind <- unlist(pbmclapply(Setting1_Region_Info$region_info,function(x) which(sapply(Setting1_TB_Regions_Raw$snps_in_gene_regions_parsed,
                                                                                             function(y) any(x$candidate_snps$ID %in% y$wgs_snps$ID))),mc.cores = 20))
  Setting1_TB_Regions <- Setting1_TB_Regions_Raw$snps_in_gene_regions_parsed[matched_ind]
  existing_tags = '../data/WGS_Host_Data/joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz'
  h3a_rsid <- system(paste0('~/Software/bcftools query -f "%ID\n" ',existing_tags),intern = T)
  
  raw_vcf <- '../data/WGS_Host_Data/training_set/joined.hg19.nodup.nomismap.training.vcf.gz'
  system(glue::glue("~/Software/plink --vcf {raw_vcf} --make-bed --out {Setting1_Out_Dir}joined.hg19.nodup.nomismap.training"))
  
  tagger_tags <- list()
  for(i in 1:length(Setting1_Region_Info$region_info)){
    print(i)
    n_tags_to_select <- length(setdiff(Setting1_Region_Info$tags_to_add[[i]],h3a_rsid))
    if(n_tags_to_select < 0){
      tagger_tags[[i]] <- c()
      next
    }
    
    cur_dir <- glue::glue("{Setting1_Out_Dir}Region{i}/")
    system(glue::glue("mkdir -p {cur_dir}"))
    
    #Existing H3A Tags
    cur_existing_tag <- Setting1_Region_Info$region_info[[i]]$existing_tag_snps$ID
    write(cur_existing_tag,glue::glue("{cur_dir}existing_tags.txt"))
    
    #Candidate Tags in region
    cur_candidate_tags <- Setting1_Region_Info$region_info[[i]]$candidate_tag_snps %>% dplyr::select(ID,Final_Score,Assay_Type) %>% dplyr::distinct(ID,.keep_all = T)
    cur_candidate_tags$Final_Score[cur_candidate_tags$Final_Score > 1] <- 1
    
    #All SNPs in the region
    if(grepl(x=Setting1_Region_Info$region_info[[i]]$region_id,pattern = 'HLA') | grepl(x=Setting1_Region_Info$region_info[[i]]$region_id,pattern = 'KIR')){
      all_snps <- Setting1_TB_Regions[[i]]$wgs_snps %>% dplyr::filter(MAF > 0.01)
    }else{
      all_snps <- Setting1_TB_Regions[[i]]$wgs_snps %>% dplyr::filter(MAF > 0.05)
    }
    all_snps <- union(union(cur_existing_tag,all_snps$ID),cur_candidate_tags$ID)
    write(all_snps,glue::glue("{cur_dir}all_snps.txt"))
    
    #Excl Tags = A SNP but not a candidate Tag
    excl_tags <- setdiff(all_snps,union(cur_candidate_tags$ID,cur_existing_tag))
    write(excl_tags,glue::glue("{cur_dir}excl_tags.txt"))
    
    
    #Candidate Tags in region
    design_scores <- cur_candidate_tags %>% dplyr::select(ID,Final_Score)
    data.table::fwrite(design_scores,file = glue::glue("{cur_dir}designscores.txt"),sep = ' ',col.names = F,row.names = F,quote = F)
    
    #Extract Ped File
    cur_chr <- unique(Setting1_Region_Info$region_info[[i]]$candidate_tag_snps$CHROM)
    if(cur_chr == 'X'){
      system(glue::glue("~/Software/plink --bfile {Setting1_Out_Dir}joined.hg19.nodup.nomismap.training --chr X XY --extract {cur_dir}all_snps.txt --recode --out {cur_dir}Region{i}"))
    }else{
      system(glue::glue("~/Software/plink --bfile {Setting1_Out_Dir}joined.hg19.nodup.nomismap.training --chr {cur_chr} --extract {cur_dir}all_snps.txt --recode --out {cur_dir}Region{i}"))
    }
    #Rewrite Info and Ped file
    ped_file <- data.table::fread(glue::glue("{cur_dir}Region{i}.ped"))
    ped_file$V1 <- ped_file$V2
    ped_file$V6 <- 0
    data.table::fwrite(ped_file,col.names = F,row.names = F,quote = F,sep = ' ',file = glue::glue("{cur_dir}Region{i}.ped"))
    map_file <- data.table::fread(glue::glue("{cur_dir}Region{i}.map"))
    data.table::fwrite(map_file %>% dplyr::select(V2,V4),col.names = F,row.names = F,quote = F,sep = ' ',file = glue::glue("{cur_dir}Region{i}.info"))
    
    #Run Tagger
    system(glue::glue("java -jar -Xmx150000m ~/Software/Haploview.jar -nogui -memory 150000 -pedfile {cur_dir}Region{i}.ped -info {cur_dir}Region{i}.info -designScores {cur_dir}designscores.txt -includeTagsFile {cur_dir}existing_tags.txt -excludeTagsFile {cur_dir}excl_tags.txt -pairwiseTagging -maxNumTags {n_tags_to_select + length(cur_existing_tag)} -minGeno 0.00001 -hwcutoff 0 -missingCutoff 1 -minMAF 0.0 -mintagdistance 0 -tagrsqcutoff 0.8 -out {cur_dir}Region{i} > {cur_dir}Region{i}.tagger.log"))
    selected_tags <- setdiff(data.table::fread(glue::glue("{cur_dir}Region{i}.TESTS"),header = F)$V1,cur_existing_tag)
    if(length(selected_tags) != n_tags_to_select){
      #If tagger cannot select enough tags (non-variable tags tagging itself), pick random tags to avoid bias.
      tagger_tags[[i]] <- c(selected_tags,sample(setdiff(cur_candidate_tags$ID,union(cur_existing_tag,selected_tags)),size = n_tags_to_select - length(selected_tags)))
      write(paste0('Tagger Failed:',i,'(',length(selected_tags),'/',n_tags_to_select,')'),file = '../results/Tagger/Setting1/tagger_tags.log',append = T)
    }else{
      tagger_tags[[i]] <- selected_tags
    }
    
  }
  saveRDS(tagger_tags,file = '../results/Tagger/Setting1/tagger_tags.rds')
  tags_to_add_merged <- unlist(tagger_tags)
  tags_to_add_merged <- tags_to_add_merged[!is.na(tags_to_add_merged)]
  system(paste0('mkdir -p ','../results/Tagger/Addtl_Tags/'))
  all_tb_snps <- unlist(lapply(Setting1_TB_Regions$snps_in_gene_regions_parsed,function(x) x$cur_region_tb_gwas))
  missing_tb_snps <- all_tb_snps[!all_tb_snps %in% union(h3a_rsid,unlist(Setting1_Region_Info$tags_to_add))]
  write(c(missing_tb_snps,union(h3a_rsid,tags_to_add_merged)),file = '../results/Tagger/Addtl_Tags/h3a.with.addlt.tags.setting1.tagger.txt',sep = '\n',append = F)
  WriteAddOnFile(lift_over_vcf = '../data/WGS_Host_Data/joined.hg19.nodup.nomismap.vcf.gz',phased_vcf = '../results/testing_vcfs/',out_prefix = 'addlt.tags.setting1.tagger',out_path = '../results/Tagger/Addtl_Tags/',random = F)
}


RunTaggerSetting2 <- function(Setting2_Out_Dir = '../results/Tagger/Setting2/'){
  load('../results/Tag_SNP_Selection/setting2.rda')
  remove(list = setdiff(ls(),c('region_info','Setting2_Out_Dir')))
  
  WriteAddOnFile <- function(lift_over_vcf,phased_vcf,out_prefix,out_path){
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
    
    addtl_tags_dir <- paste0(out_path,'/Setting2/')
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

  tags_incl <- data.table::fread('../results/WGS_Addtl_Tags/h3a.with.addlt.tags.setting2.txt',header = F)
  Setting2_Region_Info <- readRDS('../results/WGS_Addtl_Tags/addlt.tags.setting2.rds')
  regions_incl <- which(Setting2_Region_Info$tags_to_add %in% tags_incl$V1)

  Setting2_Region_Info$tags_to_add <- Setting2_Region_Info$tags_to_add[regions_incl]
  Setting2_Region_Info$tag_info <- Setting2_Region_Info$max_tag_info[regions_incl]
  Setting2_Region_Info$region_info <- region_info[sapply(Setting2_Region_Info$tag_info,function(x) x$region_id)]

  Setting2_TB_Regions_Raw <- readRDS( '../results/WGS_Addtl_Tags/Setting2/TB_gene_regions_Setting2_baseline.rds')

  existing_tags = data.table::fread('../results/Tagger/Addtl_Tags/h3a.with.addlt.tags.setting1.tagger.txt',header = F)$V1

  unique_regions <- unique(sapply(Setting2_Region_Info$tag_info,function(x) x$region_id))
  Setting2_TB_Regions <- Setting2_TB_Regions_Raw$snps_in_gene_regions_parsed[unique_regions]

  if(!file.exists(glue::glue("{Setting2_Out_Dir}joined.hg19.nodup.nomismap.training.bim"))){
    raw_vcf <- '../data/WGS_Host_Data/training_set/joined.hg19.nodup.nomismap.training.vcf.gz'
    system(glue::glue("~/Software/plink --vcf {raw_vcf} --make-bed --out {Setting2_Out_Dir}joined.hg19.nodup.nomismap.training"))
  }

  tagger_tags <- pbmclapply(1:length(unique_regions),function(i){
    cur_region <- unique_regions[i]
    n_tags_to_select <- sum(sapply(Setting2_Region_Info$tag_info,function(x) x$region_id) == cur_region)

    cur_dir <- glue::glue("{Setting2_Out_Dir}Region{i}/")
    system(glue::glue("mkdir -p {cur_dir}"))

    #Existing H3A Tags
    cur_existing_tag <- dplyr::filter(Setting2_TB_Regions[[i]]$wgs_snps,ID %in% existing_tags)$ID
    write(cur_existing_tag,glue::glue("{cur_dir}existing_tags.txt"))

    #Candidate Tags in region
    cur_candidate_tags <- Setting2_Region_Info$region_info[[which(sapply(Setting2_Region_Info$tag_info,function(x) x$region_id) == cur_region)[1]]]$candidate_tag_snps %>%
      dplyr::select(ID,Final_Score,Assay_Type) %>% dplyr::distinct(ID,.keep_all = T)
    cur_candidate_tags$Final_Score[cur_candidate_tags$Final_Score > 1] <- 1

    #All SNPs in the region
    all_snps <- Setting2_TB_Regions[[i]]$wgs_snps %>% dplyr::filter(MAF > 0.05)
    all_snps <- union(union(cur_existing_tag,all_snps$ID),cur_candidate_tags$ID)
    write(all_snps,glue::glue("{cur_dir}all_snps.txt"))

    #Excl Tags = A SNP but not a candidate Tag
    excl_tags <- setdiff(all_snps,union(cur_candidate_tags$ID,cur_existing_tag))
    write(excl_tags,glue::glue("{cur_dir}excl_tags.txt"))


    #Candidate Tags in region
    design_scores <- cur_candidate_tags %>% dplyr::select(ID,Final_Score)
    data.table::fwrite(design_scores,file = glue::glue("{cur_dir}designscores.txt"),sep = ' ',col.names = F,row.names = F,quote = F)

    #Extract Ped File
    cur_chr <- unique(Setting2_TB_Regions[[i]]$wgs_snps$CHROM)
    if(cur_chr == 'X'){
      system(glue::glue("~/Software/plink --bfile {Setting2_Out_Dir}joined.hg19.nodup.nomismap.training --chr X XY --extract {cur_dir}all_snps.txt --recode --out {cur_dir}Region{i}"))
    }else{
      system(glue::glue("~/Software/plink --bfile {Setting2_Out_Dir}joined.hg19.nodup.nomismap.training --chr {cur_chr} --extract {cur_dir}all_snps.txt --recode --out {cur_dir}Region{i}"))
    }
    #Rewrite Info and Ped file
    ped_file <- data.table::fread(glue::glue("{cur_dir}Region{i}.ped"))
    ped_file$V1 <- ped_file$V2
    ped_file$V6 <- 0
    data.table::fwrite(ped_file,col.names = F,row.names = F,quote = F,sep = ' ',file = glue::glue("{cur_dir}Region{i}.ped"))
    map_file <- data.table::fread(glue::glue("{cur_dir}Region{i}.map"))
    data.table::fwrite(map_file %>% dplyr::select(V2,V4),col.names = F,row.names = F,quote = F,sep = ' ',file = glue::glue("{cur_dir}Region{i}.info"))

    #Run Tagger
    system(glue::glue("java -jar -Xmx150000m ~/Software/Haploview.jar -nogui -memory 150000 -pedfile {cur_dir}Region{i}.ped -info {cur_dir}Region{i}.info -designScores {cur_dir}designscores.txt -includeTagsFile {cur_dir}existing_tags.txt -excludeTagsFile {cur_dir}excl_tags.txt -pairwiseTagging -maxNumTags {n_tags_to_select + length(cur_existing_tag)} -minGeno 0.00001 -hwcutoff 0 -missingCutoff 1 -minMAF 0.0 -mintagdistance 0 -tagrsqcutoff 0.8 -out {cur_dir}Region{i} > {cur_dir}Region{i}.tagger.log"))
    selected_tags <- setdiff(data.table::fread(glue::glue("{cur_dir}Region{i}.TESTS"),header = F)$V1,cur_existing_tag)
    if(length(selected_tags) != n_tags_to_select){
      #If tagger cannot select enough tags (non-variable tags tagging itself), pick random tags to avoid bias.
      return(c(selected_tags,sample(setdiff(cur_candidate_tags$ID,union(cur_existing_tag,selected_tags)),size = n_tags_to_select - length(selected_tags))))
      write(paste0('Tagger Failed:',i,'(',length(selected_tags),'/',n_tags_to_select,')'),file = '../results/Tagger/Setting2/tagger_tags.log',append = T)
    }else{
      return(selected_tags)
    }

  },mc.cores = 30)
  saveRDS(tagger_tags,file = '../results/Tagger/Setting2/tagger_tags.rds')
  
  #Extract tags to add based on tagger output
  tagger_snps <- unlist(readRDS('../results/Tagger/Setting2/tagger_tags.rds'))
  extra_snps <- data.table::fread('../results/Tag_SNP_Selection/extra_snps.csv',header = T)$Locus_Name
  mt_snps <- data.table::fread('../results/WGS_Addtl_Tags/mt.tags.setting2.txt',header = F)$V1

  addtl_tags_setting2 <- data.table::fread('../results/WGS_Addtl_Tags/h3a.with.addlt.tags.setting2.txt',header = F)$V1
  tagger_tags_setting1 <- data.table::fread('../results/Tagger/Addtl_Tags/h3a.with.addlt.tags.setting1.tagger.txt',header = F)$V1
  tagger_addtl_tags_setting1 <- setdiff(tagger_tags_setting1,
                                        data.table::fread('../data/WGS_Host_Data/joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz.pos',header = F)$V3)
  
  #Write out add-on tags
  out_prefix <- 'addlt.tags.setting2.tagger'
  write(c(tagger_tags_setting1,intersect(c(mt_snps,extra_snps),addtl_tags_setting2),tagger_snps),file = paste0('../results/Tagger/Addtl_Tags/h3a.with.',out_prefix,'.txt'))
  write(c(intersect(c(mt_snps,extra_snps),addtl_tags_setting2),tagger_snps,tagger_addtl_tags_setting1),file = paste0('../results/Tagger/Addtl_Tags/',out_prefix,'.txt'))
  
  WriteAddOnFile(lift_over_vcf = '../data/WGS_Host_Data/joined.hg19.nodup.nomismap.vcf.gz',phased_vcf = '../results/testing_vcfs/',
                 out_prefix = out_prefix,out_path = '../results/Tagger/Addtl_Tags/')
  
}


# RunTaggerSetting2Full <- function(){
#   WriteAddOnFile <- function(lift_over_vcf,phased_vcf,out_prefix,out_path){
#     addtl_tags_dir <- out_path
#     system(paste0('mkdir -p ',addtl_tags_dir))
#     
#     #Write unphased vcf for AFGR Imputation
#     system(paste0('~/Software/plink2 --vcf ',gsub(lift_over_vcf,pattern = '.gz',replacement = ''),' --extract ',out_path,'h3a.with.',out_prefix,'.txt',' --export vcf bgz --out ',
#                   paste0(addtl_tags_dir,out_prefix)))
#     system(paste0('~/Software/bcftools index -t --threads 5 ',addtl_tags_dir,out_prefix,'.vcf.gz'))
#     testing_samples <- readRDS(paste0('../data/WGS_Host_Data/testing_set/testing_split.rds'))
#     system(paste0('~/Software/bcftools view -O v -s \"',
#                   paste(testing_samples,collapse = ','),'\" ',paste0(addtl_tags_dir,out_prefix,'.vcf.gz'),' > ',
#                   paste0(addtl_tags_dir,out_prefix,'.testing.vcf')))
#     system(paste0('bgzip -c ',paste0(addtl_tags_dir,out_prefix,'.testing.vcf'),'>',paste0(addtl_tags_dir,out_prefix,'.testing.vcf.gz')))
#     system(paste0('~/Software/bcftools index -t --threads 5 ',paste0(addtl_tags_dir,out_prefix,'.testing.vcf.gz')))
#     #Write phased VCF for WGS Imputation
#     chr <- c(seq(1,22,1),'X','PAR')
# 
#     addtl_tags_dir <- paste0(out_path,'/Setting2/')
#     system(paste0('mkdir -p ',addtl_tags_dir))
#     #Write files with h3a and addtl tags for each chr (including X)
#     for(cur_chr in chr){
#       system(paste0('gunzip -c ',phased_vcf,'chr',cur_chr,'.phased.testing.vcf.gz > ',phased_vcf,'chr',cur_chr,'.phased.testing.vcf'))
#       system(paste0('~/Software/plink2 --vcf ',phased_vcf,'chr',cur_chr,'.phased.testing.vcf',' --extract ',out_path,
#                     '/h3a.with.',out_prefix,'.txt',' --export vcf bgz --out ',
#                     paste0(addtl_tags_dir,'chr',cur_chr,'.',out_prefix,'.testing')))
#       system(paste0('~/Software/bcftools index -t --threads 5 ',addtl_tags_dir,'chr',cur_chr,'.',out_prefix,'.testing.vcf.gz'))
#     }
#     
#   }
#   
#   out_prefix <- 'addlt.tags.setting2.tagger'
#   load('../results/Tag_SNP_Selection/setting2.rda')
#   remove(list = setdiff(ls(),'region_info'))
#   tags_incl <- data.table::fread('../results/WGS_Addtl_Tags/h3a.with.addlt.tags.setting2.txt',header = F)
#   Setting2_Region_Info <- readRDS('../results/WGS_Addtl_Tags/addlt.tags.setting2.rds')
#   regions_incl <- which(Setting2_Region_Info$tags_to_add %in% tags_incl$V1)
# 
#   Setting2_Region_Info$tags_to_add <- Setting2_Region_Info$tags_to_add[regions_incl]
#   Setting2_Region_Info$tag_info <- Setting2_Region_Info$max_tag_info[regions_incl]
#   Setting2_Region_Info$region_info <- region_info[regions_incl]
# 
#   TB_Regions_Setting2 <- readRDS('../results/WGS_Addtl_Tags/Setting2/TB_gene_regions_Setting2.rds')
#   TB_Regions_Setting2 <- TB_Regions_Setting2$snps_in_gene_regions_parsed[sapply(Setting2_Region_Info$tag_info,function(x) x$region_id)]
#   Setting2_Region_Info$wgs_snps <- lapply(TB_Regions_Setting2,function(x) x$wgs_snps)
# 
#   existing_tags <- list()
#   all_snps <- list()
#   candidate_tags <- list()
#   design_scores <- list()
# 
#   for(i in 1:length(Setting2_Region_Info$region_info)){
#     print(i)
#     #Existing H3A Tags
#     existing_tags[[i]] <- Setting2_Region_Info$region_info[[i]]$existing_tag_snps$ID
# 
#     #Candidate Tags in region
#     cur_candidate_tags <- Setting2_Region_Info$region_info[[i]]$candidate_tag_snps %>% dplyr::select(ID,Final_Score,Assay_Type) %>% dplyr::distinct(ID,.keep_all = T)
#     cur_candidate_tags$Final_Score[cur_candidate_tags$Final_Score > 1] <- 1
#     candidate_tags[[i]] <- cur_candidate_tags$ID
# 
#     #All SNPs in the region
#     cur_all_snps <- Setting2_Region_Info$wgs_snps[[i]] %>% dplyr::filter(MAF > 0.05)
#     all_snps[[i]] <- union(union(existing_tags[[i]],cur_all_snps$ID),cur_candidate_tags$ID)
# 
#     #Candidate Tags in region
#     design_scores[[i]] <- cur_candidate_tags %>% dplyr::select(ID,Final_Score)
# 
#   }
#   out_dir <- "../results/Tagger/Setting2/"
#   system(glue::glue("mkdir -p {out_dir}"))
#   write(unique(unlist(all_snps)),glue::glue('{out_dir}all_snps.txt'))
#   write(unlist(existing_tags),glue::glue('{out_dir}existing_tags.txt'))
# 
#   #Excl Tags = A SNP but not a candidate/existing Tag
#   excl_tags <- setdiff(unlist(all_snps),union(unlist(candidate_tags),unlist(existing_tags)))
#   write(unique(excl_tags),glue::glue("{out_dir}excl_tags.txt"))
#   data.table::fwrite(do.call(rbind,design_scores) %>%  dplyr::distinct(ID,.keep_all = T) %>% dplyr::filter(ID %in% unlist(candidate_tags)),file = glue::glue("{out_dir}designscores.txt"),col.names = F,row.names = F,quote = F,sep = ' ')
# 
#   raw_vcf <- '../data/WGS_Host_Data/training_set/joined.hg19.nodup.nomismap.training.vcf.gz'
#   system(glue::glue("~/Software/plink --vcf {raw_vcf} --extract {out_dir}all_snps.txt --recode --out {out_dir}Setting2.all.snps"))
# 
#   #Rewrite Info and Ped file
#   ped_file <- data.table::fread(glue::glue("{out_dir}Setting2.all.snps.ped"))
#   ped_file$V1 <- ped_file$V2
#   ped_file$V6 <- 0
#   data.table::fwrite(ped_file,col.names = F,row.names = F,quote = F,sep = ' ',file = glue::glue("{out_dir}Setting2.all.snps.ped"))
#   map_file <- data.table::fread(glue::glue("{out_dir}Setting2.all.snps.map"))
#   chrom_cumsum <- cumsum(as.numeric(sapply(unique(map_file$V1),function(x)max(map_file$V4[map_file$V1 == x]))))
#   for(i in 2:length(unique(map_file$V1))){
#     cur_chr <- unique(map_file$V1)[i]
#     map_file$V4[map_file$V1 == cur_chr] <- map_file$V4[map_file$V1 == cur_chr] + chrom_cumsum[i-1]
#   }
#   data.table::fwrite(map_file %>% dplyr::select(V2,V4),col.names = F,row.names = F,quote = F,sep = ' ',file = glue::glue("{out_dir}Setting2.all.snps.info"))
#   #Run Tagger
#   system(glue::glue("java -jar -Xmx150000m ~/Software/Haploview.jar -nogui -memory 150000 -pedfile {out_dir}Setting2.all.snps.ped -info {out_dir}Setting2.all.snps.info -designScores {out_dir}designscores.txt -includeTagsFile {out_dir}existing_tags.txt -excludeTagsFile {out_dir}excl_tags.txt -pairwiseTagging -dontaddtags -minGeno 0.00001 -hwcutoff 0 -missingCutoff 1 -minMAF 0.0 -mintagdistance 0 -tagrsqcutoff 0.8 -out {out_dir}Setting2.noaddtags > {out_dir}Setting2.noaddtags.tagger.log"))
#   system(glue::glue("java -jar -Xmx150000m ~/Software/Haploview.jar -nogui -memory 100000 -pedfile {out_dir}Setting2.all.snps.ped -info {out_dir}Setting2.all.snps.info -designScores {out_dir}designscores.txt -includeTagsFile {out_dir}existing_tags.txt -excludeTagsFile {out_dir}excl_tags.txt -pairwiseTagging -minGeno 0.00001 -hwcutoff 0 -missingCutoff 1 -minMAF 0.0 -mintagdistance 0 -tagrsqcutoff 0.8 -maxNumTags {36685+2504} -out {out_dir}Setting2 > {out_dir}Setting2.tagger.log"))
#   system(glue::glue("java -jar -Xmx150000m ~/Software/Haploview.jar -nogui -memory 100000 -pedfile {out_dir}Setting2.all.snps.ped -info {out_dir}Setting2.all.snps.info -designScores {out_dir}designscores.txt -includeTagsFile {out_dir}existing_tags.txt -excludeTagsFile {out_dir}excl_tags.txt -pairwiseTagging -minGeno 0.00001 -hwcutoff 0 -missingCutoff 1 -minMAF 0.0 -mintagdistance 0 -tagrsqcutoff 0.8 -maxNumTags {36685+5000} -out {out_dir}Setting2.full > {out_dir}Setting2.full.tagger.log"))
#   
#   #Select Number of SNPs matching those added under setting 2
#   extra_snps <- data.table::fread('../results/Tag_SNP_Selection/extra_snps.csv',header = T)
#   mt_snps <- data.table::fread('../results/WGS_Addtl_Tags/mt.tags.setting2.txt',header = F)
#   
#   addtl_tags_setting1 <- data.table::fread('../results/WGS_Addtl_Tags/h3a.with.addlt.tags.setting1.txt',header = F)
#   addtl_tags_setting2 <- data.table::fread('../results/WGS_Addtl_Tags/h3a.with.addlt.tags.setting2.txt',header = F)
#   N_SNPs <- length(setdiff(addtl_tags_setting2$V1,c(extra_snps$Locus_Name,mt_snps$V1,addtl_tags_setting1$V1)))
#   
#   #Extract tags to add based on tagger output
#   tagger_snps <- data.table::fread('../results/Tagger/Setting2/Setting2.full.TESTS',header = F)
#   setting_1_snps <- nrow(data.table::fread(glue::glue('{out_dir}existing_tags.txt'),header=F))
#   tagger_snps_addtl <- tagger_snps$V1[(setting_1_snps+1):nrow(tagger_snps)]
#   
#   #Remove existing tagger Setting 1 tags
#   tagger_addtl_tags_setting1 <- data.table::fread('../results/Tagger/Addtl_Tags/h3a.with.addlt.tags.setting1.tagger.txt',header = F)
#   tagger_snps_addtl <- tagger_snps_addtl[sort(which(is.na(match(tagger_snps_addtl,tagger_addtl_tags_setting1$V1))))]
#   tagger_snps_addtl <- tagger_snps_addtl[1:N_SNPs]
#   
#   #Write out add-on tags
#   write(c(tagger_addtl_tags_setting1$V1,intersect(c(mt_snps$V1,extra_snps$Locus_Name),addtl_tags_setting2$V1),tagger_snps_addtl),file = paste0('../results/Tagger/Addtl_Tags/h3a.with.',out_prefix,'.txt'))
#   write(c(intersect(c(mt_snps$V1,extra_snps$Locus_Name),addtl_tags_setting2$V1),tagger_snps_addtl),file = paste0('../results/Tagger/Addtl_Tags/',out_prefix,'.txt'))
#   
#   WriteAddOnFile(lift_over_vcf = '../data/WGS_Host_Data/joined.hg19.nodup.nomismap.vcf.gz',phased_vcf = '../results/testing_vcfs/',
#                  out_prefix = out_prefix,out_path = '../results/Tagger/Addtl_Tags/')
#   
# }
args <- commandArgs(trailingOnly = T)
if(args[[1]] == 'setting1'){
  RunTaggerSetting1()
}else if(args[[1]] == 'setting2'){
  RunTaggerSetting2()
}else{
  stop('Invalid Setting')
}
