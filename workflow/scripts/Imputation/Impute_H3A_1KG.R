library(dplyr)
library(parallel)

# #Write out Imputation Files
# system('mkdir -p ../results/1000_Genomes/')
# TB_DAR_path <- '../data/WGS_Host_Data/joined.hg19.nodup.nomismap'
# TB_DAR_pos <- data.table::fread(paste0(TB_DAR_path,'.bim')) %>% dplyr::filter(!V1 %in% c('Y'))
# 
# TB_DAR_vcf_h3a <- '../data/WGS_Host_Data/joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz'
# TB_DAR_vcf_pos_h3a <- data.table::fread('../data/WGS_Host_Data/joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz.pos') %>% dplyr::filter(!V1 %in% c('Y'))
# 
# KG_vcf <- '../data/1000_Genomes/joined.1000genomes.withchrX.vcf.gz'
# KG_vcf_pos <- data.table::fread('../data/1000_Genomes/joined.1000genomes.withchrX.vcf.gz.pos')
# 
# consensus_h3a <- dplyr::inner_join(TB_DAR_vcf_pos_h3a,KG_vcf_pos,by=c('V1'='V1','V2'='V2','V4'='V4','V5'='V5'))
# snps_to_keep_TB_DAR <- consensus_h3a$V3.x
# write(snps_to_keep_TB_DAR,'../results/1000_Genomes/consensus_snps_TBDAR.txt')
# system(glue::glue("~/Software/plink2 --vcf {TB_DAR_vcf_h3a} --extract ../results/1000_Genomes/consensus_snps_TBDAR.txt --export vcf bgz --out ../results/1000_Genomes/TB_DAR"))
# system(glue::glue("~/Software/plink2 --bfile {TB_DAR_path} --threads 22 --freq --chr 1-22 --out ../results/1000_Genomes/TBDAR"))
# 
# snps_to_keep_KG <- consensus_h3a$V3.y
# write(snps_to_keep_KG,'../results/1000_Genomes/consensus_snps_1KG.txt')
# 
# KG_id_mapping <-data.table::fread('../data/1000_Genomes/20130606_g1k.ped') %>% dplyr::filter(Population %in% c('MSL','GWD','YRI','ESN','LWK')) %>%
#   dplyr::filter(!Relationship %in% c('child','Child2','paternal brother','paternal father','paternal grandmother'))
# mclapply(1:length(unique(KG_id_mapping$Population)),function(i){
#   cur_pop <- unique(KG_id_mapping$Population)[i]
#   cur_ID <- KG_id_mapping %>% dplyr::filter(Population == cur_pop)
#   write(cur_ID$`Individual ID`,glue::glue("../results/1000_Genomes/{cur_pop}.IDs.txt"))
#   system(glue::glue("~/Software/plink2 --vcf {KG_vcf} --threads 22 --extract ../results/1000_Genomes/consensus_snps_1KG.txt --keep ../results/1000_Genomes/{cur_pop}.IDs.txt --export vcf bgz --out ../results/1000_Genomes/1KG_{cur_pop}"))
#   system(glue::glue("~/Software/plink2 --vcf {KG_vcf} --threads 22 --keep ../results/1000_Genomes/{cur_pop}.IDs.txt --freq --out ../results/1000_Genomes/1KG_{cur_pop}"))
# },mc.cores = length(unique(KG_id_mapping$Population)))


#Compare Imputation
datasets <- c('ESN','GWD','LWK','MSL','TBDAR','YRI')

info_file <- lapply(datasets,function(x){
  print(glue::glue('INFO:{x}'))
  system(glue::glue("~/Software/bcftools query -f '%CHROM %POS %ID %REF %ALT %INFO/RefPanelAF %INFO/INFO\n' ../results/1000_Genomes/AFGR_Imputed/{x}H3A.vcfs/{x}H3A.vcf.gz > ../results/1000_Genomes/AFGR_Imputed/{x}H3A.vcfs/{x}H3A.info.txt"))
  df = data.table::fread(glue::glue("awk '{{ if (($6 < 0.99) && ($6 > 0.01)) {{ print }} }}' ../results/1000_Genomes/AFGR_Imputed/{x}H3A.vcfs/{x}H3A.info.txt"),header = F,select = c(6,7))
  colnames(df) <- c('AFGR_AF','INFO')
  return(df)
})
names(info_file) <- datasets
saveRDS(info_file,'../results/1000_Genomes/info_file.rds')

info_file_X <- lapply(datasets,function(x){
  print(glue::glue('INFO_X:{x}'))
  system(glue::glue("~/Software/bcftools query -f '%CHROM %POS %ID %REF %ALT %INFO/RefPanelAF %INFO/INFO\n' ../results/1000_Genomes/AFGR_Imputed/{x}H3A.vcfs/X.vcf.gz > ../results/1000_Genomes/AFGR_Imputed/{x}H3A.vcfs/X.info.txt"))
  df = data.table::fread(glue::glue("awk '{{ if (($6 < 0.99) && ($6 > 0.01)) {{ print }} }}' ../results/1000_Genomes/AFGR_Imputed/{x}H3A.vcfs/X.info.txt"),header = F,select = c(6,7))
  colnames(df) <- c('AFGR_AF','INFO')
  return(df)
})
names(info_file_X) <- datasets
saveRDS(info_file_X,'../results/1000_Genomes/info_file_X.rds')


# maf_file <- lapply(datasets,function(x){
#   print(glue::glue('MAF:{x}'))
#   if(x != 'TBDAR'){
#     maf_file <- data.table::fread(glue::glue("awk '{{ if (($5 < 0.95) && ($5 > 0.05)) {{ print }} }}' ../results/1000_Genomes/1KG_{x}.afreq"),select = c(2,5))
#     colnames(maf_file) <- c('ID','AF')
#     maf_file <- maf_file %>% dplyr::left_join(KG_vcf_pos %>% dplyr::select(CHROM=V1,POS=V2,ID=V3,REF=V4,ALT=V5),by=c('ID'='ID'))
#   }else{
#     maf_file <- data.table::fread(glue::glue("awk '{{ if (($5 < 0.95) && ($5 > 0.05)) {{ print }} }}' ../results/1000_Genomes/TBDAR.afreq"),select = c(2,5))
#     colnames(maf_file) <- c('ID','AF')
#     maf_file <- maf_file %>% dplyr::left_join(TB_DAR_pos %>% dplyr::select(CHROM=V1,POS=V4,ID=V2,REF=V6,ALT=V5),by=c('ID'='ID'))
#   }
#   maf_file$MAF <- pmin(maf_file$AF,1-maf_file$AF)
#   maf_file <- maf_file %>% dplyr::select(-AF,-ID)
#   return(maf_file)
# })
# names(maf_file) <- datasets
# saveRDS(maf_file,'../results/1000_Genomes/maf_df.rds')
# 
# info_file <- lapply(datasets,function(x){
#   print(glue::glue('INFO:{x}'))
# 
#   df = data.table::fread(glue::glue("../results/1000_Genomes/AFGR_Imputed/{x}H3A.vcfs/{x}H3A.info.txt"),header = F,select = c(1,2,4,5,6,7))
#   colnames(df) <- c('CHROM','POS','REF','ALT','AFGR_AF','INFO')
#   return(dplyr::left_join(maf_file[[x]],df,by=c('CHROM','POS','REF','ALT')))
# })
# names(info_file) <- datasets
# saveRDS(info_file,'../results/1000_Genomes/info_df.rds')
