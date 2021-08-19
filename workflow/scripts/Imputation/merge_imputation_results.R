#Merge imputation results from AFGR and TB-DAR, based on better imputation accuracy per SNP.

library(glue)
library(dplyr)
library(data.table)

MergeImputationResults <- function(afgr_imputed_path,afgr_prefix,wgs_imputed_path,wgs_unimputed_info,out_dir,n_cores){
  system(glue('mkdir -p {out_dir}'))
  system(glue('mkdir -p {out_dir}/tmp/'))
  
  #Read AFGR Imputed info 
  afgr_info <- data.table::fread(paste0(afgr_imputed_path,afgr_prefix,'.info.txt'),header = F)
  colnames(afgr_info) <- c('CHROM','POS','ID','REF','ALT','AF','INFO')
  afgr_nodup <- data.table::fread(paste0(afgr_imputed_path,afgr_prefix,'.nodup.vcf.gz.pos'),header = F)
  colnames(afgr_nodup) <- c('CHROM','POS','ID','REF','ALT')
  afgr_info <- dplyr::inner_join(afgr_nodup,afgr_info,by=c('CHROM'='CHROM','POS'='POS','REF'='REF','ALT'='ALT')) %>% dplyr::select(CHROM,POS,ID=ID.x,REF,ALT,INFO,AF.AFGR=AF)
  remove(afgr_nodup)
  gc();
  
  #Extract WGS Imputed Info
  system(paste0("bcftools query -f '%CHROM %POS %ID %REF %ALT %INFO/R2 %INFO/AF\n' ",wgs_imputed_path,'WGS.imputed.vcf.gz > ',wgs_imputed_path,'WGS.imputed.info.txt'))
  #system(paste0("bcftools query -f '%CHROM %POS %ID %REF %ALT %INFO/R2 %INFO/AF\n' ",wgs_imputed_path,'WGS.imputed.chrX.dose.vcf.gz > ',wgs_imputed_path,'WGS.imputed.chrX.info.txt'))
  #system(paste0("bcftools query -f '%CHROM %POS %ID %REF %ALT %INFO/R2 %INFO/AF\n' ",wgs_imputed_path,'WGS.imputed.chrPAR.dose.vcf.gz > ',wgs_imputed_path,'WGS.imputed.chrXPAR.info.txt'))
  
  # if(!'WGS.imputed.chrPAR.info.txt' %in% dir(wgs_imputed_path)){
  #   wgs_info <- rbind(rbind(data.table::fread(paste0(wgs_imputed_path,'WGS.imputed.info.txt'),header = F),
  #                           data.table::fread(paste0(wgs_imputed_path,'WGS.imputed.chrX.info.txt'),header = F)))
  # }else{
  #   wgs_info <- rbind(rbind(data.table::fread(paste0(wgs_imputed_path,'WGS.imputed.info.txt'),header = F),
  #                           data.table::fread(paste0(wgs_imputed_path,'WGS.imputed.chrX.info.txt'),header = F)),
  #                     data.table::fread(paste0(wgs_imputed_path,'WGS.imputed.chrPAR.info.txt'),header = F))
  #   
  # }
  wgs_info <- data.table::fread(paste0(wgs_imputed_path,'WGS.imputed.info.txt'),header = F)
  colnames(wgs_info) <- c('CHROM','POS','ID','REF','ALT','R2','AF')
  wgs_unimputed_info <- data.table::fread(wgs_unimputed_info)
  wgs_unimputed_info <- dplyr::select(wgs_unimputed_info,'CHROM'=V1,'POS'=V2,'ID'=V3,'REF'=V4,'ALT'=V5)
  wgs_info <- dplyr::inner_join(wgs_unimputed_info,wgs_info,by=c('CHROM'='CHROM','POS'='POS','REF'='REF','ALT'='ALT')) %>% dplyr::select(CHROM,POS,ID.WGS=ID.x,ID.WGS.imp=ID.y,REF,ALT,R2,AF.WGS=AF)
  remove(wgs_unimputed_info)
  gc();
  
  #Merge WGS imputed and AFGR imputed
  wgs_afgr_imp_merged <- dplyr::full_join(afgr_info,wgs_info,by = c('CHROM'='CHROM','POS'='POS','REF'='REF','ALT'='ALT')) %>% dplyr::select(CHROM,POS,ID.AFGR=ID,REF,ALT,INFO,ID.WGS,ID.WGS.imp,R2,AF.WGS,AF.AFGR)
  afgr_snps <- which(wgs_afgr_imp_merged$INFO >= wgs_afgr_imp_merged$R2)
  afgr_snps <- union(afgr_snps,which(is.na(wgs_afgr_imp_merged$R2)))
  
  wgs_snps <- which(wgs_afgr_imp_merged$INFO < wgs_afgr_imp_merged$R2)
  wgs_snps <- union(wgs_snps,which(is.na(wgs_afgr_imp_merged$INFO)))
  
  write(wgs_afgr_imp_merged$ID.AFGR[afgr_snps],file=paste0(out_dir,'/tmp/afgr_imp_rsid.txt'))
  write(wgs_afgr_imp_merged$ID.WGS.imp[wgs_snps],file=paste0(out_dir,'/tmp/wgs_imp_rsid.txt'))
  write(wgs_afgr_imp_merged$ID.WGS[wgs_snps],file=paste0(out_dir,'/tmp/wgs_rsid.txt'))
  
  system(glue('plink2 --vcf {paste0(afgr_imputed_path,afgr_prefix,".nodup.vcf")} --extract {paste0(out_dir,"/tmp/afgr_imp_rsid.txt")} --export vcf bgz --out {paste0(out_dir,"/tmp/afgr_info_selected")}'))
  system(glue('plink2 --vcf {paste0(wgs_imputed_path,"WGS.imputed.vcf.gz")} --extract {paste0(out_dir,"/tmp/wgs_imp_rsid.txt")} --export vcf bgz --out {paste0(out_dir,"/tmp/wgs_imp_info_selected")}'))
  system(glue('bcftools index -t --threads {n_cores} {paste0(out_dir,"/tmp/wgs_imp_info_selected.vcf.gz")}'))
  system(glue('bcftools index -t --threads {n_cores} {paste0(out_dir,"/tmp/afgr_info_selected.vcf.gz")}'))
  #Keep consensus samples
  system(glue('bcftools query -l {paste0(out_dir,"/tmp/afgr_info_selected.vcf.gz")} > {paste0(out_dir,"/tmp/afgr_info_selected.vcf.gz")}.samples'))
  system(glue('bcftools query -l {paste0(out_dir,"/tmp/wgs_imp_info_selected.vcf.gz")} > {paste0(out_dir,"/tmp/wgs_imp_info_selected.vcf.gz")}.samples'))
  consensus_samples <- intersect(data.table::fread(paste0(out_dir,"/tmp/afgr_info_selected.vcf.gz.samples"),header = F)$V1,
                             data.table::fread(paste0(out_dir,"/tmp/wgs_imp_info_selected.vcf.gz.samples"),header = F)$V1)
  write(consensus_samples,paste0(out_dir,"/tmp/consensus_samples"))
  system(glue('bcftools view -S {paste0(out_dir,"/tmp/consensus_samples")} -O z -o {paste0(out_dir,"/tmp/afgr_info_selected.consensus.vcf.gz")} {paste0(out_dir,"/tmp/afgr_info_selected.vcf.gz")}'))
  system(glue('bcftools view -S {paste0(out_dir,"/tmp/consensus_samples")} -O z -o {paste0(out_dir,"/tmp/wgs_imp_info_selected.consensus.vcf.gz")} {paste0(out_dir,"/tmp/wgs_imp_info_selected.vcf.gz")}'))
  system(glue('bcftools index -t --threads {n_cores} {paste0(out_dir,"/tmp/wgs_imp_info_selected.consensus.vcf.gz")}'))
  system(glue('bcftools index -t --threads {n_cores} {paste0(out_dir,"/tmp/afgr_info_selected.consensus.vcf.gz")}'))
  
  #Merge Two reference panels
  system(glue('bcftools concat -a {paste0(out_dir,"/tmp/wgs_imp_info_selected.consensus.vcf.gz")} {paste0(out_dir,"/tmp/afgr_info_selected.consensus.vcf.gz")} -O z -o {paste0(out_dir,"/tmp/wgs.afgr.merged.vcf.gz")}'))
  system(glue('bcftools index -t --threads {n_cores} {paste0(out_dir,"/tmp/wgs.afgr.merged.vcf.gz")}'))
  system(glue('bcftools sort -T {out_dir} -O z -o {paste0(out_dir,"/tmp/wgs.afgr.merged.sorted.vcf.gz")} {paste0(out_dir,"/tmp/wgs.afgr.merged.vcf.gz")}'))
  system(glue('bcftools index -t --threads {n_cores} {paste0(out_dir,"/tmp/wgs.afgr.merged.sorted.vcf.gz")}'))
  
  #Write out merged file 
  wgs_afgr_imp_merged_afgr <- wgs_afgr_imp_merged[afgr_snps,] %>% dplyr::select(CHROM,POS,ID=ID.AFGR,REF,ALT,INFO,AF=AF.WGS,AF.AFGR)
  wgs_afgr_imp_merged_afgr$AF[is.na(wgs_afgr_imp_merged_afgr$AF)] <- wgs_afgr_imp_merged_afgr$AF.AFGR[is.na(wgs_afgr_imp_merged_afgr$AF)]
  wgs_afgr_imp_merged_afgr <- dplyr::select(wgs_afgr_imp_merged_afgr,-'AF.AFGR')
  wgs_afgr_imp_merged_wgs <- wgs_afgr_imp_merged[wgs_snps,] %>% dplyr::select(CHROM,POS,ID=ID.WGS,REF,ALT,INFO=R2,AF=AF.WGS)
  wgs_afgr_imp_merged_afgr$Source = 'AFGR'
  wgs_afgr_imp_merged_wgs$Source = 'WGS'
  wgs_afgr_imp_merged <- rbind(wgs_afgr_imp_merged_afgr,wgs_afgr_imp_merged_wgs)                                                     
  data.table::fwrite(wgs_afgr_imp_merged,file = paste0(out_dir,'WGS_AFGR_Merged_Info.txt'),sep = ' ')
  
  unique_chr <- unique(wgs_afgr_imp_merged$CHROM)
  for(chr in unique_chr){
    data.table::fwrite(dplyr::filter(wgs_afgr_imp_merged,CHROM==chr) %>% dplyr::arrange(POS),file = paste0(out_dir,'WGS_AFGR_Merged_Info_chr',chr,'.txt'),sep = ' ')
  }
}

args <- commandArgs(trailingOnly = T)
MergeImputationResults(args[[1]],args[[2]],args[[3]],args[[4]],args[[5]],args[[6]])

# afgr_imputed_path <- '/home/mxu/G2G_TB/AFGRTesting.vcfs/'
# afgr_prefix <- 'AFGR.imputed'
# wgs_imputed_path <- '/home/mxu/G2G_TB/WGS_Imputed/'
# out_dir <- '~/G2G_TB/AFGRAndWGSMergedImpFirstPass/'
# MergeImputationResults(afgr_imputed_path,afgr_prefix,wgs_imputed_path,out_dir)
# 
# afgr_imputed_path <- '/home/mxu/G2G_TB/WGS_Host_Data/WGS_Addtl_Tags/AFGRSecondPass0p9.vcfs/'
# afgr_prefix <- 'AFGRSecondPass.imputed'
# wgs_imputed_path <- '/home/mxu/G2G_TB/WGS_Host_Data/WGS_Addtl_Tags/WGSSecondPass0p9.vcfs/'
# out_dir <- '~/G2G_TB/WGS_Host_Data/WGS_Addtl_Tags/AFGRAndWGSMergedImpSecondPass/'
# MergeImputationResults(afgr_imputed_path,afgr_prefix,wgs_imputed_path,out_dir)
# 
# afgr_imputed_path <- '/home/mxu/G2G_TB/WGS_Host_Data/WGS_Addtl_Tags/AFGRThirdPassSum.vcfs//'
# afgr_prefix <- 'AFGRThirdPass.imputed'
# wgs_imputed_path <- '/home/mxu/G2G_TB/WGS_Host_Data/WGS_Addtl_Tags/WGSThirdPassSum.vcfs/'
# out_dir <- '~/G2G_TB/WGS_Host_Data/WGS_Addtl_Tags/AFGRAndWGSMergedImpThirdPassSum/'
# MergeImputationResults(afgr_imputed_path,afgr_prefix,wgs_imputed_path,out_dir)
# 
# afgr_imputed_path <- '/home/mxu/G2G_TB/WGS_Host_Data/WGS_Addtl_Tags/AFGRThirdPassAvg.vcfs/'
# afgr_prefix <- 'AFGRThirdPass.imputed'
# wgs_imputed_path <- '/home/mxu/G2G_TB/WGS_Host_Data/WGS_Addtl_Tags/WGSThirdPassAvg.vcfs/'
# out_dir <- '~/G2G_TB/WGS_Host_Data/WGS_Addtl_Tags/AFGRAndWGSMergedImpThirdPassAvg/'
# MergeImputationResults(afgr_imputed_path,afgr_prefix,wgs_imputed_path,out_dir)
