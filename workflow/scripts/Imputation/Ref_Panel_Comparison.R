library(filematrix)
library(dplyr)

#Convert raw file into dosage and filematrix
WriteFileMatrix <- function(file_path){
  system(glue::glue("awk \'{{print $2}}\' {file_path} | rs -c -C -T > {gsub(x=file_path,pattern='raw',replacement='dosage')}"))
  system(glue::glue("cut -f 9- {file_path} | rs -c -C -T >> {gsub(x=file_path,pattern='raw',replacement='dosage')}"))

  fm = fm.create.from.text.file(textfilename = gsub(x=file_path,pattern='.raw',replacement='.dosage'),
                                filenamebase = gsub(x=file_path,pattern='.raw',replacement=''),type = 'integer')
  
}

#Correlation between imputed and groud truth
FindCor <- function(mat1,mat2,out_path,n_cores = 1){
  mat1_colnames <- colnames(mat1)
  mat2_colnames <- colnames(mat2)
  
  shared_variants <- intersect(mat1_colnames,mat2_colnames)
  mat1_index <- match(shared_variants,mat1_colnames)
  mat2_index <- match(shared_variants,mat2_colnames)
  
  #rand_variants <- shared_variants[sample(1:length(shared_variants),size = 1000,replace = F)]
  cor_vect <- unlist(pbmcapply::pbmclapply(1:length(shared_variants),function(i) as.vector(cor(mat1[,mat1_index[i]],mat2[,mat2_index[i]],use = 'complete.obs')^2),mc.cores = n_cores,mc.preschedule = T))
  names(cor_vect) <- shared_variants
  saveRDS(cor_vect,out_path)
}

# system(glue::glue("~/Software/plink2 --bfile ../data/WGS_Host_Data/joined.hg19.nodup.nomismap --maf 0.01 --geno 0.5 --set-all-var-ids @:#:\\$r:\\$a --make-bed --out ../data/WGS_Host_Data/joined.hg19.nodup.nomismap.maf0p01"))
# system(glue::glue("~/Software/plink2 --freq --bfile ../data/WGS_Host_Data/joined.hg19.nodup.nomismap.maf0p01 --out ../data/WGS_Host_Data/joined.hg19.nodup.nomismap.maf0p01"))

# maf_file <- data.table::fread('../data/WGS_Host_Data/joined.hg19.nodup.nomismap.maf0p01.bim',select = c('V1','V2','V4','V5','V6')) %>%
#   dplyr::left_join(data.table::fread('../data/WGS_Host_Data/joined.hg19.nodup.nomismap.maf0p01.afreq',select = c('ID','ALT_FREQS'),header = T),by=c('V2'='ID')) %>% dplyr::select(CHROM=V1,ID=V2,POS = V4,REF=V6,ALT=V5,ALT_FREQS) %>%
#   dplyr::filter(CHROM != 'X')
# maf_file$MAF <- sapply(maf_file$ALT_FREQS,function(x) ifelse(x > 0.5,1-x,x))
# maf_file$CHROM <- as.numeric(maf_file$CHROM)
# system(glue::glue("~/Software/plink2 --bfile ../data/WGS_Host_Data/joined.hg19.nodup.nomismap.maf0p01 --export A --out ../data/WGS_Host_Data/joined.hg19.nodup.nomismap.maf0p01"))

# CAAPA_Info <- data.table::fread('../data/CAAPA.vcfs/CAAPA.info',header = F) %>%
#   select(CHROM = V1,POS=V2,ID=V3,REF=V4,ALT=V5,Rsq=V7,Imputed=V6) %>%
# dplyr::inner_join(maf_file,by=c('CHROM'='CHROM','POS'='POS','REF'='REF','ALT'='ALT'))
# data.table::fwrite(cbind(CAAPA_Info$CHROM,CAAPA_Info$POS,CAAPA_Info$POS),'../data/CAAPA.vcfs/CAAPA.maf0p01.SNPs',col.names = F,row.names = F,sep = ' ')
# system(glue::glue("~/Software/plink2 --vcf ../data/CAAPA.vcfs/CAAPA.vcf.gz --extract range ../data/CAAPA.vcfs/CAAPA.maf0p01.SNPs --set-all-var-ids @:#:\\$r:\\$a --make-bed --out ../data/CAAPA.vcfs/CAAPA.maf0p01"))
# system(glue::glue("~/Software/plink2 --bfile ../data/CAAPA.vcfs/CAAPA.maf0p01 --export A --out ../data/CAAPA.vcfs/CAAPA.maf0p01"))

# AFGR_Info <- data.table::fread('../data/AFGR.vcfs/AFGR.info',header = T,select = c('CHROM','POS','ID','REF','ALT','INFO')) %>%
#   dplyr::inner_join(maf_file,by=c('CHROM'='CHROM','POS'='POS','REF'='REF','ALT'='ALT'))
# data.table::fwrite(cbind(AFGR_Info$CHROM,AFGR_Info$POS,AFGR_Info$POS),'../data/AFGR.vcfs/AFGR.maf0p01.SNPs',col.names = F,row.names = F,sep = ' ')
# system(glue::glue("~/Software/plink2 --vcf ../data/AFGR.vcfs/AFGR.vcf.gz --extract range ../data/AFGR.vcfs/AFGR.maf0p01.SNPs --set-all-var-ids @:#:\\$r:\\$a --make-bed --out ../data/AFGR.vcfs/AFGR.maf0p01"))
# system(glue::glue("~/Software/plink2 --bfile ../data/AFGR.vcfs/AFGR.maf0p01 --export A --out ../data/AFGR.vcfs/AFGR.maf0p01"))

# HRC_Info <- data.table::fread('../data/HRC.vcfs/HRC.info',header = T,select = c('CHROM','POS','REF','ALT','INFO')) %>%
#   dplyr::inner_join(maf_file,by=c('CHROM'='CHROM','POS'='POS','REF'='REF','ALT'='ALT'))
# data.table::fwrite(cbind(HRC_Info$CHROM,HRC_Info$POS,HRC_Info$POS),'../data/HRC.vcfs/HRC.maf0p01.SNPs',col.names = F,row.names = F,sep = ' ')
# system(glue::glue("~/Software/plink2 --vcf ../data/HRC.vcfs/HRC.vcf.gz --extract range ../data/HRC.vcfs/HRC.maf0p01.SNPs --set-all-var-ids @:#:\\$r:\\$a --make-bed --out ../data/HRC.vcfs/HRC.maf0p01"))
# system(glue::glue("~/Software/plink2 --bfile ../data/HRC.vcfs/HRC.maf0p01 --export A --out ../data/HRC.vcfs/HRC.maf0p01"))


# TopMed_hg19 <- data.table::fread('../data/TopMed.vcfs/TopMed.hg19.info',header = F)
# TopMed_hg19 <- TopMed_hg19 %>% select(CHROM = V1,POS=V2,ID=V3,REF=V4,ALT=V5,Rsq=V7,Imputed=V6) %>%
#   dplyr::inner_join(maf_file %>% dplyr::mutate(CHROM_TopMed=paste0('chr',CHROM)),by=c('CHROM'='CHROM_TopMed','POS'='POS','REF'='REF','ALT'='ALT'))
# data.table::fwrite(cbind(TopMed_hg19$CHROM,TopMed_hg19$POS,TopMed_hg19$POS),'../data/TopMed.vcfs/TopMed.maf0p01.SNPs',col.names = F,row.names = F,sep = ' ')
# system(glue::glue("~/Software/plink2 --vcf ../data/TopMed.vcfs/TopMed.hg19.vcf.gz --allow-extra-chr --extract range ../data/TopMed.vcfs/TopMed.maf0p01.SNPs --set-all-var-ids @:#:\\$r:\\$a --new-id-max-allele-len 100 --chr 1-22 --make-bed --out ../data/TopMed.vcfs/TopMed.hg19.maf0p01"))
# system(glue::glue("~/Software/plink2 --bfile ../data/TopMed.vcfs/TopMed.hg19.maf0p01 --export A --out ../data/TopMed.vcfs/TopMed.hg19.maf0p01"))


# WriteFileMatrix("../data/CAAPA.vcfs/CAAPA.maf0p01.raw")
# WriteFileMatrix("../data/WGS_Host_Data/joined.hg19.nodup.nomismap.maf0p01.raw")
# WriteFileMatrix("../data/AFGR.vcfs/AFGR.maf0p01.raw")
# WriteFileMatrix("../data/HRC.vcfs/HRC.maf0p01.raw")
# WriteFileMatrix("../data/TopMed.vcfs/TopMed.hg19.maf0p01.raw")

# WGS_Mat <- fm.open("../data/WGS_Host_Data/joined.hg19.nodup.nomismap.maf0p01")
# AFGR_Mat <- fm.open("../data/AFGR.vcfs/AFGR.maf0p01")
# CAAPA_Mat <- fm.open("../data/CAAPA.vcfs/CAAPA.maf0p01")
# HRC_Mat <- fm.open("../data/HRC.vcfs/HRC.maf0p01")
# TopMed_Mat <- fm.open("../data/TopMed.vcfs/TopMed.hg19.maf0p01")
# 
# FindCor(WGS_Mat,AFGR_Mat,"../data/AFGR.vcfs/AFGR.cor.rds")
# FindCor(WGS_Mat,CAAPA_Mat,"../data/CAAPA.vcfs/CAAPA.cor.rds")
# FindCor(WGS_Mat,HRC_Mat,"../data/HRC.vcfs/HRC.cor.rds")
# FindCor(WGS_Mat,TopMed_Mat,"../data/TopMed.vcfs/TopMed.cor.rds")

MAF_Bins <- list(c(0.01,0.05),c(0.05,0.1),c(0.1,0.5))
AFGR_Info_Bins_Mean <- sapply(MAF_Bins,function(x) mean(AFGR_Info$INFO[AFGR_Info$MAF > x[1] & AFGR_Info$MAF <= x[2]],na.rm = T))
AFGR_Info_Bins_Perc <- sapply(MAF_Bins,function(x) sum(AFGR_Info$INFO[AFGR_Info$MAF > x[1] & AFGR_Info$MAF <= x[2]] > 0.8,na.rm = T) / sum(maf_file$MAF > x[1] & maf_file$MAF <= x[2],na.rm = T))
AFGR_Info_Cor <- readRDS("../data/AFGR.vcfs/AFGR.cor.rds")
AFGR_Info_Cor <- sapply(AFGR_Info_Cor,function(x) ifelse(is.numeric(x),x,NA))
AFGR_Info_Bins_Cor <- sapply(MAF_Bins,function(x) mean(AFGR_Info_Cor[paste0(AFGR_Info$ID.y[AFGR_Info$MAF > x[1] & AFGR_Info$MAF <= x[2]],'_',
                                                                            AFGR_Info$REF[AFGR_Info$MAF > x[1] & AFGR_Info$MAF <= x[2]])],na.rm = T))

CAAPA_Info_Bins_Mean <- sapply(MAF_Bins,function(x) mean(as.numeric(CAAPA_Info$Rsq)[CAAPA_Info$MAF > x[1] & CAAPA_Info$MAF <= x[2]],na.rm = T))
CAAPA_Info_Bins_Perc <- sapply(MAF_Bins,function(x) sum(as.numeric(CAAPA_Info$Rsq)[CAAPA_Info$MAF > x[1] & CAAPA_Info$MAF <= x[2]] > 0.8,na.rm = T) / sum(maf_file$MAF > x[1] & maf_file$MAF <= x[2],na.rm = T))
CAAPA_Info_Cor <- readRDS("../data/CAAPA.vcfs/CAAPA.cor.rds")
CAAPA_Info_Cor <- sapply(CAAPA_Info_Cor,function(x) ifelse(is.numeric(x),x,NA))
CAAPA_Info_Bins_Cor <- sapply(MAF_Bins,function(x) mean(setNames(as.numeric(CAAPA_Info_Cor), names(CAAPA_Info_Cor))[paste0(CAAPA_Info$ID.y[CAAPA_Info$MAF > x[1] & CAAPA_Info$MAF <= x[2]],'_',
                                                                            CAAPA_Info$REF[CAAPA_Info$MAF > x[1] & CAAPA_Info$MAF <= x[2]])],na.rm = T))

HRC_Info_Bins_Mean <- sapply(MAF_Bins,function(x) mean(HRC_Info$INFO[HRC_Info$MAF > x[1] & HRC_Info$MAF <= x[2]],na.rm = T))
HRC_Info_Bins_Perc <- sapply(MAF_Bins,function(x) sum(HRC_Info$INFO[HRC_Info$MAF > x[1] & HRC_Info$MAF <= x[2]] > 0.8,na.rm = T) / sum(maf_file$MAF > x[1] & maf_file$MAF <= x[2],na.rm = T))
HRC_Info_Cor <- readRDS("../data/HRC.vcfs/HRC.cor.rds")
HRC_Info_Cor <- sapply(HRC_Info_Cor,function(x) ifelse(is.numeric(x),x,NA))
HRC_Info_Bins_Cor <- sapply(MAF_Bins,function(x) mean(HRC_Info_Cor[paste0(HRC_Info$ID[HRC_Info$MAF > x[1] & HRC_Info$MAF <= x[2]],'_',
                                                                            HRC_Info$REF[HRC_Info$MAF > x[1] & HRC_Info$MAF <= x[2]])],na.rm = T))
# TopMed_Info_Bins_Mean <- sapply(MAF_Bins,function(x) mean(as.numeric(TopMed_hg19$Rsq)[TopMed_hg19$MAF > x[1] & TopMed_hg19$MAF <= x[2]],na.rm = T))
# TopMed_Info_Bins_Perc <- sapply(MAF_Bins,function(x) sum(as.numeric(TopMed_hg19$Rsq)[TopMed_hg19$MAF > x[1] & TopMed_hg19$MAF <= x[2]] > 0.8,na.rm = T) / sum(maf_file$MAF > x[1] & maf_file$MAF <= x[2],na.rm = T))
# TopMed_Info_Cor <- readRDS("../data/TopMed.vcfs/TopMed.cor.rds")
# TopMed_Info_Cor <- sapply(TopMed_Info_Cor,function(x) ifelse(is.numeric(x),x,NA))
# TopMed_Info_Bins_Cor <- sapply(MAF_Bins,function(x) mean(setNames(as.numeric(TopMed_Info_Cor), names(TopMed_Info_Cor))[paste0(TopMed_hg19$ID.y[TopMed_hg19$MAF > x[1] & TopMed_hg19$MAF <= x[2]],'_',
#                                                                                                                               TopMed_hg19$REF[TopMed_hg19$MAF > x[1] & TopMed_hg19$MAF <= x[2]])],na.rm = T))
