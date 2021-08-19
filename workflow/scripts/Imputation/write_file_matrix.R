#Calculate r2 between Imputed and WGS ground truth genotypes

library(snpStats)
library(pbmcapply)
library(filematrix)
library(ggplot2)
library(dplyr)
library(glue)

LoadDosage <- function(chr=NULL,pos=NULL,REF=NULL,ALT=NULL,vcf_path,incl_samples = NULL,chr_prefix='',names=F){
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

WriteChunks <- function(file_matrix_pointer,path,chr){
  start_index <- 1
  rsid <- list()
  for(i in 1:length(chr)) {
    print(paste0('Writing Chr:',chr[i]))
    print(paste0('Start:',start_index))
    cur_matrix <- t(LoadDosage(chr = chr[i],vcf_path = path,names = T))
    file_matrix_pointer[,start_index:(start_index + ncol(cur_matrix) - 1)] <- cur_matrix
    start_index <- start_index + ncol(cur_matrix)
    rsid[[i]] <- colnames(cur_matrix)
  }
  return(rsid)
}

SaveMatrix <- function(imp_path,prefix,out_path){
  #Get SNP Positions
  pos_table <- data.table::fread(paste0(imp_path,'.pos'),header = F)
  #Number of SNPs to include
  n_snps <- nrow(pos_table)
  #Chrom ID of the vcf file
  unique_chr <- unique(pos_table$V1)
  samples_to_incl <- system(paste0('bcftools query -l ',imp_path),intern = T)
  #Create File Matrix
  mat_obj <- fm.create(filenamebase = paste0(out_path,prefix,'_mat'),
                           nrow = length(samples_to_incl),
                           ncol = n_snps,
                           type = 'integer')
  print(paste0('Writing',prefix,' Matrix'))
  #Write row name and col name to newly created file matrix.
  rsid <- WriteChunks(mat_obj,imp_path,unique_chr)
  colnames(mat_obj) <- unlist(rsid)
  rownames(mat_obj) <- samples_to_incl
}


args <- commandArgs(trailingOnly = T)
input_path <- args[[1]]
out_prefix <- args[[2]]
out_path <- args[[3]]
SaveMatrix(input_path,out_prefix,out_path)

#Path to Imputation Results
# hrc_imp_path <- '~/G2G_TB/HRC.vcfs/HRC.imputed.nodup'
# afgr_imp_path <- '~/G2G_TB/AFGRTesting.vcfs/AFGR.imputed.nodup'
# afgr_imp_second_pass_path <- '~/G2G_TB/WGS_Host_Data/WGS_Addtl_Tags/AFGRSecondPass0p9.vcfs/AFGRSecondPass.imputed.nodup'
# afgr_imp_third_pass_sum_path <- '~/G2G_TB/WGS_Host_Data/WGS_Addtl_Tags/AFGRThirdPassSum.vcfs/AFGRThirdPass.imputed.nodup'
# afgr_imp_third_pass_avg_path <- '~/G2G_TB/WGS_Host_Data/WGS_Addtl_Tags/AFGRThirdPassAvg.vcfs/AFGRThirdPass.imputed.nodup'
# wgs_path <- '~/G2G_TB/WGS_Host_Data/joined.hg19.nodup.nomismap'
# wgs_imp_path <- '~/G2G_TB/WGS_Imputed/WGS.imputed'
# wgs_imp_second_pass_path <- '~/G2G_TB/WGS_Host_Data/WGS_Addtl_Tags/WGSSecondPass0p9.vcfs/WGS.imputed'
# wgs_imp_third_pass_sum_path <- '~/G2G_TB/WGS_Host_Data/WGS_Addtl_Tags/WGSThirdPassSum.vcfs/WGS.imputed'
# wgs_imp_third_pass_avg_path <- '~/G2G_TB/WGS_Host_Data/WGS_Addtl_Tags/WGSThirdPassAvg.vcfs/WGS.imputed'

#Save Imputation Results as filematrix
# SaveMatrix(hrc_imp_path,'hrc','~/G2G_TB/Imputation_Eval/')
# SaveMatrix(afgr_imp_path,'afgr','~/G2G_TB/Imputation_Eval/')
# SaveMatrix(afgr_imp_second_pass_path,'afgr_second_pass','~/G2G_TB/Imputation_Eval/')
# SaveMatrix(afgr_imp_third_pass_sum_path,'afgr_third_pass_sum','~/G2G_TB/Imputation_Eval/')
# SaveMatrix(afgr_imp_third_pass_avg_path,'afgr_third_pass_avg','~/G2G_TB/Imputation_Eval/')
# SaveMatrix(wgs_path,'wgs','~/G2G_TB/Imputation_Eval/')
# SaveMatrix(wgs_imp_path,'wgs_imp','~/G2G_TB/Imputation_Eval/')
# SaveMatrix(wgs_imp_second_pass_path,'wgs_imp_second_pass','~/G2G_TB/Imputation_Eval/')
# SaveMatrix(wgs_imp_third_pass_sum_path,'wgs_imp_third_pass_sum','~/G2G_TB/Imputation_Eval/')
# SaveMatrix(wgs_imp_third_pass_avg_path,software_dir,'wgs_imp_third_pass_avg','~/G2G_TB/Imputation_Eval/')