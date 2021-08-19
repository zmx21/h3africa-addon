library(filematrix)
library(pbmcapply)
#Evaluate Impuation Performance (r2 and Correlation)
RunCorComparison <- function(h3a_snps,mat1,mat2,mat1_colnames,mat2_colnames,snps_to_excl,samples_to_incl,out_name,non_h3a){
  mat1_rownames <- rownames(mat1)
  mat2_rownames <- rownames(mat2)
  samples_to_incl <- samples_to_incl[samples_to_incl %in% mat1_rownames & samples_to_incl %in% mat2_rownames]
  mat1_samples_to_incl <- match(samples_to_incl,mat1_rownames)
  mat2_samples_to_incl <- match(samples_to_incl,mat2_rownames)
  
  if(non_h3a){
    non_h3a_snps <- intersect(mat1_colnames,setdiff(mat2_colnames,h3a_snps))
    non_h3a_snps <- setdiff(non_h3a_snps,snps_to_excl)
    non_h3a_index_mat1 <- match(non_h3a_snps,mat1_colnames)
    non_h3a_index_mat2 <- match(non_h3a_snps,mat2_colnames)
    remove(mat1_colnames);
    remove(mat2_colnames);
    gc();
    
    non_h3a_cor <- unlist(pbmclapply(1:length(non_h3a_snps),function(i){ 
      cor_result <- tryCatch(
        expr = {cor(as.vector(mat1[mat1_samples_to_incl,non_h3a_index_mat1[i]]),as.vector(mat2[mat2_samples_to_incl,non_h3a_index_mat2[i]]),use = 'complete.obs')},
        error = function(e) {return(NA)},
        warning = function(w){ return(NA)}
      )
      return(cor_result)
    }
    ,mc.cores = 1,ignore.interactive=T))
    saveRDS(list(snps = non_h3a_snps,cor = non_h3a_cor),file = paste0('../results/Imputation_Eval/',out_name,'_non_h3a_cor.rds'))
  }else{
    consensus_h3a_snps <- intersect(h3a_snps,intersect(mat1_colnames,mat2_colnames))
    consensus_h3a_snps <- setdiff(consensus_h3a_snps,snps_to_excl)
    h3a_index_mat1 <- match(consensus_h3a_snps,mat1_colnames)
    h3a_index_mat2 <- match(consensus_h3a_snps,mat2_colnames)
    remove(mat1_colnames);
    remove(mat2_colnames);
    gc();
    
    h3a_cor <- unlist(pbmclapply(1:length(consensus_h3a_snps),function(i){ 
      cor_result <- tryCatch(
        expr = {cor(as.vector(mat1[mat1_samples_to_incl,h3a_index_mat1[i]]),as.vector(mat2[mat2_samples_to_incl,h3a_index_mat2[i]]),use = 'complete.obs')},
        error = function(e) {return(NA)},
        warning = function(w){ return(NA)}
      )
      return(cor_result)
    }
    ,mc.cores = 1,ignore.interactive=T))
    saveRDS(list(snps = consensus_h3a_snps,cor = h3a_cor),file = paste0('../results/Imputation_Eval/',out_name,'_h3a_cor.rds'))
    
  }
  
}

RunAFGRComp <- function(wgs_vcf_path,h3a_snps_path,afgr_imp_path,out_prefix,pass_num,random = T){
  h3a_snps<- system(paste0('~/Software/bcftools query -f "%ID\n" ',h3a_snps_path),intern = T)
  wgs_bim <- data.table::fread(paste0(gsub(wgs_vcf_path,pattern = '.vcf.gz',replacement = '.bim')))
  afgr_bim <- data.table::fread(paste0(afgr_imp_path,'.bim'))
  wgs_afgr <- dplyr::inner_join(wgs_bim,afgr_bim,by = c('V2' = 'V2'))
  wgs_afgr_excl_biallelic <- dplyr::filter(wgs_afgr,(V5.x != V5.y | V6.x != V6.y) & (V4.x == V4.y))
  wgs_afgr_excl_mapping_error <- dplyr::filter(wgs_afgr,V4.x != V4.y)
  wgs_afgr_excl <- union(wgs_afgr_excl_biallelic$V2,wgs_afgr_excl_mapping_error$V2)
  remove(wgs_afgr)
  remove(afgr_bim)
  gc();
  wgs_mat <- fm.open('../results/Imputation_Eval/wgs_mat')
  if(pass_num == 'baseline'){
    afgr_mat <- fm.open('../results/Imputation_Eval/afgr_mat')
  }else if(pass_num == 'setting1'){
    if(random){
      afgr_mat <- fm.open('../results/Imputation_Eval/afgr_setting1_random_mat')
    }else{
      afgr_mat <- fm.open('../results/Imputation_Eval/afgr_setting1_mat')
    }
  }else if (pass_num == 'setting2'){
    if(random){
      afgr_mat <- fm.open('../results/Imputation_Eval/afgr_setting2_random_mat')
    }else{
      afgr_mat <- fm.open('../results/Imputation_Eval/afgr_setting2_mat')
    }
  }
  testing_set <- readRDS('../data/WGS_Host_Data/testing_set/testing_split.rds')
  RunCorComparison(h3a_snps,wgs_mat,afgr_mat,colnames(wgs_mat),colnames(afgr_mat),wgs_afgr_excl,testing_set,out_prefix,T)
}
RunWGSComp <- function(wgs_path,h3a_snps_path,out_prefix,pass_num,random = T){
  h3a_snps<- system(paste0('~/Software/bcftools query -f "%ID\n" ',h3a_snps_path),intern = T)
  wgs_bim <- data.table::fread(paste0(wgs_path,'.bim'))
  wgs_mat <- fm.open('../results/Imputation_Eval/wgs_mat')
  if(pass_num == 'baseline'){
    wgs_imp_mat <- fm.open('../results/Imputation_Eval/wgs_imp_baseline_mat')
    wgs_imp_bim <- data.table::fread('../results/Baseline/WGS_Imputed/WGS.imputed.vcf.gz.pos')
  }else if(pass_num == 'setting1'){
    if(random){
      wgs_imp_mat <- fm.open('../results/Imputation_Eval/wgs_imp_setting1_random_mat')
      wgs_imp_bim <- data.table::fread('../results/WGS_Addtl_Tags/Setting1Random/WGS_Imputed/WGS.imputed.vcf.gz.pos')
    }else{
      wgs_imp_mat <- fm.open('../results/Imputation_Eval/wgs_imp_setting1_mat')
      wgs_imp_bim <- data.table::fread('../results/WGS_Addtl_Tags/Setting1/WGS_Imputed/WGS.imputed.vcf.gz.pos')
    }
  }else if(pass_num == 'setting2'){
    if(random){
      wgs_imp_mat <- fm.open('../results/Imputation_Eval/wgs_imp_setting2_random_mat')
      wgs_imp_bim <- data.table::fread('../results/WGS_Addtl_Tags/Setting2Random/WGS_Imputed/WGS.imputed.vcf.gz.pos')
    }else{
      wgs_imp_mat <- fm.open('../results/Imputation_Eval/wgs_imp_setting2_mat')
      wgs_imp_bim <- data.table::fread('../results/WGS_Addtl_Tags/Setting2/WGS_Imputed/WGS.imputed.vcf.gz.pos')
    }
  }
  #Excl Y chromosome from analysis
  wgs_imp_excl <- dplyr::filter(wgs_bim,!V1 %in% c('X',as.character(seq(1,22,1))))
  testing_set <- readRDS('../data/WGS_Host_Data/testing_set/testing_split.rds') #readRDS('~/G2G_TB/WGS_Host_Data/testing_set/testing_split.rds')
  gc();
  
  wgs_bim_X <- dplyr::filter(wgs_bim,V1 == 'X')
  wgs_imp_bim_X <- dplyr::filter(wgs_imp_bim,V1 == 'X')
  wgs_imp_bim_X_joined <-  dplyr::left_join(wgs_imp_bim_X,wgs_bim_X,by=c('V1'='V1','V2'='V4'))
  wgs_bim_Non_X <- dplyr::filter(wgs_bim,V1 %in% as.character(seq(1,22,1)))
  wgs_imp_colnames <- c(wgs_bim_Non_X$V2,wgs_imp_bim_X_joined$V2.y)
  
  RunCorComparison(h3a_snps,wgs_mat,wgs_imp_mat,wgs_bim$V2,wgs_imp_colnames,wgs_imp_excl$V2,testing_set,out_prefix,T)
}

args <- commandArgs(trailingOnly = T)
if(args[[1]] == 'baseline'){
  if(args[[2]] == 'AFGR'){
    RunAFGRComp('../data/WGS_Host_Data/joined.hg19.nodup.nomismap.vcf.gz',
                '../data/WGS_Host_Data/joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz',
                '../data/AFGRTesting.vcfs/AFGR.imputed.nodup','afgr_','baseline',random = F)
    
  }else if(args[[2]] == 'WGS'){
    RunWGSComp('../data/WGS_Host_Data/joined.hg19.nodup.nomismap.vcf.gz',
               '../data/WGS_Host_Data/joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz','wgs_imp_','baseline',random = F)
  }

}else if(args[[1]] == 'setting1'){
  if(as.logical(args[[3]])){
    if(args[[2]] == 'AFGR'){
      RunAFGRComp('../data/WGS_Host_Data/joined.hg19.nodup.nomismap.vcf.gz',
                  '../data/WGS_Host_Data/joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz',
                  '../results/WGS_Addtl_Tags/Setting1Random/AFGRSetting1Random.vcfs/AFGRSetting1Random.nodup','afgr_setting1_random','setting1',random = T)
      
    }else if(args[[2]] == 'WGS'){
      RunWGSComp('../data/WGS_Host_Data/joined.hg19.nodup.nomismap',
                 '../data/WGS_Host_Data/joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz',
                 'wgs_imp_setting1_random','setting1',random = T)
      
    }
    
  }else{
    if(args[[2]] == 'AFGR'){
      RunAFGRComp('../data/WGS_Host_Data/joined.hg19.nodup.nomismap.vcf.gz',
                  '../data/WGS_Host_Data/joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz',
                  '../results/WGS_Addtl_Tags/Setting1Random/AFGRSetting1.vcfs/AFGRSetting1.nodup','afgr_setting1','setting1',random = F)
      
    }else if(args[[2]] == 'WGS'){
      RunWGSComp('../data/WGS_Host_Data/joined.hg19.nodup.nomismap',
                 '../data/WGS_Host_Data/joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz',
                 'wgs_imp_setting1','setting1',random = F)
      
    }
    
  }
}else if(args[[1]] == 'setting2'){
  if(as.logical(args[[3]])){
    if(args[[2]] == 'AFGR'){
      RunAFGRComp('../data/WGS_Host_Data/joined.hg19.nodup.nomismap.vcf.gz',
                  '../data/WGS_Host_Data/joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz',
                  '../results/WGS_Addtl_Tags/Setting2Random/AFGRSetting2Random.vcfs/AFGRSetting2Random.nodup','afgr_setting2_random','setting2',random = T)
      
    }else if(args[[2]] == 'WGS'){
      RunWGSComp('../data/WGS_Host_Data/joined.hg19.nodup.nomismap',
                 '../data/WGS_Host_Data/joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz',
                 'wgs_imp_setting2_random','setting2',random = T)
      
    }
    
  }else{
    if(args[[2]] == 'AFGR'){
      RunAFGRComp('../data/WGS_Host_Data/joined.hg19.nodup.nomismap.vcf.gz',
                  '../data/WGS_Host_Data/joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz',
                  '../results/WGS_Addtl_Tags/Setting2/AFGRSetting2.vcfs/AFGRSetting2.nodup','afgr_setting2','setting2',random = F)
      
    }else if(args[[2]] == 'WGS'){
      RunWGSComp('../data/WGS_Host_Data/joined.hg19.nodup.nomismap',
                 '../data/WGS_Host_Data/joined.hg19.nodup.nomismap.bypos.h3achip.vcf.gz',
                 'wgs_imp_setting2','setting2',random = F)
      
    }
    
  }
}