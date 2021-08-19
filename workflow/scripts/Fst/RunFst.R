library(dplyr)
out_path <- '../results/Fst/'
TB_DAR_Path <- '../data/WGS_Host_Data/joined.hg19'
KG_Path <- '../data/1000_Genomes/joined.1000genomes'
n_cores <- 40
snp_missingness_thresh <- 0.1
hwe_thresh <- 1e-5

#Extract superpopulation
pop_map <- data.table::fread('../data/1000_Genomes/20130606_g1k.ped')
super_pop <- data.table::fread('../data/1000_Genomes/20131219.populations.tsv')
pop_map <- dplyr::left_join(pop_map,super_pop,by=c('Population' = 'Population Code'))

pop_to_keep= 'AFR'
pop_map_afr <- dplyr::filter(pop_map,`Super Population` == 'AFR') %>% dplyr::filter(!Relationship %in% c('child','Child2','paternal brother','paternal father','paternal grandmother'))
# pop_map_afr$`Family ID` <- pop_map_afr$`Individual ID`
# 
# data.table::fwrite(pop_map_afr %>% dplyr::select(`Family ID`,`Individual ID`),file = paste0(out_path,pop_to_keep,'.pop'),row.names = F,col.names = F,quote = F,sep = ' ')
# 
# for(i in 1:length(unique(pop_map_afr$Population))){
#   cur_pop <- unique(pop_map_afr$Population)[i]
#   data.table::fwrite(pop_map_afr %>% dplyr::filter(Population == cur_pop) %>% dplyr::select(`Family ID`,`Individual ID`),file = paste0(out_path,cur_pop,'.pop'),row.names = F,col.names = F,quote = F,sep = ' ')
# }
# 
# #Convert ID of 1KG and keep AFR samples
# system(
#   glue::glue("~/Software/plink2 --bfile {KG_Path} --threads {n_cores} --set-all-var-ids @:#[b37]\\$r,\\$a --keep {out_path}{pop_to_keep}.pop  --make-bed --autosome --out {out_path}1KG.tmp1"
#   )
# )
# 
# #Convert ID of TB-DAR
# system(
#   glue::glue("~/Software/plink2 --bfile {TB_DAR_Path} --threads {n_cores} --set-all-var-ids @:#[b37]\\$r,\\$a --make-bed --autosome --out {out_path}TBDAR.tmp1"
#   )
# )
# 
# #1KG Filtered
# system(
#   glue::glue("~/Software/plink2 --bfile {out_path}1KG.tmp1 --threads {n_cores} --geno {snp_missingness_thresh} --hwe {hwe_thresh} --make-bed --autosome --out {out_path}1KG.tmp2"
#   )
# )
# 
# #TB-DAR Filtered
# system(
#   glue::glue("~/Software/plink2 --bfile {out_path}TBDAR.tmp1 --threads {n_cores} --geno {snp_missingness_thresh} --hwe {hwe_thresh} --make-bed --autosome --out {out_path}TBDAR.tmp2"
#   )
# )
# 
# 
# #1KG and TBDAR consensus
# KG_bim <- data.table::fread(glue::glue("{out_path}1KG.tmp2.bim"))
# TBDAR_bim <- data.table::fread(glue::glue("{out_path}TBDAR.tmp2.bim"))
# 
# consensus_snps <- dplyr::inner_join(KG_bim,TBDAR_bim,by=c('V1'='V1','V4'='V4','V5'='V5','V6'='V6'))
# write(unique(consensus_snps$V2.y),file = file(glue::glue("{out_path}consensus")))
# 
# system(
#   glue::glue("~/Software/plink2 --bfile {out_path}1KG.tmp2 --threads {n_cores} --extract {out_path}consensus --make-bed --out {out_path}1KG.consensus.tmp"
#   )
# )
# system(
#   glue::glue("~/Software/plink2 --bfile {out_path}TBDAR.tmp2 --threads {n_cores} --extract {out_path}consensus --make-bed --out {out_path}TBDAR.consensus.tmp"
#   )
# )
# 
# #Merged consensus
# system(
#   glue::glue(
#     "~/Software/plink --bfile {out_path}1KG.consensus.tmp --bmerge {out_path}TBDAR.consensus.tmp --threads {n_cores} --keep-allele-order --make-bed --out {out_path}TBDAR.1KG.merged"
#   )
# )
# fam_file <- data.table::fread(glue::glue("{out_path}TBDAR.1KG.merged.fam"))
# fam_file$V1 <- fam_file$V2
# data.table::fwrite(fam_file,file = glue::glue("{out_path}TBDAR.1KG.merged.fam"),col.names = F,row.names = F,sep = '\t',quote = F)
# 
# 
# fam_file_jned <- dplyr::left_join(fam_file,pop_map_afr,by=c('V1'='Individual ID')) %>% dplyr::select(V1,V2,V3,V4,V5,V6=Population)
# fam_file_jned$V6[is.na(fam_file_jned$V6)] <- 'TBDAR'
# data.table::fwrite(fam_file_jned,file = glue::glue("{out_path}TBDAR.1KG.merged.pedind"),col.names = F,row.names = F,sep = '\t',quote = F)
# 
# #MAF > 0.05
# system(
#   glue::glue(
#     "~/Software/plink --bfile {out_path}TBDAR.1KG.merged --maf 0.05 --keep-allele-order --make-bed --out {out_path}TBDAR.1KG.merged.maf0p05"
#   )
# )
# 
# fam_file <- data.table::fread(glue::glue("{out_path}TBDAR.1KG.merged.maf0p05.fam"))
# fam_file$V1 <- fam_file$V2
# data.table::fwrite(fam_file,file = glue::glue("{out_path}TBDAR.1KG.merged.maf0p05.fam"),col.names = F,row.names = F,sep = '\t',quote = F)
# 
# 
# fam_file_jned <- dplyr::left_join(fam_file,pop_map_afr,by=c('V1'='Individual ID')) %>% dplyr::select(V1,V2,V3,V4,V5,V6=Population)
# fam_file_jned$V6[is.na(fam_file_jned$V6)] <- 'TBDAR'
# data.table::fwrite(fam_file_jned,file = glue::glue("{out_path}TBDAR.1KG.merged.maf0p05.pedind"),col.names = F,row.names = F,sep = '\t',quote = F)
# 
# 
# #Convert to map and ped
# system(
#   glue::glue(
#     "~/Software/plink --bfile {out_path}TBDAR.1KG.merged --threads {n_cores} --export ped --out {out_path}TBDAR.1KG.merged"
#   )
# )
# 
# system(
#   glue::glue(
#     "awk '{{print $1\" \"$2\" \"$3\" \"$4\" \"$5\" \"$6}}' {out_path}TBDAR.1KG.merged.ped > {out_path}TBDAR.1KG.merged.pedind"
#   )
# )
# system(
#   glue::glue(
#     "mv {out_path}TBDAR.1KG.merged.ped {out_path}TBDAR.1KG.merged.ped.raw"
#   )
# )
# system(
#   glue::glue(
#     "awk '{{print substr($0, index($0, $7))}}' {out_path}TBDAR.1KG.merged.ped.raw > {out_path}TBDAR.1KG.merged.ped"
#   )
# )
# 
# pedind <- data.table::fread(glue::glue("{out_path}TBDAR.1KG.merged.pedind"))
# pedind_jned <- dplyr::left_join(pedind,pop_map_afr,by=c('V1'='Individual ID')) %>% dplyr::select(V1,V2,V3,V4,V5,V6=Population)
# pedind_jned$V6[is.na(pedind_jned$V6)] <- 'TBDAR'
# data.table::fwrite(pedind_jned,glue::glue("{out_path}TBDAR.1KG.merged.pedind"),col.names = F,row.names = F,sep = ' ',quote = F)
# 
# 
# system(
#   glue::glue(
#     "~/Software/plink --bfile {out_path}TBDAR.1KG.merged --threads {n_cores} --export ped --maf 0.05 --out {out_path}TBDAR.1KG.merged.maf0p05"
#   )
# )
# system(
#   glue::glue(
#     "awk '{{print $1\" \"$2\" \"$3\" \"$4\" \"$5\" \"$6}}' {out_path}TBDAR.1KG.merged.maf0p05.ped > {out_path}TBDAR.1KG.merged.maf0p05.pedind"
#   )
# )
# 
# system(
#   glue::glue(
#     "mv {out_path}TBDAR.1KG.merged.maf0p05.ped {out_path}TBDAR.1KG.merged.maf0p05.ped.raw"
#   )
# )
# 
# system(
#   glue::glue(
#     "awk '{{print substr($0, index($0, $7))}}' {out_path}TBDAR.1KG.merged.maf0p05.ped.raw > {out_path}TBDAR.1KG.merged.maf0p05.ped"
#   )
# )
# pedind <- data.table::fread(glue::glue("{out_path}TBDAR.1KG.merged.maf0p05.pedind"))
# pedind_jned <- dplyr::left_join(pedind,pop_map_afr,by=c('V1'='Individual ID')) %>% dplyr::select(V1,V2,V3,V4,V5,V6=Population)
# pedind_jned$V6[is.na(pedind_jned$V6)] <- 'TBDAR'
# data.table::fwrite(pedind_jned,glue::glue("{out_path}TBDAR.1KG.merged.maf0p05.pedind"),col.names = F,row.names = F,sep = ' ',quote = F)
# 
# cur_wd <- getwd()
# setwd(out_path)
# system(glue::glue("./run_smartpca"))
# setwd(cur_wd)


