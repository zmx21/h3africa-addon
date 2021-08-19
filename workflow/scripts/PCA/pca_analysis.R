library(bigsnpr)
library(ggplot2)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg19)

#### Specify Path of Inputs ####
software_dir <- '~/Software/'
script_dir <- './scripts/PCA/'
thousand_genome_path <- '../data/1000_Genomes/'
out_path <- '../results/PCA/'
data_dir_path <- '../data/WGS_Host_Data/'
lift_over_vcf <- paste0(data_dir_path,'joined.hg19.vcf.gz')

#### Clumping and Long Range LD ####
#Download 1000 genome data
system(paste0(script_dir,'download_1000genomes_vcf.sh',' ',thousand_genome_path))
system(paste0(script_dir,'join_1000genomes_vcf.sh',' ',thousand_genome_path))

#Compare to 1000 genome data
#Download from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped
thousand_genome_data <- data.table::fread(paste0(thousand_genome_path,'20130606_g1k.ped'))
#Download from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv
pop_tbl <- data.table::fread(paste0(thousand_genome_path,'20131219.populations.tsv'))
thousand_genome_data <- thousand_genome_data %>% dplyr::left_join(pop_tbl,c('Population'='Population Code'))
afr_tbl <- dplyr::filter(thousand_genome_data,`Super Population` == 'AFR') %>%
  dplyr::filter(!Relationship %in% c('child','Child2','paternal brother','paternal father','paternal grandmother')) %>% dplyr::select('FID'='Individual ID','IID'='Individual ID','Super_Population'='Super Population','Population')
fam_file <- data.table::fread(gsub(x = lift_over_vcf,pattern = '.vcf.gz',replacement = '.fam'))
afr_tbl <- rbind(afr_tbl,cbind(fam_file %>% dplyr::select('FID'=V1,'IID'=V2),data.frame(`Super_Population`=rep('TB-DAR',nrow(fam_file)),`Population`=rep('TB-DAR',nrow(fam_file)))))
data.table::fwrite(afr_tbl %>% dplyr::select(FID,IID),sep = ' ',col.names = F,file = paste0(out_path,'AFR.IDs.txt'))
tbdar <- data.table::fread(paste0(gsub(lift_over_vcf,pattern = '.vcf.gz',replacement = ''),'.bim'))
system(paste0(software_dir,'plink2 --vcf ',thousand_genome_path,'joined.1000genomes.vcf.gz --make-bed --out ',thousand_genome_path,'joined.1000genomes'))
thousand_genome <- data.table::fread(paste0(thousand_genome_path,'joined.1000genomes.bim'))
consensus_snps <- intersect(tbdar$V2,thousand_genome$V2)
consensus_snps <- consensus_snps[which(consensus_snps != '.')]
tbdar <- dplyr::filter(tbdar,V2 %in% consensus_snps) %>% dplyr::left_join(thousand_genome,c('V2'='V2'))
filt_consensus_snps <- dplyr::filter(tbdar,!V2 %in% tbdar$V2[which(duplicated(tbdar$V2))])%>%
  dplyr::filter_('V1.x == V1.y') %>% dplyr::filter_('V4.x == V4.y') %>%
  dplyr::filter_('V5.x == V5.y')  %>%  dplyr::filter_('V6.x == V6.y')
consensus_rsid_file <- paste0(out_path,'1000genomes.tbdar.consensus.rsid.txt')
write(filt_consensus_snps$V2,consensus_rsid_file)
merged_bed <- 'merged.1000genomes.TBDAR'
#Merge TBDAR file with 1000 Genomes
system(glue::glue("{script_dir}merge_tbdar_thousand_genome.sh {gsub(x=lift_over_vcf,pattern = '.vcf.gz',replacement = '')} {thousand_genome_path}joined.1000genomes {out_path} {consensus_rsid_file} {merged_bed}"))

#Convert to BigSNPR
#Detection of LRLD
#Calcuated PCA using GCTA
CalculatePCA <- function(out_path,prefix,ncores = 40,afr_only = F){
  if(file.exists(glue::glue("{out_path}{prefix}.rds"))){
    system(glue::glue("rm {out_path}{prefix}.rds"))
  }
  if(file.exists(glue::glue("{out_path}{prefix}.bk"))){
    system(glue::glue("rm {out_path}{prefix}.bk"))
  }
  rds <- snp_readBed(glue::glue("{out_path}{prefix}.bed"), backingfile = glue::glue("{out_path}{prefix}"))
  merged_1KG_TBDAR <- snp_attach(glue::glue("{out_path}{prefix}.rds"))
  G <- merged_1KG_TBDAR$genotypes
  CHR <- merged_1KG_TBDAR$map$chromosome
  POS <- merged_1KG_TBDAR$map$physical.pos

  svd_LRLD <- snp_autoSVD(G,CHR,POS,ncores = ncores)
  
  remove(merged_1KG_TBDAR)
  remove(rds)
  gc()
  
  LRLD <- attr(svd_LRLD,'lrldr')
  LRLD$ID <- seq(1:nrow(LRLD))
  data.table::fwrite(LRLD,glue::glue("{out_path}{prefix}.LRLD.txt"),col.names = F,row.names = F,quote = F,sep = ' ')
  if(afr_only){
    system(glue::glue("~/Software/plink --bfile {out_path}merged.1000genomes.TBDAR.pcafilt --extract {out_path}merged.1000genomes.TBDAR.pcafilt.prune.in  --exclude range {out_path}{prefix}.LRLD.txt --keep {out_path}AFR.IDs.txt --make-bed --out {out_path}{prefix}.pruned"))
    
  }else{
    system(glue::glue("~/Software/plink --bfile {out_path}{prefix} --indep-pairwise 50 5 0.5 --out {out_path}{prefix}"))
    system(glue::glue("~/Software/plink --bfile {out_path}{prefix} --extract {out_path}{prefix}.prune.in --exclude range {out_path}{prefix}.LRLD.txt --make-bed --out {out_path}{prefix}.pruned"))
    
  }
  system(glue::glue("~/Software/gcta --bfile {out_path}{prefix}.pruned --make-grm --threads {ncores} --out {out_path}{prefix}.pruned"))
  system(glue::glue("~/Software/gcta --grm {out_path}{prefix}.pruned --pca 20 --threads {ncores} --out {out_path}{prefix}.pruned"))
  system(glue::glue("~/Software/gcta --bfile {out_path}{prefix}.pruned --pc-loading {out_path}{prefix}.pruned --threads {ncores} --out {out_path}{prefix}.pruned"))
  
  return(svd_LRLD)
  
}
Clump_LRLD <- CalculatePCA(out_path,"merged.1000genomes.TBDAR.pcafilt")
Clump_LRLD_AFR <- CalculatePCA(out_path,"merged.1000genomes.TBDAR.pcafilt.AFR",afr_only = T)

#Read in the eigenvectors
eigen_vect_df_thousand_genome <- data.table::fread(paste0(out_path,merged_bed,'.pcafilt.pruned.eigenvec')) %>% dplyr::select(-V1) %>% dplyr::rename(IID=V2) %>%
  dplyr::left_join(thousand_genome_data,c('IID'='Individual ID'))
colnames(eigen_vect_df_thousand_genome)[colnames(eigen_vect_df_thousand_genome) %in% paste0('V',seq(3,22))] <- paste0('PC',seq(1,20))
eigen_vect_df_thousand_genome$`Data Set` <- factor(ifelse(is.na(eigen_vect_df_thousand_genome$Population),'TB DAR','1000 Genomes'),level = c('1000 Genomes','TB DAR'))
eigen_vect_df_thousand_genome$`Super Population` <- as.factor(as.character(eigen_vect_df_thousand_genome$`Super Population`))
eigen_vect_df_thousand_genome$`Super Population` <-factor(eigen_vect_df_thousand_genome$`Super Population` ,
                                                          levels = rev(levels(eigen_vect_df_thousand_genome$`Super Population`)))

eigen_val <- as.vector(data.table::fread(paste0(out_path,merged_bed,'.pcafilt.pruned.eigenval'))$V1)
var_explained <- eigen_val / sum(eigen_val) * 100

#Plots (Fig S2)

pc_plot_thousand_genome_PC1_PC2 <- ggplot(data = eigen_vect_df_thousand_genome %>% dplyr::arrange(`Data Set`)) +
  aes(x=PC1,y=PC2,color=`Super Population`,shape=`Data Set`) +
  geom_point() + scale_shape_manual(values = c(20,4)) +
  xlab(paste0('PC1 (',signif(var_explained[1],2),'%)'))  +
  ylab(paste0('PC2 (',signif(var_explained[2],2),'%)')) + 
  scale_color_discrete(breaks = levels(eigen_vect_df_thousand_genome$`Super Population`))

eigen_vect_df_afr <- data.table::fread(paste0(out_path,merged_bed,'.pcafilt.AFR.pruned.eigenvec')) %>% dplyr::select(-V1) %>% dplyr::rename(IID=V2) %>%
  dplyr::left_join(thousand_genome_data,c('IID'='Individual ID')) %>% dplyr::mutate(POP_ID = ifelse(is.na(`Population Description`),NA,paste0(`Population Description`,' (',`Population`,')')))
colnames(eigen_vect_df_afr)[colnames(eigen_vect_df_afr) %in% paste0('V',seq(3,22))] <- paste0('PC',seq(1,20))

eigen_val_afr <- as.vector(data.table::fread(paste0(out_path,merged_bed,'.pcafilt.AFR.pruned.eigenval'))$V1)
var_explained_afr <- eigen_val_afr / sum(eigen_val_afr) * 100
eigen_vect_df_afr$`POP_ID` <- as.factor(eigen_vect_df_afr$`POP_ID`)
eigen_vect_df_afr$`Data Set` <- ifelse(is.na(eigen_vect_df_afr$Population),'TB DAR','1000 Genomes')

pc_plot_afr_pc1_pc2 <- ggplot(data = eigen_vect_df_afr %>% dplyr::arrange(`Data Set`)) +
  aes(x=PC1,y=PC2,color=`POP_ID`,shape=`Data Set`) +
  geom_point() + scale_shape_manual(values = c(20,4)) +
  xlab(paste0('PC1 (',signif(var_explained_afr[1],2),'%)'))  +
  ylab(paste0('PC2 (',signif(var_explained_afr[2],2),'%)')) + 
  scale_color_discrete(breaks = levels(eigen_vect_df_afr$`POP_ID`)) + labs(color = "Population")

pc_plot <- ggpubr::ggarrange(plotlist = list(pc_plot_thousand_genome_PC1_PC2 + theme(plot.margin = unit(c(1,10.8,1,1), "lines")),pc_plot_afr_pc1_pc2 + theme(plot.margin = unit(c(1,1,1,1), "lines"))),nrow = 2,labels = c('A)','B)'))


snp_loadings <- data.table::fread(glue::glue("{out_path}merged.1000genomes.TBDAR.pcafilt.AFR.pruned.pcl")) %>% dplyr::select(SNP,contains("loading"))
snp_coord <- data.table::fread(glue::glue("{out_path}merged.1000genomes.TBDAR.pcafilt.AFR.pruned.bim")) %>% dplyr::select(SNP=V2,CHROM=V1,POS=V4)
snp_loadings <- dplyr::left_join(snp_loadings,snp_coord,by=c('SNP'='SNP'))
snp_loadings$CUM_POS <- snp_loadings$POS
chr.lengths = seqlengths(Hsapiens)[1:22]
chr_cum_sum <- cumsum(as.numeric(chr.lengths))
for(i in 2:length(unique(snp_loadings$CHROM))){
  cur_chr <- unique(snp_loadings$CHROM)[i]
  snp_loadings$CUM_POS[snp_loadings$CHROM == cur_chr] <- snp_loadings$POS[snp_loadings$CHROM == cur_chr]  + chr_cum_sum[i-1]
}
loading_plot <- ggpubr::ggarrange(plotlist = list(ggplot2::ggplot(data= snp_loadings,aes(x=CUM_POS,y=pc1_loading)) + geom_bin2d(bins = 30) + scale_fill_continuous(type = "viridis") + xlab('SNP_Pos (Cumulative)') + ggtitle('PC1'),
                       ggplot2::ggplot(data= snp_loadings,aes(x=CUM_POS,y=pc2_loading)) + geom_bin2d(bins = 30) + scale_fill_continuous(type = "viridis") + xlab('SNP_Pos (Cumulative)')+ ggtitle('PC2'),
                       ggplot2::ggplot(data= snp_loadings,aes(x=CUM_POS,y=pc3_loading)) + geom_bin2d(bins = 30) + scale_fill_continuous(type = "viridis") + xlab('SNP_Pos (Cumulative)')+ ggtitle('PC3'),
                       ggplot2::ggplot(data= snp_loadings,aes(x=CUM_POS,y=pc4_loading)) + geom_bin2d(bins = 30) + scale_fill_continuous(type = "viridis") + xlab('SNP_Pos (Cumulative)')+ ggtitle('PC4')),common.legend = T)
