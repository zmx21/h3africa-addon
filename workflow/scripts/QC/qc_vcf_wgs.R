#Run QC, liftover, and extract H3A SNPs

library(dplyr)
library(data.table)
library(here)
library(ggplot2)
library(ggrepel)
library(grid)
library(gridExtra)
#### Specify Path of Inputs ####
software_dir <- paste0(here(),'/software/')
script_dir <- paste0(here(),'/src/')
data_dir_path <- '~/G2G_TB/WGS_Host_Data/'
thousand_genome_path <- '~/G2G_TB/1000_Genomes/'
chip_info_file <- '~/G2G_TB/Chip_Data/h3achip_dbsnp150.tsv'
vcf_file <- 'WGS_Fellay.hg38.joint.118h-1947437863.genotyped.vcf.gz'
clin_file <- 'metadata_shipped_Lausanne_31072019_MZW.txt'
genome_build_dir_path <- '~/G2G_TB/genome_builds/'
sent_sample_file <- 'list_blood_samplesTBDAR_shipped_Lausanne_2018.csv'

#### Re-write VCF with consensus samples and correct names####
#Get Sample IDs in the VCF File
geno_samples <- system(paste0(software_dir,'bcftools query -l ',
                              data_dir_path,vcf_file),intern = T)
geno_samples_rename <- geno_samples
geno_samples_rename[which(geno_samples_rename == 'WGS_Fellay.80499_Fellay_Labs')] <- 'WGS_Fellay.80499'
geno_samples_rename[which(geno_samples_rename == 'WGS_Fellay.80499_TB_DAR')] <- 'WGS_Fellay.80449'

#Rename VCF with correct IDs
id_to_rename <- geno_samples[which(geno_samples_rename != geno_samples)]
id_renamed_to <- geno_samples_rename[which(geno_samples_rename != geno_samples)]
write(paste0(id_to_rename,' ',id_renamed_to),file = paste0(data_dir_path,'name_replacement.txt'),sep = '\n')
system(paste0(software_dir,'bcftools reheader -s ',data_dir_path,'name_replacement.txt ',data_dir_path,vcf_file,' > ',
             data_dir_path,gsub(vcf_file,pattern = 'vcf.gz',replacement = 'renamed.vcf.gz')))
renamed_vcf_file <- gsub(vcf_file,pattern = 'vcf.gz',replacement = 'renamed.vcf.gz')
system(paste0(software_dir,'bcftools index -t --threads 5 ',data_dir_path,renamed_vcf_file))

geno_samples_IDs <- sapply(geno_samples_rename,function(x) gsub(x=x,pattern = 'WGS_Fellay.',replacement = ''))
#Get Sample IDs of the Clinical Data File
clin_samples <- data.table::fread(paste0(data_dir_path,clin_file))$patient_id
#Remove extra id from labels
clin_samples <- sapply(clin_samples,function(x) gsub(x=x,pattern = '-31',replacement = ''))
#Get Sample IDs of samples which have been sent out. 
sent_samples <- as.character(data.table::fread(paste0(data_dir_path,sent_sample_file))$PATIENT_ID)
#Get set of samples which existed in all two files
consensus_samples <- geno_samples_IDs #intersect(geno_samples_IDs,clin_samples) For chip design, keep all samples!
#Get sample missingness, and remove those with >0.5
system(paste0(software_dir,'plink2 --vcf ',data_dir_path,renamed_vcf_file,' --missing sample-only --out ',data_dir_path,renamed_vcf_file))
missingness_tbl <- data.table::fread(paste0(data_dir_path,renamed_vcf_file,'.smiss'))
high_missingness_sample <- dplyr::filter(missingness_tbl,F_MISS > 0.5)$'#IID'
consensus_samples <- setdiff(consensus_samples,gsub(high_missingness_sample,pattern = 'WGS_Fellay.',replacement = ''))

#### Filter VCF File and Only Keeop Consensus####
#Filter vcf files for:
#                     1. GATK PASS
#                     2. Biallelic SNPs only
filter_pipe <- paste0(software_dir,'bcftools view --threads 5 -O u --max-alleles 2 --types \"snps\" ',data_dir_path,renamed_vcf_file,
                      " | ",software_dir,"bcftools filter --threads 5 -O u -i \"FILTER==\'PASS\'\" | ")
#Rewrite vcf file by including only samples which are in the consensus and low missingness
geno_samples_to_include <- unlist(sapply(consensus_samples,
                                         function(x) paste0('WGS_Fellay.',x)))
names(geno_samples_to_include) <- NULL
write(geno_samples_to_include,file = paste0(data_dir_path,'consensus_samples.txt'))

consensus_vcf <- paste0(data_dir_path,gsub(x=renamed_vcf_file,pattern = '.vcf.gz',replacement = '.consensus.filt.vcf.gz'))
system(paste0(filter_pipe,software_dir,'bcftools view --threads 5 -O z -s \"',
             paste(geno_samples_to_include,collapse = ','),'\" ',' > ',consensus_vcf))
system(paste0(software_dir,'bcftools index --threads 5 -t ',consensus_vcf))

#### Re-write Phenotype File ####
#Rewrite the phenotype file such that sample IDs match, in the same order as vcf file
clin_data <- data.table::fread(paste0(data_dir_path,clin_file))
clin_data$patient_id <- sapply(clin_data$patient_id,function(x) gsub(x=x,pattern = '-31',replacement = ''))
clin_data <- clin_data[clin_data$patient_id %in% consensus_samples,]
clin_data$patient_id <- sapply(clin_data$patient_id,function(x) paste0('WGS_Fellay.',x))
#Consolidate one repeating row
dup_index <- which(clin_data$patient_id == clin_data$patient_id[duplicated(clin_data$patient_id)])
clin_data[dup_index[1],]$LINEAGE = paste(clin_data[dup_index[1],]$LINEAGE,clin_data[dup_index[2],]$LINEAGE,sep = '|')
clin_data <- clin_data[-dup_index[2],]
data.table::fwrite(clin_data,paste0(data_dir_path,'consensus_metadata_shipped_Lausanne_31072019_MZW.txt'))

#### Liftover to hg19/GRCh37 from hg38/GRCh38 ####
#Call bash script: liftover.sh
lift_over_vcf <- 'joined.hg19.vcf.gz'
system(paste0(script_dir,'liftover.sh ',genome_build_dir_path,' ',consensus_vcf,' ',data_dir_path,' ',lift_over_vcf))

#### Sex Check ####
lift_over_vcf <- paste0(data_dir_path,'joined.hg19.vcf.gz')
fam_file <- data.table::fread(paste0(gsub(lift_over_vcf,pattern = '.vcf.gz',replacement = ''),'.fam'),fill=T)
clin_sex_info <- dplyr::left_join(fam_file,clin_data,by = c('V2' = 'patient_id'))
fam_file$V5[which(clin_sex_info$`patient:sex` == 'male')] <- 1
fam_file$V5[which(clin_sex_info$`patient:sex` == 'female')] <- 2
fam_file$V5[which(is.na(clin_sex_info$`patient:sex`))]<-0
fam_file$V6 <- -9
data.table::fwrite(fam_file,paste0(gsub(lift_over_vcf,pattern = '.vcf.gz',replacement = '.fam')),col.names = F,sep = ' ')
system(paste0(script_dir,'sex_check.sh ',gsub(lift_over_vcf,pattern = '.vcf.gz',replacement = '')))
sex_check_tbl <- data.table::fread(paste0(gsub(lift_over_vcf,pattern = '.vcf.gz',replacement = ''),'.sexcheck'))
sex_check_tbl <- dplyr::left_join(sex_check_tbl,clin_data,by = c('IID' = 'patient_id'))
sex_check_plot <- ggplot(sex_check_tbl, aes(x=`F`, y=`YCOUNT`,color=`patient:sex`)) +
  geom_point() + 
  geom_label_repel(data= dplyr::filter(sex_check_tbl,`F` < -0 & `YCOUNT` > 15000),aes(label=IID)) + 
  guides(colour = guide_legend(override.aes = list(shape = 16)))

#### Kingship Analysis ####
system(paste0('cp ',gsub(x=lift_over_vcf,pattern = '.vcf.gz',replacement = '.bim'),' ',gsub(x=lift_over_vcf,pattern = '.vcf.gz',replacement = '.diff.fam.bim')))
system(paste0('cp ',gsub(x=lift_over_vcf,pattern = '.vcf.gz',replacement = '.bed'),' ',gsub(x=lift_over_vcf,pattern = '.vcf.gz',replacement = '.diff.fam.bed')))
system(paste0('cp ',gsub(x=lift_over_vcf,pattern = '.vcf.gz',replacement = '.fam'),' ',gsub(x=lift_over_vcf,pattern = '.vcf.gz',replacement = '.diff.fam.fam')))
fam_file <- data.table::fread(gsub(x=lift_over_vcf,pattern = '.vcf.gz',replacement = '.diff.fam.fam'),header = F)
fam_file$V1 <- fam_file$V2
data.table::fwrite(fam_file,gsub(x=lift_over_vcf,pattern = '.vcf.gz',replacement = '.diff.fam.fam'),sep = ' ',quote = F,col.names = F)

system(paste0(software_dir,'KING/king -b ',gsub(x=lift_over_vcf,pattern = '.vcf.gz',replacement = '.diff.fam'),'.bed --related --prefix ',gsub(x=lift_over_vcf,pattern = '.vcf.gz',replacement = '.diff.fam')))
kingship_tbl <- data.table::fread(file = paste0(gsub(x=lift_over_vcf,pattern = '.vcf.gz',replacement = '.diff.fam'),'.kin0'))
kingship_tbl <- dplyr::filter(kingship_tbl,InfType == 'PO' | InfType == 'Dup/MZ')
mono_twins <- dplyr::filter(kingship_tbl,InfType=='Dup/MZ')
write(c(paste(paste0('0 ',mono_twins$ID1),collapse = '\n'),paste(paste0('0 ',mono_twins$ID2),collapse = '\n')),sep = '\n',file = paste0(paste0(gsub(x=lift_over_vcf,pattern = '.vcf.gz',replacement = ''),'.excl.samples.txt')))
kingship_tbl_img <- tableGrob(dplyr::select(kingship_tbl,ID1,ID2,KinshipCoeff = Kinship,PropIBD,InfType))
grid.arrange(kingship_tbl_img)

#### Extract SNPs within array and write VCF File ####
#Remove mismapped SNPs
wgs_file <- 'joined.hg19.nodup.vcf.gz'
system(paste0(software_dir,'bcftools query -f "%CHROM %POS %ID %REF %ALT\n" ',data_dir_path,wgs_file,' > ',data_dir_path,wgs_file,'.pos'))
hg19_pos <- data.table::fread(paste0(data_dir_path,wgs_file,'.pos'),header = F)
dup_rows <- which(duplicated(dplyr::select(hg19_pos,V1,V2,V4,V5)))

rows_to_remove <- c()
for(i in 1:length(dup_rows)){
  cur_chr <- hg19_pos$V1[dup_rows[i]]
  cur_pos <- hg19_pos$V2[dup_rows[i]]
  
  hg19_snp_id <- system(paste0(software_dir,'bcftools view -r ',paste0(cur_chr,':',cur_pos),
                               ' ',data_dir_path,wgs_file,' | ',software_dir,'bcftools query -f "%ID\n" '),
                        intern = T)
  hg19_ref <- system(paste0(software_dir,'bcftools view -r ',paste0(cur_chr,':',cur_pos),
                            ' ',data_dir_path,wgs_file,' | ',software_dir,'bcftools query -f "%REF\n" '),
                     intern = T)
  hg19_alt <- system(paste0(software_dir,'bcftools view -r ',paste0(cur_chr,':',cur_pos),
                            ' ',data_dir_path,wgs_file,' | ',software_dir,'bcftools query -f "%ALT\n" '),
                     intern = T)
  
  dbsnp_snp_id <-  system(paste0(software_dir,'bcftools view -r ',paste0(cur_chr,':',cur_pos),
                                 ' ','~/G2G_TB/genome_builds/human_9606_b150_GRCh37p13/00-All.vcf.gz',' | ',
                                 software_dir,'bcftools query -f "%ID\n" '),
                          intern = T)
  dbsnp_ref <-  system(paste0(software_dir,'bcftools view -r ',paste0(cur_chr,':',cur_pos),
                              ' ','~/G2G_TB/genome_builds/human_9606_b150_GRCh37p13/00-All.vcf.gz',' | ',
                              software_dir,'bcftools query -f "%REF\n" '),
                       intern = T)
  dbsnp_alt<-  system(paste0(software_dir,'bcftools view -r ',paste0(cur_chr,':',cur_pos),
                             ' ','~/G2G_TB/genome_builds/human_9606_b150_GRCh37p13/00-All.vcf.gz',' | ',
                             software_dir,'bcftools query -f "%ALT\n" '),
                      intern = T)
  
  wrong_snp <- setdiff(hg19_snp_id,dbsnp_snp_id)
  correct_snp <- setdiff(hg19_snp_id,wrong_snp)
  if(length(wrong_snp) != 1 | length(hg19_snp_id) != 2 ){
    wrong_snp <- hg19_snp_id
  }else if(!((hg19_ref[which(hg19_snp_id == correct_snp)] == dbsnp_ref[which(dbsnp_snp_id == correct_snp)]) & 
             (hg19_ref[which(hg19_snp_id == correct_snp)] == dbsnp_ref[which(dbsnp_snp_id == correct_snp)]))){
    wrong_snp <- hg19_snp_id
  }
  print(paste0('Liftover:',paste(hg19_snp_id,collapse = ','),
               ' dbsnp:',paste(dbsnp_snp_id,collapse = ','),
               ' removed:',paste(wrong_snp,collapse = ',')))
  rows_to_remove <- c(rows_to_remove,wrong_snp)
}
write(rows_to_remove,file = paste0(data_dir_path,'mismapped_snps.txt'))
#Write no dup and no mismapped vcf file
system(paste0(software_dir,'plink2 --const-fid --bfile ',
              data_dir_path,gsub(x=wgs_file,pattern = '.vcf.gz',replacement = ''),' --exclude ',
              data_dir_path,'mismapped_snps.txt --make-bed --out ',data_dir_path,gsub(x=wgs_file,pattern = '.vcf.gz',replacement = '.nomismap')))
system(paste0(software_dir,'plink2 --const-fid --bfile ',
              data_dir_path,gsub(x=wgs_file,pattern = '.vcf.gz',replacement = ''),' --exclude ',
              data_dir_path,'mismapped_snps.txt --export vcf --out ',data_dir_path,gsub(x=wgs_file,pattern = '.vcf.gz',replacement = '.nomismap')))
system(paste0('bgzip -c ',data_dir_path,gsub(x=wgs_file,pattern = '.vcf.gz',replacement = '.nomismap.vcf'),' > ',data_dir_path,gsub(x=wgs_file,pattern = '.vcf.gz',replacement = '.nomismap.vcf.gz')))
system(paste0('bcftools index -t --threads 5 ',data_dir_path,gsub(x=wgs_file,pattern = '.vcf.gz',replacement = '.nomismap.vcf.gz')))
system(paste0('bcftools query -f "%CHROM %POS %ID %REF %ALT\n" ',data_dir_path,gsub(x=wgs_file,pattern = '.vcf.gz',replacement = '.nomismap.vcf.gz'),' > ',data_dir_path,gsub(x=wgs_file,pattern = '.vcf.gz',replacement = '.nomismap.vcf.gz.pos')))
system(paste0("bcftools view ",data_dir_path,gsub(x=wgs_file,pattern = '.vcf.gz',replacement = '.nomismap.vcf.gz')," | bcftools +fill-tags -- -t AF | bcftools query -f '%CHROM %POS %ID %REF %ALT %INFO/AF\n' ",' > ',data_dir_path,gsub(x=wgs_file,pattern = '.vcf.gz',replacement = '.nomismap.vcf.gz.info')))

#Extract H3A snps 
source(paste0(script_dir,'QC/intersect_h3a_wgs.R'))
system(paste0(script_dir,'QC/extract_h3a_snps.sh ',data_dir_path,' ',gsub(lift_over_vcf,pattern = '.vcf.gz',replacement = '.nodup.nomismap')))

#### Generate Training and Testing Set ####
split_perc <- 3/4
geno_samples_to_include <- setdiff(geno_samples_to_include,c(mono_twins$ID1,mono_twins$ID2))
if(!'training_set' %in% dir(data_dir_path)){
  if(!'testing_split.txt' %in% dir(paste0(data_dir_path,'training_set/'))){
  training_samples <- sample(geno_samples_to_include,size = floor(length(geno_samples_to_include) * split_perc),replace = F)
  testing_samples <- setdiff(geno_samples_to_include,training_samples)
  # saveRDS(training_samples,paste0(data_dir_path,'training_set/training_split.rds'))
  # saveRDS(testing_samples,paste0(data_dir_path,'testing_set/testing_split.rds'))
  write(training_samples,file = paste0(data_dir_path,'training_set/training_split.txt'),sep = '\t',append = F) 
  write(testing_samples,file = paste0(data_dir_path,'testing_set/testing_split.txt'),sep = '\t',append = F)
  }else{
    training_samples <- readRDS(paste0(data_dir_path,'training_set/training_split.rds'))
    testing_samples <- readRDS(paste0(data_dir_path,'testing_set/testing_split.rds'))
  }
}else{
  training_samples <- readRDS(paste0(data_dir_path,'training_set/training_split.rds'))
  testing_samples <- readRDS(paste0(data_dir_path,'testing_set/testing_split.rds'))
}
# h3a_chip_file <- paste0(gsub(lift_over_vcf,pattern = '.vcf.gz',replacement = '.nodup.nomismap'),'.bypos.h3achip.v1.vcf.gz')
h3a_chip_file <- paste0(gsub(lift_over_vcf,pattern = '.vcf.gz',replacement = '.nodup.nomismap'),'.bypos.h3achip.vcf.gz')
h3a_chip_file <- gsub(h3a_chip_file,pattern = data_dir_path,replacement = '')
wgs_file_no_mismap <- 'joined.hg19.nodup.nomismap.vcf.gz'

system(paste0('bgzip -c ',data_dir_path,gsub(x=h3a_chip_file,pattern = '.vcf.gz',replacement = '.vcf'),' > ',data_dir_path,h3a_chip_file))
system(paste0('bcftools index -t --threads 5 ',data_dir_path,h3a_chip_file))

#H3A Training
system(paste0(software_dir,'bcftools view --threads 5 -O z -s \"',
                            paste(training_samples,collapse = ','),'\" ',data_dir_path,h3a_chip_file,' > ',paste0(data_dir_path,'training_set/',gsub(h3a_chip_file,pattern = '.vcf.gz',replacement = '.training.vcf.gz'))))
system(paste0(software_dir,'bcftools index -t --threads 5 ',paste0(data_dir_path,'training_set/',gsub(h3a_chip_file,pattern = '.vcf.gz',replacement = '.training.vcf.gz'))))
system(paste0('gunzip -c ',paste0(data_dir_path,'training_set/',gsub(h3a_chip_file,pattern = '.vcf.gz',replacement = '.training.vcf.gz')),' > ',
              paste0(data_dir_path,'training_set/',gsub(h3a_chip_file,pattern = '.vcf.gz',replacement = '.training.vcf'))))

#WGS Training
system(paste0(software_dir,'bcftools view --threads 5 -O z -s \"',
              paste(training_samples,collapse = ','),'\" ',data_dir_path,wgs_file_no_mismap,' > ',paste0(data_dir_path,'training_set/',gsub(wgs_file_no_mismap,pattern = '.vcf.gz',replacement = '.training.vcf.gz'))))
system(paste0(software_dir,'bcftools index -t --threads 5 ',paste0(data_dir_path,'training_set/',gsub(wgs_file_no_mismap,pattern = '.vcf.gz',replacement = '.training.vcf.gz'))))
system(paste0('gunzip -c ',paste0(data_dir_path,'training_set/',gsub(wgs_file_no_mismap,pattern = '.vcf.gz',replacement = '.training.vcf.gz')),' > ',
              paste0(data_dir_path,'training_set/',gsub(wgs_file_no_mismap,pattern = '.vcf.gz',replacement = '.training.vcf'))))

#H3A Testing
system(paste0(software_dir,'bcftools view --threads 5 -O z -s \"',
              paste(testing_samples,collapse = ','),'\" ',data_dir_path,h3a_chip_file,' > ',paste0(data_dir_path,'testing_set/',gsub(h3a_chip_file,pattern = '.vcf.gz',replacement = '.testing.vcf.gz'))))
system(paste0(software_dir,'bcftools index -t --threads 5 ',paste0(data_dir_path,'testing_set/',gsub(h3a_chip_file,pattern = '.vcf.gz',replacement = '.testing.vcf.gz'))))
system(paste0('gunzip -c ',paste0(data_dir_path,'testing_set/',gsub(h3a_chip_file,pattern = '.vcf.gz',replacement = '.testing.vcf.gz')),' > ',
              paste0(data_dir_path,'testing_set/',gsub(h3a_chip_file,pattern = '.vcf.gz',replacement = '.testing.vcf'))))

#WGS Testing
system(paste0(software_dir,'bcftools view --threads 5 -O z -s \"',
              paste(testing_samples,collapse = ','),'\" ',data_dir_path,wgs_file_no_mismap,' > ',paste0(data_dir_path,'testing_set/',gsub(wgs_file_no_mismap,pattern = '.vcf.gz',replacement = '.testing.vcf.gz'))))
system(paste0(software_dir,'bcftools index -t --threads 5 ',paste0(data_dir_path,'testing_set/',gsub(wgs_file_no_mismap,pattern = '.vcf.gz',replacement = '.testing.vcf.gz'))))
system(paste0('gunzip -c ',paste0(data_dir_path,'testing_set/',gsub(wgs_file_no_mismap,pattern = '.vcf.gz',replacement = '.testing.vcf.gz')),' > ',
              paste0(data_dir_path,'testing_set/',gsub(wgs_file_no_mismap,pattern = '.vcf.gz',replacement = '.testing.vcf'))))
