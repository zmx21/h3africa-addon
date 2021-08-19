#Extract H3A SNPs from WGS data, based on manifest.

library(pbmcapply)
library(dplyr)
library(biomaRt)
library(GWASTools)
h3a_v2_path <- '../data/H3Africa/v2/H3Africa_2019_20037295_A1.csv'
wgs_pos_path <- '../data/WGS_Host_Data/joined.hg19.nodup.nomismap.vcf.gz.pos'

#Load h3a snps and sequenced snps
h3a_v2 <- data.table::fread(h3a_v2_path,skip=5) %>% dplyr::select(CHR = Chr,POS = MapInfo,SNP,IlmnID,Name,IlmnStrand)
wgs_pos <- data.table::fread(wgs_pos_path)
colnames(wgs_pos) <- c('CHR','POS','ID','REF','ALT')
wgs_pos_by_chr <- lapply(c(seq(1,22,1),'X','Y','MT'),function(x) dplyr::filter(wgs_pos,CHR == x))
names(wgs_pos_by_chr) <- c(seq(1,22,1),'X','Y','MT')
#Parse SNP
h3a_v2$AlelleA <- unlist(pbmclapply(h3a_v2$SNP,function(x) gsub(pattern = '\\[',replacement = '',x=strsplit(x,split = '/')[[1]][1]),mc.cores = 10))
h3a_v2$AlelleB <- unlist(pbmclapply(h3a_v2$SNP,function(x) gsub(pattern = '\\]',replacement = '',x=strsplit(x,split = '/')[[1]][2]),mc.cores = 10))

#Processed H3A PAR Regions
data("pseudoautosomal.hg19")
kgp_to_rs <- data.table::fread('../data/InfiniumOmni2-5-8v1-5_A1_b151_rsids.txt')
h3a_v2_PAR <- dplyr::filter(h3a_v2,CHR %in% 'XY') 
h3a_v2_PAR1 <- dplyr::filter(h3a_v2_PAR,pseudoautosomal.hg19$start.base[1] <= POS & POS <= pseudoautosomal.hg19$end.base[1])
h3a_v2_PAR2 <- dplyr::filter(h3a_v2_PAR,pseudoautosomal.hg19$start.base[2] <= POS & POS <= pseudoautosomal.hg19$end.base[2])
h3a_v2_PAR3 <- dplyr::filter(h3a_v2_PAR,pseudoautosomal.hg19$start.base[3] <= POS & POS <= pseudoautosomal.hg19$end.base[3])

h3a_v2_PAR3_Y <- dplyr::filter(h3a_v2_PAR,pseudoautosomal.hg19$start.base[6] <= POS & POS <= pseudoautosomal.hg19$end.base[6])


#Illumina provided probes which have been mapped to a chr
h3a_v2_mapped <- dplyr::filter(h3a_v2,CHR %in% c(seq(1,22,1),'X','Y','MT'))
h3a_v2_PAR_X <- rbind(h3a_v2_PAR1,h3a_v2_PAR2,h3a_v2_PAR3)
h3a_v2_PAR_X$CHR <- 'X'
h3a_v2_mapped <- rbind(h3a_v2_mapped,h3a_v2_PAR_X)

#Reverse complement of a nucleotide.
GetReverseComplement <- function(nucleotide){
  if(nucleotide == 'A'){
    return('T')
  }else if (nucleotide == 'T'){
    return('A')
  }else if (nucleotide == 'C'){
    return('G')
  }else if (nucleotide == 'G'){
    return('C')
  }
}

#For all illumina provided probes, check if they're sequenced and if the nucleotide match
wgs_pos_to_incl <- pbmclapply(1:nrow(h3a_v2_mapped),function(i) {
  cur_wgs_snp <- dplyr::filter(wgs_pos_by_chr[[as.character(h3a_v2_mapped$CHR[i])]],POS == h3a_v2_mapped$POS[i])
  if(nrow(cur_wgs_snp) == 0){
    return(F)
  }
  for(j in 1:nrow(cur_wgs_snp)){
    #Check if SNP is reverse complement
    if((cur_wgs_snp$REF[j] == GetReverseComplement(h3a_v2_mapped$AlelleA[i]) & cur_wgs_snp$ALT[j] == GetReverseComplement(h3a_v2_mapped$AlelleB[i])) |
       (cur_wgs_snp$REF[j] == GetReverseComplement(h3a_v2_mapped$AlelleB[i]) & cur_wgs_snp$ALT[j] == GetReverseComplement(h3a_v2_mapped$AlelleA[i]))){
      return(T)
    }
    #Check if SNP is reverse
    else if((cur_wgs_snp$REF[j] == h3a_v2_mapped$AlelleA[i] & cur_wgs_snp$ALT[j] == h3a_v2_mapped$AlelleB[i]) |
            (cur_wgs_snp$REF[j] == h3a_v2_mapped$AlelleB[i] & cur_wgs_snp$ALT[j] == h3a_v2_mapped$AlelleA[i])){
      return(T)
    }
  }
  return(F)
},mc.cores = 15)
saveRDS(wgs_pos_to_incl,'../data/WGS_Host_Data/wgs_pos_to_incl.rds')
wgs_pos_to_incl <- readRDS('../data/WGS_Host_Data/wgs_pos_to_incl.rds')
h3a_v2_mapped_to_incl <- h3a_v2_mapped[unlist(wgs_pos_to_incl),] %>% dplyr::select(CHR,START = POS) %>% dplyr::mutate(END = START)

#Look at unmapped Illumina probes,which has rsid
h3a_v2_unmapped <- dplyr::filter(h3a_v2,!(IlmnID %in% h3a_v2_mapped$IlmnID) & CHR == 0)
snp_mart = useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_snp")

h3a_v2_unmapped_to_incl <- pbmclapply(1:nrow(h3a_v2_unmapped),function(i) {
  cur_name <- h3a_v2_unmapped$Name[i]
  if(grepl(pattern='rs',x = cur_name)){
    if(grepl(pattern='h3a-req',x = cur_name)){
      cur_rs <- gsub(pattern = 'h3a-req-',replacement = '',x = cur_name)
      cur_rs <- gsub(pattern = '_ver2',replacement = '',x = cur_rs)
      cur_rs <- gsub(pattern = '_ver1',replacement = '',x = cur_rs)
    }else if(grepl(pattern='seq-',x = cur_name)){
      cur_rs <- gsub(pattern = 'seq-',replacement = '',x = cur_name)
    }else if(grepl(pattern='h3a_37',x = cur_name)){
      cur_rs <- strsplit(x = cur_name,split = '_')[[1]][5]
    }else{
      cur_rs <- cur_name
    }
    cur_h3a_snp <- getBM(attributes = c('chr_name','chrom_start','chrom_end','allele'), filters = 'snp_filter', values = cur_rs, mart = snp_mart) %>%
    dplyr::filter(chr_name %in% c(seq(1,22,1),'X','Y','MT'))
    cur_h3a_snp$REF <- sapply(cur_h3a_snp$allele,function(x) strsplit(x,split = '/')[[1]][1])
    cur_h3a_snp$ALT <- sapply(cur_h3a_snp$allele,function(x) strsplit(x,split = '/')[[1]][2])
    cur_h3a_snp <-  cur_h3a_snp %>% dplyr::select(CHR=chr_name,POS = chrom_start,REF,ALT)
  }else if(grepl(pattern='h3a_37',x = cur_name) & !grepl(pattern='rs',x = cur_name)){
    cur_chr <- strsplit(x = cur_name,split = '_')[[1]][3]
    cur_pos <- strsplit(x = cur_name,split = '_')[[1]][4]
    cur_ref <- strsplit(x = cur_name,split = '_')[[1]][5]
    cur_alt <- strsplit(x = cur_name,split = '_')[[1]][6]
    cur_h3a_snp <- data.frame(CHR = cur_chr,POS = cur_pos,REF = cur_ref,ALT = cur_alt)
  }else{
    return(data.frame(CHR = character(),START = numeric(),END = numeric()))
  }
  if(nrow(cur_h3a_snp) == 0){
    return(data.frame(CHR = character(),START = numeric(),END = numeric()))
  }
  cur_wgs_snp <- dplyr::filter(wgs_pos_by_chr[[as.character(cur_h3a_snp$CHR)]],POS == cur_h3a_snp$POS)
  if(nrow(cur_wgs_snp) == 0){
    return(data.frame(CHR = character(),START = numeric(),END = numeric()))
  }
  for(j in 1:nrow(cur_wgs_snp)){
    #Check if SNP is reverse complement
    if((cur_wgs_snp$REF[j] == GetReverseComplement(cur_h3a_snp$REF) & cur_wgs_snp$ALT[j] == GetReverseComplement(cur_h3a_snp$ALT)) | 
       (cur_wgs_snp$REF[j] == GetReverseComplement(cur_h3a_snp$ALT) & cur_wgs_snp$ALT[j] == GetReverseComplement(cur_h3a_snp$REF))){
      return(data.frame(CHR = cur_wgs_snp$CHR[j],START = cur_wgs_snp$POS[j],END = cur_wgs_snp$POS[j]))
    }
    #Check if SNP is reverse
    else if((cur_wgs_snp$REF[j] == cur_h3a_snp$REF & cur_wgs_snp$ALT[j] == cur_h3a_snp$ALT) | 
            (cur_wgs_snp$REF[j] == cur_h3a_snp$ALT & cur_wgs_snp$ALT[j] == cur_h3a_snp$REF)){
      return(data.frame(CHR = cur_wgs_snp$CHR[j],START = cur_wgs_snp$POS[j],END = cur_wgs_snp$POS[j]))
    }
  }
  return(data.frame(CHR = character(),START = numeric(),END = numeric()))
},mc.cores = 5)
saveRDS(h3a_v2_unmapped_to_incl,'../data/WGS_Host_Data/h3a_v2_unmapped_to_incl.rds')

h3a_v2_unmapped_to_incl <- do.call(rbind,h3a_v2_unmapped_to_incl)

h3a_v2_to_incl <- rbind(h3a_v2_mapped_to_incl,h3a_v2_unmapped_to_incl) %>% dplyr::distinct() %>% dplyr::arrange(CHR,START)
data.table::fwrite(h3a_v2_to_incl,'../data/WGS_Host_Data/h3a_pos.txt',col.names = F,sep = ' ',row.names = F)

