#download from: https://github.com/svohr/mixemt/blob/master/mixemt/phylotree/mtDNA_tree_Build_17.csv
MtTreeCsv <- data.table::fread('../data/MT/mtDNA_tree_Build_17.csv',header = F)

haplogroup_snp_df <- data.frame(POS=integer(),AncestralAllele=character(),DerivedAllele = character(),Haplogroup = character())

for(i in 1:nrow(MtTreeCsv)){
  cur_row <- MtTreeCsv[i,]
  #Find empty cells (equivalent to depth of tree)
  cur_row_empty <- apply(cur_row,1,function(x) x=="")
  #Get cur subtree haplogroup
  cur_sub_tree_haplo_ind <- min(which(!cur_row_empty))
  cur_sub_tree_haplo <- cur_row[,..cur_sub_tree_haplo_ind]
  cur_sub_tree_haplo <- as.vector(t(cur_sub_tree_haplo))
  #Get Mutations present for the current haplogroup
  cur_sub_tree_mut_ind <- min(which(!cur_row_empty)) + 1
  cur_sub_tree_mut <- cur_row[,..cur_sub_tree_mut_ind]
  if(cur_sub_tree_mut == ''){
    next
  }
  cur_sub_tree_mut_split <- unlist(strsplit(as.vector(t(cur_sub_tree_mut)),split = ' '))
  cur_sub_tree_mut_split <- cur_sub_tree_mut_split[cur_sub_tree_mut_split != ""]

  #Remove exclamation mark from mutations (https://www.phylotree.org/rCRS-oriented_version.htm)
  cur_sub_tree_mut_split <- sapply(cur_sub_tree_mut_split,function(x)gsub(pattern = "!|\\)|\\(",replacement = '',x=x))
  mut_character_length <- sapply(cur_sub_tree_mut_split,nchar)
  #Parse Mutations based on REF_POS_ALT
  cur_ref <- sapply(1:length(cur_sub_tree_mut_split),function(i)strsplit(x=cur_sub_tree_mut_split[i],split='')[[1]][1])
  cur_alt <- sapply(1:length(cur_sub_tree_mut_split),function(i)strsplit(x=cur_sub_tree_mut_split[i],split='')[[1]][mut_character_length[i]])
  #Find if SNP format
  snp_format <- mapply(function(x,y) grepl(pattern = 'A|T|C|G',x = x) & grepl(pattern = 'A|T|C|G|a|t|c|g',x = y),cur_ref,cur_alt)
  cur_pos <- sapply(1:length(cur_sub_tree_mut_split),function(i) paste0(strsplit(x=cur_sub_tree_mut_split[i],split='')[[1]][2:(mut_character_length[i]-1)],collapse = ''))
  cur_pos <- as.numeric(cur_pos)
  haplogroup_snp_df <- rbind(haplogroup_snp_df,data.frame(POS=cur_pos[snp_format],
                                                          AncestralAllele=cur_ref[snp_format],
                                                          DerivedAllele=cur_alt[snp_format],
                                                          Haplogroup = rep(cur_sub_tree_haplo,sum(snp_format))))

  #Find Deletions
  if(any(cur_alt=='d')){
    del_format <- cur_alt=='d'
    for(k in which(del_format)){
      if(!grepl(cur_ref[k],pattern = 'A|T|C|G')){
        cur_pos[k] <- sapply(k,function(i) paste0(strsplit(x=cur_sub_tree_mut_split[i],split='d')[[1]][1],collapse = ''))
        cur_ref[k] <- ''
      }
    }
    haplogroup_snp_df <- rbind(haplogroup_snp_df,data.frame(POS=cur_pos[del_format],
                                                            AncestralAllele=cur_ref[del_format],
                                                            DerivedAllele=cur_alt[del_format],
                                                            Haplogroup = rep(cur_sub_tree_haplo,sum(del_format))))

  }

  #Find Insertions
  cur_ins <- sapply(cur_sub_tree_mut_split,function(x)grepl(pattern = "\\.",x=x))
  if(any(cur_ins)){
    cur_pos[cur_ins] <- sapply(which(cur_ins),function(i) paste0(strsplit(x=cur_sub_tree_mut_split[i],split='\\.')[[1]][1],collapse = ''))
    cur_alt[cur_ins] <- sapply(which(cur_ins),function(i) paste0(strsplit(x=cur_sub_tree_mut_split[i],split='\\.')[[1]][2],collapse = ''))
    haplogroup_snp_df <- rbind(haplogroup_snp_df,data.frame(POS=cur_pos[cur_ins],
                                                            AncestralAllele=rep("",sum(cur_ins)),
                                                            DerivedAllele=cur_alt[cur_ins],
                                                            Haplogroup = rep(cur_sub_tree_haplo,sum(cur_ins))))

  }


}
GRCh38_dbsnp_MT <- data.table::fread(
  '../data/genome_builds/human_9606_b151_GRCh38p7/chrMT.GRCh38.dbsnp151')
colnames(GRCh38_dbsnp_MT) <- c('CHR','ID','POS','REF','ALT','RV')
GRCh38_dbsnp_MT_SNPs <- dplyr::filter(GRCh38_dbsnp_MT,grepl(',',ALT) | (nchar(ALT) == nchar(REF)))

haplogroup_snp_df_snps <- haplogroup_snp_df %>% dplyr::filter(AncestralAllele!='' & DerivedAllele !='d')
haplogroup_snp_df_snps$POS <- as.numeric(haplogroup_snp_df_snps$POS)
haplogroup_snp_df_mapped_snps <- dplyr::left_join(haplogroup_snp_df_snps,GRCh38_dbsnp_MT_SNPs %>%
                                               dplyr::select(POS,ID,REF.dbSNP=REF,ALT.dbSNP=ALT),by=c('POS'='POS')) %>%
  dplyr::filter(!is.na(ID))
data.table::fwrite(haplogroup_snp_df_mapped_snps,file = '../data/MT/mtDNA_tree_Build_17_parsed_ZMX.csv')

#Now we manually take care of indels
GRCh38_dbsnp_MT_indels <- dplyr::filter(GRCh38_dbsnp_MT,!grepl(',',ALT) & (nchar(ALT) != nchar(REF)))
haplogroup_snp_df_indels <- dplyr::filter(haplogroup_snp_df %>% dplyr::filter(AncestralAllele == '' | DerivedAllele == '')) %>% dplyr::distinct(POS,.keep_all=T)
matching_list <- vector(mode = 'list',length = nrow(haplogroup_snp_df_indels))
for(i in 1:nrow(haplogroup_snp_df_indels)){
  if(is.na(as.numeric(haplogroup_snp_df_indels$POS[i]))){
    min_range <- as.numeric(strsplit(x=haplogroup_snp_df_indels$POS[i],split = '-')[[1]][1])
    max_range <- as.numeric(strsplit(x=haplogroup_snp_df_indels$POS[i],split = '-')[[1]][2])
  }else{
    min_range <- as.numeric(haplogroup_snp_df_indels$POS[i])
    max_range <- as.numeric(haplogroup_snp_df_indels$POS[i])
    
  }
  possible_match <- dplyr::filter(GRCh38_dbsnp_MT_indels,POS > min_range - 20 & POS < max_range + 20)
  if(nrow(possible_match) > 0){
    matching_list[[i]] <- cbind(haplogroup_snp_df_indels[i,],possible_match)
  }else{
    matching_list[[i]] <- NULL
  }
}
#Manual Check Passed: 1,17,34,39,43,47,54
matching_list[c(1,17,34,39,43,47,54)]
rsids_indels <- c('rs369704279','rs375488999','rs377245343','rs369786048','rs66492218')


