#!/bin/bash
TBDAR_joined=$1
thousand_genomes_dir=$3
thousand_genomes_joined=$2
consensus_rsid_file=$4
output=$5
#Write consense file of tbdar and 1000 genomes
~/Software/plink2 --bfile $TBDAR_joined --extract $consensus_rsid_file --const-fid --keep $thousand_genomes_dir"tanz.pop" --make-bed --out $thousand_genomes_dir"TBDAR.consensus"
~/Software/plink2 --bfile $thousand_genomes_joined --extract $consensus_rsid_file --const-fid --make-bed --out $thousand_genomes_dir"1000genomes.consensus"
#Join 1000 genomes and tbdar vcf file
~/Software/plink --bfile $thousand_genomes_dir"1000genomes.consensus" -bmerge $thousand_genomes_dir"TBDAR.consensus" --make-bed --out $thousand_genomes_dir$output
~/Software/plink2 --bfile $thousand_genomes_dir$output --geno 0 --maf 0.1 --not-chr 23-26 --make-bed --out $thousand_genomes_dir$output".pcafilt"
~/Software/plink2 --bfile $thousand_genomes_dir$output --geno 0 --maf 0.1 --keep $thousand_genomes_dir"AFR.IDs.txt" --make-bed --out $thousand_genomes_dir$output".pcafilt.AFR"
 