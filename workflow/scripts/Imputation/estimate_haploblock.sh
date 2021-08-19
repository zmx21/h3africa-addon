#!/bin/bash
INPUT_DIR="/mnt/data2/xu/G2G_TB/WGS_Host_Data/"
VCF_PREFIX="joined.hg19.nodup.nomismap"
CHR=$1
../software/plink --bfile $INPUT_DIR$VCF_PREFIX --blocks no-pheno-req --chr $CHR --out $INPUT_DIR$VCF_PREFIX"."$CHR
#../software/SMIGPP/smigpp --vcf $INPUT_DIR$VCF_PREFIX"chr"$CHR".phased.vcf" --maf 0.05 --out $INPUT_DIR$VCF_PREFIX"chr"$CHR > $INPUT_DIR"haplo.block.est.log.chr"$CHR".out"
