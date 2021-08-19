#!/bin/bash
tbdar_dir=$1
joined_tbdar=$2
../software/plink2 --bfile $joined_tbdar --extract bed1 $tbdar_dir"h3a_pos.txt" --export vcf --out $joined_tbdar".bypos.h3achip"
../software/plink2 --bfile $joined_tbdar --extract bed1 $tbdar_dir"h3a_pos.txt" --make-bed --out $joined_tbdar".bypos.h3achip"
# ../software/plink2 --bfile $joined_tbdar --extract bed1 $tbdar_dir"h3a_pos_v1.txt" --export vcf --out $joined_tbdar".bypos.h3achip.v1"
# ../software/plink2 --bfile $joined_tbdar --extract bed1 $tbdar_dir"h3a_pos_v1.txt" --make-bed --out $joined_tbdar".bypos.h3achip.v1"