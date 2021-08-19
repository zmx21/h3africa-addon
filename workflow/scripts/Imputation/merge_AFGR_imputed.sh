#!/bin/bash
OUTPUT=$1
bcftools concat {1..22}.vcf.gz X.vcf.gz -O z -o $OUTPUT".vcf.gz"
bcftools index -t --threads 5 $OUTPUT".vcf.gz"
bcftools query -f "%CHROM %POS %ID %REF %ALT %INFO/RefPanelAF %INFO/INFO\n" $OUTPUT".vcf.gz" > $OUTPUT".info.txt"
plink2 --vcf $OUTPUT".vcf.gz" --max-alleles 2 --set-missing-var-ids "@:#[b37]\$r,\$a" --const-fid --make-bed --out $OUTPUT
cut -f 2 $OUTPUT".bim" | sort | uniq -d > $OUTPUT".duprs"
plink2 --bfile $OUTPUT --exclude $OUTPUT".duprs" --make-bed --out $OUTPUT".nodup"
plink2 --bfile $OUTPUT".nodup" --export vcf --out $OUTPUT".nodup"
bgzip -c $OUTPUT".nodup.vcf" > $OUTPUT".nodup.vcf.gz"
bcftools index -t --threads 5 $OUTPUT".nodup.vcf.gz"
bcftools query -f "%CHROM %POS %ID %REF %ALT\n" $OUTPUT".nodup.vcf.gz" > $OUTPUT".nodup.vcf.gz.pos"
#plink2 --bfile $OUTPUT".nodup" --freq --out $OUTPUT".nodup"
