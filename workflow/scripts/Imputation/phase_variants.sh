#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 28
#SBATCH --nodes 1
#SBATCH --time 3-00:00:00
#SBATCH --account gr-fe
DATA_DIR='/home/zmxu/G2G_TB/WGS_Host_Data/'
INPUT_PREFIX=$1
OUT_DIR=$2
CHR=$3
GENO_FILT=$4
N_THREAD=28
NO_MISSING_OUT=$OUT_DIR$INPUT_PREFIX".nomissing"
if [ $CHR == 'PAR1' ] || [ $CHR == 'PAR2' ]
then
	#Extract CHR and Remove Samples with gender assignment issues (and SNPs with geno missing if specified)
	if [ $GENO_FILT == 'F' ]
	then
		plink2 --bfile $NO_MISSING_OUT --chr $CHR --remove $OUT_DIR"sex_mismatch.txt" --export vcf --out $NO_MISSING_OUT".chr"$CHR
	else
		plink2 --bfile $NO_MISSING_OUT --chr $CHR --geno 0 --remove $OUT_DIR"sex_mismatch.txt" --export vcf --out $NO_MISSING_OUT".chr"$CHR
	fi
	#Run Phasing
	shapeit -V $NO_MISSING_OUT".chr"$CHR".vcf" \
		-M "/home/zmxu/G2G_TB/Ref_Panel/1000GP_Phase3/genetic_map_chrX_"$CHR"_combined_b37.txt" \
		-O $NO_MISSING_OUT".chr"$CHR".phased" \
		--thread $N_THREAD \
		--force
	#Convert to vcfs
	shapeit -convert \
		--input-haps $NO_MISSING_OUT".chr"$CHR".phased" \
		--output-vcf $OUT_DIR"phased_vcfs/"$INPUT_PREFIX".nomissing.chr"$CHR".phased.vcf"
	

elif [ $CHR == 'X' ]
then
	#Extract CHR and Remove Samples with gender assignment issues (and SNPs with geno missing if specified)
	if [ $GENO_FILT == 'F' ]
	then
		plink2 --bfile $NO_MISSING_OUT --chr $CHR --remove $OUT_DIR"sex_mismatch.txt" --export vcf --out $NO_MISSING_OUT".chr"$CHR
	else
		plink2 --bfile $NO_MISSING_OUT --chr $CHR --geno 0 --remove $OUT_DIR"sex_mismatch.txt" --export vcf --out $NO_MISSING_OUT".chr"$CHR
	fi
	#Remove SNPs which are haploid heterozygous
	awk '{print $3}' $NO_MISSING_OUT".hh" | uniq > $OUT_DIR"excl_hap_het_id.txt"
	plink2 --vcf $NO_MISSING_OUT".chr"$CHR".vcf" --exclude $OUT_DIR"excl_hap_het_id.txt" --export vcf --out $NO_MISSING_OUT".chr"$CHR".nohet"
	#Run Phasing
	shapeit -V $NO_MISSING_OUT".chr"$CHR".nohet.vcf" \
		-M "/home/zmxu/G2G_TB/Ref_Panel/1000GP_Phase3/genetic_map_chrX_nonPAR_combined_b37.txt" \
		-O $NO_MISSING_OUT".chr"$CHR".phased" \
		--chrX \
		--thread $N_THREAD \
		--force
	shapeit -convert \
		--input-haps $NO_MISSING_OUT".chr"$CHR".phased" \
		--output-vcf $OUT_DIR"phased_vcfs/"$INPUT_PREFIX".nomissing.chr"$CHR".phased.vcf"
else
	#Extract CHR (and SNPs with geno missing if specified)
	if [ $GENO_FILT == 'F' ]
	then
		plink2 --vcf $NO_MISSING_OUT".vcf" --chr $CHR --export vcf --out $NO_MISSING_OUT".chr"$CHR
	else
		plink2 --vcf $NO_MISSING_OUT".vcf" --chr $CHR --geno 0 --export vcf --out $NO_MISSING_OUT".chr"$CHR
	fi
	#Run Phasing
	shapeit -V $NO_MISSING_OUT".chr"$CHR".vcf" \
		-M "/home/zmxu/G2G_TB/Ref_Panel/1000GP_Phase3/genetic_map_chr"$CHR"_combined_b37.txt" \
		-O $NO_MISSING_OUT".chr"$CHR".phased" \
		--thread $N_THREAD \
		--force
	shapeit -convert \
		--input-haps $NO_MISSING_OUT".chr"$CHR".phased" \
		--output-vcf $OUT_DIR"phased_vcfs/"$INPUT_PREFIX".nomissing.chr"$CHR".phased.vcf"
fi
