#!/bin/bash
#Liftover from hg38/GRCh38 to hg19/GRCh37 using picard LiftoverVcf
GENOME_BUILD_DIR=$1
INPUT=$2
OUTPUT_PATH=$3
OUTPUT_FILE=$4
for chr in {{1..22},"X","Y"}; 
do
	#Extract current chr from input file
	bcftools view -r "chr"$chr -O v --threads 5 $INPUT > $OUTPUT_PATH"temp_chr"$chr".vcf"
	#Invoke liftover
	java -jar /mnt/data2/xu/G2G_TB/software/picard.jar LiftoverVcf \
		I=$OUTPUT_PATH"temp_chr"$chr".vcf" \
		O=$OUTPUT_PATH"raw.chr"$chr"."$OUTPUT_FILE \
		CHAIN=$GENOME_BUILD_DIR"hg38ToHg19.over.chain.gz" \
		REJECT=$OUTPUT_PATH"chr"$chr".rejected.liftover.vcf" \
		WRITE_ORIGINAL_POSITION=true \
		R=$GENOME_BUILD_DIR"hg19.fa.gz" > $OUTPUT_PATH"chr"$chr".liftover.log" 2>&1
	#Extract variants which were lifted over to the same chr
	bcftools view -r "chr"$chr -O z --threads 5 $OUTPUT_PATH"raw.chr"$chr"."$OUTPUT_FILE > $OUTPUT_PATH"chr"$chr"."$OUTPUT_FILE
	#Index the file containing variants which stayed in the same chr
	bcftools index --tbi --threads 5 $OUTPUT_PATH"chr"$chr"."$OUTPUT_FILE
	#Remove the temporary file
	rm $OUTPUT_PATH"temp_chr"$chr".vcf"
done
#Merge all files, with variants which stayed in the same chr
bcftools concat --threads 5 -O z -o $OUTPUT_PATH$OUTPUT_FILE $OUTPUT_PATH"chr"{1..22}"."$OUTPUT_FILE $OUTPUT_PATH"chrX."$OUTPUT_FILE $OUTPUT_PATH"chrY."$OUTPUT_FILE
rm $OUTPUT_PATH"chr"*"vcf.gz"*
bcftools index --tbi --threads 5 $OUTPUT_PATH"joined."$OUTPUT_FILE
#Delete all raw files 
rm $OUTPUT_PATH"raw.chr"*".vcf.gz"*
#Convert to PLINK format
plink2 --max-alleles 2 --set-missing-var-ids "@:#[b37]\$r,\$a" --const-fid --vcf $OUTPUT_PATH$OUTPUT_FILE --make-bed --out $OUTPUT_PATH${OUTPUT_FILE/'.vcf.gz'/''}
cut -f 2 $OUTPUT_PATH${OUTPUT_FILE/'.vcf.gz'/'.bim'} | sort | uniq -d > $OUTPUT_PATH${OUTPUT_FILE/'.vcf.gz'/'.dupvars.txt'}
../software/plink2 --const-fid --bfile $OUTPUT_PATH${OUTPUT_FILE/'.vcf.gz'/''} --exclude $OUTPUT_PATH${OUTPUT_FILE/'.vcf.gz'/'.dupvars.txt'} --max-alleles 2 --snps-only --make-bed --out $OUTPUT_PATH${OUTPUT_FILE/'.vcf.gz'/'.nodup'}
../software/plink2 --bfile $OUTPUT_PATH${OUTPUT_FILE/'.vcf.gz'/'.nodup'} --export vcf --out $OUTPUT_PATH${OUTPUT_FILE/'.vcf.gz'/'.nodup'}
bgzip -c $OUTPUT_PATH${OUTPUT_FILE/'.vcf.gz'/'.nodup.vcf'} > $OUTPUT_PATH${OUTPUT_FILE/'.vcf.gz'/'.nodup.vcf.gz'}
bcftools index -t --threads 5 $OUTPUT_PATH${OUTPUT_FILE/'.vcf.gz'/'.nodup.vcf.gz'}
