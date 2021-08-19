# Using population-specific add-on polymorphisms to improve genotype imputation in underrepresented populations
https://www.biorxiv.org/content/10.1101/2021.02.03.429542v1

## Data  
#### Add-on SNPs for the H3Africa Array based on the TB-DAR cohort: [369154_H3A.AddOn_08_05_2020.score.csv](./results/369154_H3A.AddOn_08_05_2020.score.csv)
#### Prioritized Gene Regions: [Prioritized_regions.csv](./data/Prioritized_regions.csv)

## Code
### 1. QC
#### Main Script: [qc_vcf_wgs.R](./src/QC/qc_vcf_wgs.R)
   - Sample/SNP based QC (Quality, Sex Check, Kinship).
   - Liftover from GRCh38 to GRCh37.
   - Extract SNPs on the H3Africa array.
   - Split data into training and testing set.

### 2. Imputation and Add-on Tag Selection 
#### A Snakemake file has been created to run both Imputation, Add-on Tag Selection, and Evaluation of Imputation Results: [Snakefile](./workflow/Snakefile)
  - Input files and settings should be specified in the config file: [config.yaml](./config/config.yaml)
  - Snakemake should be run for Setting 1 first, then Setting 2. The working directory should be [workflow](./workflow). 
  - All outputs will be stored in the [results](./results) directory. 

### 3. Figures
  - Figures for add-on selection are constructed using [plots.R](./workflow/scripts/plots/plots.R)
  - Figures for PCA are constructed using [pca_analysis.R](./workflow/scripts/PCA/pca_analysis.R)
  - Figures for Fst are constructed using [AFRMap.R](./workflow/scripts/Fst/AFRMap.R)

