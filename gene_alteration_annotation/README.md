# Annotation of gene alteration by copy-number
These scripts have been modified from Zhou et al. 2022, available here [https://github.com/GavinHaLab/crpc-sv-pattern-study]

## How to run

### For gene annotation
1. Create a matrix using makeMatrixFromTITAN-ICHOR_segBased_hg38_fromCNA.R. This should create .Rdata file that should be used for the next step of gene annotation run.
This step can be run as follows:
Rscript /path/to/makeMatrixFromTITAN-ICHOR_segBased_hg38_fromCNA.R /path/to/Optimal_ploidy_manual /path/to/metadata/mart_export-2.txt /path/to/metadata/tRCC_allSamples.txt common Gene 1000 output_name

2. Run generate_comut_file_for_CN.R using .Rdata file generated from step 1

### For histogram plot

1. Run runTitanMatrix.sh - This step is similar to runFreq.sh (This takes 100kb_geneMats.RData from step 1)
2. Run runFreq.sh (This takes 100kb_geneMats.RData from step 1)
