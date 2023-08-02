# Annotation of gene alteration by copy-number
These scripts have been modified from Zhou et al. 2022 [https://insight.jci.org/articles/view/161370], available here [https://github.com/GavinHaLab/crpc-sv-pattern-study]

## How to run

### For gene annotation
1. Create a matrix using makeMatrixFromTITAN-ICHOR_segBased_hg38_fromCNA.R. This should create .Rdata file that should be used for the next step of gene annotation run.

  This step can be run as follows:

  Rscript /path/to/makeMatrixFromTITAN-ICHOR_segBased_hg38_fromCNA.R /path/to/Optimal_ploidy_manual /path/to/metadata/mart_export-2.txt 
  /path/to/metadata/tRCC_allSamples.txt common Gene 1000 output_name

2. Run generate_comut_file_for_CN.R using .Rdata file generated from step 1

### For histogram plot

1. Run runTitanMatrix.sh -  This step generates a file.RData to be used in the next step. It can be run as follows:

  Rscript /path/to/makeMatrixFromTITAN-ICHOR_segBased_hg38.R path/to/titan_cluster_solution path/to/metadata/100kb.bins.bed path/to/all_samples/     
  common chrPosn 10000 path/to/output_run

2. Run runFreq.sh - This takes 100kb_geneMats.RData from step 1 and can be run as follows

  Rscript getGeneFrequency_filterLength_groupTest_hg38_tRCC.R /path/to/file.RData 0 /path/to/metadata/100kb.bins.bed 0 chrPosn 1e7 overall fdr 
  /path/to/metadata/cgc_tier1.txt output_geneFreq.txt

3. Run CN_Boxplot.R - This takes the data generated above and creates plots of CN information across the genome. It also requires the files Census_allSat May 13 21_30_39 2023.csv and centromere_locations. 


