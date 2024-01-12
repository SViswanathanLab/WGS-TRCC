## For copy number alterations and gene frequency in the tRCC WGS cohort to generate figS2B

The scripts have been modified from Zhou et al. 2022 [https://insight.jci.org/articles/view/161370], available here [https://github.com/GavinHaLab/crpc-sv-pattern-study]

## How to run

### For gene frequency

1. Run runTitanMatrix.sh -  This step generates a matrix file of genes with file.RData to be used in the next step. It can be run as follows:

 ```
  Rscript /path/to/makeMatrixFromTITAN-ICHOR_segBased_hg38.R path/to/titan_cluster_solution path/to/metadata/100kb.bins.bed path/to/all_samples/     
  common chrPosn 10000 path/to/output_run
```

2. Run runFreq.sh - This takes 100kb_geneMats.RData from step 1 and can be run as follows

```
  Rscript getGeneFrequency_filterLength_groupTest_hg38_tRCC.R /path/to/file.RData 0 /path/to/metadata/100kb.bins.bed 0 chrPosn 1e7 overall fdr 
  /path/to/metadata/cgc_tier1.txt output_geneFreq.txt
```

3. Run CN_boxplots.R - This takes the data generated above and creates plots of CN information across the genome. It also requires the files Census_allSat May 13 21_30_39 2023.csv and centromere_locations. This generates fig S2B histogram

Figure S2B:

tRCC_recurrrent_CN_opt_sol_manual_geneFreq_100kb.txt (copy number loss and gain across chromosomes, format includes columns for chromosome, location start, location end, CN loss, and CN gain)

centromere_locations (contains information on the starting and ending positions of the centromeres of each chromosome, downloaded as BED file from the UCSC Genome Browser)

Census_allSat May 13 21_30_39 2023.csv (list of genes that have some evidence of being associated with cancer and corresponding information regarding location, confidence of association, gene effect, etc. More information can be found on cancer.sanger.ac.uk/census)

CN_boxplots.R plots the copy number loss and gain across chromosomes and identifies cancer-related genes that are enriched or deleted in specific areas.
