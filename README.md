# WGS-TRCC

## Whole genome sequence analyses of translocation renal cell carcinoma (tRCC)

This repository contains codes used in the manuscript "A genetic basis for sex differences revealed in Xp11 translocation renal cell carcinoma" https://www.biorxiv.org/content/10.1101/2023.08.04.552029v1 to analyze 29 tRCC samples (in-house datasets) from 15 patients and other tRCCs samples from publicly available datasets as previously described in the following studies-

TCGA, Durinck et al., 2014 [https://www.nature.com/articles/ng.3146]

Malouf et al., 2014 [https://aacrjournals.org/clincancerres/article/20/15/4129/206245/Next-Generation-Sequencing-of-Translocation-Renal] 

MSK Panel, MSK WES, Oncopanel

Sato et al., 2013 [https://www.nature.com/articles/ng.2699]

Sun et al., 2021 [https://www.nature.com/articles/s41467-021-25618-z]

Qu et al., 2022 [https://www.nature.com/articles/s41467-022-34460-w]

## chrX rearrangements analysis

For rearrangement analysis on chrXi we used data from-

Sturm et al., 2016 [https://www.cell.com/fulltext/S0092-8674(16)00055-6]

Gao et al., 2018 for TCGA RNA-Seq data and the Cancer Cell Line Encyclopedia (CCLE) [https://doi.org/10.1016/j.celrep.2018.03.050]


## ASE

This directory contains ASEReadCounter run to determine allele-specific expression for RNA-seq data at germline heterozygous sites.

## CNA

This directory contains codes for copy number analysis for somatic DNA copy number analysis and rearrangement analysis. As described in the method section, for linked reads samples with more variable total sequence coverage at 10kb intervals, the code used is 
```
TRCC_LR.bin_len.bsh
```
For FFPE and other standard WGS samples, the code used is 
```
TRCC.bsh
```

## SV_10X_analysis
To run SVABA on all TRCC samples with matched-normals, we used snakemake workflow modified from Gavin Ha's lab [https://github.com/GavinHaLab/TitanCNA_SV_WGS] 

## TitanCNA_SV_WGS

This directory contains code for copy number analysis by TITAN on published dataset of Qu et al.'s with their baitset design. samples_qu.yaml are the tumor-normal samples from their data while samples.yaml is an example run on our institutional samples used for running Titan. 

## TitanCNA_10X_snakemake

This directory contains code and snakemake worklow to run Titan on linked-read samples. We ran Titan for copy no. analyses on our four linked read samples (TRCC5, TRCC6, TRCC8 and 
TRCC10). This Titan run differs from the standard WGS Titan run because it starts from the BAM files aligned using Long Ranger software.

## NOTE

Upon release, our study website and data for the url to reference will be: https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs003008.v1.p1





