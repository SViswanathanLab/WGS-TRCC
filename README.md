# WGS-TRCC

## Whole genome sequence analyses of translocation renal cell carcinoma (tRCC)

This repository contains codes used in the manuscript "A genetic basis for cancer sex differences revealed in Xp11 translocation renal cell carcinoma" https://www.biorxiv.org/content/10.1101/2023.08.04.552029v1 to analyze 29 tRCC samples (in-house datasets) from 15 patients and published tRCC whole-exome sequences.

The published data used in this study are from- 

TCGA, Durinck et al., 2014 [https://www.nature.com/articles/ng.3146]

Malouf et al., 2014 [https://aacrjournals.org/clincancerres/article/20/15/4129/206245/Next-Generation-Sequencing-of-Translocation-Renal] MSK Panel, MSK WES, Oncopanel

Sato et al., 2013 [https://www.nature.com/articles/ng.2699]

Sun et al., 2021 [https://www.nature.com/articles/s41467-021-25618-z]

Qu et al., 2022 [https://www.nature.com/articles/s41467-022-34460-w]

Sturm et al., 2016 [https://www.cell.com/fulltext/S0092-8674(16)00055-6]

## ASE

This directory contains ASEReadCounter run to determine allele-specific expression for RNA-seq data at germline heterozygous sites.

## CNA

This directory contains codes for copy number analysis for somatic DNA copy number analysis and rearrangement analysis. As described in the method section, for linked reads samples with more variable total sequence coverage at 10kb intervals, the code used is 
```
TRCC_LR.bin_len.bsh
```
## SVABA

## TitanCNA_SV_WGS

This directory contains code copy number analysis by TITAN on a published dataset of Qu et al.'s baitset design.







