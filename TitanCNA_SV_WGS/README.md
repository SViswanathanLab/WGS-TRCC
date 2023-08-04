## *Snakemake workflow for TitanCNA analysis of tRCC FFPE and standard WGS samples* ##

## Description

This workflow uses the original Titan Snakemake workflow from [https://github.com/gavinha/TitanCNA] starting from analysis-ready bam files for tRCC samples TRCC3, TRCC4, TRCC7, TRCC9, TRCC11, TRCC12, TRCC13, TRCC14, TRCC15, TRCC17, TRCC18 and for publsihed datasets. 

## How to run

```
cd /path/to/TitanCNA_SV_WGS
```
## Invoking Sankemake run

```
path/to/bin/snakemake -s path/to/TitanCNA_10X_snakemake/TitanCNA.snakefile --cluster-config /path/to/config/cluster_qsub.yaml --cluster "qsub -l h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -pe {cluster.pe} -binding {cluster.binding}" --jobs 30 --verbose --rerun-incomplete
```

