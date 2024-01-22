#!/bin/bash

cd /home/ma1111/tools/TitanCNA_SV_WGS

/mnt/storage/apps/miniconda3/bin/snakemake -s /home/ma1111/tools/TitanCNA_SV_WGS/TitanCNA.snakefile --cluster-config /home/ma1111/tools/TitanCNA_SV_WGS/config/cluster_qsub.yaml --cluster "qsub -l h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -pe {cluster.pe} -binding {cluster.binding}" --jobs 30 --verbose --rerun-incomplete
