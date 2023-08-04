#!/bin/bash

#module load use.own
#source /mnt/storage/apps/anaconda3/etc/profile.d/conda.sh
#conda activate base

#module load /mnt/storage/spack/share/spack/modules/linux-ubuntu16.04-x86_64/samtools-1.9-gcc-5.4.0-jjq5nua
#module load r-3.6.1-gcc-5.4.0-xtlkk5k

cd /home/ma1111/tools/TitanCNA_SV_WGS

/mnt/storage/apps/anaconda3/bin/snakemake -s /home/ma1111/tools/TitanCNA_SV_WGS/combineSvabaTitan.snakefile --cluster-config /home/ma1111/tools/TitanCNA_SV_WGS/config/cluster_qsub.yaml --cluster "qsub -l h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -pe {cluster.pe} -binding {cluster.binding}" --jobs 30 --verbose --rerun-incomplete 
