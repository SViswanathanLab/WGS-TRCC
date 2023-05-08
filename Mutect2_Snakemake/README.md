# Mutect2

mutect2 run for multiple input bams and splits by chromosome.

## Run
```
/mnt/storage/apps/miniconda3/bin/snakemake -s /home/mi724/Tools/Snakemake/Mutect2/mutect2.snakefile --cluster-config /home/mi724/Tools/Snakemake/Mutect2/config/cluster_qsub.yaml --cluster "qsub -l h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -pe {cluster.pe} -binding {cluster.binding}" --jobs 30 --rerun-incomplete
```
