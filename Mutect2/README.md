# Somatic mutation calls with Mutect2 from GATK

The scripts on each subfolder were used for somatic variant calls of tRCC samples for FFPE-WGS, Frozen-LR and Frozen-WGS.

## Reference files used throughout the analyses are-

Homo_sapiens_assembly38.fasta

Homo_sapiens_assembly38.dict

af-only-gnomad.hg38.vcf.gz

M2PoN_4.0_WGS_for_public.vcf

mutect2-contamination-variants.hg38_primary.vcf.gz

wgs_metrics_intervals.interval_list

Homo_sapiens_assembly38.index_bundle

cytoBand_acen.hg38.interval_list

### How to run

```
path/to/bin/snakemake -s path/to/Snakemake/Mutect2/gather_mutect_calls.snakefile --cluster-config path/to/Mutect2/config/cluster_qsub.yaml --cluster "qsub -l h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -pe {cluster.pe} -binding {cluster.binding}" --jobs 30 --rerun-incomplete

```