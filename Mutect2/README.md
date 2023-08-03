# Somatic mutation calls with Mutect2 from GATK

This repository contains Mutect2 run with Snakemake pipeline following GATK best practices on analysis-ready bams. 

The scripts on each subfolder were used for somatic variant calls of tRCC samples for FFPE-WGS, Frozen-LR and Frozen-WGS. Sample14, Sample17, and Samples12_16 (joint calling mode) folders represent example Mutect2 calls with Snakemake pipeline for frozen-WGS. Config files are found in each of the sample folders. 

This repository contains all the folders needed to run the snakemake. config/ in each of the sample sub-folders contains a configuration file, with logs/ and results/ generated.

For future runs, config/samples.yaml and config/config.yaml can be updated as required with instructions in mutect2.snakefile.

The scripts have been modified from [https://github.com/GavinHaLab/mutect2_snakemake]

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
