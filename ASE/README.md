# ASE Pipeline

1. Run ASE_example.py on your sample. In brief, it will perform the following steps:

Convert input RNA bam to fastq
Re-align fastq using TopHat2
Preprocess the resulting bam
Convert input DNA bam to fastq
Re-align fastq using bwa
Preprocess the resulting bam
Run HaplotypeCaller on the bam
Run GenotypeGVCFs on the g.vcf
Filter the resulting VCF to SNVs
Run ASE using this vcf and the realigned RNA bam

The command to submit ASE_example.py is: qsub -t 1:1 ASE_submit.qsub

2. Plot ASE output across chrX genes (escapees excluded) using WGS_ASE_plotting_example_chrX.py. Ensure you have a .seg file generated through copy number pipelines in addition to the ASE output (from step 1).
