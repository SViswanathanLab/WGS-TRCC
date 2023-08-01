configfile: "config/samples_male.yaml"
configfile: "config/config_male.yaml" 

import glob
import re
def getFullPathToFile(base, filepath):
	print(glob.glob(''.join([filepath, base, "/", base, ".table"])))
	return glob.glob(''.join([filepath, base, "/", base, ".table"]))


rule all:
	input:
		expand("results/GetPileupSummaries/{tumor}/pileup_summaries_{chromosomes}.table",tumor=config["base_file_name"],chromosomes=config["chromosomes"]),
		expand("results/GatherPileupSummaries/{tumor}/{tumor}.table",tumor=config["base_file_name"]),
		expand("results/CalculateContamination/{tumor}/{tumor}_contamination.table",tumor=config["normals"]),
		expand("results/CalculateContamination/{tumor}/{tumor}.segments.table",tumor=config["normals"]),
		expand("results/FilterMutectCalls/{tumor}/filtered_all.vcf.gz",tumor=config["normals"]),
		expand("results/FilterMutectCalls/{tumor}/filtering_stats.tsv",tumor=config["normals"]),
		expand("results/FilterMutectCalls/{tumor}/filtered.vcf.gz",tumor=config["normals"]),
		expand("results/MergeBamOuts/{tumor}/unsorted.out.bam",tumor=config["normals"]),
		expand("results/MergeBamOuts/{tumor}/bamout.bam",tumor=config["normals"]),
		expand("results/MergeBamOuts/{tumor}/bamout.bai",tumor=config["normals"]),
		expand("results/FilterAlignmentArtifacts/{tumor}/{tumor}_filtered.realigned.vcf.gz",tumor=config["normals"]),
		expand("results/FilterAlignmentArtifacts/{tumor}/{tumor}.final.vcf.gz",tumor=config["normals"])


rule GetPileupSummaries:
	input:
		filepaths = lambda wildcards: config["base_file_name"][wildcards.tumor]
	output:
		"results/GetPileupSummaries/{tumor}/pileup_summaries_{chromosomes}.table"
	params:
		reference_genome = config["reference_genome"],
		gatk = config["gatk_path"],
		variants_for_contamination = config["variants_for_contamination"]
	log:
		"logs/GetPileupSummaries/{tumor}_get_pileup_summaries_{chromosomes}.txt"
	shell:
		"({params.gatk} GetPileupSummaries \
		-R {params.reference_genome} \
		-I {input.filepaths} \
		-V {params.variants_for_contamination} \
		-L {wildcards.chromosomes} \
		-O {output}) 2> {log}"

rule GatherPileupSummaries:
	output:
		"results/GatherPileupSummaries/{tumor}/{tumor}.table"
	params:
		reference_dict = config["reference_dict"],
		gatk = config["gatk_path"],
		chromosomes=config["chromosomes"]
	log:
		"logs/GatherPileupSummaries/{tumor}.log"
	shell:
		"""
		all_pileup_inputs=`for chrom in {params.chromosomes}; do
		printf -- "-I results/GetPileupSummaries/{wildcards.tumor}/pileup_summaries_${{chrom}}.table "; done`
		
		({params.gatk} GatherPileupSummaries \
		--sequence-dictionary {params.reference_dict} \
		$all_pileup_inputs \
		-O {output}) 2> {log}
		"""
	
rule CalculateContamination:
	input:
		tumor_pileup=lambda wildcards: getFullPathToFile(wildcards.tumor, "results/GatherPileupSummaries/"),
		normal_pileup=lambda wildcards: getFullPathToFile(config["normals"][wildcards.tumor], "results/GatherPileupSummaries/")
	output:
		contamination_table="results/GatherPileupSummaries/{tumor}/{tumor}_contamination.table",
		tumor_segmentation="results/GatherPileupSummaries/{tumor}/{tumor}.segments.table"
		
	params:
		gatk = config["gatk_path"]
	log:
		"logs/CalculateContamination/{tumor}/{tumor}_contamination.log"
	shell:
		"({params.gatk} CalculateContamination \
   		-I {input.tumor_pileup} \
   		-matched {input.normal_pileup} \
		--tumor-segmentation {output.tumor_segmentation} \
   		-O {output.contamination_table}) 2> {log}"

rule FilterMutectCalls:
	input:
		unfiltered_vcf = "results/GatherVcfs/{tumor}/gathered_unfiltered.vcf.gz",
		vcf_index = "results/GatherVcfs/{tumor}/gathered_unfiltered.vcf.gz.tbi",
		segments_table = "results/CalculateContamination/{tumor}/{tumor}.segments.table",
		contamination_table = "results/CalculateContamination/{tumor}/{tumor}_contamination.table",
		read_orientation_model = "results/LearnReadOrientationModel/{tumor}/read_orientation_model.tar.gz",
		mutect_stats = "results/MergeMutectStats/{tumor}/mutect_merged.stats"
	output:
		filtered_vcf = "results/FilterMutectCalls/{tumor}/filtered_all.vcf.gz",
		filtering_stats = "results/FilterMutectCalls/{tumor}/filtering_stats.tsv"
	params:
		gatk = config["gatk_path"],
		reference_genome = config["reference_genome"]
	log:
		"logs/FilterMutectCalls/{tumor}/{tumor}_filter_mutect_calls.txt"
	shell:
		"({params.gatk} FilterMutectCalls \
		-R {params.reference_genome} \
		-V {input.unfiltered_vcf} \
		--tumor-segmentation {input.segments_table} \
		--contamination-table {input.contamination_table} \
		--ob-priors {input.read_orientation_model} \
		--stats {input.mutect_stats} \
		--filtering-stats {output.filtering_stats} \
		-O {output.filtered_vcf}) 2> {log}"

rule SelectVariantsForFilterMutectCalls:
	input:
		filtered_all="results/FilterMutectCalls/{tumor}/filtered_all.vcf.gz"
	output:
		filtered_vcf="results/FilterMutectCalls/{tumor}/filtered.vcf.gz"
	params:
		gatk = config["gatk_path"],
		reference_genome = config["reference_genome"],
		interval_list = config["interval_list"]
	log:
		"logs/FilterMutectCalls/{tumor}/SelectVariantsForFilterMutectCalls.txt"
	shell:
		"({params.gatk} SelectVariants \
		-R {params.reference_genome} \
		-L {params.interval_list}\
		-V {input.filtered_all} \
		-O {output.filtered_vcf} \
		--exclude-filtered) 2> {log}"


rule MergeBamOuts:
	output:
		unsorted_output = "results/MergeBamOuts/{tumor}/unsorted.out.bam",
		bam_out = "results/MergeBamOuts/{tumor}/bamout.bam",
		bam_bai = "results/MergeBamOuts/{tumor}/bamout.bai"
	params:
		reference_genome = config["reference_genome"],
		java = config["java"],
		picard_jar = config["picard_jar"],
		chromosomes=config["chromosomes"]
	log:
		"logs/MergeBamOuts/{tumor}/MergeBamOuts.log"
	shell:
		"""
		all_bamout_inputs=`for chrom in {params.chromosomes}; do
		printf -- "I=results/mutect2/{wildcards.tumor}/{wildcards.tumor}_${{chrom}}_bamout.bam "; done`
		
		({params.java} -jar {params.picard_jar} GatherBamFiles \
		R={params.reference_genome} \
		$all_bamout_inputs \
		O={output.unsorted_output}
		
		{params.java} -jar {params.picard_jar} SortSam \
		I={output.unsorted_output} \
		O={output.bam_out} \
		SORT_ORDER=coordinate \
		VALIDATION_STRINGENCY=LENIENT
		
		{params.java} -jar {params.picard_jar} BuildBamIndex \
		I={output.bam_out} \
		O={output.bam_bai} \
		VALIDATION_STRINGENCY=LENIENT) 2> {log}"""


rule FilterAlignmentArtifacts:
	input:
		filtered_vcf="results/FilterMutectCalls/{tumor}/filtered.vcf.gz",
		mergedBamout="results/MergeBamOuts/{tumor}/bamout.bam"
	output:
		filtered_vcf="results/FilterAlignmentArtifacts/{tumor}/{tumor}_filtered.realigned.vcf.gz"
	params:
		gatk = config["gatk_path"],
		reference_genome = config["reference_genome"],
		realignment_index_bundle=config["realignment_index_bundle"]
	log:
		"logs/FilterAlignmentArtifacts/{tumor}/FilterAlignmentArtifacts.txt"
	shell:
		"({params.gatk} FilterAlignmentArtifacts \
		-R {params.reference_genome} \
		-V {input.filtered_vcf} \
		-I {input.mergedBamout} \
		--bwa-mem-index-image {params.realignment_index_bundle} \
		-O {output.filtered_vcf}) 2> {log}"

rule SelectVariantsForFilterAlignmentArtifacts:
	input:
		filtered_vcf="results/FilterAlignmentArtifacts/{tumor}/{tumor}_filtered.realigned.vcf.gz"
	output:
		final_vcf="results/FilterAlignmentArtifacts/{tumor}/{tumor}.final.vcf.gz"
	params:
		gatk = config["gatk_path"],
		reference_genome = config["reference_genome"],
		cytoBand_list=config["cytoBand_list"]
	log:
		"logs/SelectVariantsForFilterAlignmentArtifacts/{tumor}/SelectVariantsForFilterAlignmentArtifacts.txt"
	shell:
		"({params.gatk} SelectVariants \
		-R {params.reference_genome} \
		-XL {params.cytoBand_list}\
		-V {input.filtered_vcf} \
		-O {output.final_vcf} \
		--exclude-filtered) 2> {log}"

