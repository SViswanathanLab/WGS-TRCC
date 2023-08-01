configfile: "config/samples.yaml"
configfile: "config/config.yaml" 

rule all:
  input:
    expand("results/GetPileupSummaries/{tumor}/pileup_summaries_{chromosomes}.table",tumor=config["base_file_name"],chromosomes=config["chromosomes"]),
    expand("results/MergeBamOuts/{tumor}/unsorted.out.bam",tumor=config["normals"]),
    expand("results/MergeBamOuts/{tumor}/bamout.bam",tumor=config["normals"]),
    expand("results/MergeBamOuts/{tumor}/bamout.bai",tumor=config["normals"]),
    expand("results/GatherPileupSummaries/{tumor}/{tumor}.table",tumor=config["base_file_name"])


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
		
		params.java} -jar {params.picard_jar} SortSam \
		I={output.unsorted_output} \
		O={output.bam_out} \
		SORT_ORDER=coordinate \
		VALIDATION_STRINGENCY=LENIENT
		
		{params.java} -jar {params.picard_jar} BuildBamIndex \
		I={output.bam_out} \
		O={output.bam_bai} \
		VALIDATION_STRINGENCY=LENIENT) 2> {log}
		"""

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
