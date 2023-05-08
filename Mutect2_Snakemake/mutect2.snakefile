configfile: "config/samples.yaml"
configfile: "config/config.yaml" 

rule all:
    input:
        expand("results/{base_file_name}/unfiltered_{chromosomes}.vcf.gz",base_file_name=config["base_file_name"],chromosomes=config["chromosomes"]),
        expand("results/{base_file_name}/unfiltered_{chromosomes}.vcf.gz.tbi",base_file_name=config["base_file_name"],chromosomes=config["chromosomes"]),
        expand("results/{base_file_name}/unfiltered_{chromosomes}_f1r2.tar.gz",base_file_name=config["base_file_name"],chromosomes=config["chromosomes"]),
        expand("results/{base_file_name}/unfiltered_{chromosomes}.vcf.gz.stats",base_file_name=config["base_file_name"],chromosomes=config["chromosomes"]),
        expand("results/{base_file_name}/mutect_merged.stats", base_file_name = config["base_file_name"]),
        expand("results/MergeMutectStats/{base_file_name}/mutect_merged.stats",base_file_name=config["base_file_name"]),
	expand("results/GatherVcfs/{base_file_name}/gathered_unfiltered.vcf.gz",base_file_name=config["base_file_name"]),
	expand("results/LearnReadOrientationModel/{base_file_name}/read_orientation_model.tar.gz", base_file_name = config["base_file_name"]),
	expand("results/GatherVcfs/{base_file_name}/gathered_unfiltered.vcf.gz.tbi", base_file_name = config["base_file_name"])

rule mutect2:
    input:
        tumor_filepath = config["samples"]
    output:
        vcf = temp("results/mutect2/{base_file_name}/unfiltered_{chromosomes}.vcf.gz"),
        tbi = temp("results/mutect2/{base_file_name}/unfiltered_{chromosomes}.vcf.gz.tbi"),
        tar = temp("results/mutect2/{base_file_name}/unfiltered_{chromosomes}_f1r2.tar.gz"),
        stats = temp("results/mutect2/{base_file_name}/unfiltered_{chromosomes}.vcf.gz.stats")
    params:
        reference_genome = config["reference_genome"],
        germline_resource = config["germline_resource"],
        gatk = config["gatk_path"],
        panel_of_normals = config["panel_of_normals"],
        normals = config["normals"]
    log:
        "logs/mutect2/{base_file_name}_{chromosomes}_mutect2.txt"
    shell:
        "({params.gatk} Mutect2 \
        -reference {params.reference_genome} \
        -input {input.tumor_filepath} \
        -normal {params.normals} \
        -intervals {wildcards.chromosomes} \
        --germline-resource {params.germline_resource} \
        --f1r2-tar-gz {output.tar} \
        --panel-of-normals {params.panel_of_normals} \
        -output {output.vcf}) 2> {log}"


rule MergeMutectStats:
	output:
		"results/MergeMutectStats/{base_file_name}/mutect_merged.stats"
	params:
		gatk = config["gatk_path"],
		chromosomes=config["chromosomes"]
	log:
		"logs/MergeMutectStats/{base_file_name}_merge_mutect_stats.txt"
	shell:
		"""
		all_stat_inputs=`for chrom in {params.chromosomes}; do
		printf -- "-stats results/{wildcards.base_file_name}/unfiltered_$chrom.vcf.gz.stats "; done`
		
		({params.gatk} MergeMutectStats \
		$all_stat_inputs \
		-O {output}) 2> {log}"""

rule GatherVcfs:
	output:
		"results/GatherVcfs/{base_file_name}/gathered_unfiltered.vcf.gz"
	params:
		java = config["java"],
		picard_jar = config["picard_jar"],
		chromosomes=config["chromosomes"]
	log:
		"logs/GatherVcfs/{base_file_name}_gather_mutect_calls.txt"
	shell:
		"""
		all_vcf_inputs=`for chrom in {params.chromosomes}; do
		printf -- "I=results/{wildcards.base_file_name}/unfiltered_$chrom.vcf.gz "; done`
	
		({params.java} -jar {params.picard_jar} GatherVcfs \
		$all_vcf_inputs \
		O={output}) 2> {log}"""


rule LearnReadOrientationModel:
	input:
		chr1_tar = "results/{base_file_name}/unfiltered_chr1_f1r2.tar.gz",
		chr2_tar = "results/{base_file_name}/unfiltered_chr2_f1r2.tar.gz",
		chr3_tar = "results/{base_file_name}/unfiltered_chr3_f1r2.tar.gz",
		chr4_tar = "results/{base_file_name}/unfiltered_chr4_f1r2.tar.gz",
		chr5_tar = "results/{base_file_name}/unfiltered_chr5_f1r2.tar.gz",
		chr6_tar = "results/{base_file_name}/unfiltered_chr6_f1r2.tar.gz",
		chr7_tar = "results/{base_file_name}/unfiltered_chr7_f1r2.tar.gz",
		chr8_tar = "results/{base_file_name}/unfiltered_chr8_f1r2.tar.gz",
		chr9_tar = "results/{base_file_name}/unfiltered_chr9_f1r2.tar.gz",
		chr10_tar = "results/{base_file_name}/unfiltered_chr10_f1r2.tar.gz",
		chr11_tar = "results/{base_file_name}/unfiltered_chr11_f1r2.tar.gz",
		chr12_tar = "results/{base_file_name}/unfiltered_chr12_f1r2.tar.gz",
		chr13_tar = "results/{base_file_name}/unfiltered_chr13_f1r2.tar.gz",
		chr14_tar = "results/{base_file_name}/unfiltered_chr14_f1r2.tar.gz",
		chr15_tar = "results/{base_file_name}/unfiltered_chr15_f1r2.tar.gz",
		chr16_tar = "results/{base_file_name}/unfiltered_chr16_f1r2.tar.gz",
		chr17_tar = "results/{base_file_name}/unfiltered_chr17_f1r2.tar.gz",
		chr18_tar = "results/{base_file_name}/unfiltered_chr18_f1r2.tar.gz",
		chr19_tar = "results/{base_file_name}/unfiltered_chr19_f1r2.tar.gz",
		chr20_tar = "results/{base_file_name}/unfiltered_chr20_f1r2.tar.gz",
		chr21_tar = "results/{base_file_name}/unfiltered_chr21_f1r2.tar.gz",
		chr22_tar = "results/{base_file_name}/unfiltered_chr22_f1r2.tar.gz",
		chrX_tar = "results/{base_file_name}/unfiltered_chrX_f1r2.tar.gz",
		chrY_tar = "results/{base_file_name}/unfiltered_chrY_f1r2.tar.gz"
	output:
		"results/LearnReadOrientationModel/{base_file_name}/read_orientation_model.tar.gz"
	params:
		gatk = config["gatk_path"]
	log:
		"logs/learn_read_orientation_model/{base_file_name}_learn_read_orientation_model.txt"
	shell:
		"({params.gatk} LearnReadOrientationModel \
		-I {input.chr1_tar} \
		-I {input.chr2_tar} \
		-I {input.chr3_tar} \
		-I {input.chr4_tar} \
		-I {input.chr5_tar} \
		-I {input.chr6_tar} \
		-I {input.chr7_tar} \
		-I {input.chr8_tar} \
		-I {input.chr9_tar} \
		-I {input.chr10_tar} \
		-I {input.chr11_tar} \
		-I {input.chr12_tar} \
		-I {input.chr13_tar} \
		-I {input.chr14_tar} \
		-I {input.chr15_tar} \
		-I {input.chr16_tar} \
		-I {input.chr17_tar} \
		-I {input.chr18_tar} \
		-I {input.chr19_tar} \
		-I {input.chr20_tar} \
		-I {input.chr21_tar} \
		-I {input.chr22_tar} \
		-I {input.chrX_tar} \
		-I {input.chrY_tar} \
		-O {output}) 2> {log}"

rule IndexFeatureFile:
	input:
		vcf = expand("results/GatherVcfs/{base_file_name}/gathered_unfiltered.vcf.gz",base_file_name=config["base_file_name"])
	output:
		"results/GatherVcfs/{base_file_name}/gathered_unfiltered.vcf.gz.tbi"
	params:
		gatk = config["gatk_path"]
	log:
		"logs/IndexFeatureFile/{base_file_name}.log"
	shell:
		"({params.gatk} IndexFeatureFile \
		-I {input.vcf} \
		-O {output}) 2> {log}"
