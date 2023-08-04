configfile: "config/config.yaml"
configfile: "config/samples.yaml"

rule correctDepth:
	input:
		expand("results/ichorCNA/{tumor}/{tumor}.cna.seg", tumor=config["pairings"]),
		expand("results/ichorCNA/{tumor}/{tumor}.seg.txt", tumor=config["pairings"]),
		expand("results/ichorCNA/{tumor}/{tumor}.params.txt", tumor=config["pairings"]),
		expand("results/ichorCNA/{tumor}/{tumor}.correctedDepth.txt", tumor=config["pairings"]),
		expand("results/readDepth/{samples}.bin{binSize}.wig", samples=config["samples"], binSize=str(config["binSize"]))

rule read_counter:
	input:
		lambda wildcards: config["samples"][wildcards.samples]
	output:
		"results/readDepth/{samples}.bin{binSize}.wig"		
	params:
		readCounter=config["readCounterScript"],
		binSize=config["binSize"],
		qual="20",
		chrs=config["chrs"]
	resources:
		mem=4
	log:
		"logs/readDepth/{samples}.bin{binSize}.log"
	shell:
		"{params.readCounter} {input} -c {params.chrs} -w {params.binSize} -q {params.qual} > {output} 2> {log}"

rule ichorCNA:
	input:
		tum="results/readDepth/{tumor}.bin" + str(config["binSize"]) + ".wig",
		norm=lambda wildcards: "results/readDepth/" + config["pairings"][wildcards.tumor] + ".bin" + str(config["binSize"]) + ".wig"
	output:
		corrDepth="results/ichorCNA/{tumor}/{tumor}.correctedDepth.txt",
		param="results/ichorCNA/{tumor}/{tumor}.params.txt",
		cna="results/ichorCNA/{tumor}/{tumor}.cna.seg",
		segTxt="results/ichorCNA/{tumor}/{tumor}.seg.txt",
		#seg="results/ichorCNA/{tumor}/{tumor}.seg",
		#rdata="results/ichorCNA/{tumor}/{tumor}.RData",
	params:
		outDir="results/ichorCNA/{tumor}/",
		rscript=config["ichorCNA_rscript"],
		libdir=config["ichorCNA_libdir"],
		id="{tumor}",
		ploidy=config["ichorCNA_ploidy"],
		normal=config["ichorCNA_normal"],
		sex=config["sex"],
		genomeBuild=config["genomeBuild"],
		genomeStyle=config["genomeStyle"],
		gcwig=config["ichorCNA_gcWig"],
		mapwig=config["ichorCNA_mapWig"],
		repTimeWig=config["ichorCNA_repTimeWig"],
		estimateNormal=config["ichorCNA_estimateNormal"],
		estimatePloidy=config["ichorCNA_estimatePloidy"],
		estimateClonality=config["ichorCNA_estimateClonality"],
		scStates=config["ichorCNA_scStates"],
		maxCN=config["ichorCNA_maxCN"],
		includeHOMD=config["ichorCNA_includeHOMD"],
		chrs=config["ichorCNA_chrs"],
		#chrTrain=config["ichorCNA_chrTrain"],
		likModel=config["ichorCNA_likModel"],
		centromere=config["centromere"],
		exons=config["ichorCNA_exons"],
		txnE=config["ichorCNA_txnE"],
		txnStrength=config["ichorCNA_txnStrength"],
		fracReadsChrYMale="0.001",
		minMapScore=config["ichorCNA_minMapScore"],
		maxFracGenomeSubclone=config["ichorCNA_maxFracGenomeSubclone"],
		maxFracCNASubclone=config["ichorCNA_maxFracCNASubclone"],
		normal2IgnoreSC=config["ichorCNA_normal2IgnoreSC"],
		scPenalty=config["ichorCNA_scPenalty"],
		plotFileType=config["ichorCNA_plotFileType"],
		plotYlim=config["ichorCNA_plotYlim"]
	resources:
		mem=4
	log:
		"logs/ichorCNA/{tumor}.log"	
	shell:
		"Rscript {params.rscript} --libdir {params.libdir} --id {params.id} --WIG {input.tum} --repTimeWig {params.repTimeWig} --sex {params.sex} --gcWig {params.gcwig} --mapWig {params.mapwig} --NORMWIG {input.norm} --ploidy \"{params.ploidy}\" --normal \"{params.normal}\" --maxCN {params.maxCN} --includeHOMD {params.includeHOMD} --genomeStyle {params.genomeStyle} --genomeBuild {params.genomeBuild} --chrs \"{params.chrs}\" --estimateNormal {params.estimateNormal} --estimatePloidy {params.estimatePloidy} --estimateScPrevalence {params.estimateClonality} --scStates \"{params.scStates}\" --likModel {params.likModel} --minMapScore {params.minMapScore} --maxFracGenomeSubclone {params.maxFracGenomeSubclone} --maxFracCNASubclone {params.maxFracCNASubclone} --normal2IgnoreSC {params.normal2IgnoreSC} --scPenalty {params.scPenalty} --centromere {params.centromere} --exons.bed {params.exons} --txnE {params.txnE} --txnStrength {params.txnStrength} --fracReadsInChrYForMale {params.fracReadsChrYMale} --plotFileType {params.plotFileType} --plotYLim \"{params.plotYlim}\" --outDir {params.outDir} > {log} 2> {log}"
