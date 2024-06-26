configfile: "config/config.yaml"
configfile: "config/samples.yaml"

CHRS = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX', 'chrY']
#CHRS = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']

import glob
def getLRFullPath(base, filename):
  return glob.glob(''.join([base, filename]))


rule correctMolCov:
  input: 
  	expand("results/bxTile/{samples}/{samples}.bxTile.{chr}.bed", samples=config["samples"], chr=CHRS),
  	#expand("results/bxMol/{samples}/{samples}.bxMol.{chr}.bed", samples=config["samples"], chr=CHRS),
  	#expand("results/moleculeCoverage/{tumor}/{tumor}.cna.seg", tumor=config["pairings"]),
  	#expand("results/moleculeCoverage/{tumor}/{tumor}.seg.txt", tumor=config["pairings"]),
  	#expand("results/moleculeCoverage/{tumor}/{tumor}.params.txt", tumor=config["pairings"]),
  	expand("results/moleculeCoverage/{tumor}/{tumor}.BXcounts.txt", tumor=config["pairings"]),
  	expand("results/moleculeCoverage/{tumor}/{tumor}.cna.seg", tumor=config["pairings"])
  	#expand("results/bxTile/{samples}/", samples=config["samples"])
  	#expand("results/moleculeCoverageBXmol/{tumor}/{tumor}.BXcounts.txt", tumor=config["pairings"])

# rule mkdir:
# 	output:
# 		directory("results/bxTile/{samples}/")
# 	params:
# 		mem=config["std_mem"],
# 		runtime=config["std_runtime"],
# 		pe=config["std_numCores"]
# 	shell:
# 		"mkdir -p {output}"


rule bxTile:
	input:
		lambda wildcards: getLRFullPath(config["samples"][wildcards.samples], config["bamFileName"])
	output:		
		"results/bxTile/{samples}/{samples}.bxTile.{chr}.bed"
	params:
		bxTools=config["bxTools"],
		samTools=config["samTools"],
		mapQual=config["bx_mapQual"],
		bedFile=config["bx_bedFileRoot"]
	log:
		"logs/bxTile/{samples}/{samples}.bxTile.{chr}.log"
	shell:
		"{params.samTools} view -h -F 0x4 -q {params.mapQual} {input} {wildcards.chr} | {params.bxTools} tile - -b {params.bedFile}.{wildcards.chr}.bed > {output} 2> {log}"


rule moleculeCoverage:
	input:
		expand("results/bxTile/{samples}/{samples}.bxTile.{chr}.bed", samples=config["samples"], chr=CHRS),
		#tumDir="results/bxTile/{tumor}/",
		#normDir=lambda wildcards: "results/bxTile/" + config["pairings"][wildcards.tumor] + "/"
	output:
		corrDepth="results/moleculeCoverage/{tumor}/{tumor}.BXcounts.txt",
		cna="results/moleculeCoverage/{tumor}/{tumor}.cna.seg",
		segTxt="results/moleculeCoverage/{tumor}/{tumor}.seg.txt",
		paramTxt="results/moleculeCoverage/{tumor}/{tumor}.params.txt",
		#outDir="results/moleculeCoverage/{tumor}/",
	params:
		tumDir="results/bxTile/{tumor}/",
		normDir=lambda wildcards: "results/bxTile/" + config["pairings"][wildcards.tumor] + "/",
		outDir="results/moleculeCoverage/{tumor}/",
		molCovScript=config["molCov_script"],
		id="{tumor}",
		minReadsPerBX=config["molCov_minReadsPerBX"],
		genomeStyle=config["genomeStyle"],
		chrs=config["molCov_chrs"],
		maxCN=config["molCov_maxCN"],
		likModel=config["molCov_likModel"],
		gcwig=config["molCov_gcWig"],
		mapwig=config["molCov_mapWig"],
		#reptimewig=config["molCov_repTimeWig"],
		titanLibDir=config["TitanCNA_libdir"],
		ichorLibDir=config["ichorCNA_libdir"],
		centromere=config["centromere"]
	log:
		"logs/moleculeCoverage/{tumor}.molCov.log"	
	shell:		
		"Rscript {params.molCovScript} --id {params.id} --tumorBXDir {params.tumDir} --normalBXDir {params.normDir} --minReadsPerBX {params.minReadsPerBX} --genomeStyle {params.genomeStyle} --chrs \"{params.chrs}\" --maxCN {params.maxCN} --likModel {params.likModel} --gcWig {params.gcwig} --mapWig {params.mapwig} --libdirTitanCNA {params.titanLibDir} --libdirIchorCNA {params.ichorLibDir} --outDir {params.outDir} --centromere {params.centromere} > {log} 2> {log}"
