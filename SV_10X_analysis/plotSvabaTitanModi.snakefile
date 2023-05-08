configfile: "config/configPlot.yaml"
configfile: "config/samplesPlot.yaml"

import glob
def getTITANpath(base, id, ext):
  return glob.glob(''.join([base, "optimalClusterSolution/", id, "_cluster*", ext]))

CHRS = ['X']

rule all:
  input: 
  	expand("results/plotSvabaTitanModi/{tumor}/{tumor}_CNA-SV-BX_titan_Modi_chr{chr}.{format}", tumor=config["pairings"], chr=CHRS, format=config["plot_format"]),
	
		
rule plotSvabaTitanModi:
	input:
		svabaVCF="results/combineSvabaTitan/{tumor}/{tumor}.svabaTitan.sv.PoNToolFilter.bedpe",
		titanBinFile=lambda wildcards: getTITANpath(config["titan_results"], wildcards.tumor, ".titan.ichor.cna.txt"),
		titanSegFile=lambda wildcards: getTITANpath(config["titan_results"], wildcards.tumor, ".titan.ichor.seg.noSNPs.txt"),
		titanParamFile=lambda wildcards: getTITANpath(config["titan_results"], wildcards.tumor, ".params.txt")
	output:
		"results/plotSvabaTitanModi/{tumor}/{tumor}_CNA-SV-BX_titan_Modi_chr{chr}.{format}"
	params:
		plotSVCNscript=config["plotSVCN_script"],
		tenXfuncs=config["tenX_funcs"],
		svabafuncs=config["svaba_funcs"],
		plotfuncs=config["plot_funcs"],
		libdir=config["titan_libdir"],
		genomeBuild=config["genomeBuild"],
		genomeStyle=config["genomeStyle"],
		cytobandFile=config["cytobandFile"],
		zoom=config["plot_zoom"],
		chrs=config["plot_chrs"],
		start=config["plot_startPos"],
		end=config["plot_endPos"],
		ylim=config["plot_ylim"],
		geneFile=config["plot_geneFile"],
		size=config["plot_size"],
		format=config["plot_format"]	
	log:
		"logs/plotSvabaTitan/{tumor}/{tumor}_CNA-SV-BX_titan_Modi_chr{chr}.{format}.log"
	shell:
		"Rscript {params.plotSVCNscript} --id {wildcards.tumor} --tenX_funcs {params.tenXfuncs} --svaba_funcs {params.svabafuncs} --plot_funcs {params.plotfuncs} --titan_libdir {params.libdir} --svFile {input.svabaVCF} --titanBinFile {input.titanBinFile} --titanSegFile {input.titanSegFile} --titanParamFile {input.titanParamFile} --chrs {wildcards.chr} --genomeBuild {params.genomeBuild} --genomeStyle {params.genomeStyle} --cytobandFile {params.cytobandFile} --start {params.start} --end {params.end} --zoom {params.zoom} --plotYlim \"{params.ylim}\" --geneFile {params.geneFile} --plotCNAtype \"titan\" --plotSize \"{params.size}\" --outPlotFile {output} > {log} 2> {log}" 

