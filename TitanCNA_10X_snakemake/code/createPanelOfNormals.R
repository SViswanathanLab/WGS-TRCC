# file:   createPanelOfNormals.R
# author: Gavin Ha, Ph.D.
# institution: Fred Hutchinson Cancer Research Center#
# contact: <gha@fredhutch.org>

# website: https://github.com/gavinha/TitanCNA_10X_snakemake
# date:   November 12, 2018
# description: Create a Panel of Normals using "bxTools tile" data from matched normals


library(GenomeInfoDb)
library(data.table)
library(HMMcopy)
library(optparse)

options(stringsAsFactors=FALSE, scipen=0)
options(bitmapType='cairo')

option_list <- list(
	make_option(c("--gcWig"), type = "character", help = "GC Wig file for reference genome"),
	make_option(c("--mapWig"), type = "character", default=NULL, help = "Mappabiliy Wig file for reference genome"),
	make_option(c("-f", "--filelist"), type = "character", help = "List of of wig files."),
	make_option(c("-o", "--outfile"), type = "character", help = "Output file."),
	make_option(c("--chrs"), type="character", default="c(1:22,\"X\")", help = "Specify chromosomes to analyze."),
	make_option(c("--chrNormalize"), type="character", default="c(1:22)", help = "Specify chromosomes to normalize GC/mappability biases"),
	make_option(c("--maleChrXLogRThres"), type="numeric", default=-0.80, help = "ChrX Log ratio threshold to confirm as male gender."),
	make_option(c("--minReadsPerBX"), type="integer", default=2, help="Minimum number of reads per barcode. [Default: %default]"),
	make_option(c("-e", "--exons.bed"), type = "character", default=NULL, help = "Path to bed file containing exon regions."),
	make_option(c("--centromere"), type="character", default=NULL, help = "File containing Centromere locations; if not provided then will use hg19 version from ichorCNA package. Default: [%default]"),
	make_option(c("--rmCentromereFlankLength"), type="numeric", default=1e5, help="Length of region flanking centromere to remove. Default: [%default]"),
	make_option(c("--method"), type = "character", default="median", help="Median or Mean."),
	make_option(c("--ylim"), type = "character", default="c(-2,2)", help="Y-limits for plotting of mean/median log ratios"),
	make_option(c("--genomeStyle"), type = "character", default="NCBI", help = "Chr naming convention. NCBI (e.g. 1) or UCSC (e.g. chr1). Default: [%default]"),
	make_option(c("--libdirTitanCNA"), type = "character", default=NULL, help = "TitanCNA library path. Usually exclude this argument unless custom modifications have been made to the TitanCNA R package code and the user would like to source those R files. Default: [%default]"),
	make_option(c("--libdirIchorCNA"), type = "character", default=NULL, help = "ichorCNA library path. Usually exclude this argument unless custom modifications have been made to the ichorCNA R package code and the user would like to source those R files. Default: [%default]")
)
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)

#id <- opt$id
gcWig <- opt$gcWig
mapWig <- opt$mapWig
filelist <- opt$filelist
exons.bed <- opt$exons.bed  # "0" if none specified
centromere <- opt$centromere
flankLength <- opt$rmCentromereFlankLength
genomeStyle <- opt$genomeStyle
method <- opt$method
minReadsPerBX <- opt$minReadsPerBX
outfile <- opt$outfile
chrs <- as.character(eval(parse(text = opt$chrs)))
ylim <- eval(parse(text = opt$ylim))
chrNormalize <- as.character(eval(parse(text=opt$chrNormalize)))
maleChrXLogRThres <- opt$maleChrXLogRThres
outImageFile <- paste0(outfile, "_", method, ".rds")

libdirTitanCNA <- opt$libdirTitanCNA
libdirIchorCNA <- opt$libdirIchorCNA
if (!is.null(libdirTitanCNA) && libdirTitanCNA != "None"){
	source(paste0(libdirTitanCNA, "/R/haplotype.R"))
}
if (!is.null(libdirIchorCNA) && libdirIchorCNA != "None"){
	source(paste0(libdirIchorCNA, "/R/utils.R"))
	source(paste0(libdirIchorCNA, "/R/EM.R"))
	source(paste0(libdirIchorCNA, "/R/output.R"))
	source(paste0(libdirIchorCNA, "/R/segmentation.R"))
	source(paste0(libdirIchorCNA, "/R/plotting.R"))
}else{
	library(ichorCNA)
}


if (!is.null(centromere)){
	centromere <- read.delim(centromere,header=T,stringsAsFactors=F,sep="\t")
}
if (!is.null(exons.bed) && exons.bed != "None"){
	targetedSequences <- read.delim(exons.bed, header=F, sep="\t", skip=86)
}else{
	targetedSequences <- NULL
}


## set genome style for chromosome names
chrsAll <- c(chrs, "Y")
seqlevelsStyle(chrs) <- genomeStyle
seqlevelsStyle(chrsAll) <- genomeStyle
seqlevelsStyle(chrNormalize) <- genomeStyle

## load list of bxTile directories
files <- read.delim(filelist, header = FALSE, stringsAsFactors=FALSE, sep ="\t")[, 1]

## load GC and mappability files 
message("Reading GC and mappability files")
if (is.null(gcWig)) {
  stop("GC Wig file required")
} else {
  gc <- wigToGRanges(gcWig)
}
if (is.null(mapWig)) {
  message("Normalizing without mappability Wig file.")
  map <- NULL
} else {
  map <- wigToGRanges(mapWig)
}

normalGR <- NULL
## Analyze the normal samples ##
for (i in 1:length(files)){
	### LOAD NORMAL FILES ###
	sid <- basename(files[i])
	message("Loading normal file:", files[i])
	normal_reads <- loadBXcountsFromBEDDir(files[i], chrs = chrsAll, minReads = minReadsPerBX)
  	normal_reads$BX.medianNorm <- log2(normal_reads$BX.count / median(normal_reads$BX.count, na.rm=T))
		
	### CORRECT TUMOUR DATA FOR GC CONTENT AND MAPPABILITY BIASES ###
	message("Correcting ", sid)
	## add gc and map to RangedData object ##
	normal_counts <- loadReadCountsFromWig(normal_reads, chrs = chrs, genomeStyle = genomeStyle,
					   gc = gc, map = map, centromere = centromere, flankLength = flankLength, 
					   targetedSequences = targetedSequences, 
					   chrNormalize = chrNormalize, mapScoreThres = 0.9)
	gender <- normal_counts$gender

	## add to 
	if (is.null(normalGR)){
		normalGR <- normal_counts$counts
		colnames(normalGR)[which(colnames(normalGR) %in% "copy")] <- sid
	}else{
		normalGR[[sid]] <- normal_counts$counts$copy
	}
	
	chrXMedian <- gender$chrXMedian
	chrXInd <- grepl("X", normalGR$space)
	## Normalize chrX ##
	normalGR[[sid]][chrXInd] <- normalGR[[sid]][chrXInd] - chrXMedian
	
}
saveRDS(normalGR, file = outImageFile)

## compute median ##
mat <- as.data.frame(normalGR)
mat <- mat[, 11:ncol(mat)]
if (method == "median"){
   medianVal <- apply(mat, 1, median, na.rm = TRUE)
}else if (method == "mean"){
   medianVal <- apply(mat, 1, mean, na.rm = TRUE)
}else{
  stop("method is not specified as median or mean.")
}
normalGR[["Median"]] <- medianVal

## output to text file ##
write.table(as.data.frame(normalGR[,"Median"]), file=paste0(outfile, "_", method, ".txt"), col.names=TRUE, row.names=F, quote=F, sep="\t")

## save GR object ##
saveRDS(normalGR, file = outImageFile)
