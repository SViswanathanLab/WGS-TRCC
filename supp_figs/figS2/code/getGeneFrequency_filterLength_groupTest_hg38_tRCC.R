
library(data.table)
library(GenomicRanges)
library(plyr)
library(ggplot2)
library(plyranges)
library(dplyr)
library(ggrepel)


options(stringsAsFactors=FALSE)

args <- commandArgs(TRUE)

matFile <- args[1] #R.data matrix file
sampleList <- args[2] #tRCC sample names
annotFile <- args[3] # ensembl annotation file mart-export2
groupType <- args[4] # not sure
headerType <- args[5] # Gene or chrPosn
lengthThreshold <- as.numeric(args[6]) # not sure
cnaType <- args[7] #{'extreme','overall'}
p.adjustMethod <- args[8] # 1e7 ?
#gisticFile <- args[9] #MA changed # gistic file all-lesions.txt
geneFile <- args[10] # tier1_genes file
outRoot <- args[11]
outImage <- paste0(outRoot, ".RData")
#save.image(outImage)

genomeStyle <- "UCSC"
genomeBuild <- "hg38"
seqinfo <- Seqinfo(genome=genomeBuild)
numCores <- 10
chrs <- c(1:22, "X")
clonalThres <- 0.85
signifLevel <- 0.05
ylim <- c(0,6)

load(matFile)# loads all gene matrices
geneCN <- geneCNmat
geneLen <- geneLenmat

#print(head(geneLen))

genes <- fread(annotFile)
if (headerType == "Gene"){
	setnames(genes, c("cdsStart", "cdsEnd"), c("Start", "End"))
	genes <- genes[, .(Chr = Chr[1], Start = min(Start), End = max(End), Karyotype_band=Karyotype_band[1], strand=strand[1]), by=Gene]
	genes[, Length := End - Start]
	genes <- genes[Length > 1]
	genes[, chrPosn := paste0(Chr,":", Start, "-", End)]
}else{ #if (!headerType %in% colnames(genes)){
   colnames(genes)[1:3] <- c("Chr", "Start", "End")
}
seqlevelsStyle(genes$Chr) <- genomeStyle
seqlevelsStyle(chrs) <- genomeStyle
seqinfo <- seqinfo[chrs]
genes <- genes[Chr %in% chrs]
genes$Chr <- factor(genes$Chr, levels = chrs)
genes <- genes[order(Chr, Start)]
numGenes <- nrow(genes)

# order the genes in the matrices
geneCN <- geneCN[, genes[[headerType]]]
geneLen <- geneLen[, genes[[headerType]]]

numGenes <- ncol(geneCN)
numSamples <- nrow(geneCN)

if (nrow(genes) != numGenes){
	stop("Number of genes in matrices don't match gene list.")
}

## load sample list if given ##
if (sampleList != "0"){
	groups <- fread(sampleList)
}else{
	groups <- data.table(Sample = rownames(geneCN), 1)
	colnames(groups)[2] <- groupType
}
groups <- groups[, c(1,which(names(groups)==groupType)), with=F]
types <- unique(groups[, get(groupType)])
numGroups <- length(types)
samples <- groups[[1]]

geneCN <- geneCN[rownames(geneCN) %in% samples, ]


#save.image(outImage)
####################################################
######## CORRECT COPY NUMBER BY MEDIAN CN ##########
####################################################
auto.chr.genes <- genes[!grep("X", Chr), get(headerType)]
cn.autoChrGenes.median <- t(apply(geneCN[, auto.chr.genes], 1, function(x){ x - median(x, na.rm=T) }))
chrX.genes <- genes[grep("X", Chr), get(headerType)]
cn.chrXGenes.median <- t(apply(geneCN[, chrX.genes], 1, function(x){ x - median(x, na.rm=T) }))
cn.median.norm <- cbind(cn.autoChrGenes.median, cn.chrXGenes.median)
geneCN <- cn.median.norm

####################################################
###### FUNCTION TO COMPUTE FREQUENCY ###############
####################################################
getFrequency <- function(mat, neut = NULL, loss=1, gain=3, na.is.neut = FALSE, numSamples = NULL, ind = TRUE){
	if (is.null(numSamples)){
	  numSamples <- nrow(mat)
	}

	if (is.null(neut)){
		indGain <- (mat >= gain) & ind
		indLoss <- (mat <= loss) & ind
		neut <- (gain - loss)
	}else{
		indGain <- (mat >= neut) & ind
		indLoss <- (mat <= neut) & ind
	}
	if (na.is.neut){
	  #mat[is.na(mat)] <- neut
	  indGain[is.na(indGain)] <- FALSE
	  indLoss[is.na(indLoss)] <- FALSE
	}

	gainFreq <- colSums(indGain, na.rm=TRUE) / numSamples	
	lossFreq <- colSums(indLoss, na.rm=TRUE) / numSamples
	return(list(gainFreq=gainFreq, lossFreq=lossFreq))
}

#gain <- 3
#loss <- 1
# relative copy number to neutral = 0
# the above correction subtracts the median CN from all bins/genes
gain <- 1 
loss <- -1
neut <- NULL
if (cnaType == "extreme"){
	gain <- 3
	loss <- -2	
}

## Overall (clonal + subclonal) copy number frequency ##
freqAll <- getFrequency(geneCN, neut = neut, loss=loss, gain=gain, numSamples = nrow(groups), na.is.neut=T)

## copy number filtered by length threshold ##
freqLen <- getFrequency(geneCN, neut = neut, loss=loss, gain=gain, numSamples = nrow(groups), na.is.neut=T, 
    ind = geneLen[samples, ] <= lengthThreshold)



outMat <- data.table(genes[, 1:3], All_loss=freqAll$lossFreq, All_gain=freqAll$gainFreq)
outFile <- paste(outRoot, "_geneFreq.txt", sep="")
## output to text file ##
write.table(format(outMat,digits=4), file=outFile, col.names=T, row.names=F, quote=F, sep="\t")

outMat$pos <- rowMeans(cbind(outMat$Start, outMat$End))

genes.df <- as.data.frame(genes)
genes.df$Chr <- as.character(genes.df$Chr)
genes.gr <- makeGRangesFromDataFrame(genes.df, keep.extra.columns=TRUE)


#################################################
################ PLOT GENE FREQUENCIES ##########

## not required for tRCC small cohort samples- MA changed here

##load GISTIC data to overlay hits over frequency plot
#gistic <- read.table(gisticFile, sep="\t", head=T, stringsAsFactors=F)
#gistic_peaks_only <- gistic[!grepl("CN",gistic[,1]),]
#cn_type <- unlist(lapply(strsplit(gistic_peaks_only$Unique.Name,' '),'[[',1)) #Amplification or Deletion

#Gain_idx <- which(cn_type == "Amplification")
#Loss_idx <- which(cn_type == "Deletion")
#Peak_label <- gsub(" ","",paste(cn_type, gistic_peaks_only$Descriptor,sep="_"))

## load CRPC gene list to annotate to GISTIC peak
#CRPC_gene_list <- read.table(geneFile, sep="\t", head=T)
#regions <- makeGRangesFromDataFrame(CRPC_gene_list, TRUE)

#peak_type <- "Wide.Peak.Limits"


#Peak <- gsub("\\(probes.*$","",gsub(" ", "", gistic_peaks_only[,peak_type]))
#Peak_chr <- unlist(lapply(strsplit(Peak,":"),'[[',1))
#Peak_start <- as.numeric(unlist(lapply(strsplit(sub(".*:", "", Peak),"-"),'[[',1)))
#Peak_end <- as.numeric(unlist(lapply(strsplit(sub(".*:", "", Peak),"-"),'[[',2)))

#Gain_regions.gr <- GRanges(seqnames=Peak_chr[Gain_idx], ranges=IRanges(start = Peak_start[Gain_idx], end =Peak_end[Gain_idx]),type=Peak_label[Gain_idx])
#Loss_regions.gr <- GRanges(seqnames=Peak_chr[Loss_idx], ranges=IRanges(start = Peak_start[Loss_idx], end =Peak_end[Loss_idx]),type=Peak_label[Loss_idx])

#Gain_intersectRegion <- join_overlap_inner(Gain_regions.gr, regions)
#Loss_intersectRegion <- join_overlap_inner(Loss_regions.gr, regions)
#gene_annotation <- rbind(as.data.frame(Gain_intersectRegion), as.data.frame(Loss_intersectRegion))

## output to text file ##
#write.table(gene_annotation,sprintf("GISTIC_%s_with_gene_annotation.txt", peak_type), sep='\t', col.names=T, row.names=F, quote=F)

## save image ##
#save.image(file=outImage)

############################
#find overlap between amp peak & bins
## again not using this for tRCC- MA changed here

#gistic_intersectGainRegion <- as.data.frame(subsetByOverlaps(genes.gr, Gain_regions.gr))
#colnames(gistic_intersectGainRegion)[1:3] <- c("Chr","Start","End")

#find overlap between amp peak related genes and bins
#gistic.gain.gr <- makeGRangesFromDataFrame(gistic_intersectGainRegion, TRUE)

#Gain.overlap <- findOverlaps(gistic.gain.gr, Gain_intersectRegion)
#Gain.matched <- gistic.gain.gr[queryHits(Gain.overlap)]
# Add the metadata 

