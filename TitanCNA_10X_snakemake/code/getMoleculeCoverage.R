# file:   getMoleculeCoverage.R
# author: Gavin Ha, Ph.D.
# institution: Fred Hutchinson Cancer Research Center#
# contact: <gha@fredhutch.org>

# ichorCNA: https://github.com/broadinstitute/ichorCNA
# HMMcopy website: http://compbio.bccrc.ca/software/hmmcopy/ and https://www.bioconductor.org/packages/release/bioc/html/HMMcopy.html
# date:   July 19, 2018
# description: 

library(optparse)

option_list <- list(
  make_option(c("-t", "--tumorBXDir"), type = "character", help = "Path to directory containing tumor bed files for each chromosome containing BX tags."),
  make_option(c("-n", "--normalBXDir"), type = "character", default=NULL, help = "Path to directory containing normal bed files for each chromosome containing BX tags."),
  make_option(c("--minReadsPerBX"), type="integer", default=2, help="Minimum number of reads per barcode. [Default: %default]"),
  make_option(c("--gcWig"), type = "character", help = "Path to GC-content WIG file; Required"),
  make_option(c("--mapWig"), type = "character", default=NULL, help = "Path to mappability score WIG file. Default: [%default]"),
  #make_option(c("--repTimeWig"), type = "character", default=NULL, help ="Path to replication timing WIG file. Default: [%default]"),
  make_option(c("--normalPanel"), type="character", default=NULL, help="Median corrected depth from panel of normals. Default: [%default]"),
  make_option(c("--exons.bed"), type = "character", default=NULL, help = "Path to bed file containing exon regions. Default: [%default]"),
  make_option(c("--likModel"), type="character", default="t", help="Likelihood model to use. \"t\" or \"gaussian\". Use \"gaussian\" for faster runtimes. Default: [%default]"),
  make_option(c("--id"), type = "character", default="test", help = "Patient ID. Default: [%default]"),
  make_option(c("--centromere"), type="character", default=NULL, help = "File containing Centromere locations; if not provided then will use hg19 version from ichorCNA package. Default: [%default]"),
  make_option(c("--rmCentromereFlankLength"), type="numeric", default=1e5, help="Length of region flanking centromere to remove. Default: [%default]"),
  make_option(c("--normal"), type="character", default="0.5", help = "Initial normal contamination; can be more than one value if additional normal initializations are desired. Default: [%default]"),
  make_option(c("--scStates"), type="character", default="NULL", help = "Subclonal states to consider. Default: [%default]"),
  make_option(c("--coverage"), type="numeric", default=NULL, help = "PICARD sequencing coverage. Default: [%default]"),
  make_option(c("--lambda"), type="character", default="NULL", help="Initial Student's t precision; must contain 4 values (e.g. c(1500,1500,1500,1500)); if not provided then will automatically use based on variance of data. Default: [%default]"),
  make_option(c("--lambdaScaleHyperParam"), type="numeric", default=3, help="Hyperparameter (scale) for Gamma prior on Student's-t precision. Default: [%default]"),
  #	make_option(c("--kappa"), type="character", default=50, help="Initial state distribution"),
  make_option(c("--ploidy"), type="character", default="2", help = "Initial tumour ploidy; can be more than one value if additional ploidy initializations are desired. Default: [%default]"),
  make_option(c("--maxCN"), type="numeric", default=7, help = "Total clonal CN states. Default: [%default]"),
  make_option(c("--estimateNormal"), type="logical", default=TRUE, help = "Estimate normal. Default: [%default]"),
  make_option(c("--estimateScPrevalence"), type="logical", default=TRUE, help = "Estimate subclonal prevalence. Default: [%default]"),
  make_option(c("--estimatePloidy"), type="logical", default=TRUE, help = "Estimate tumour ploidy. Default: [%default]"),
  make_option(c("--maxFracCNASubclone"), type="numeric", default=0.7, help="Exclude solutions with fraction of subclonal events greater than this value. Default: [%default]"),
  make_option(c("--maxFracGenomeSubclone"), type="numeric", default=0.5, help="Exclude solutions with subclonal genome fraction greater than this value. Default: [%default]"),
  make_option(c("--minSegmentBins"), type="numeric", default=50, help="Minimum number of bins for largest segment threshold required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction."),
  make_option(c("--altFracThreshold"), type="numeric", default=0.05, help="Minimum proportion of bins altered required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction. Default: [%default]"),
  make_option(c("--genomeBuild"), type="character", default="hg38", help="Geome build. Default: [%default]"),
  make_option(c("--genomeStyle"), type = "character", default="NCBI", help = "Chr naming convention. NCBI (e.g. 1) or UCSC (e.g. chr1). Default: [%default]"),
  make_option(c("--chrNormalize"), type="character", default="c(1:22)", help = "Specify chromosomes to normalize GC/mappability biases. Default: [%default]"),
  make_option(c("--chrTrain"), type="character", default="c(1:22)", help = "Specify chromosomes to estimate params. Default: [%default]"),
  make_option(c("--chrs"), type="character", default="c(1:22,\"X\")", help = "Specify chromosomes to analyze. Default: [%default]"),
  make_option(c("--normalizeMaleX"), type="logical", default=TRUE, help = "If male, then normalize chrX by median. Default: [%default]"),
  make_option(c("--fracReadsInChrYForMale"), type="numeric", default=0.001, help = "Threshold for fraction of reads in chrY to assign as male. Default: [%default]"),
  make_option(c("--includeHOMD"), type="logical", default=TRUE, help="If FALSE, then exclude HOMD state. Useful when using large bins (e.g. 1Mb). Default: [%default]"),
  make_option(c("--txnE"), type="numeric", default=0.9999999, help = "Self-transition probability. Increase to decrease number of segments. Default: [%default]"),
  make_option(c("--txnStrength"), type="numeric", default=1e7, help = "Transition pseudo-counts. Exponent should be the same as the number of decimal places of --txnE. Default: [%default]"),
  make_option(c("--plotFileType"), type="character", default="pdf", help = "File format for output plots. Default: [%default]"),
	make_option(c("--plotYLim"), type="character", default="c(-2,2)", help = "ylim to use for chromosome plots. Default: [%default]"),
  make_option(c("--outDir"), type="character", default="./", help = "Output Directory. Default: [%default]"),
  make_option(c("--libdirTitanCNA"), type = "character", default=NULL, help = "TitanCNA library path. Usually exclude this argument unless custom modifications have been made to the TitanCNA R package code and the user would like to source those R files. Default: [%default]"),
  make_option(c("--libdirIchorCNA"), type = "character", default=NULL, help = "ichorCNA library path. Usually exclude this argument unless custom modifications have been made to the ichorCNA R package code and the user would like to source those R files. Default: [%default]"),
  make_option(c("--cores"), type="numeric", default = 1, help = "Number of cores to use for EM. Default: [%default]")
)
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)
options(scipen=0, stringsAsFactors=F, bitmapType='cairo')

library(GenomeInfoDb)
library(GenomicRanges)
library(data.table)
library(HMMcopy)
library(TitanCNA)
library(foreach)
library(doMC)

registerDoMC()
options(cores = opt$cores)
message("getMoleculeCoverage: Using ", getDoParWorkers(), " cores.")

patientID <- opt$id
tumour_file <- opt$tumorBXDir
normal_file <- opt$normalBXDir
minReadsPerBX <- opt$minReadsPerBX
gcWig <- opt$gcWig
mapWig <- opt$mapWig
#repTimeWig <- opt$repTimeWig
normal_panel <- opt$normalPanel
exons.bed <- opt$exons.bed  # "0" if none specified
likModel <- opt$likModel
centromere <- opt$centromere
flankLength <- opt$rmCentromereFlankLength
normal <- eval(parse(text = opt$normal))
scStates <- eval(parse(text = opt$scStates))
lambda <- eval(parse(text = opt$lambda))
lambdaScaleHyperParam <- opt$lambdaScaleHyperParam
estimateNormal <- opt$estimateNormal
estimatePloidy <- opt$estimatePloidy
estimateScPrevalence <- opt$estimateScPrevalence
maxFracCNASubclone <- opt$maxFracCNASubclone
maxFracGenomeSubclone <- opt$maxFracGenomeSubclone
minSegmentBins <- opt$minSegmentBins
altFracThreshold <- opt$altFracThreshold
ploidy <- eval(parse(text = opt$ploidy))
coverage <- opt$coverage
maxCN <- opt$maxCN
txnE <- opt$txnE
txnStrength <- opt$txnStrength
normalizeMaleX <- as.logical(opt$normalizeMaleX)
includeHOMD <- as.logical(opt$includeHOMD)
fracReadsInChrYForMale <- opt$fracReadsInChrYForMale
outDir <- opt$outDir
outPlotDir <- paste0(outDir, "/", patientID) 
libdirTitanCNA <- opt$libdirTitanCNA
libdirIchorCNA <- opt$libdirIchorCNA
plotFileType <- opt$plotFileType
plotYLim <- eval(parse(text=opt$plotYLim))
gender <- NULL
outImage <- paste0(outDir,"/", patientID,".RData")
## create output directories for each sample ##
dir.create(paste0(outDir, "/"), recursive = TRUE)
dir.create(paste0(outPlotDir, "/"), recursive = TRUE)
  
## set genome style for chromosome names
genomeStyle <- opt$genomeStyle
genomeBuild <- opt$genomeBuild
chrs <- eval(parse(text = opt$chrs));
chrsAll <- c(chrs, "Y")
chrTrain <- as.character(eval(parse(text=opt$chrTrain))); 
chrNormalize <- as.character(eval(parse(text=opt$chrNormalize))); 
seqlevelsStyle(chrs) <- genomeStyle
seqlevelsStyle(chrsAll) <- genomeStyle
seqlevelsStyle(chrNormalize) <- genomeStyle
seqlevelsStyle(chrTrain) <- genomeStyle

maxiter <- 50

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

## load seqinfo 
seqinfo <- getSeqInfo(genomeBuild, genomeStyle, chrs = chrs)

## FILTER BY EXONS IF PROVIDED ##
## add gc and map to RangedData object ##
if (is.null(exons.bed) || exons.bed == "None" || exons.bed == "NULL"){
  targetedSequences <- NULL
}else{
  targetedSequences <- read.delim(exons.bed, header=T, sep="\t")  
}

## load PoN
if (is.null(normal_panel) || normal_panel == "None" || normal_panel == "NULL"){
	normal_panel <- NULL
}

if (is.null(centromere) || centromere == "None" || centromere == "NULL"){ # no centromere file provided
	centromere <- system.file("extdata", "GRCh37.p13_centromere_UCSC-gapTable.txt", 
			package = "ichorCNA")
}
centromere <- read.delim(centromere,header=T,stringsAsFactors=F,sep="\t")

save.image(outImage)

## LOAD GC/MAP/REPTIME WIG FILES ###
message("Reading GC and mappability files")
gc <- wigToGRanges(gcWig)
if (is.null(gc)){
    stop("GC wig file not provided but is required")
}
map <- wigToGRanges(mapWig)
if (is.null(map)){
  message("No mappability wig file input, excluding from correction")
}

## LOAD IN WIG FILES ##
numSamples <- 1
counts <- list()
tumour_copy <- list()
for (i in 1:numSamples) {
  id <- patientID

  ### LOAD TUMOUR AND NORMAL FILES ###
  message("Loading tumour files from ", tumour_file)
  tumour_doc <- loadBXcountsFromBEDDir(tumour_file, chrs = chrsAll, minReads = minReadsPerBX)
  tumour_doc$BX.medianNorm <- log2(tumour_doc$BXcounts / median(tumour_doc$BXcounts, na.rm=T))
  
  ## LOAD GC/MAP WIG FILES ###
  message("Correcting Tumour")
  counts[[id]] <- loadReadCountsFromWig(tumour_doc, chrs = chrs, genomeStyle = genomeStyle,
  									   gc = gc, map = map, repTime = NULL,
                                       centromere = centromere, flankLength = flankLength, 
                                       targetedSequences = targetedSequences, 
                                       chrNormalize = chrNormalize, mapScoreThres = 0.9)
  gender <- counts[[id]]$gender
 	## load in normal file if provided 
 	if (!is.null(normal_file)){
		message("Loading normal files from ", normal_file)
		normal_doc <- loadBXcountsFromBEDDir(normal_file, chrs = chrsAll, minReads = minReadsPerBX)
		normal_doc$BX.medianNorm <- log2(normal_doc$BXcounts / median(normal_doc$BXcounts, na.rm=T))
		message("Correcting Normal")
		counts.normal <- loadReadCountsFromWig(normal_doc, chrs=chrs, gc=gc, map=map, repTime = NULL, 
        genomeStyle = genomeStyle, centromere=centromere, flankLength = flankLength, 
        targetedSequences=targetedSequences, chrNormalize = chrNormalize, mapScoreThres = 0.9)
		normal_copy <- counts.normal$counts #as(counts$counts, "GRanges") 
    counts[[id]]$counts$cor.gc.normal <- counts.normal$counts$cor.gc
    counts[[id]]$counts$cor.map.normal <- counts.normal$counts$cor.map
    counts[[id]]$counts$cor.rep.normal <- counts.normal$counts$cor.rep
		gender.normal <- counts.normal$gender
	}else{
	  normal_copy <- NULL
	}

	### DETERMINE GENDER ###
	## if normal file not given, use chrY, else use chrX
	message("Determining gender...", appendLF = FALSE)
	gender.mismatch <- FALSE
	if (!is.null(normal_copy)){
	  if (gender$gender != gender.normal$gender){ #use tumour # use normal if given
		# check if normal is same gender as tumour
		  gender.mismatch <- TRUE
		}
	}
	message("Gender ", gender$gender)

  ## NORMALIZE GENOME-WIDE BY MATCHED NORMAL ##
  tumour_copy[[id]] <- counts[[id]]$counts 
	tumour_copy[[id]]$tumour.copy <- tumour_copy[[id]]$copy
	tumour_copy[[id]]$tumor.BX.count <- tumour_copy[[id]]$reads
	tumour_copy[[id]]$normal.copy <- normal_copy$copy
	tumour_copy[[id]]$normal.BX.count <- normal_copy$reads
	tumour_copy[[id]]$normal.BX.medianNorm <- normal_copy$BX.medianNorm
	tumour_copy[[id]]$copy <- tumour_copy[[id]]$tumour.copy - normal_copy$copy

	### OUTPUT FILE ###
	### PUTTING TOGETHER THE COLUMNS IN THE OUTPUT ###
	outMat <- as.data.frame(tumour_copy[[id]])
	outMat <- outMat[, c("seqnames","start","end","tumour.copy")]
	colnames(outMat)[1:3] <- c("chr","start","end")
	outFile <- paste0(outDir,"/",id,".BXcounts.txt")
	message(paste("Outputting to:", outFile))
	write.table(outMat, file=outFile, row.names=F, col.names=T, quote=F, sep="\t")

} ## end of for each sample
save.image(outImage)

chrInd <- as.character(seqnames(tumour_copy[[1]])) %in% chrTrain
## get positions that are valid
if (grepl("gauss", likModel, ignore.case = TRUE)) {
  valid <- tumour_copy[[1]]$valid & !is.na(tumour_copy[[1]]$copy)
  if (length(tumour_copy) >= 2){
    for (i in 2:length(tumour_copy)){ 
      valid <- valid & tumour_copy[[i]]$valid & !is.na(tumour_copy[[i]]$copy)
    } 
  }
}else{
  valid <- tumour_copy[[1]]$valid
}


### RUN HMM ###
## store the results for different normal and ploidy solutions ##
ptmTotalSolutions <- proc.time() # start total timer
results <- list()
loglik <- as.data.frame(matrix(NA, nrow = length(normal) * length(ploidy), ncol = 7, 
                 dimnames = list(c(), c("init", "n_est", "phi_est", "BIC", 
                 												"Frac_genome_subclonal", "Frac_CNA_subclonal", "loglik"))))
counter <- 1
compNames <- rep(NA, nrow(loglik))
mainName <- rep(NA, length(normal) * length(ploidy))
#### restart for purity and ploidy values ####
for (n in normal){
  for (p in ploidy){
    if (n == 0.95 & p != 2) {
        next
    }
    logR <- as.data.frame(lapply(tumour_copy, function(x) { x$copy })) # NEED TO EXCLUDE CHR X #
    param <- getDefaultParameters(logR[valid & chrInd, , drop=F], maxCN = maxCN, includeHOMD = includeHOMD, 
                ct.sc=scStates, ploidy = floor(p), e=txnE, e.same = 50, strength=txnStrength, likModel = likModel)
    param$phi_0 <- rep(p, numSamples)
    param$n_0 <- rep(n, numSamples)
    
    ############################################
    ######## CUSTOM PARAMETER SETTINGS #########
    ############################################
  	K <- length(param$ct)
  	logR.var <- 1 / ((apply(logR, 2, sd, na.rm = TRUE) / sqrt(length(param$ct))) ^ 2)
	param$lambda <- matrix(logR.var, length(param$ct))
	param$lambda[param$ct >= 4] <- logR.var
	param$lambda[param$ct == max(param$ct)] <- logR.var / 10
	param$lambda[param$ct %in% c(1,2,3)] <- logR.var
	param$lambda[param$ct %in% c(0)] <- logR.var / 2 		
				
	#############################################
	################ RUN HMM ####################
	#############################################
    hmmResults.cor <- HMMsegment(tumour_copy, valid, dataType = "copy", 
                                 param = param, chrTrain = chrTrain, maxiter = maxiter,
                                 estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
                                 estimateSubclone = estimateScPrevalence, verbose = TRUE)
    
    for (s in 1:numSamples){
			iter <- hmmResults.cor$results$iter
			id <- names(hmmResults.cor$cna)[s]
	
      ## correct integer copy number based on estimated purity and ploidy
      correctedResults <- correctIntegerCN(cn = hmmResults.cor$cna[[s]],
            segs = hmmResults.cor$results$segs[[s]],
            purity = 1 - hmmResults.cor$results$n[s, iter], ploidy = hmmResults.cor$results$phi[s, iter],
            cellPrev = 1 - hmmResults.cor$results$sp[s, iter],
            maxCNtoCorrect.autosomes = maxCN, maxCNtoCorrect.X = maxCN, minPurityToCorrect = 0.03,
            gender = gender$gender, chrs = chrs, correctHOMD = includeHOMD)
      hmmResults.cor$results$segs[[s]] <- correctedResults$segs
      hmmResults.cor$cna[[s]] <- correctedResults$cn

			## convert full diploid solution (of chrs to train) to have 1.0 normal or 0.0 purity
			## check if there is an altered segment that has at least a minimum # of bins
			segsS <- hmmResults.cor$results$segs[[s]]
			segsS <- segsS[segsS$chr %in% chrTrain, ]
			segAltInd <- which(segsS$event != "NEUT")
			maxBinLength = -Inf
			if (sum(segAltInd) > 0){
        maxInd <- which.max(segsS$end[segAltInd] - segsS$start[segAltInd] + 1)
        maxSegRD <- GRanges(seqnames=segsS$chr[segAltInd[maxInd]], 
                  ranges=IRanges(start=segsS$start[segAltInd[maxInd]], end=segsS$end[segAltInd[maxInd]]))
        hits <- findOverlaps(query=maxSegRD, subject=tumour_copy[[s]][valid, ])
        maxBinLength <- length(subjectHits(hits))
      }
			## check if there are proportion of total bins altered 
			# if segment size smaller than minSegmentBins, but altFrac > altFracThreshold, then still estimate TF
			altInd <- hmmResults.cor$cna[[s]][hmmResults.cor$cna[[s]]$chr %in% chrTrain, "event"] == "NEUT"
			altFrac <- sum(!altInd, na.rm=TRUE) / length(altInd)
			if ((maxBinLength <= minSegmentBins) & (altFrac <= altFracThreshold)){
			   hmmResults.cor$results$n[s, iter] <- 1.0
    	}
      ## plot solution ##
      outPlotFile <- paste0(outPlotDir, "/", id, "_genomeWide_", "n", n, "-p", p)
      mainName[counter] <- paste0(id, ", n: ", n, ", p: ", p, ", log likelihood: ", signif(hmmResults.cor$results$loglik[hmmResults.cor$results$iter], digits = 4))
      plotGWSolution(hmmResults.cor, s=s, outPlotFile=outPlotFile, plotFileType=plotFileType, seqinfo = seqinfo,
                     plotYLim=plotYLim, estimateScPrevalence=estimateScPrevalence, main=mainName[counter])
    }
    iter <- hmmResults.cor$results$iter
    results[[counter]] <- hmmResults.cor
    loglik[counter, "loglik"] <- signif(hmmResults.cor$results$loglik[iter], digits = 4)
    subClonalBinCount <- unlist(lapply(hmmResults.cor$cna, function(x){ sum(x$subclone.status) }))
    fracGenomeSub <- subClonalBinCount / unlist(lapply(hmmResults.cor$cna, function(x){ nrow(x) }))
    fracAltSub <- subClonalBinCount / unlist(lapply(hmmResults.cor$cna, function(x){ sum(x$copy.number != 2) }))
    fracAltSub <- lapply(fracAltSub, function(x){if (is.na(x)){0}else{x}})
    loglik[counter, "Frac_genome_subclonal"] <- paste0(signif(fracGenomeSub, digits=2), collapse=",")
    loglik[counter, "Frac_CNA_subclonal"] <- paste0(signif(as.numeric(fracAltSub), digits=2), collapse=",")
    loglik[counter, "init"] <- paste0("n", n, "-p", p)
    loglik[counter, "n_est"] <- paste(signif(hmmResults.cor$results$n[, iter], digits = 2), collapse = ",")
    loglik[counter, "phi_est"] <- paste(signif(hmmResults.cor$results$phi[, iter], digits = 4), collapse = ",")

    counter <- counter + 1
  }
}
## get total time for all solutions ##
elapsedTimeSolutions <- proc.time() - ptmTotalSolutions
message("Total ichorCNA Runtime: ", format(elapsedTimeSolutions[3] / 60, digits = 2), " min.")

### SAVE R IMAGE ###
save.image(outImage)

### SELECT SOLUTION WITH LARGEST LIKELIHOOD ###
if (estimateScPrevalence){ ## sort but excluding solutions with too large % subclonal 
	fracInd <- which(loglik[, "Frac_CNA_subclonal"] <= maxFracCNASubclone & 
						 		   loglik[, "Frac_genome_subclonal"] <= maxFracGenomeSubclone)
	if (length(fracInd) > 0){ ## if there is a solution satisfying % subclonal
		ind <- fracInd[order(loglik[fracInd, "loglik"], decreasing=TRUE)]
	}else{ # otherwise just take largest likelihood
		ind <- order(as.numeric(loglik[, "loglik"]), decreasing=TRUE) 
	}
}else{#sort by likelihood only
  ind <- order(as.numeric(loglik[, "loglik"]), decreasing=TRUE) 
}

#new loop by order of solutions (ind)
outPlotFile <- paste0(outDir, "/", id, "/", id, "_genomeWide_all_sols")
for(i in 1:length(ind)) {
  hmmResults.cor <- results[[ind[i]]]
  turnDevOff <- FALSE
  turnDevOn <- FALSE
  if (i == 1){
  	turnDevOn <- TRUE
  }
  if (i == length(ind)){
  	turnDevOff <- TRUE
  }
  plotGWSolution(hmmResults.cor, s=s, outPlotFile=outPlotFile, plotFileType="pdf", 
                     plotYLim=plotYLim, estimateScPrevalence=estimateScPrevalence, 
                     turnDevOn = turnDevOn, turnDevOff = turnDevOff, main=mainName[ind[i]])
}

hmmResults.cor <- results[[ind[1]]]
hmmResults.cor$results$loglik <- as.data.frame(loglik)
hmmResults.cor$results$gender <- gender$gender
hmmResults.cor$results$chrYCov <- gender$chrYCovRatio
hmmResults.cor$results$chrXMedian <- gender$chrXMedian
hmmResults.cor$results$coverage <- coverage

outputHMM(cna = hmmResults.cor$cna, segs = hmmResults.cor$results$segs, 
                      results = hmmResults.cor$results, patientID = patientID, outDir=outDir)
outFile <- paste0(outDir, "/", patientID, ".params.txt")
outputParametersToFile(hmmResults.cor, file = outFile)

## plot solutions for all samples 
plotSolutions(hmmResults.cor, tumour_copy, chrs, outDir, counts, numSamples=numSamples,
              plotFileType=plotFileType, plotYLim=plotYLim, plotSegs = FALSE,
              estimateScPrevalence=estimateScPrevalence, maxCN=maxCN)

