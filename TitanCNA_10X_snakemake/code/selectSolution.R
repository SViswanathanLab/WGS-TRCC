#' selectSolutions.R
#' author: Gavin Ha 
#' Fred Hutchinson Cancer Research Center
#' contact: <gha@fredhutch.org
#' date:  July 24, 2018

library(optparse)

option_list <- list(
  make_option(c("--ploidyRun2"), type = "character", help = "Directory containing TitanCNA results for initialized ploidy=2. [Required]"),
  make_option(c("--ploidyRun3"), type = "character", default = 0, help = "Directory containing TitanCNA results for initialized ploidy=3. [Required if --ploidyRun4 not given]"),
  make_option(c("--ploidyRun4"), type = "character", default = 0, help = "Directory containing TitanCNA results for initialized ploidy=4. [Required if --ploidyRun3 not given]"),
  make_option(c("--threshold"), type = "numeric", default=0.05, help = "Proportion ploidyRun2 likelihood greater than than ploidyRun3/4 by at least this value. [Default %default]"),
  make_option(c("--outFile"), type = "character", help = "Output file containing a list of solutions chosen for all samples. [Required]")
  #make_option(c("--outDir"), type = "character", help = "Output directory containing a list of solutions chosen for all samples. [Required]")
)

parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
print(opt)

library(stringr)
library(data.table)

phi2Dir <- opt$ploidyRun2
phi3Dir <- opt$ploidyRun3
phi4Dir <- opt$ploidyRun4
threshold <- opt$threshold
outFile <- opt$outFile


homd.snp.thres <- 1500
homd.snp.prop <- 0.01

phi2Files <- list.files(phi2Dir, pattern="params.txt", full.names=T)
if (!is.null(phi3Dir) && phi3Dir != "NULL"){
  phi3Files <- list.files(phi3Dir, pattern="params.txt", full.names=T)
}else{
  phi3Files <- NULL
}
if (!is.null(phi4Dir) && phi4Dir != "NULL"){
  phi4Files <- list.files(phi4Dir, pattern="params.txt", full.names=T)
}else{
  phi4Files <- NULL
}
if (is.null(phi3Files) && is.null(phi4Files)){
  stop("phi3Dir or phi4Dir must be provided.")
}
samples <-  gsub(".params.txt", "", basename(phi2Files))
numSamples <- length(samples)
patients <- unique(str_match(samples, "(.+)_cluster")[,2])
numPatients <- length(patients)

########################################################
##### FUNCTION TO FORMAT PARAMS FOR OUTPUT #############
########################################################
formatParams <- function(params){
	id <- colnames(params)
	barcode <- strsplit(id, "_cluster")[[1]][1]
	cellPrev <- strsplit(params[grepl("Clonal cluster cellular prevalence", 
	  rownames(params)), 1], " ")[[1]]
	numClust <- length(cellPrev)
	cellPrev <- paste0(format(as.numeric(cellPrev), digits=4, nsmall=4), collapse=",")
	norm <- as.numeric(params[grepl("Normal contamination estimate", rownames(params)), 1])
	purity <- 1 - norm
	ploidy <- as.numeric(params[grepl("Average tumour ploidy estimate", rownames(params)), 1])
	loglik <- as.numeric(params[grepl("likelihood", rownames(params)), 1])
	sdbw <- as.numeric(params[grepl("S_Dbw validity index \\(Both\\)", rownames(params)), 1])
	return(list(id=id, barcode=barcode, numClust=numClust, cellPrev=cellPrev, 
	  purity=purity, norm=norm, ploidy=ploidy, loglik=loglik, sdbw=sdbw))
}
#save.image()

getParamAllClusters <- function(phiSamples, phiStr = "2"){
	phiParams <- NULL#data.frame(id=NA, barcode=NA, numClust=NA, 
	  #cellPrev=NA, purity=NA, norm=NA, ploidy=NA, loglik=NA, sdbw=NA)
	for (j in 1:length(phiSamples)){
	  phi <- read.delim(phiSamples[j], header=F, row.names=1, stringsAsFactors=F, sep="\t")
	  colnames(phi) <- gsub(".params.txt", "", basename(phiSamples[j]))	
	  param.df <- data.frame(formatParams(phi), stringsAsFactors=F)  
	  phisegs <- fread(gsub(".params.txt", ".segs.txt", phiSamples[j]))
	  stateDist <- phisegs[, sum(Length.snp.), by = Corrected_Copy_Number]
	  # HOMD total proportion greater than homd.snp.prop
	  if (0 %in% stateDist$Corrected_Copy_Number){
		  ind <- stateDist[Corrected_Copy_Number == 0, V1] / sum(stateDist[, V1]) > homd.snp.prop
		  # HOMD seg larger than homd.snp.thres
		  ind <- ind | sum(phisegs[Corrected_Copy_Number == 0, Length.snp.] > homd.snp.thres) > 0
		  if (ind){
			param.df[, "loglik"] <- -Inf
			param.df[, "sdbw"] <- Inf
		  }
    	}
    	phiParams <- rbind(phiParams, param.df )
    }
	#phiParams <- na.omit(phiParams)
	phiParams <- cbind(Phi=phiStr, phiParams, path = gsub(".params.txt", "", phiSamples))
	return(phiParams)
}

########################################################
################### MAIN FUNCTION ######################
########################################################
#optSolution <- matrix(NA, ncol = 9, nrow = numSamples, dimnames = list(c(), c("SampleID_cluster", "SampleID_Barcode", "NumClonalClusters", "CellularPrevalence", "Purity", "NormalContam", "Ploidy", "LogLik", "Selected.Solution")))
#outLink <- paste0(basename(outDir), ".sh")
#fc <- file(outLink, "w+")
optSolutionAll <- data.table()
for (i in 1:numPatients){
  id <- patients[i]
  phi2Samples <- grep(paste0(id, "_cluster"), phi2Files, value=T)
  if (length(phi2Samples) > 0){
    phi2Params <- getParamAllClusters(phi2Samples, "2")
  }
	
  phi3Params <- NULL
  phi3Params$loglik <- NA
  if (!is.null(phi3Files)){
    phi3Samples <- grep(paste0(id, "_cluster"), phi3Files, value=T)	
    if (length(phi3Samples) > 0){
      phi3Params <- getParamAllClusters(phi3Samples, "3")
    }
  }
	
	phi4Params <- NULL
	phi4Params$loglik <- NA
  if (!is.null(phi4Files)){
    phi4Samples <- grep(paste0(id, "_cluster"), phi4Files, value=T)
    if (length(phi4Samples) > 0){
	    phi4Params <- getParamAllClusters(phi4Samples, "4")
    }
  }

  ## select ploidy based on which set of ploidy solutions has consistently lower loglik
  maxInd <- apply(cbind(phi2Params$loglik + (abs(phi2Params$loglik)*threshold), phi3Params$loglik, phi4Params$loglik), 1, which.max)
  ploidySolInd <- as.numeric(names(which.max(table(maxInd))))
  optPloidy <- switch(ploidySolInd, phi2Params, phi3Params, phi4Params)
  optSolution <- optPloidy[which.min(optPloidy$sdbw), ]
  optSolutionAll <- rbind(optSolutionAll, optSolution)
  
}
write.table(optSolutionAll, file = outFile, col.names = T, row.names = F, quote = F, sep = "\t")



