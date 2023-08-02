title: "Gene Breakpoint Map"
author: "Fiona McBride"
date: "2023-02-15"

library(tidyverse)

## adapted from: https://github.com/greymonroe/genemodel/blob/master/R/genemodel.plot.R

# function to plot gene features, customize the colors of utrs, exons, and borders 

genemodel_plot <- function(model, start, bpstop, orientation, xaxis=TRUE, border_color, exon_color, UTR_color){
  par(mar=c(1,1,3,1), cex=1)
  model<-cbind(model[,1], as.data.frame(stringr::str_split_fixed(model$coordinates, "-", 2)))
  colnames(model)<-c("feature", "start", "bpstop")
  model$start<-as.numeric(as.character(model$start));model$bpstop<-as.numeric(as.character(model$bpstop))
  
  length<-bpstop-start
  
  if (orientation=="reverse"){
    model$newstart<-bpstop-model$bpstop+1
    model$newstop<-bpstop-model$start+1
    model<-model[which(model$feature!="exon"),]
    model<-model[which(model$feature!="ORF"),]
    model<-model[order(model$newstart),]
    
    plot(1, type="l",axes=F,ann=FALSE, xlim=c(start-.03*length, bpstop), ylim=c(-1, .5))
    
    ## modified to make colors flexible
    for (i in 2:nrow(model)){
      type<-model$feature[i]
      if (type=="coding_region"){
        rect(model$newstart[i], 0, model$newstop[i], .2, col = exon_color, border=exon_color, lwd=1)
      } 
      
      if (type=="intron"){
        mid<-mean(c(model$newstart[i], model$newstop[i]))
        segments(x0=model$newstart[i],y0=.1,x1=mid,y1=.2, lwd=1, col=border_color)
        segments(x0=model$newstop[i],y0=.1,x1=mid,y1=.2, lwd=1, col=border_color)
      }
      
      if (type=="utr one"){
        rect(model$newstart[i], 0.08, model$newstop[i], 0.12, col = UTR_color, border=border_color, lwd=1)
      }
      
      if (type=="utr two"){
        rect(model$newstart[i], 0.08, model$newstop[i], 0.12, col = UTR_color, border=border_color, lwd=1)
      }
    }
    
    ## removed x and y assignments
    
    endtype<-model$feature[1]
    if (endtype=="coding_region") {rect(model$newstart[1], 0.08, model$newstop[1], 0.12, col = UTR_color, border=border_color, lwd=1)
    }
    else {rect(model$newstart[1], 0.08, model$newstop[1], 0.12, col = UTR_color, border=border_color, lwd=1)
    }
  }
  
  
  if (orientation=="forward"){
    model$newstart<-start+model$start-1
    model$newstop<-start+model$bpstop-1
    model<-model[which(model$feature!="exon"),]
    model<-model[which(model$feature!="ORF"),]
    model<-model[order(model$newstop, decreasing = T),]
    
    plot(1, type="l",axes=F,ann=FALSE, xlim=c(start, bpstop+.03*length), ylim=c(-1, .5))
    
    ## modified to make colors flexible
    for (i in 2:nrow(model)){
      type<-model$feature[i]
      
      if (type=="coding_region"){
        rect(model$newstart[i], 0, model$newstop[i], .2, col = exon_color, border = exon_color, lwd=1)
      }
      
      if (type=="intron"){
        mid<-mean(c(model$newstart[i], model$newstop[i]))
        segments(x0=model$newstart[i],y0=.1,x1=mid,y1=.2, lwd=1, col=border_color)
        segments(x0=model$newstop[i],y0=.1,x1=mid,y1=.2, lwd=1, col=border_color)
      }
      
      if (type=="utr two"){
        rect(model$newstart[i], 0.08, model$newstop[i], 0.12, col = UTR_color, border=border_color, lwd=1)
      }
      
      if (type=="utr one"){
        rect(model$newstart[i], 0.08, model$newstop[i], 0.12, col = UTR_color, border=border_color, lwd=1)
      }
    }
    
    ## removed x and y assignments
    
    endtype<-model$feature[1]
    if (endtype=="coding_region") {rect(model$newstart[1], 0.08, model$newstop[1], 0.12, col = UTR_color, border=border_color, lwd=1)
    }
    else {rect(model$newstart[1], 0.08, model$newstop[1], 0.12, col = UTR_color, border=border_color, lwd=1)
    }
  }
  
  if (xaxis==T){Axis(side=3, labels=T, cex.axis=0.7)}
}

## adapted from: https://github.com/greymonroe/genemodel/blob/master/R/mutation.plot.R

# function to add markers for breakpoint locations
mutation_plot <- function(start, stop, text="", drop=-0.15, col="red", haplotypes=NULL){
  rect(start, .2, stop, drop+.01*length(haplotypes), col=col, border=col)
  text( stop, drop, text, cex=0.7, col=col, pos=4, offset=0.1+.1*length(haplotypes))
  
  ## modified points drop value
  for (i in 1:length(haplotypes)) points(stop, drop-(i-1)*0.1, col=haplotypes[i], pch=20)
}

# function to format data from a BED file for the plotting scripts above

gene_plot_df <- function(gene_file, gene_ID, DNA_breaks, RNA_breaks){
  
  # read in data
  gene <- read.csv(gene_file, sep='\t', header=FALSE)
  
  # relabel columns 
  gene <- gene %>% 
    rename(
      chromosome = V1,
      chromStart = V2,
      chromEnd = V3,
      name = V4,
      score = V5, 
      strand = V6,
      thickStart = V7,
      thickEnd = V8,
      itemRbg = V9,
      blockCount = V10,
      blockSizes = V11,
      blockStarts = V12
    )
  
  # subset the BED file to only the gene_ID of interest (the most complete gene was selected for these purposes)
  gene <- gene[gene$name == gene_ID, ]
  
  # format the subsetted BED file to extract relevant information (info in the BED file is all in one line, here they're put into lists to more easily perform calculations later)
  exon_size <- as.numeric(gene$blockSizes %>% str_split_1(pattern = ","))
  exon_size <- exon_size[-(length(exon_size))]
  
  exon_start <- as.numeric(gene$blockStarts %>% str_split_1(pattern = ","))
  exon_start <- exon_start[-(length(exon_start))]
  
  ## calculate the start and end positions of the introns and exons based on the exon start positions and sizes
  start_pos <- c(0)
  end_pos <- c()
  
  # calculate start and end of utr 1 (start at 0, end at 0 + size of utr1)
  first_utr <- gene$thickStart - gene$chromStart 
  exon_one = (exon_size[1] - first_utr) # first_utr_end (to) first_utr_end + exon_size[1] - first_utr
  exon_one_end = first_utr + exon_one
  
  end_pos <- append(end_pos, first_utr)
  
  # start and end of exon 1
  start_pos <- append(start_pos, first_utr)
  end_pos <- append(end_pos, exon_one_end)
  
  # calculate the start and end of all introns and exons but the last exon
  for (i in 2:(length(exon_size) - 1)){
    
    # intron
    start_pos <- append(start_pos, end_pos[length(end_pos)])
    end_pos <- append(end_pos, exon_start[i])
    
    # exon
    exon_two_end = exon_start[i] + exon_size[i]
    start_pos <- append(start_pos, exon_start[i])
    end_pos <- append(end_pos, exon_two_end)
  }
  
  # calculate the start and end of the last exon and utr2
  start_pos <- append(start_pos, end_pos[length(end_pos)])
  end_pos <- append(end_pos, exon_start[length(exon_start)])
  
  second_utr <- gene$chromEnd - gene$thickEnd
  last_exon <- exon_size[(length(exon_size))]
  last_exon_length <- last_exon - second_utr
  
  start_pos <- append(start_pos, end_pos[length(end_pos)])
  last_exon_end <- end_pos[length(end_pos)] + last_exon_length
  end_pos <- append(end_pos, last_exon_end)
  
  start_pos <- append(start_pos, end_pos[length(end_pos)])
  second_utr_end <- end_pos[length(end_pos)] + second_utr
  end_pos <- append(end_pos, second_utr_end)
  
  # concatenate the start and end positions with a dash to create coordinates
  coords <- paste(start_pos, end_pos, sep = "-")
  
  # label the start and end coordinates with the feature name
  ## assumes the first object is utr one and the last object is utr two
  regions = c()
  
  for (i in 1:length(coords)){
    if (i == 1){
      regions <- append(regions, "utr one")
    }
    if (i != 1 & i != length(coords) & i %% 2 == 0){
      regions <- append(regions, "coding_region")
    }
    if (i != 1 & i != length(coords) & i %% 2 == 1){
      regions <- append(regions, "intron")
    }
    if (i == length(coords)){
      regions <- append(regions, "utr two")
    }
  }
  
  # create a dataframe of coordinates and feature type for plotting
  plot_gene <- data.frame(type=regions, coordinates=coords)
  
  ## other gene information
  # print if the gene is on the forward or reverse strand (important for understanding the order of the exons on the plot)
  if (gene["strand"] == "+"){
    print("orientation is forward")
  }
  
  if (gene["strand"] == "-"){
    print("orientation is reverse")
  }
  
  # if DNA and RNA breakpoint locations are provided, calculate their position relative to the gene 
  ## gene features had to be plotted starting from 0 as a result of how the plotting function is built, so all of the coordinates have to be adjusted accordingly
  if(missing(DNA_breaks)){
    DNA_break = "No DNA breaks supplied"} 
  
  else{
    DNA_break = c()
    
    for (i in DNA_breaks){
      new_DNA_break = (end_pos[length(end_pos)] - (gene$chromEnd - i)) 
      DNA_break <- append(DNA_break, new_DNA_break)}
  }
  
  if(missing(RNA_breaks)){
    RNA_break = "No RNA breaks supplied"} 
  
  else{
    RNA_break = c()
    
    for (i in RNA_breaks){
      new_RNA_break = end_pos[length(end_pos)] - (gene$chromEnd - i)
      RNA_break <- append(RNA_break, new_RNA_break)}
  }
  
  # print the last position to set the range of the plot
  print(paste("bpstop = ", end_pos[length(end_pos)]))
  
  # return DNA and RNA breaks if applicable, and the df to be fed into the plotting function
  returns <- list("DNA_breaks" = DNA_break, "RNA_breaks" = RNA_break, "gene_df" = plot_gene)
  return(returns)
}


## DTRCC6
# run formatting function to create the df and get information
TFE3_6 <- gene_plot_df("TFE3_whole", "ENST00000315869.8", DNA_breaks = c(49034295, 49034299), RNA_breaks = c(49034251))

# can remove hashtag to export a pdf or svg of the gene plot
#pdf(file = 'TFE3_6.pdf')
#svg(file = 'TFE3_6.svg')

# plot the gene
## bpstop should be a value a little larger than the end value of the gene (bpstop) in order to have enough room to see the whole gene 
## exon in genes that are on the reverse strand but plotted forward will be in descending order
genemodel_plot(TFE3_6$gene_df, start=1, bpstop=16200, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "steelblue3", UTR_color = "lightsteelblue1")
# exon order 10 - 1

# add DNA and RNA lollipop markers based on the DNA and RNA breakpoints calculated in the formatting function
mutation_plot(5526, 5526, text="", col="black", drop=.35, haplotypes=c("red")) # RNA

mutation_plot(5570, 5570, text="", col="black", drop=-.10, haplotypes=c("blue")) # DNA
mutation_plot(5574, 5574, text="", col="black", drop=-.13, haplotypes=c("blue"))

## will need to manually change the axis numbers since the plots were adjusted to start at 0

SFPQ_6 <- gene_plot_df("SFPQ_whole", "ENST00000357214.6", DNA_breaks = c(35181558, 35181562), RNA_breaks = c(35187001))

# because of where the DNA breakpoint was relative to this gene, extra nucleotides needed to be added manually to the df that was generated in the formatting function

SFPQ_increased <- SFPQ_6$gene_df
coords_list <- SFPQ_increased$coordinates

new_coord_list = c()
for (i in 1:length(coords_list)){
  coords <- coords_list[i] %>% str_split_1(pattern = "-") %>% as.numeric(.)
  coords_shifted <- coords + 5000
  
  new_coord_list <- append(new_coord_list, coords_shifted)
}

even <- c()
odd <- c()
for (i in 1:length(new_coord_list)){
  if (i %% 2 == 0){
    even <- append(even, new_coord_list[i])
  }
  if (i %% 2 == 1){
    odd <- append(odd, new_coord_list[i])
  }
}

odd <- c(0, odd)
even <- c(5000, even)
new_coords <- paste(odd, even, sep = "-")

new_type <- c("intron", SFPQ_increased$type)

SFPQ_increased <- data.frame(type=new_type, coordinates=new_coords)

#pdf(file = 'SFPQ_6_manual.pdf')
#svg(file = 'SFPQ_6_manual.svg')

genemodel_plot(model=SFPQ_increased, start=1, bpstop=16000, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "sienna4", UTR_color = "sienna3")

mutation_plot(9066, 9066, text="", col="black", drop=.35, haplotypes=c("red"))

mutation_plot(3623, 3623, text="", col="black", drop=-.10, haplotypes=c("blue"))
mutation_plot(3627, 3627, text="", col="black", drop=-.13, haplotypes=c("blue"))

## DTRCC8

TFE3_8 <- gene_plot_df("TFE3_whole", "ENST00000315869.8", DNA_breaks = c(49036303, 49036296), RNA_breaks = c(49034251))

#pdf(file = 'TFE3_8.pdf')
#svg(file = 'TFE3_8.svg')

genemodel_plot(TFE3_8$gene_df, start=1, bpstop=16200, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "steelblue3", UTR_color = "lightsteelblue1")

mutation_plot(5526, 5526, text="", col="black", drop=.35, haplotypes=c("red"))

mutation_plot(7578, 7578, text="", col="black", drop=-.10, haplotypes=c("blue")) #DNA
mutation_plot(7571, 7571, text="", col="black", drop=-.13, haplotypes=c("blue"))

ASPSCR1_8 <- gene_plot_df("ASPSCR1_whole", "ENST00000306739.9", DNA_breaks=c(82002033, 82001961), RNA_breaks=c(81996846))

#pdf(file = 'ASPSCR1_8.pdf')
#svg(file = 'ASPSCR1_8.svg')

genemodel_plot(model=ASPSCR1_8$gene_df, start=1, bpstop=40000, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "hotpink3", UTR_color = "hotpink1")

mutation_plot(24405, 24405, text="", col="black", drop=-.10, haplotypes=c("blue")) #DNA
mutation_plot(24333, 24333, text="", col="black", drop=-.13, haplotypes=c("blue"))

mutation_plot(19218, 19218, text="", col="black", drop=.35, haplotypes=c("red"))

## DTRCC10

TFE3_10 <- gene_plot_df("TFE3_whole", "ENST00000315869.8", DNA_breaks = c(49037093, 49037298), RNA_breaks = c(49038010, 49034251))

#pdf(file = 'TFE3_10.pdf')
#svg(file = 'TFE3_10.svg')

genemodel_plot(TFE3_10$gene_df, start=1, bpstop=16200, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "steelblue3", UTR_color = "lightsteelblue1")

mutation_plot(5526, 5526, text="", col="black", drop=.35, haplotypes=c("red"))

mutation_plot(8368, 8368, text="", col="black", drop=-.10, haplotypes=c("blue")) #DNA
mutation_plot(8573, 8573, text="", col="black", drop=-.10, haplotypes=c("blue"))

MED15_10 <- gene_plot_df("MED15_whole", "ENST00000263205.11", DNA_breaks=c(20571325, 20571496), RNA_breaks = c(20568631, 20575113))

#pdf(file = 'MED15_10.pdf')
#svg(file = 'MED15_10.svg')

genemodel_plot(model=MED15_10$gene_df, start=1, bpstop=81000, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "firebrick3", UTR_color = "firebrick1")

mutation_plot(61022, 61022, text="", col="black", drop=.35, haplotypes=c("red")) #RNA

mutation_plot(63716, 63716, text="", col="black", drop=-.10, haplotypes=c("blue")) #DNA
mutation_plot(63887, 63887, text="", col="black", drop=-.13, haplotypes=c("blue"))
```

## DTRCC11 

TFE3_11 <- gene_plot_df("TFE3_whole", "ENST00000315869.8", DNA_breaks = c(49039399, 49039398), RNA_breaks = c(49040455))

#pdf(file = 'TFE3_11.pdf')
#svg(file = 'TFE3_11.svg')

genemodel_plot(TFE3_11$gene_df, start=1, bpstop=16200, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "steelblue3", UTR_color = "lightsteelblue1")

mutation_plot(11730, 11730, text="", col="black", drop=.35, haplotypes=c("red")) #RNA

mutation_plot(10674, 10674, text="", col="black", drop=-.10, haplotypes=c("blue")) #DNA
mutation_plot(10673, 10673, text="", col="black", drop=-.15, haplotypes=c("blue"))

PRCC_11 <- gene_plot_df("PRCC_whole", "ENST00000271526.9", DNA_breaks = c(156788672, 156788667), RNA_breaks = c(156791697))

PRCC_11$RNA_breaks

#pdf(file = 'PRCC_11.pdf')
#svg(file = 'PRCC_11.svg')

genemodel_plot(model=PRCC_11$gene_df, start=1, bpstop=33900, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "palegreen4", UTR_color = "palegreen3")

mutation_plot(21138, 21138, text="", col="black", drop=-.10, haplotypes=c("blue")) #DNA
mutation_plot(21133, 21133, text="", col="black", drop=-.15, haplotypes=c("blue"))

mutation_plot(24163, 24163, text="", col="black", drop=.35, haplotypes=c("red"))


## DTRCC12_2018

TFE3_2018 <- gene_plot_df("TFE3_whole", "ENST00000315869.8", DNA_breaks = c(49038227, 49038228), RNA_breaks = c(49038114))

#pdf(file = 'TFE3_2018.pdf')
#svg(file = 'TFE3_2018.svg')

genemodel_plot(TFE3_2018$gene_df, start=1, bpstop=16200, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "steelblue3", UTR_color = "lightsteelblue1")

mutation_plot(9389, 9389, text="", col="black", drop=.35, haplotypes=c("red")) #RNA

mutation_plot(9502, 9502, text="", col="black", drop=-.10, haplotypes=c("blue")) #DNA
mutation_plot(9503, 9503, text="", col="black", drop=-.13, haplotypes=c("blue"))

PRCC_2018 <- gene_plot_df("PRCC_whole", "ENST00000271526.9", DNA_breaks = c(156774403, 156774440), RNA_breaks = c(156768239))

#pdf(file = 'PRCC_2018.pdf')
#svg(file = 'PRCC_2018.svg')

genemodel_plot(model=PRCC_2018$gene_df, start=1, bpstop=33900, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "palegreen4", UTR_color = "palegreen3")

mutation_plot(6869, 6869, text="", col="black", drop=-.10, haplotypes=c("blue")) #DNA
mutation_plot(6906, 6906, text="", col="black", drop=-.13, haplotypes=c("blue"))

mutation_plot(705, 705, text="", col="black", drop=.35, haplotypes=c("red"))


## DTRCC3

TFE3_3 <- gene_plot_df("TFE3_whole", "ENST00000315869.8", DNA_breaks = c(49038126, 49038126))

#pdf(file = "TFE3_3.pdf")
#svg(file = "TFE3_3.svg")

genemodel_plot(TFE3_3$gene_df, start=1, bpstop=16200, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "steelblue3", UTR_color = "lightsteelblue1")

mutation_plot(9401, 9401, text="", col="black", drop=-.13, haplotypes=c("blue"))
mutation_plot(9401, 9401, text="", col="black", drop=-.10, haplotypes=c("blue"))

PRCC_3 <- gene_plot_df("PRCC_whole", "ENST00000271526.9", DNA_breaks = c(156788343, 156788331))
PRCC_3$DNA_breaks

#pdf(file = "PRCC_3.pdf")
#svg(file = "PRCC_3.svg")

genemodel_plot(model=PRCC_3$gene_df, start=1, bpstop=33900, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "palegreen4", UTR_color = "palegreen3")

mutation_plot(20809, 20809, text="", col="black", drop=-.10, haplotypes=c("blue")) #DNA
mutation_plot(20797, 20797, text="", col="black", drop=-.13, haplotypes=c("blue"))



## DTRCC4

TFE3_4 <- gene_plot_df("TFE3_whole", "ENST00000315869.8", DNA_breaks = c(49037882, 49037901))

#pdf(file = "TFE3_4.pdf")
#svg(file = "TFE3_4.svg")

genemodel_plot(TFE3_4$gene_df, start=1, bpstop=16200, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "steelblue3", UTR_color = "lightsteelblue1")

mutation_plot(9157, 9157, text="", col="black", drop=-.13, haplotypes=c("blue"))
mutation_plot(9176, 9176, text="", col="black", drop=-.10, haplotypes=c("blue"))

NONO_4 <- gene_plot_df("NONO_whole", "ENST00000373841.5", DNA_breaks = c(71299453, 71299455))

NONO_4$gene_df <- NONO_4$gene_df[-c(1, 2),]
NONO_4$gene_df <- rbind(data.frame(type = "utr one", coordinates = "0-89"), NONO_4$gene_df)

#pdf(file = 'NONO_4.pdf')
#svg(file = 'NONO_4.svg')

genemodel_plot(model=NONO_4$gene_df, start=1, bpstop=20000, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "mediumpurple4", UTR_color = "mediumpurple1")

mutation_plot(15821, 15821, text="", col="black", drop=-.13, haplotypes=c("blue"))
mutation_plot(15823, 15823, text="", col="black", drop=-.10, haplotypes=c("blue"))


## DTRCC5

TFE3_5 <- gene_plot_df("TFE3_whole", "ENST00000315869.8", DNA_breaks = c(49038152, 49038151), RNA_breaks = c(49038197))

#pdf(file = "TFE3_5.pdf")
#svg(file = "TFE3_5.svg")

genemodel_plot(TFE3_5$gene_df, start=1, bpstop=16200, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "steelblue3", UTR_color = "lightsteelblue1")

mutation_plot(9472, 9472, text="", col="black", drop=.35, haplotypes=c("red"))

mutation_plot(9427, 9427, text="", col="black", drop=-.13, haplotypes=c("blue"))
mutation_plot(9426, 9426, text="", col="black", drop=-.10, haplotypes=c("blue"))

RBM10_5 <- gene_plot_df("RBM10_whole", "ENST00000329236.8", DNA_breaks = c(47183719, 47183719), RNA_breaks = c(47185055))

#pdf(file = 'RBM10_5.pdf')
#svg(file = 'RBM10_5.svg')

genemodel_plot(model=RBM10_5$gene_df, start=1, bpstop=45000, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "goldenrod", UTR_color = "khaki")

mutation_plot(39609, 39609, text="", col="black", drop=.35, haplotypes=c("red"))

mutation_plot(38273, 38273, text="", col="black", drop=-.13, haplotypes=c("blue"))
mutation_plot(38273, 38273, text="", col="black", drop=-.10, haplotypes=c("blue"))


## DTRCC7

TFE3_7 <- gene_plot_df("TFE3_whole", "ENST00000315869.8", DNA_breaks = c(49038179, 49038180))

#pdf(file = "TFE3_7.pdf")
#svg(file = "TFE3_7.svg")

genemodel_plot(TFE3_7$gene_df, start=1, bpstop=16200, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "steelblue3", UTR_color = "lightsteelblue1")

mutation_plot(9454, 9454, text="", col="black", drop=-.13, haplotypes=c("blue"))
mutation_plot(9455, 9455, text="", col="black", drop=-.10, haplotypes=c("blue"))

RBM10_7 <- gene_plot_df("RBM10_whole", "ENST00000329236.8", DNA_breaks = c(47183129, 47183133))

#pdf(file = 'RBM10_7.pdf')
#svg(file = 'RBM10_7.svg')

genemodel_plot(model=RBM10_7$gene_df, start=1, bpstop=43000, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "goldenrod", UTR_color = "khaki")

mutation_plot(37683, 37683, text="", col="black", drop=-.13, haplotypes=c("blue"))
mutation_plot(37687, 37687, text="", col="black", drop=-.10, haplotypes=c("blue"))

## DTRCC9

TFE3_9 <- gene_plot_df("TFE3_whole", "ENST00000315869.8", DNA_breaks = c(49043616, 49043612))

#pdf(file = "TFE3_9.pdf")
#svg(file = "TFE3_9.svg")

genemodel_plot(TFE3_9$gene_df, start=1, bpstop=16200, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "steelblue3", UTR_color = "lightsteelblue1")

mutation_plot(14891, 14891, text="", col="black", drop=-.13, haplotypes=c("blue"))
mutation_plot(14887, 14887, text="", col="black", drop=-.10, haplotypes=c("blue"))

SFPQ_9 <- gene_plot_df("SFPQ_whole", "ENST00000357214.6", DNA_breaks = c(35188059, 35188070))

#pdf(file = 'SFPQ_9.pdf')
#svg(file = 'SFPQ_9.svg')

genemodel_plot(model=SFPQ_9$gene_df, start=1, bpstop=11000, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "sienna4", UTR_color = "sienna3")

mutation_plot(5124, 5124, text="", col="black", drop=-.13, haplotypes=c("blue"))
mutation_plot(5135, 5135, text="", col="black", drop=-.10, haplotypes=c("blue"))

## DTRCC13

TFE3_13 <- gene_plot_df("TFE3_whole", "ENST00000315869.8", DNA_breaks = c(49035276, 49035276))

#pdf(file = "TFE3_13.pdf")
#svg(file = "TFE3_13.svg")

genemodel_plot(TFE3_13$gene_df, start=1, bpstop=16200, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "steelblue3", UTR_color = "lightsteelblue1")

mutation_plot(6551, 6551, text="", col="black", drop=-.13, haplotypes=c("blue"))
mutation_plot(6551, 6551, text="", col="black", drop=-.10, haplotypes=c("blue"))

MED15_13 <- gene_plot_df("MED15_whole", "ENST00000263205.11", DNA_breaks = c(20572549, 20572567))

#pdf(file = 'MED15_13.pdf')
#svg(file = 'MED15_13.svg')

genemodel_plot(model=MED15_13$gene_df, start=1, bpstop=81000, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "firebrick3", UTR_color = "firebrick1")

mutation_plot(64940, 64940, text="", col="black", drop=-.13, haplotypes=c("blue"))
mutation_plot(64958, 64958, text="", col="black", drop=-.10, haplotypes=c("blue")) 

## DTRCC14

TFE3_14 <- gene_plot_df("TFE3_whole", "ENST00000315869.8", DNA_breaks = c(49039387, 49039395))

#pdf(file = "TFE3_14.pdf")
#svg(file = "TFE3_14.svg")

genemodel_plot(TFE3_14$gene_df, start=1, bpstop=16200, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "steelblue3", UTR_color = "lightsteelblue1")

mutation_plot(10662, 10662, text="", col="black", drop=-.13, haplotypes=c("blue"))
mutation_plot(10670, 10670, text="", col="black", drop=-.10, haplotypes=c("blue"))

LUC7L3_14 <- gene_plot_df("LUC7L3_whole", "ENST00000505658.6", DNA_breaks = c(50746816, 50746816))

#pdf(file = "LUC7L3_14.pdf")
#svg(file = "LUC7L3_14.svg")

genemodel_plot(LUC7L3_14$gene_df, start=1, bpstop=40000, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "darkcyan", UTR_color = "paleturquoise")

mutation_plot(27214, 27214, text="", col="black", drop=-.13, haplotypes=c("blue"))
mutation_plot(27214, 27214, text="", col="black", drop=-.10, haplotypes=c("blue"))


## DTRCC15

TFE3_15 <- gene_plot_df("TFE3_whole", "ENST00000315869.8", DNA_breaks = c(49038924, 49038924))

#pdf(file = "TFE3_15.pdf")
#svg(file = "TFE3_15.svg")

genemodel_plot(TFE3_15$gene_df, start=1, bpstop=16200, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "steelblue3", UTR_color = "lightsteelblue1")

mutation_plot(10199, 10199, text="", col="black", drop=-.13, haplotypes=c("blue"))
mutation_plot(10199, 10199, text="", col="black", drop=-.10, haplotypes=c("blue"))

PRCC_15 <- gene_plot_df("PRCC_whole", "ENST00000271526.9", DNA_breaks = c(156792390, 156792396))

#pdf(file = "PRCC_15.pdf")
#svg(file = "PRCC_15.svg")

genemodel_plot(model=PRCC_15$gene_df, start=1, bpstop=33900, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "palegreen4", UTR_color = "palegreen3")

mutation_plot(24856, 24856, text="", col="black", drop=-.10, haplotypes=c("blue")) #DNA
mutation_plot(24862, 24862, text="", col="black", drop=-.13, haplotypes=c("blue"))


## DTRCC17

TFE3_17 <- gene_plot_df("TFE3_whole", "ENST00000315869.8", DNA_breaks = c(49038453, 49038501))

#pdf(file = "TFE3_17.pdf")
#svg(file = "TFE3_17.svg")

genemodel_plot(TFE3_17$gene_df, start=1, bpstop=16200, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "steelblue3", UTR_color = "lightsteelblue1")

mutation_plot(9728, 9728, text="", col="black", drop=-.13, haplotypes=c("blue"))
mutation_plot(9776, 9776, text="", col="black", drop=-.10, haplotypes=c("blue"))

PRCC_17 <- gene_plot_df("PRCC_whole", "ENST00000271526.9", DNA_breaks = c(156773675, 156773675))

#pdf(file = "PRCC_17.pdf")
#svg(file = "PRCC_17.svg")

genemodel_plot(model=PRCC_17$gene_df, start=1, bpstop=33900, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "palegreen4", UTR_color = "palegreen3")

mutation_plot(6141, 6141, text="", col="black", drop=-.10, haplotypes=c("blue")) #DNA
mutation_plot(6141, 6141, text="", col="black", drop=-.13, haplotypes=c("blue"))


## DTRCC18

TFE3_18 <- gene_plot_df("TFE3_whole", "ENST00000315869.8", DNA_breaks = c(49051451, 49051452))

#pdf(file = "TFE3_18.pdf")
#svg(file = "TFE3_18.svg")

genemodel_plot(TFE3_18$gene_df, start=1, bpstop=25000, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "steelblue3", UTR_color = "lightsteelblue1")

mutation_plot(22726, 22726, text="", col="black", drop=-.13, haplotypes=c("blue"))
mutation_plot(22727, 22727, text="", col="black", drop=-.10, haplotypes=c("blue"))

SFPQ_18 <- gene_plot_df("SFPQ_whole", "ENST00000357214.6", DNA_breaks = c(35188009, 35188032))

#pdf(file = 'SFPQ_18.pdf')
#svg(file = 'SFPQ_18.svg')

genemodel_plot(model=SFPQ_18$gene_df, start=1, bpstop=11000, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "sienna4", UTR_color = "sienna3")

mutation_plot(5074, 5074, text="", col="black", drop=-.13, haplotypes=c("blue"))
mutation_plot(5097, 5097, text="", col="black", drop=-.10, haplotypes=c("blue"))


## joint TFE3

TFE3_joint <- gene_plot_df("TFE3_whole", "ENST00000315869.8", DNA_breaks = c(49038126, 49037882, 49038152, 49034295, 49038179, 49036296, 49043612, 49037065, 49039398, 49038228, 49035276, 49039387, 49038924, 49038453, 49051451))

#pdf(file = "joint_TFE3_labeled.pdf")

# labeled with the sample each marker came from

genemodel_plot(TFE3_joint$gene_df, start=1, bpstop=25000, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "steelblue3", UTR_color = "lightsteelblue1")

mutation_plot(9401, 9401, text="3", col="darkred", drop=-.18, haplotypes=c("white"))

mutation_plot(9157, 9157, text="4", col="darkred", drop=-.05, haplotypes=c("white"))

mutation_plot(9427, 9427, text="5", col="darkred", drop=.30, haplotypes=c("white"))

mutation_plot(5570, 5570, text="6", col="darkred", drop=-.05, haplotypes=c("white")) 

mutation_plot(9454, 9454, text="7", col="darkred", drop=-.12, haplotypes=c("white"))

mutation_plot(7571, 7571, text="8", col="darkred", drop=-.05, haplotypes=c("white")) 

mutation_plot(14887, 14887, text="9", col="darkred", drop=-.10, haplotypes=c("white"))

mutation_plot(8340, 8340, text="10", col="darkred", drop=-.10, haplotypes=c("white"))

mutation_plot(10673, 10673, text="11", col="darkred", drop=-.05, haplotypes=c("white"))

mutation_plot(9503, 9503, text="12", col="darkred", drop=.35, haplotypes=c("white")) 

mutation_plot(6551, 6551, text="13", col="darkred", drop=-.10, haplotypes=c("white")) 

mutation_plot(10662, 10662, text="14", col="darkred", drop=-.10, haplotypes=c("white")) 

mutation_plot(10199, 10199, text="15", col="darkred", drop=-.13, haplotypes=c("white")) 

mutation_plot(9728, 9728, text="17", col="darkred", drop=-.05, haplotypes=c("white")) 

mutation_plot(22726, 22726, text="18", col="darkred", drop=-.10, haplotypes=c("white")) 


#pdf(file = "joint_TFE3.pdf")
#svg(file = "joint_TFE3.svg")

# plot with markers but no labels 

genemodel_plot(TFE3_joint$gene_df, start=1, bpstop=25000, orientation="forward", xaxis=T, border_color = "grey40", exon_color = "steelblue3", UTR_color = "lightsteelblue1")

mutation_plot(9401, 9401, text="", col="darkred", drop=-.16, haplotypes=c("white"))

mutation_plot(9157, 9157, text="", col="darkred", drop=-.05, haplotypes=c("white"))

mutation_plot(9427, 9427, text="", col="darkred", drop=-.16, haplotypes=c("white"))

mutation_plot(5570, 5570, text="", col="darkred", drop=-.05, haplotypes=c("white")) 

mutation_plot(9454, 9454, text="", col="darkred", drop=-.16, haplotypes=c("white"))

mutation_plot(7571, 7571, text="", col="darkred", drop=-.05, haplotypes=c("white")) 

mutation_plot(14887, 14887, text="", col="darkred", drop=-.10, haplotypes=c("white"))

mutation_plot(8340, 8340, text="", col="darkred", drop=-.10, haplotypes=c("white"))

mutation_plot(10673, 10673, text="", col="darkred", drop=-.07, haplotypes=c("white"))

mutation_plot(9503, 9503, text="", col="darkred", drop=-.16, haplotypes=c("white")) 

mutation_plot(6551, 6551, text="", col="darkred", drop=-.10, haplotypes=c("white")) 

mutation_plot(10662, 10662, text="", col="darkred", drop=-.07, haplotypes=c("white")) 

mutation_plot(10199, 10199, text="", col="darkred", drop=-.13, haplotypes=c("white")) 

mutation_plot(9728, 9728, text="", col="darkred", drop=-.05, haplotypes=c("white")) 

mutation_plot(22726, 22726, text="", col="darkred", drop=-.10, haplotypes=c("white")) 
