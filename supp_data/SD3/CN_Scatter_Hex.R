title: "CN_scatter_hex"
author: "Fiona McBride"
date: "2023-05-22"


# load libraries
library(ggplot2)
library(tidyverse)
library(dplyr)
library(he)
library(ggpubr)
library(gtools)
library(gridExtra)

# calculate allelic copy number from the corrected ratio

clean_file_corr <- function(file_path){
  file_df <- read.csv(file_path, header = T, sep = "\t")
  
  file_df$CN_a <- (file_df$Corrected_Ratio * file_df$logR_Copy_Number)
  file_df$CN_b <- (1 - file_df$Corrected_Ratio) * file_df$logR_Copy_Number
  
  subset_df <- file_df[, c("Chr", "Position", "Start", "End", "CN_a", "CN_b")]
  clean_df <- na.omit(subset_df)
  
  return(clean_df)
}

# calculate allelic copy number from the raw allelic ratios

clean_file_noCorr <- function(file_path, ploidy_num){
  file_df <- read.csv(file_path, header = T, sep = "\t")
  
  file_df$CN_a <- (file_df$AllelicRatio * (2^file_df$LogRatio) * ploidy_num)
  file_df$CN_b <- (1 - file_df$AllelicRatio) * (2^file_df$LogRatio) * ploidy_num
  
  subset_df <- file_df[, c("Chr", "Position", "Start", "End", "CN_a", "CN_b")]
  clean_df <- na.omit(subset_df)
  
  return(clean_df)
}

# make hexbin plots of CN_a and CN_b per chromosome

hexbin_plot <- function(chrom_sub, chromosome, upper_axis_limit){
  
  plot <- ggplot(chrom_sub, aes(CN_a, CN_b)) +  
    scale_fill_gradient(low = "white", high = "black") +
    stat_binhex(show.legend = F, bins = 100) + theme_classic() +
    ggtitle(chromosome) + xlab("CN_a") + ylab("CN_b") +
    xlim(0, upper_axis_limit) + ylim(0, upper_axis_limit)
  
  return(plot)
}

# make scatter plots of CN_a and CN_b per chromosome

scatter_plot <- function(chrom_sub, chromosome, upper_axis_limit){
  
  plot <- ggplot(chrom_sub,aes(CN_a, CN_b)) +  
    geom_point(size = 0.5) + theme_classic() +
    ggtitle(chromosome) + xlab("CN_a") + ylab("CN_b") +
    xlim(0, upper_axis_limit) + ylim(0, upper_axis_limit)
  
  
  return(plot)
}

# arrange the plots in a grid formation by chromosome

plot_arrangement <- function(df, plot_function, upper_axis_limit){
  
  plots <- list()
  
  plots[[1]] <- plot_function(df, "Genome", upper_axis_limit)
  
  chroms <- mixedsort(unique(df$Chr))
  
  for(chr in chroms){
    chrom_sub <- df[df$Chr == chr, ]
    plots[[which(chroms == chr) + 1]] <- plot_function(chrom_sub, chr, upper_axis_limit)
  }
  
  return(do.call(grid.arrange, plots))
}

# function to average allelic copy number data in 100kb bins

bin_100kb <- function(df){
  a_10 <- c()
  b_10 <- c()
  chrom <- c()
  start <- c()
  end <- c()
  
  for(chromosome in unique(df$Chr)){
    current_loc <- c(0)
    
    subset <- df[df$Chr == chromosome, ]
    
    while (sum(current_loc) < nrow(subset)){
      kb_100 <- subset[subset$Start >= subset[1+ sum(current_loc), ]$Start & 
                         subset$Start <= subset[1 + sum(current_loc), ]$Start + 100000, ]
      
      avg_a <- sum(kb_100$CN_a)/nrow(kb_100)
      avg_b <- sum(kb_100$CN_b)/nrow(kb_100)
      
      a_10 <- append(a_10, avg_a)
      b_10 <- append(b_10, avg_b)
      
      start <- append(start, subset[1+ sum(current_loc), ]$Start)
      end <- append(end, subset[1 + sum(current_loc), ]$Start + 100000)
      
      current_loc <- append(current_loc, nrow(kb_100))
    }
    chrom <- append(chrom, rep(chromosome, length(current_loc) -1))
  }
  
  bin_df <- data.frame("CN_a" = a_10, "CN_b" = b_10, "Chr" = chrom, "Start" = start, "End" = end)
  
  return(bin_df)
}




## Example with DTRCC11_Rkidney

### ploidy 2 corrected

# path to titan cna file
## metastatic sites in DTRCC18 are in a separate subfolder

TRCC13_2_corr <- clean_file_corr("~/Ploidy 2/DTRCC_13_LLN_cluster1.titan.ichor.cna.txt")

# path where you want the image to go
#pdf("CN_scatter_hex/Hex Plots/DTRCC_13_LLN/DTRCC13_corr2_hex.pdf", height = 15, width = 15)

## upper axis limit set to 6 for ploidy 4
## plotting options are scatter_plot and hexbin_plot
DTRCC13_corr2_hex <- plot_arrangement(TRCC13_2_corr, hexbin_plot, 4)


### ploidy 2 corrected binned
# bin the data
TRCC13_2_corr_bin <- bin_100kb(TRCC13_2_corr)
head(TRCC)

# path for the image
#pdf("CN_scatter_hex/Hex Plots/DTRCC_13_LLN/DTRCC13_corr2_bin10_hex.pdf", height = 15, width = 15)

## upper axis limit set to 6 for ploidy 4
## plotting options are scatter_plot and hexbin_plot
DTRCC13_corr2_bin10_hex <- plot_arrangement(TRCC13_2_corr_bin, hexbin_plot, 4)


### ploidy 2 not corrected
TRCC13_2_noCorr <- clean_file_noCorr("~/Ploidy 2/DTRCC_13_LLN_cluster1.titan.ichor.cna.txt", 2)

#pdf("CN_scatter_hex/Hex Plots/DTRCC_13_LLN/DTRCC13_noCorr2_hex.pdf", height = 15, width = 15)

DTRCC13_noCorr2_hex <- plot_arrangement(TRCC13_2_noCorr, hexbin_plot, 4)


## ploidy 2 not corrected binned
TRCC13_2_noCorr_bin <- bin_100kb(TRCC13_2_noCorr)

#pdf("CN_scatter_hex/Hex Plots/DTRCC13_LLN/DTRCC13_noCorr4_bin10_hex.pdf")
DTRCC13_noCorr2_bin10_hex <- plot_arrangement(TRCC13_2_noCorr_bin, hexbin_plot, 4)



### ploidy 4 corrected
TRCC13_4_corr <- clean_file_corr("~/Ploidy 4/DTRCC_13_LLN_cluster1.titan.ichor.cna.txt")

#pdf("CN_scatter_hex/Hex Plots/DTRCC13_LLN/DTRCC13_corr4_hex.pdf")
DTRCC13_corr4_hex <- plot_arrangement(TRCC13_4_corr, hexbin_plot, 6)


### ploidy 4 corrected binned
TRCC13_corr4_bin10 <- bin_100kb(TRCC13_4_corr)

#pdf("CN_scatter_hex/Hex Plots/DTRCC13_LLN/DTRCC13_corr4_bin10_hex.pdf")
TRCC13_corr4_bin10_hex <- plot_arrangement(TRCC13_corr4_bin10, hexbin_plot, 6)


### ploidy 4 not corrected
TRCC13_4_noCorr <- clean_file_noCorr("~/Ploidy 4/DTRCC_13_LLN_cluster1.titan.ichor.cna.txt", 4)

#pdf("CN_scatter_hex/Hex Plots/DTRCC13_LLN/DTRCC13_noCorr4_hex.pdf")
DTRCC13_noCorr4_hex <- plot_arrangement(TRCC13_4_noCorr, hexbin_plot, 6)


### ploidy 4 not corrected binned
TRCC13_4_noCorr_bin10 <- bin_100kb(TRCC13_4_noCorr)

#pdf("CN_scatter_hex/Hex Plots/DTRCC13_LLN/DTRCC13_noCorr4_bin10_hex.pdf")
DTRCC13_noCorr4_bin10_hex <- plot_arrangement(TRCC13_4_noCorr_bin10, hexbin_plot, 6)


