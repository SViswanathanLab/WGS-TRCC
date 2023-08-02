title: "CN Boxplot"
author: "Fiona McBride"
date: "2023-04-28"

library(ggplot2)
library(tidyverse)

# load the data file
freq_data <- read.table("tRCC_recurrent_CN_hist_final_ploidy_no_sample16_geneFreq.txt", header=TRUE, sep="\t")

# multiply the column of CN loss by -1 to make all the values negative so they go below the x axis
freq_data$loss_neg <- lapply(freq_data$All_loss, "*", -1)

# calculate the bin that is at the very edge of each chromosome to be able to add dividing lines

chr_list <- unique(freq_data$Chr)

chrom_lines <- c(0)

for (chrom in chr_list) {
  chr_subset <- freq_data[freq_data$Chr == chrom, ]
  chrom_lines <- append(chrom_lines, (nrow(chr_subset) + chrom_lines[length(chrom_lines)]))
}

# add bins across all the data to be able to continuously plot the whole genome
freq_data$bin <- c(1:nrow(freq_data))

# read in a data file containing the centromere location of each chromosome

centromeres_data <- read.table("centromere_locations")
centromeres_data <- centromeres_data %>% rename("chromosome" = "V1", "start" = "V2", "end" = "V3", "sample" = "V4")

# based on the centromere data above, find the first bin that starts the centromeric gap for each chromosome
# these bins were manually identified from the file read in above; un-comment the file import and column renaming to check  
centromeres <- c(1236, 3429, 5822, 7397, 9288, 
                 11221, 12925, 14370, 15810, 17152, 
                 18625, 19799, 20966, 22126, 23213, 
                 24379, 25165, 25932, 26813, 27419, 
                 27924, 28427, 29370)

# plot the data across the entire genome with chromosomes separated by sold black lines and the p arm of each chromosome overlayed with a grey box
titan_manual <- ggplot(freq_data, aes(x = bin)) +
  annotate("rect", xmin = chrom_lines[-24], xmax = centromeres, ymin = -1, ymax = 1,
           alpha = .4,fill = "grey") + 
  geom_bar(aes(y = All_gain), stat = "identity", fill = "red") + 
  geom_bar(aes(y = as.numeric(loss_neg)), stat = "identity", fill = "blue") + 
  theme_classic() + 
  geom_vline(xintercept = chrom_lines[-1], linetype="solid") +
  ylim(-1, 1) +
  xlab("Chromosome") + ylab("CN Loss/Gain") + theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent"))

#ggsave("titan_manual.pdf", titan_manual, height = 4, width = 8)

# print chromosomes that have gain greater than 0.45
gain_threshold <- freq_data[freq_data$All_gain >= 0.45, ]
unique(gain_threshold$Chr)

# print chromosomes that have loss greater than 0.45
loss_threshold <- freq_data[freq_data$loss_neg <= -0.45, ]
unique(loss_threshold$Chr)

# load the cancer gene census file
cancer_genes <- read.csv("Census_allSat May 13 21_30_39 2023.csv")

# split the "genome location" column so that the chromosome is in a separate column
chrom_number <- strsplit(cancer_genes$Genome.Location, split = ":")

chromosome = c()
location = c()

for (i in (1:length(chrom_number))){
  chromosome <- append(chromosome, (chrom_number[[i]][1]))
  location <- append(location, chrom_number[[i]][2])
}

cancer_genes$chrom <- chromosome
cancer_genes$bp_range <- location

# chromosome X shows up as "NA" in the dataframe, replace NA with X in the chromosome column
cancer_genes$chrom <- replace_na(cancer_genes$chrom, "X")

# split the "bp range" column so the start and end positions are in separate columns
range_split <- strsplit(cancer_genes$bp_range, split = "-")

start = c()
end = c()

for (i in (1:length(range_split))){
  start <- append(start, as.numeric(range_split[[i]][1]))
  end <- append(end, as.numeric(range_split[[i]][2]))
}

cancer_genes$start <- start 
cancer_genes$end <- end


# Oncogenes
# filter the data file to the chromosomes of interest, genes that are classified as Tier 1 (known evidence of cancer contribution), and are labeled as oncogenes

chrom1_onco <- cancer_genes[cancer_genes$chrom == "1" & cancer_genes$Tier == 1 & cancer_genes$Role.in.Cancer %in% c("oncogene, fusion", "oncogene"), ]

chrom17_onco <- cancer_genes[cancer_genes$chrom == "17" & cancer_genes$Tier == 1 & cancer_genes$Role.in.Cancer %in% c("oncogene, fusion", "oncogene"), ]

## TSGs
# filter the data file to the chromosomes of interest, genes that are classified as Tier 1 (known evidence of cancer contribution), and are labeled as TSGs

chrom1_TSG <- cancer_genes[cancer_genes$chrom == "1" & cancer_genes$Tier == 1 & cancer_genes$Role.in.Cancer %in% c("TSG, fusion", "TSG") , ]

chrom6_TSG <- cancer_genes[cancer_genes$chrom == "6" & cancer_genes$Tier == 1 & cancer_genes$Role.in.Cancer %in% c("TSG, fusion", "TSG") , ]

chrom9_TSG <- cancer_genes[cancer_genes$chrom == "9" & cancer_genes$Tier == 1 & cancer_genes$Role.in.Cancer %in% c("TSG, fusion", "TSG") , ]

chrom14_TSG <- cancer_genes[cancer_genes$chrom == "14" & cancer_genes$Tier == 1 & cancer_genes$Role.in.Cancer %in% c("TSG, fusion", "TSG") , ]

chrom17_TSG <- cancer_genes[cancer_genes$chrom == "17" & cancer_genes$Tier == 1 & cancer_genes$Role.in.Cancer %in% c("TSG, fusion", "TSG") , ]

## Filter genes on p and q arms
### split chromosome by centromere location in "centromere_data" df above
chrom1_q <- chrom1_onco[as.numeric(chrom1_onco$start) >= 122026459, ]

bin1q <- c()

for (i in c(1:nrow(chrom1_q))){
  condition <- freq_data[as.numeric(freq_data$Start) < as.numeric(chrom1_q$start[i]) & as.numeric(freq_data$End) > as.numeric(chrom1_q$start[i]) & freq_data$Chr == "chr1", ]
  bin1q <- append(bin1q, condition$bin)
}

chrom1_q$bin <- bin1q

key_genes_1q <- freq_data[freq_data$bin %in% bin1q, ]
key_genes_1q$gene <- chrom1_q$Gene.Symbol

chrom17_q <- chrom17_onco[as.numeric(chrom17_onco$start) >= 22813679, ]

bin17q <- c()

for (i in c(1:nrow(chrom17_q))){
  condition <- freq_data[as.numeric(freq_data$Start) < as.numeric(chrom17_q$start[i]) & as.numeric(freq_data$End) > as.numeric(chrom17_q$start[i]) & freq_data$Chr == "chr17", ]
  bin17q <- append(bin17q, condition$bin)
}

chrom1_p <- chrom1_TSG[as.numeric(chrom1_TSG$start) <= 122026459, ]

bin1p <- c()

for (i in c(1:nrow(chrom1_p))){
  condition <- freq_data[as.numeric(freq_data$Start) < as.numeric(chrom1_p$start[i]) & as.numeric(freq_data$End) > as.numeric(chrom1_p$start[i]) & freq_data$Chr == "chr1", ]
  bin1p <- append(bin1p, condition$bin)
}

chrom6_q <- chrom6_TSG[as.numeric(chrom6_TSG$start) >= 58553888, ]

bin6q <- c()

for (i in c(1:nrow(chrom6_q))){
  condition <- freq_data[as.numeric(freq_data$Start) < as.numeric(chrom6_q$start[i]) & as.numeric(freq_data$End) > as.numeric(chrom6_q$start[i]) & freq_data$Chr == "chr6", ]
  bin6q <- append(bin6q, condition$bin)
}

chrom9_p <- chrom9_TSG[as.numeric(chrom9_TSG$start) <= 43389635, ]

bin9p <- c()

for (i in c(1:nrow(chrom9_p))){
  condition <- freq_data[as.numeric(freq_data$Start) < as.numeric(chrom9_p$start[i]) & as.numeric(freq_data$End) > as.numeric(chrom9_p$start[i]) & freq_data$Chr == "chr9", ]
  bin9p <- append(bin9p, condition$bin)
}

# NO GENES IDENTIFIED IN CHROM 14 REGION

chrom17_p <- chrom17_TSG[as.numeric(chrom17_TSG$start) <= 22813679, ]

bin17p <- c()

for (i in c(1:nrow(chrom17_p))){
  condition <- freq_data[as.numeric(freq_data$Start) < as.numeric(chrom17_p$start[i]) & as.numeric(freq_data$End) > as.numeric(chrom17_p$start[i]) & freq_data$Chr == "chr17", ]
  bin17p <- append(bin17p, condition$bin)
}



## Genome Plot
# overview of where all the identified genes fall on the chromosome wide CN plot
## green lines for genes on the q arm, orange lines for genes on the p arm 

overview_plot <- ggplot(freq_data, aes(x = bin)) + 
  geom_bar(aes(y = All_gain), stat = "identity", fill = "red") + 
  geom_bar(aes(y = as.numeric(loss_neg)), stat = "identity", fill = "blue") + theme_classic() + 
  geom_vline(xintercept = chrom_lines[-1], linetype="solid") +
  geom_vline(xintercept = centromeres, linetype = "dotted") + 
  geom_vline(xintercept = bin1q, linetype = "solid", color = "green", size = 0.15) +
  geom_vline(xintercept = bin17q, linetype = "solid", color = "green", size = 0.15) +
  geom_vline(xintercept = bin1p, linetype = "solid", color = "orange", size = 0.15) +
  geom_vline(xintercept = bin6q, linetype = "solid", color = "orange", size = 0.15) +
  geom_vline(xintercept = bin9p, linetype = "solid", color = "orange", size = 0.15) +
  geom_vline(xintercept = bin17p, linetype = "solid", color = "orange", size = 0.15) +
  ylim(-1, 1) +
  xlab("Chromosome") + ylab("CN Loss/Gain")


# plot chromosome 1 with all of the labeled genes from the cancer gene census
chrom1 <- ggplot(freq_data[freq_data$Chr == "chr1", ], aes(x = bin)) + 
  geom_bar(aes(y = All_gain), stat = "identity", fill = "red") + 
  geom_bar(aes(y = as.numeric(loss_neg)), stat = "identity", fill = "blue") + theme_classic() + 
  geom_vline(xintercept = bin1q, linetype = "solid", color = "green", size = 0.15) +
  geom_vline(xintercept = bin1p, linetype = "solid", color = "orange", size = 0.15) +
  annotate("text", x =26, y = -0.9, label = "TNFRSF14", angle = 90, size = 3) +
  annotate("text", x =62, y = -0.65, label = "RPL22", angle = 90, size = 3) +
  annotate("text", x =90, y = -0.9, label = "CAMTA1", angle = 90, size = 3) +
  annotate("text", x =159, y = -0.9, label = "SPEN", angle = 90, size = 3) +
  annotate("text", x =171, y = -0.7, label = "SDHB", angle = 90, size = 3) +
  annotate("text", x =267, y = -0.9, label = "ARID1A", angle = 90, size = 3) +
  annotate("text", x =352, y = -0.9, label = "SFPQ", angle = 90, size = 3) +
  annotate("text", x =454, y = -0.9, label = "MUTYH", angle = 90, size = 3) +
  annotate("text", x =510, y = -0.6, label = "CDKN2C", angle = 90, size = 3) +
  annotate("text", x =530, y = -0.9, label = "EPS15", angle = 90, size = 3) +
  annotate("text", x =853, y = -0.9, label = "BCL10", angle = 90, size = 3) +
  annotate("text", x =929, y = -0.9, label = "RPL5", angle = 90, size = 3) +
  annotate("text", x =1144, y = -0.9, label = "TRIM33", angle = 90, size = 3) +
  annotate("text", x =1177, y = -0.6, label = "TENT5C", angle = 90, size = 3) +
  
  annotate("text", x =1476, y = 0.9, label = "BCL9", angle = 90, size = 3) +
  annotate("text", x =1576, y = 0.9, label = "FCRL4", angle = 90, size = 3) +
  annotate("text", x =1617, y = 0.65, label = "FCGR2B", angle = 90, size = 3) +
  annotate("text", x =1627, y = 0.9, label = "DDR2", angle = 90, size = 3) +
  annotate("text", x =1660, y = 0.8, label = "PBX1", angle = 90, size = 3) +
  annotate("text", x =1791, y = 0.9, label = "ABL2", angle = 90, size = 3) +
  annotate("text", x =2046, y = 0.9, label = "MDM4", angle = 90, size = 3) +
  annotate("text", x =2070, y = 0.65, label = "ELK4", angle = 90, size = 3) +
  annotate("text", x =2261, y = 0.9, label = "H3F3A", angle = 90, size = 3) +
  ylim(-1, 1) +
  xlab("Chromosome") + ylab("CN Loss/Gain")

# final plot with filtered/selected genes
chrom1_labeled <- ggplot(freq_data[freq_data$Chr == "chr1", ], aes(x = bin)) + 
  annotate("rect", xmin = chrom_lines[1], xmax = centromeres[1], ymin = -1, ymax = 1,
           alpha = .4,fill = "grey") +
  geom_bar(aes(y = All_gain), stat = "identity", fill = "red") + 
  geom_bar(aes(y = as.numeric(loss_neg)), stat = "identity", fill = "blue") + theme_classic() + 
  geom_vline(xintercept = c(171, 267, 510, 2046), linetype = "solid", color = "black", size = 0.2) +
  annotate("text", x =151, y = -0.9, label = "SDHB", angle = 90, size = 3) +
  annotate("text", x =247, y = -0.9, label = "ARID1A", angle = 90, size = 3) +
  annotate("text", x =490, y = -0.9, label = "CDKN2C", angle = 90, size = 3) +
  annotate("text", x =2026, y = 0.9, label = "MDM4", angle = 90, size = 3) +
  ylim(-1, 1) +
  xlab("Chromosome") + ylab("CN Loss/Gain") + theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent"))

ggsave("chrom1_labeled_manual.pdf", chrom1_labeled, width = 2, height = 2)


# put genes and gene locations in order
gene_loc_1p <- data.frame("location" = bin1p, gene = chrom1_p$Gene.Symbol)
gene_loc_1p[order(gene_loc_1p$location), ]

# put genes and gene locations in order 
gene_loc_1q <- data.frame("location" = bin1q, gene = chrom1_q$Gene.Symbol)
gene_loc_1q[order(gene_loc_1q$location), ]

# label all genes identified on chromosome 6

chrom6 <- ggplot(freq_data[freq_data$Chr == "chr6", ], aes(x = bin)) + 
  geom_bar(aes(y = All_gain), stat = "identity", fill = "red") + 
  geom_bar(aes(y = as.numeric(loss_neg)), stat = "identity", fill = "blue") + theme_classic() + 
  geom_vline(xintercept = bin6q, linetype = "solid", color = "orange", size = 0.15) +
  annotate("text", x =11675, y = -0.9, label = "PRDM1", angle = 90, size = 3.5) +
  annotate("text", x =11894, y = -0.9, label = "PTPRK", angle = 90, size = 3.5) +
  annotate("text", x =11993, y = -0.9, label = "TNFAIP3", angle = 90, size = 3.5) +
  annotate("text", x =12111, y = -0.9, label = "LATS1", angle = 90, size = 3.5) +
  annotate("text", x =12182, y = -0.9, label = "ARID1B", angle = 90, size = 3.5) +
  ylim(-1, 1) +
  xlab("Chromosome") + ylab("CN Loss/Gain")

chrom6_labeled <- ggplot(freq_data[freq_data$Chr == "chr6", ], aes(x = bin)) +
  annotate("rect", xmin = chrom_lines[6], xmax = centromeres[6], ymin = -1, ymax = 1,
           alpha = .4,fill = "grey") + 
  geom_bar(aes(y = All_gain), stat = "identity", fill = "red") + 
  geom_bar(aes(y = as.numeric(loss_neg)), stat = "identity", fill = "blue") + theme_classic() + 
  geom_vline(xintercept = c(11675, 12111, 12182), linetype = "solid", color = "black", size = 0.2) +
  annotate("text", x =11655, y = -0.9, label = "PRDM1", angle = 90, size = 3.5) +
  annotate("text", x =12091, y = -0.9, label = "LATS1", angle = 90, size = 3.5) +
  annotate("text", x =12162, y = -0.9, label = "ARID1B", angle = 90, size = 3.5) +
  ylim(-1, 1) +
  xlab("Chromosome") + ylab("CN Loss/Gain") 


# put genes and gene locations in order
gene_loc_6 <- data.frame("location" = bin6q, gene = chrom6_q$Gene.Symbol)
gene_loc_6[order(gene_loc_6$location), ]

# label all genes identified from the cancer gene atlas on chromosome 9

chrom9 <- ggplot(freq_data[freq_data$Chr == "chr9", ], aes(x = bin)) + 
  geom_bar(aes(y = All_gain), stat = "identity", fill = "red") + 
  geom_bar(aes(y = as.numeric(loss_neg)), stat = "identity", fill = "blue") + theme_classic() + 
  geom_vline(xintercept = bin9p, linetype = "solid", color = "orange", size = 0.15) +
  annotate("text", x =15424, y = -0.9, label = "CD274", angle = 90, size = 3.5) +
  annotate("text", x =15589, y = -0.9, label = "CDKN2A", angle = 90, size = 3.5) +
  annotate("text", x =15720, y = -0.9, label = "FANCG", angle = 90, size = 3.5) +
  ylim(-1, 1) +
  xlab("Chromosome") + ylab("CN Loss/Gain") 

# final plot with filtered/selected genes
chrom9_labeled <- ggplot(freq_data[freq_data$Chr == "chr9", ], aes(x = bin)) + 
  annotate("rect", xmin = chrom_lines[9], xmax = centromeres[9], ymin = -1, ymax = 1,
           alpha = .4,fill = "grey") +
  geom_bar(aes(y = All_gain), stat = "identity", fill = "red") + 
  geom_bar(aes(y = as.numeric(loss_neg)), stat = "identity", fill = "blue") + theme_classic() + 
  geom_vline(xintercept = 15589, linetype = "solid", color = "black", size = 0.2) +
  annotate("text", x =15569, y = -0.9, label = "CDKN2A", angle = 90, size = 3.5) +
  ylim(-1, 1) +
  xlab("Chromosome") + ylab("CN Loss/Gain") + theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent"))

ggsave("chrom9_labeled_manual.pdf", chrom9_labeled, width = 2, height = 2)

# put genes and gene locations in order
gene_loc_9p <- data.frame("location" = bin9p, gene = chrom9_p$Gene.Symbol)
gene_loc_9p[order(gene_loc_9p$location), ]

# label all genes identified from the cancer gene atlas on chromosome 17

chrom17 <- ggplot(freq_data[freq_data$Chr == "chr17", ], aes(x = bin)) + 
  geom_bar(aes(y = All_gain), stat = "identity", fill = "red") + 
  geom_bar(aes(y = as.numeric(loss_neg)), stat = "identity", fill = "blue") + theme_classic() + 
  geom_vline(xintercept = bin17q, linetype = "solid", color = "green", size = 0.15) +
  geom_vline(xintercept = bin17p, linetype = "solid", color = "orange", size = 0.15) +
  
  annotate("text", x =24928, y = -0.9, label = "YWHAE", angle = 90, size = 3.5) +
  annotate("text", x =24996, y = -0.9, label = "PER1", angle = 90, size = 3.5) +
  annotate("text", x =25075, y = -0.9, label = "NCOR1", angle = 90, size = 3.5) +
  annotate("text", x =25092, y = -0.9, label = "FLCN", angle = 90, size = 3.5) +
  
  annotate("text", x =25273, y = 0.9, label = "TAF15", angle = 90, size = 3.5) +
  annotate("text", x =25312, y = 0.9, label = "ERBB2", angle = 90, size = 3.5) +
  annotate("text", x =25325, y = 0.5, label = "RARA", angle = 90, size = 3.5) +
  annotate("text", x =25338, y = 0.9, label = "STAT3", angle = 90, size = 3.5) +
  annotate("text", x =25350, y = 0.7, label = "ETV4", angle = 90, size = 3.5) +
  annotate("text", x =25467, y = 0.9, label = "HLF", angle = 90, size = 3.5) +
  annotate("text", x =25487, y = 0.9, label = "MSI2", angle = 90, size = 3.5) +
  annotate("text", x =25521, y = 0.9, label = "PPM1D", angle = 90, size = 3.5) +
  annotate("text", x =25554, y = 0.9, label = "CD79B", angle = 90, size = 3.5) +
  annotate("text", x =25565, y = 0.6, label = "DDX5", angle = 90, size = 3.5) +
  annotate("text", x =25672, y = 0.9, label = "H3F3B", angle = 90, size = 3.5) +
  annotate("text", x =25690, y = 0.7, label = "SRSF2", angle = 90, size = 3.5) +
  ylim(-1, 1) +
  xlab("Chromosome") + ylab("CN Loss/Gain")


# final plot with filtered/selected genes
chrom17_labeled <- ggplot(freq_data[freq_data$Chr == "chr17", ], aes(x = bin)) + 
  annotate("rect", xmin = chrom_lines[17], xmax = centromeres[17], ymin = -1, ymax = 1,
           alpha = .4,fill = "grey") + 
  geom_bar(aes(y = All_gain), stat = "identity", fill = "red") + 
  geom_bar(aes(y = as.numeric(loss_neg)), stat = "identity", fill = "blue") + theme_classic() + 
  geom_vline(xintercept = c(25092, 25312, 25690), linetype = "solid", color = "black", size = 0.2) +
  annotate("text", x =25082, y = -0.9, label = "FLCN", angle = 90, size = 3.5) +
  annotate("text", x =25302, y = 0.9, label = "ERBB2", angle = 90, size = 3.5) +
  annotate("text", x =25680, y = 0.9, label = "SRSF2", angle = 90, size = 3.5) +
  ylim(-1, 1) +
  xlab("Chromosome") + ylab("CN Loss/Gain")

ggsave("chrom17_labeled.pdf", chrom17_labeled)

# put genes and gene locations in order
gene_loc_17p <- data.frame("location" = bin17p, gene = chrom17_p$Gene.Symbol)
gene_loc_17p[order(gene_loc_17p$location), ]

# put genes and gene locations in order
gene_loc_17q <- data.frame("location" = bin17q, gene = chrom17_q$Gene.Symbol)
gene_loc_17q[order(gene_loc_17q$location), ]

# final labeled plot of all chromosomes 
whole_genome_labeled <- ggplot(freq_data, aes(x = bin)) + 
  annotate("rect", xmin = chrom_lines[-24], xmax = centromeres, ymin = -1, ymax = 1,
           alpha = .4,fill = "grey") + 
  geom_bar(aes(y = All_gain), stat = "identity", fill = "red") + 
  geom_bar(aes(y = as.numeric(loss_neg)), stat = "identity", fill = "blue") + theme_classic() + 
  geom_vline(xintercept = chrom_lines[-1], linetype="solid") +
  geom_vline(xintercept = c(171, 267, 510, 2046, 11675, 12111, 12182, 15589, 25092, 25312, 25690), 
             linetype = "solid", color = "black", size = 0.2) +
  annotate("text", x =151, y = -0.9, label = "SDHB", angle = 90, size = 2.5) +
  annotate("text", x =247, y = -0.9, label = "ARID1A", angle = 90, size = 2.5) +
  annotate("text", x =490, y = -0.9, label = "CDKN2C", angle = 90, size = 2.5) +
  annotate("text", x =2026, y = 0.9, label = "MDM4", angle = 90, size = 2.5) +
  annotate("text", x =11655, y = -0.9, label = "PRDM1", angle = 90, size = 2.5) +
  annotate("text", x =12091, y = -0.9, label = "LATS1", angle = 90, size = 2.5) +
  annotate("text", x =12162, y = -0.9, label = "ARID1B", angle = 90, size = 2.5) +
  annotate("text", x =15569, y = -0.9, label = "CDKN2A", angle = 90, size = 2.5) +
  annotate("text", x =25082, y = -0.9, label = "FLCN", angle = 90, size = 2.5) +
  annotate("text", x =25302, y = 0.9, label = "ERBB2", angle = 90, size = 2.5) +
  annotate("text", x =25680, y = 0.9, label = "SRSF2", angle = 90, size = 2.5) +  
  ylim(-1, 1) +
  xlab("Chromosome") + ylab("CN Loss/Gain") + theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent"))

ggsave("whole_genome_labeled_manual.pdf", whole_genome_labele, width = 11, height = 2.5)
