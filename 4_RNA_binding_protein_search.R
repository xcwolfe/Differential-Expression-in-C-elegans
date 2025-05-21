## 4th ##

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
library("data.table")
library("tidyverse")
library("pheatmap")
library("rex")
library("MatrixGenerics")
library("modeest")
library("dplyr")
library("plyr")
library("grid")
library("gridBase")
library("stringi")
library("RColorBrewer")
library("ggplot2")
library("purrr")
```

```{r}
save.image("C:/Your/directory/.RData")
load("C:/Your/directory/.RData")
```


# Convert RBP names to WBGene:

# First I will create a file with the RNA binding proteins of interest:
```{r}
RBP_test_genes_file <- read.csv(file = "HUGHES_RBP_genenames.csv", header = T, row.names = 1)
RBP_test_genes <- rownames(RBP_test_genes_file)

mined_RBPs_file <- read.csv(file = "simplemine_results.csv", header = T) # Duplicates of nog-1: used the row that contained "public name" of nog-1
mined_RBPs <- mined_RBPs_file[,2:3]
RBP_list <- c(mined_RBPs[,1])
```

# Then, I would like to search the .csv files I made in type vs type 1 x 1 for these RBPs:
```{r}
setwd("D:/Zach Wolfe's DESeq analysis")

overlap_data <- fread(file = "resAWBvsVC.csv", header = TRUE)[,c(1,3,7)]
match(RBP_list, overlap_data$V1, nomatch = NA)
print(length(which(match(RBP_list, overlap_data$V1) < 486 ))) # We know our RNABP list is 484 genes long. I need to add 2 to 484 (row 1 is header (+1), less than function adds +1) I am searching the list of ALL genes and printing the number of RNABPs expressed in the top 484 entries of each comparison. Not statistically accurate but we are getting somewhere.
```

# Then, I want to display these differential expression measures in a new table:
```{r}
setwd("D:/Zach Wolfe's DESeq analysis")

RBP_results <- matrix(nrow = 46, ncol = 46, dimnames = list(singletypes,singletypes))
RBP_results[is.na(RBP_results)] <- 0

for(i in singletypes){
  for(j in singletypes){
    if(i != j){

overlap_data <- fread(file = paste("res_", i, "_vs_", j, "_.csv", sep=""))[,c(1,3,7)]
match(RBP_list, overlap_data$V1, nomatch = NA)

RBP_results[i,j] <- print(length(which(match(RBP_list, overlap_data$V1) < 312))) # 312 is top 1% significant of all genes (31173)

# RBP_results[i,j] <- write.csv(overlap_data, file=paste("RBPs", i, "vs", j, ".csv", sep="")) # Write csv of results

    }
    if(i == j) next
    }
}
```

# Print names of RNABPs from chart above:
# This will be a practice chunk to ensure everything works:
```{r, include=FALSE}
setwd("D:/Zach Wolfe's DESeq analysis")

RBP_cutoff_results <- matrix(nrow = 46, ncol = 46, dimnames = list(singletypes,singletypes))
RBP_cutoff_results[is.na(RBP_cutoff_results)] <- 0

for(i in c("ADL")){
   for(j in c("AFD")){
    if(i != j){

overlap_data <- fread(file = paste("res", i, "vs", j, ".csv", sep = ""))[,c(1,3,7)] %>% as.data.frame()

rownames(overlap_data) <- overlap_data$V1
RBP_list %in% rownames(overlap_data)
intersection <- rownames(overlap_data)
# print(overlap_data[intersection,])

comparison <- print(paste(i, "vs", j))
new_row <- print(paste(comparison, "#1"))
new_row2 <- print(paste(comparison, "#2"))
new_row3 <- print(paste(comparison, "#3"))

totaldata <- data.frame(overlap_data[intersection,])
totaldata <- na.omit(totaldata)
topDEgene <- as.data.frame(totaldata[1,], row.names = comparison, col.names = c("V1", "log2FoldChange", "padj"))
topDEgene[new_row,] <- rbind(totaldata[1,], new_row)
topDEgene[new_row2,] <- rbind(totaldata[2,], new_row2)
topDEgene[new_row3,] <- rbind(totaldata[3,], new_row3)
cutoffgenes <- as.data.frame(totaldata[which(abs(totaldata$log2FoldChange)>2 & totaldata$padj < .01),]) # Change cutoffs and abs() if desired
RBP_cutoff_results[i,j] <- print(nrow(cutoffgenes))

write.csv(cutoffgenes, file=paste("DiffExpGenes", i, "vs", j, ".csv", sep=""))
     }
    if (i == j) next
   }
}

RBP_cutoff_results[is.na(RBP_cutoff_results)] <- 0
```

# Get rid of the first row in topDEgene (created in line 105):
```{r}
topDEgene <- topDEgene[-1,] # run this AFTER the chunk above JUST ONCE, then "#" out this row
```

# Now let's run the RNABP loop in full:
```{r, include=FALSE}
setwd("D:/Zach Wolfe's DESeq analysis")

RBP_cutoff_results <- matrix(nrow = 46, ncol = 46, dimnames = list(singletypes,singletypes))
RBP_cutoff_results[is.na(RBP_results)] <- 0

for(i in singletypes){
   for(j in singletypes){
    if(i != j){

overlap_data <- fread(file = paste("res_", i, "_vs_", j, "_.csv", sep = ""))[,c(1,3,7)] %>% as.data.frame()

rownames(overlap_data) <- overlap_data$V1
RBP_list %in% rownames(overlap_data)
intersection <- rownames(overlap_data) %in% RBP_list
# print(overlap_data[intersection,])

comparison <- print(paste(i, "vs", j))
new_row <- print(paste(comparison, "#1"))
new_row2 <- print(paste(comparison, "#2"))
new_row3 <- print(paste(comparison, "#3"))

totaldata <- data.frame(overlap_data[intersection,])
totaldata <- na.omit(totaldata)
topDEgene[new_row,] <- rbind(totaldata[1,], new_row)
topDEgene[new_row2,] <- rbind(totaldata[2,], new_row2)
topDEgene[new_row3,] <- rbind(totaldata[3,], new_row3)
cutoffgenes <- as.data.frame(totaldata[which(abs(totaldata$log2FoldChange)>2 & totaldata$padj < .01),]) # Change cutoffs and abs() if desired
RBP_cutoff_results[i,j] <- print(nrow(cutoffgenes))
write.csv(cutoffgenes, file=paste("DiffExpRBPs", i, "vs", j, ".csv", sep=""))

# generate histograms using cutoffgenes:

#      if (nrow(cutoffgenes) > 0){
# cutoff_plot <- plot(cutoffgenes$log2FoldChange, y = cutoffgenes$padj, type = "p", xlim = c(-8,10), ylim = c(cutoffgenes[nrow(cutoffgenes),3],cutoffgenes[1,3]), log = "y", pch = 16, col = rainbow(15), main = paste("Differentially Experessed RBPs (WormBase ID) between", i, "and", j), xlab = "log2FoldChange value", ylab = "p-value (adj.)")
# text(cutoffgenes$log2FoldChange, cutoffgenes$padj, stri_sub(cutoffgenes$V1, 9, 14), cex = 1.5, pos = 4, col = rainbow(15))

# cutoff_barplot <- barplot(cutoffgenes$log2FoldChange, width = 0.8, col = rainbow(10), main = paste("Differentially Expressed RBPs between", i, "and", j), xlab = "log2FoldChange of Gene", ylab = "", horiz = TRUE)
# text(0, cutoff_barplot, cutoffgenes$V1, cex = 0.8, pos = 4)
#      }
#    if (nrow(cutoffgenes) == 0) next
     }
    if (i == j) next
   }
}

RBP_cutoff_results[is.na(RBP_cutoff_results)] <- 0
```

# Save as .csv:
```{r}
write.csv(topDEgene, file = "topDifferentiallyExpressedRBPsbyComparison_.csv")
write.csv(RBP_cutoff_results, file = "RBP_cutoff_results_.csv")
```

# Backup in case something happens to the matrix above:
```{r}
#RBP_cutoff_results <- read.csv(file = "RBP_cutoff_results_.csv", header = TRUE)
#rownames(RBP_cutoff_results) <- colnames(RBP_cutoff_results[2:47])
#RBP_cutoff_results <- RBP_cutoff_results[2:47]
```

# Plot heatmaps:
```{r}
plot.new()
pheatmap(RBP_cutoff_results,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         fontsize = 7.5,
         col = rev(colors),
         main = "Number of Differentially Expressed RBPs by Cell Type (log2FoldChange > 2 and adj. p-value < 0.01)")

plot.new()
pheatmap(RBP_cutoff_results,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = T,
         number_format = "%.0f",
         number_color = "black",
         fontsize = 7.5,
         fontsize_number = 6.5,
         col = colors,
         main = "Number of Differentially Expressed RBPs by Cell Type (log2FoldChange > 2 and adj. p-value < 0.01)")
```

# Calculate total log2foldchange for each gene:
```{r}
setwd("D:/Zach Wolfe's DESeq analysis")

log2foldchange_per_RBP <-  matrix(nrow = 46, ncol = 484, dimnames = list(singletypes, mined_RBPs$WormBase.Gene.ID))
log2foldchange_per_RBP[is.na(log2foldchange_per_RBP)] <- 0

for(j in singletypes){
 for(k in singletypes){
  if(j != k){

DiffExpRBPs <- fread(file = paste("DiffExpRBPs", j, "vs", k, ".csv", sep = ""))[,2:3] %>% as.data.frame()
log2foldchange_per_RBP[j,DiffExpRBPs$V1] <- print(DiffExpRBPs$log2FoldChange) + log2foldchange_per_RBP[j, DiffExpRBPs$V1]
     }
      if (j == k) next
     }
   }

log2foldchange_per_RBP_zeroes <- (colSums(log2foldchange_per_RBP) == 0)
log2foldchange_per_RBP <- log2foldchange_per_RBP[, !log2foldchange_per_RBP_zeroes]
```

# Which cell types have the most upregulation and downregulation?
```{r}
rownames(log2foldchange_per_RBP)[rowSums(log2foldchange_per_RBP) == max(rowSums(log2foldchange_per_RBP))]
rownames(log2foldchange_per_RBP)[rowSums(log2foldchange_per_RBP) == min(rowSums(log2foldchange_per_RBP))]

upreg_cell_type <- rownames(log2foldchange_per_RBP)[rowSums(log2foldchange_per_RBP) == max(rowSums(log2foldchange_per_RBP))]
downreg_cell_type <- rownames(log2foldchange_per_RBP)[rowSums(log2foldchange_per_RBP) == min(rowSums(log2foldchange_per_RBP))]

log2foldchange_per_RBP[upreg_cell_type,]
log2foldchange_per_RBP[downreg_cell_type,]
```

# Which genes have the greatest difference in regulation?
```{r}
log2foldchange_per_RBP_colnames <- colnames(log2foldchange_per_RBP)

for(k in 1:ncol(log2foldchange_per_RBP)){
print(log2foldchange_per_RBP_colnames[k])
print(max(log2foldchange_per_RBP[,k]-min(log2foldchange_per_RBP[,k])))
}
```

# Plot heatmaps:
```{r}
colors <- colorRampPalette(brewer.pal(6, "Greens"))(250)
pheatmap(log2foldchange_per_RBP,
         border_color = "black",
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         fontsize = 6,
         angle_col = 45,
         col = colors,
         main = "Total log2FoldChange Accross All Comparisons")

pheatmap(log2foldchange_per_RBP[,1:50],
         border_color = "black",
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         fontsize = 7,
         angle_col = 90,
         col = colors,
         main = "Total log2FoldChange Accross All Comparisons (first quintile)")

pheatmap(log2foldchange_per_RBP[,51:100],
         border_color = "black",
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         fontsize = 7,
         angle_col = 90,
         col = colors,
         main = "Total log2FoldChange Accross All Comparisons (second quintile)")

pheatmap(log2foldchange_per_RBP[,101:150],
         border_color = "black",
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         fontsize = 7,
         angle_col = 90,
         col = colors,
         main = "Total log2FoldChange Accross All Comparisons (third quintile)")

pheatmap(log2foldchange_per_RBP[,151:200],
         border_color = "black",
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         fontsize = 7,
         angle_col = 90,
         col = colors,
         main = "Total log2FoldChange Accross All Comparisons (fourth quintile)")

pheatmap(log2foldchange_per_RBP[,201:257],
         border_color = "black",
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         fontsize = 7,
         angle_col = 90,
         col = colors,
         main = "Total log2FoldChange Accross All Comparisons (fifth quintile)")
```

# write.csv for sum of log2foldchange per RBP (in which abs(log2FoldChange) > 2 & padj < .01) in each cell type) - colSums(0) columns removed
```{r}
write.csv(log2foldchange_per_RBP, file = "sum of log2foldchange per RBP.csv")
```

# Backup read.csv:
```{r}
log2foldchange_per_RBP <- read.csv("sum of log2foldchange per RBP.csv", row.names = 1)
```


# Combine every excel file for each cell type:

## Note: If there are 0 Diff. Exp. RBPs in a comparison, dplyr will assume that the rows in the .csv file are logical (we want them to be character/numeric). This will prevent us from merging rows in different .csv files together in the upcoming chunk. You may have to examine each Diff. Exp. RBP file and manually move them outside of your working directory if there are 0 rows:
```{r}
setwd("D:/Zach Wolfe's DESeq analysis")

ADLDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsADLvs*")
AFDDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsAFDvs*")
AIMDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsAIMvs*")
AINDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsAINvs*")
AIYDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsAIYvs*")
AINDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsAINvs*")
ASELDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsASELvs*")
ASERDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsASERvs*")
ASGDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsASGvs*")
ASIDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsASIvs*")
ASKDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsASKvs*")
AVADiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsAVAvs*")
AVEDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsAVEvs*")
AVGDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsAVGvs*")
AVHDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsAVHvs*")
AVKDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsAVKvs*")
AVLDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsAVLvs*")
AVMDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsAVMvs*")
AWADiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsAWAvs*")
AWBDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsAWBvs*")
AWCDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsAWCvs*")
BAGDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsBAGvs*")
CANDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsCANvs*")
DADiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsDAvs*")
DDDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsDDvs*")
DVCDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsDVCvs*")
I5DiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsI5vs*")
IL1DiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsIL1vs*")
IL2DiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsIL2vs*")
LUADiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsLUAvs*")
NSMDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsNSMvs*")
OLLDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsOLLvs*")
OLQDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsOLQvs*")
PHADiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsPHAvs*")
PVCDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsPVCvs*")
PVDDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsPVDvs*")
PVMDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsPVMvs*")
RIADiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsRIAvs*")
RICDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsRICvs*")
RIMDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsRIMvs*")
RISDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsRISvs*")
RMDDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsRMDvs*")
SMBDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsSMBvs*")
SMDDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsSMDvs*")
VBDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsVBvs*")
VCDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsVCvs*")
VDDiffExpRBPfiles <- list.files(pattern = "DiffExpRBPsVDvs*")
```

# Make a list of all DiffExpRBPfiles for each cell type:
```{r}
DiffExpRBPfileList <- rep(paste0(singletypes, "DiffExpRBPfiles"), len = 46)
```

# Log2FoldChange values for unique RBP names in individual "cell type vs" files:
```{r}
setwd("D:/Zach Wolfe's DESeq analysis")

LargeRBPdf <- SMBDiffExpRBPfiles %>% map_dfr(read.csv) # change cell type here
LargeRBPdf <- LargeRBPdf[,-c(1)]
uniqueRBPs <- unique(LargeRBPdf$V1)

logSums <- matrix(nrow = length(uniqueRBPs), ncol = 1, dimnames = list(uniqueRBPs, "log2foldchange"))
logSums[is.na(logSums)] <- 0

logMeans <- matrix(nrow = length(uniqueRBPs), ncol = 1, dimnames = list(uniqueRBPs, "log2foldchange"))
logMeans[is.na(logMeans)] <- 0

for(i in uniqueRBPs){
logSums[i,1] <- sum(LargeRBPdf[LargeRBPdf$V1 == i, 2]) # perhaps mean would be more appropriate than sum for this type of comparison?
logMeans[i,1] <- mean(LargeRBPdf[LargeRBPdf$V1 == i, 2])
}
```


# Save/change cell type here:
```{r}
write.csv(logSums, file = "logSumsSMB.csv")
write.csv(logMeans, file = "logMeansSMB.csv") # change cell type here
#write.csv(LargeRBPdf, file = "Most interesting RBPs in cell type AFD.csv")
```

# Make a heatmap for the logSums value in each gene for a given cell type:
```{r}
setwd("D:/Zach Wolfe's DESeq analysis")

logSums_RBP_table <- matrix(nrow = 484, ncol = 46, dimnames = list(RBP_list,singletypes))
logSums_RBP_table[is.na(logSums_RBP_table)] <- 0

for (i in singletypes){
logSums_for_graph <- fread(file = paste0("logSums", i, ".csv", sep = "")) %>% as.data.frame()
rownames(logSums_for_graph) <- logSums_for_graph[,1]
  for (j in RBP_list){
logSums_RBP_table[j,i] <- print(logSums_for_graph[j,2])
logSums_RBP_table[is.na(logSums_RBP_table)] <- 0
  }
}

colors <- colorRampPalette(brewer.pal(6, "PuRd"))(250)
pheatmap(logSums_RBP_table,
         border_color = "black",
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         number_format = "%.2f",
         number_color = "black",
         fontsize = 6,
         fontsize_number = 6,
         angle_col = 45,
         col = colors,
         main = "Sum of All log2FoldChange Values for Unique RBPs in One Cell Type vs All Other Cell Types")

logSums_RBP_table <- logSums_RBP_table[rowSums(logSums_RBP_table) != 0,] # if you want to keep rows with 0 data, "#" this row out
maxDEcelltypes <- as.data.frame(table(colnames(logSums_RBP_table))[apply(logSums_RBP_table, 1, which.max)])
maxDEcelltypes <- maxDEcelltypes[order(maxDEcelltypes$Freq, decreasing = TRUE),]
minDEcelltypes <- as.data.frame(table(colnames(logSums_RBP_table))[apply(logSums_RBP_table, 1, which.min)])
minDEcelltypes <- minDEcelltypes[order(minDEcelltypes$Freq, decreasing = TRUE),]
maxDERBPs <- as.data.frame(table(colnames(t(logSums_RBP_table)))[apply(t(logSums_RBP_table), 1, which.max)])
maxDERBPs <- count(maxDERBPs$Var1)
maxDERBPs <- maxDERBPs[order(maxDERBPs$freq, decreasing = TRUE),]
minDERBPs <- as.data.frame(table(colnames(t(logSums_RBP_table)))[apply(t(logSums_RBP_table), 1, which.min)])
minDERBPs <- count(minDERBPs$Var1)
minDERBPs <- minDERBPs[order(minDERBPs$freq, decreasing = TRUE),]
```

# Make a log2foldchange list just for each cell type here (useful for gene ontology tools):
```{r}
for (i in 1:length(singletypes)){

RBP_log2foldchange_list <- logSums_RBP_table
RBP_log2foldchange_list <- rownames_to_column(data.frame(RBP_log2foldchange_list), "gene_names")
RBP_log2foldchange_list <- RBP_log2foldchange_list[,c(1,i+1)] # column 1 contains the row_names, column i+1 contains my cell type's log2foldchange values
rownames(RBP_log2foldchange_list) <- RBP_log2foldchange_list[,1]
RBP_log2foldchange_list <- RBP_log2foldchange_list[rev(order(RBP_log2foldchange_list[,2])),] # reorder by cell type's log2foldchange value column
RBP_log2foldchange_list <- as.data.frame(RBP_log2foldchange_list[which(RBP_log2foldchange_list[,2] > 2),]) # Change log2foldchange cutoffs if desired

write.table(RBP_log2foldchange_list[,1], file = paste("RBP log2foldchange list", singletypes[i], ".tsv"), row.names = FALSE, quote = FALSE, col.names = FALSE) # this is the list we will input into a gene enrichment tool (FuncAssociate, g:Profiler, etc.)
}
```

```{r}
write.csv(logSums_RBP_table, file = "logSums_RBP_table.csv")
write.csv(RBP_log2foldchange_list[,1], file = "RBP log2foldchange list (AFD).csv", row.names = FALSE) # this is the list we will input into a gene enrichment tool (WormCat, g:Profiler, etc.)
write.csv(maxDEcelltypes, file = "maxDEcelltypes.csv")
write.csv(minDEcelltypes, file = "minDEcelltypes.csv")
write.csv(maxDERBPs, file = "maxDERBPs.csv")
write.csv(minDERBPs, file = "minDERBPs.csv")
```

# Backup:
```{r}
logSums_RBP_table <- read.csv("logSums_RBP_table.csv", row.names = 1)
```

# max Diff. Exp. RBP heatmap:
```{r}
maxs <- rowMaxs(as.matrix(logSums_RBP_table))
logSums_RBP_table <- cbind(logSums_RBP_table, maxs)
logSums_RBP_table <- logSums_RBP_table[order(unlist(logSums_RBP_table[,47]), decreasing = TRUE),]
logSums_RBP_table <- logSums_RBP_table[,-47]

# Now may be a good time to change gene names from WBGeneXXXXXXXX to common gene names using WormBase's simple mine tool:
RBP_common_names <- read.csv("RBP_common_names.csv") # see "RBP_common_names.csv" for a list of the 329 RBPs used in this portion of the analysis

rownames(logSums_RBP_table) <- RBP_common_names$Public.Name

col_pallete <- colorRampPalette(c("white", "white", "white", "gold", "darkgoldenrod1", "darkorange1"))
colors <- col_pallete(250)
pheatmap(logSums_RBP_table[1:30,],
         border_color = "black",
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         number_format = "%.2f",
         number_color = "black",
         fontsize = 6.75,
         fontsize_number = 6,
         angle_col = 45,
         col = colors,
         main = "Sum of All log2FoldChange Values for Most Dysregulated RBPs in One Cell Type vs All Other Cell Types")
```

# Check relative expression of AFD compared to other cell types:
```{r}
colSums(logSums_RBP_table)["AFD"]/nrow(logSums_RBP_table)
expression_value_RBP <- matrix(nrow = ncol(logSums_RBP_table), ncol = 2)
colnames(expression_value_RBP) <- c("Cell_type", "Normalized_logSums_value")
expression_value_RBP[,1] <- colnames(logSums_RBP_table)

for (i in 1:ncol(logSums_RBP_table)){
  expression_value_RBP[i,2] <- format(colSums(logSums_RBP_table)[i]/nrow(logSums_RBP_table), scientific = F)
}

ggplot(data = data.frame(expression_value_RBP), aes(x = expression_value_RBP[,1], y = as.numeric(expression_value_RBP[,2]))) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(title = "Normalized logSums value of RBPs per cell type", x = "Cell type", y = "Normalized logSums value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

# Which RBPs have the highest value overall in a single cell?
# Report: which RBP, what is the value, and which cell?
```{r}
log2foldchange_summary_df <- data.frame()

for (i in singletypes){
cell_type <- print(paste(i))
new_row <- print(paste(cell_type, "#1"))
new_row2 <- print(paste(cell_type, "#2"))
new_row3 <- print(paste(cell_type, "#3"))
new_row4 <- print(paste(cell_type, "#4"))
new_row5 <- print(paste(cell_type, "#5"))
new_row6 <- print(paste(cell_type, "#6"))
new_row7 <- print(paste(cell_type, "#7"))
new_row8 <- print(paste(cell_type, "#8"))
new_row9 <- print(paste(cell_type, "#9"))
new_row10 <- print(paste(cell_type, "#10"))

logSums_table_new <- logSums_RBP_table[order(logSums_RBP_table[,i], decreasing = TRUE),]
log2foldchange_summary_df[new_row, 1] <- print(rownames(logSums_table_new)[1])
log2foldchange_summary_df[new_row2, 1] <- print(rownames(logSums_table_new)[2])
log2foldchange_summary_df[new_row3, 1] <- print(rownames(logSums_table_new)[3])
log2foldchange_summary_df[new_row4, 1] <- print(rownames(logSums_table_new)[4])
log2foldchange_summary_df[new_row5, 1] <- print(rownames(logSums_table_new)[5])
log2foldchange_summary_df[new_row6, 1] <- print(rownames(logSums_table_new)[6])
log2foldchange_summary_df[new_row7, 1] <- print(rownames(logSums_table_new)[7])
log2foldchange_summary_df[new_row8, 1] <- print(rownames(logSums_table_new)[8])
log2foldchange_summary_df[new_row9, 1] <- print(rownames(logSums_table_new)[9])
log2foldchange_summary_df[new_row10, 1] <- print(rownames(logSums_table_new)[10])

log2foldchange_summary_df[new_row, 2] <- print(logSums_table_new[1,i])
log2foldchange_summary_df[new_row2, 2] <- print(logSums_table_new[2,i])
log2foldchange_summary_df[new_row3, 2] <- print(logSums_table_new[3,i])
log2foldchange_summary_df[new_row4, 2] <- print(logSums_table_new[4,i])
log2foldchange_summary_df[new_row5, 2] <- print(logSums_table_new[5,i])
log2foldchange_summary_df[new_row6, 2] <- print(logSums_table_new[6,i])
log2foldchange_summary_df[new_row7, 2] <- print(logSums_table_new[7,i])
log2foldchange_summary_df[new_row8, 2] <- print(logSums_table_new[8,i])
log2foldchange_summary_df[new_row9, 2] <- print(logSums_table_new[9,i])
log2foldchange_summary_df[new_row10, 2] <- print(logSums_table_new[10,i])
}

colnames(log2foldchange_summary_df) = c("RBP", "sum of log2foldchange")
```

```{r}
write.csv(log2foldchange_summary_df, file = "which_RBP_max_per_cell_type.csv")
```

# Now let's find which RBPs are the most differentially expressed accross ALL cell types, not just in individual cells or comparisons:
```{r}
log2foldchange_per_RBP_colsums <- colSums(log2foldchange_per_RBP)
log2foldchange_per_RBP_colsums <- log2foldchange_per_RBP_colsums[order(log2foldchange_per_RBP_colsums, decreasing = TRUE)]
log2foldchange_per_RBP_colsums <- log2foldchange_per_RBP_colsums[order(log2foldchange_per_RBP_colsums, decreasing = FALSE)]
```

```{r}
write.csv(log2foldchange_per_RBP_colsums, file = "log2foldchange of all 484 RBPs.csv")
```

# Looks like WBGene00010046 and WBGene00016260 are the "most intersting" for now. Let's examine them each more closely:
```{r}
colors <- colorRampPalette(brewer.pal(20, "PuRd"))(250)
interesting_RBP_table <- logSums_RBP_table[c("F54D1.1", "rege-1"),]

pheatmap(interesting_RBP_table,
         border_color = "black",
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = T,
         display_numbers = F,
         fontsize = 7.5,
         angle_col = 45,
         col = colors,
         main = "Sum of All log2FoldChange Values for Interesting RBPs in One Cell Type vs All Other Cell Types")

interesting_RBP_table <- logSums_RBP_table[c("F54D1.1", "rege-1", "ZK616.1", "mbl-1", "sup-12", "Y73B6BL.29"),]

pheatmap(interesting_RBP_table,
         border_color = "black",
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = T,
         display_numbers = F,
         fontsize = 7.5,
         angle_col = 45,
         col = colors,
         main = "Sum of All log2FoldChange Values for Interesting RBPs in One Cell Type vs All Other Cell Types")

pheatmap(logSums_RBP_table[,c("AFD", "ASEL", "ASER")],
         border_color = "black",
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         show_colnames = T,
         display_numbers = F,
         fontsize = 7,
         angle_col = 45,
         col = colors,
         main = "Sum of All log2FoldChange Values for AFD, ASEL, and ASER vs All Other Cell Types")
```

Copyright 2024 The Regents of the University of California

All Rights Reserved

Created by Zachery Wolfe

Department of Biochemistry

This file is part of Differential Expression in C. elegans. \
Differential Expression in C. elegans is a free pipeline: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. \
Differential Expression in C. elegans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. \
You should have received a copy of the GNU General Public License along with Differential Expression in C. elegans. If not, see <https://www.gnu.org/licenses/>.
