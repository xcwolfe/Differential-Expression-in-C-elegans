## 5th ##

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
library("data.table")
library("tidyverse")
library("pheatmap")
library("rex")
library("modeest")
library("dplyr")
library("plyr")
library("grid")
library("gridBase")
library("stringi")
library("RColorBrewer")
library("purrr")
```

```{r}
save.image("C:/Your/directory/.RData")
load("C:/Your/directory/.RData")
```


# I want to display differential expression measures in a new table:
```{r}
setwd("D:/Zach Wolfe's DESeq analysis")

untrimmed_results <- matrix(nrow = 46, ncol = 46, dimnames = list(singletypes,singletypes))
untrimmed_results[is.na(untrimmed_results)] <- 0

for(i in singletypes){
  for(j in singletypes){
    if(i != j){

untrimmed_data <- fread(file = paste("res_", i, "_vs_", j, "_.csv", sep=""))[,c(1,3,7)]
#match(RBP_list, overlap_data$V1, nomatch = NA)

untrimmed_results[i,j] <- print(length(untrimmed_data) < 312) # 312 is top 1% significant of all genes (31173)

untrimmed_results[i,j] <- write.csv(untrimmed_data, file=paste("genes", i, "vs", j, ".csv", sep="")) # Write csv of results

    }
    if(i == j) next
    }
}
```

```{r}
comparison <- print(paste(ADL, "vs", VD))
# gene <- print(overlap_data$V1[which(match(RBP_list, overlap_data$V1) < 312)])
# log2foldchange <- print(overlap_data$log2FoldChange[which(match(RBP_list, overlap_data$V1) < 312)])
# pvalue <- print(overlap_data$padj[which(match(RBP_list, overlap_data$V1) < 312)])
```

# Print names of genes from chart above:
# This will be a practice chunk to ensure everything works:
```{r, include=FALSE}
setwd("D:/Zach Wolfe's DESeq analysis")

untrimmed_cutoff_results <- matrix(nrow = 46, ncol = 46, dimnames = list(singletypes,singletypes))
untrimmed_cutoff_results[is.na(untrimmed_cutoff_results)] <- 0

for(i in c("ADL")){
   for(j in c("AFD")){
    if(i != j){

overlap_data <- fread(file = paste("res_", i, "_vs_", j, "_.csv", sep = ""))[,c(1,3,7)] %>% as.data.frame()

rownames(overlap_data) <- overlap_data$V1
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
untrimmed_cutoff_results[i,j] <- print(nrow(cutoffgenes))

write.csv(cutoffgenes, file=paste("DiffExpGenes", i, "vs", j, ".csv", sep=""))
     }
    if (i == j) next
   }
}

untrimmed_cutoff_results[is.na(untrimmed_cutoff_results)] <- 0
```

# Get rid of the first row in topDEgene (created in line 115):
```{r}
topDEgene <- topDEgene[-1,] # run this AFTER the chunk above JUST ONCE, then "#" out this row
```

# Now let's run the loop in full:
```{r}
setwd("D:/Zach Wolfe's DESeq analysis")

untrimmed_cutoff_results <- matrix(nrow = 46, ncol = 46, dimnames = list(singletypes,singletypes))
untrimmed_cutoff_results[is.na(untrimmed_cutoff_results)] <- 0

for(i in singletypes){
   for(j in singletypes){
    if(i != j){

overlap_data <- fread(file = paste("res_", i, "_vs_", j, "_.csv", sep = ""))[,c(1,3,7)] %>% as.data.frame()

rownames(overlap_data) <- overlap_data$V1
intersection <- rownames(overlap_data)
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
untrimmed_cutoff_results[i,j] <- print(nrow(cutoffgenes))

write.csv(cutoffgenes, file=paste("DiffExpGenes", i, "vs", j, ".csv", sep=""))

     }
    if (i == j) next
   }
}

untrimmed_cutoff_results[is.na(untrimmed_cutoff_results)] <- 0
```

# Save/Change topDEgene object:
```{r}
write.csv(topDEgene, file = "topDifferentiallyExpressedGenesbyComparison.csv")
write.csv(untrimmed_cutoff_results, file = "DE_genes_results.csv")
```

# Backup in case something happens to the matrix above:
```{r}
#untrimmed_cutoff_results <- read.csv(file = "untrimmed_cutoff_results.csv", header = TRUE)
#rownames(untrimmed_cutoff_results) <- colnames(untrimmed_cutoff_results[2:42])
#untrimmed_cutoff_results <- untrimmed_cutoff_results[2:42]
```

```{r}
colors <- colorRampPalette(brewer.pal(6, "Reds"))(250)
plot.new()
pheatmap(untrimmed_cutoff_results,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         fontsize = 7.5,
         col = rev(colors),
         main = "Number of Differentially Expressed Genes by Cell Type (log2FoldChange > 2 and adj. p-value < 0.01)")

plot.new()
pheatmap(untrimmed_cutoff_results,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = T,
         number_format = "%.0f",
         number_color = "black",
         fontsize = 7.5,
         fontsize_number = 6.75,
         col = colors,
         main = "Number of Differentially Expressed Genes by Cell Type (log2FoldChange > 2 and adj. p-value < 0.01)")
```

## log2foldchange per gene:
```{r include=FALSE, results='hide'}
setwd("D:/Zach Wolfe's DESeq analysis")

log2foldchange_per_gene <-  matrix(nrow = 46, ncol = 46934, dimnames = list(singletypes, genenames$Gene)) # genenames object is made in "CenGen bulk RNAseq Data PCA and DESeq"
log2foldchange_per_gene[is.na(log2foldchange_per_gene)] <- 0

for(j in singletypes){
 for(k in singletypes){
  if(j != k){

DiffExpGenes <- fread(file = paste("DiffExpGenes", j, "vs", k, ".csv", sep = ""))[,2:3] %>% as.data.frame()
log2foldchange_per_gene[j,DiffExpGenes$V1] <- print(DiffExpGenes$log2FoldChange) + log2foldchange_per_gene[j, DiffExpGenes$V1]
     }
      if (j == k) next
     }
   }

log2foldchange_per_gene_zeroes <- (colSums(log2foldchange_per_gene) == 0)
log2foldchange_per_gene <- log2foldchange_per_gene[, !log2foldchange_per_gene_zeroes]
```

# Which cell types have the most upregulation and downregulation?
```{r}
rownames(log2foldchange_per_gene)[rowSums(log2foldchange_per_gene) == max(rowSums(log2foldchange_per_gene))]
rownames(log2foldchange_per_gene)[rowSums(log2foldchange_per_gene) == min(rowSums(log2foldchange_per_gene))]

upreg_cell_type <- rownames(log2foldchange_per_gene)[rowSums(log2foldchange_per_gene) == max(rowSums(log2foldchange_per_gene))]
downreg_cell_type <- rownames(log2foldchange_per_gene)[rowSums(log2foldchange_per_gene) == min(rowSums(log2foldchange_per_gene))]

log2foldchange_per_gene[upreg_cell_type,]
log2foldchange_per_gene[downreg_cell_type,]
```

# Which genes have the greatest difference in regulation?
```{r}
log2foldchange_per_gene_colnames <- colnames(log2foldchange_per_gene)

for(k in 1:ncol(log2foldchange_per_gene)){
print(log2foldchange_per_gene_colnames[k])
print(max(log2foldchange_per_gene[,k]-min(log2foldchange_per_gene[,k])))
}
```

```{r}
colors <- colorRampPalette(brewer.pal(6, "Reds"))(250)
pheatmap(log2foldchange_per_gene,
         border_color = "black",
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         fontsize = 6,
         angle_col = 45,
         col = colors,
         main = "Total log2FoldChange Accross All Comparisons")

pheatmap(log2foldchange_per_gene[,1:50],
         border_color = "black",
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         fontsize = 7,
         angle_col = 90,
         col = colors,
         main = "Total log2FoldChange Accross All Comparisons (first 50)")

pheatmap(log2foldchange_per_gene[,51:100],
         border_color = "black",
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         fontsize = 7,
         angle_col = 90,
         col = colors,
         main = "Total log2FoldChange Accross All Comparisons (second 50)")
```

# write.csv for sum of log2foldchange per gene (in which abs(log2FoldChange) > 2 & padj < .01) in each cell type) - colSums(0) columns removed
```{r}
write.csv(log2foldchange_per_gene, file = "sum of log2foldchange per gene.csv")
```

# Backup read.csv:
```{r}
log2foldchange_per_gene <- read.csv("sum of log2foldchange per gene.csv", row.names = 1)
```


# Combine every excel file for each cell type:

# Note: If there are 0 Diff. Exp. Genes in a comparison, dplyr will assume that the rows in the .csv file are logical (we want them to be character/numeric). This will prevent us from merging rows in different .csv files together in the upcoming chunk. You may have to examine each Diff. Exp. Gene file and manually move them outside of your working directory if there are 0 rows:

## cell types with 0 Diff. Exp. Genes in my example can be found in D:\Zach Wolfe's DESeq analysis\Zero_Diff_Exp_Genes_comparisons 
```{r}
setwd("D:/Zach Wolfe's DESeq analysis")

ADLDiffExpGenefiles <- list.files(pattern = "DiffExpGenesADLvs*")
AFDDiffExpGenefiles <- list.files(pattern = "DiffExpGenesAFDvs*")
AIMDiffExpGenefiles <- list.files(pattern = "DiffExpGenesAIMvs*")
AINDiffExpGenefiles <- list.files(pattern = "DiffExpGenesAINvs*")
AIYDiffExpGenefiles <- list.files(pattern = "DiffExpGenesAIYvs*")
AINDiffExpGenefiles <- list.files(pattern = "DiffExpGenesAINvs*")
ASELDiffExpGenefiles <- list.files(pattern = "DiffExpGenesASELvs*")
ASERDiffExpGenefiles <- list.files(pattern = "DiffExpGenesASERvs*")
ASGDiffExpGenefiles <- list.files(pattern = "DiffExpGenesASGvs*")
ASIDiffExpGenefiles <- list.files(pattern = "DiffExpGenesASIvs*")
ASKDiffExpGenefiles <- list.files(pattern = "DiffExpGenesASKvs*")
AVADiffExpGenefiles <- list.files(pattern = "DiffExpGenesAVAvs*")
AVEDiffExpGenefiles <- list.files(pattern = "DiffExpGenesAVEvs*")
AVGDiffExpGenefiles <- list.files(pattern = "DiffExpGenesAVGvs*")
AVHDiffExpGenefiles <- list.files(pattern = "DiffExpGenesAVHvs*")
AVKDiffExpGenefiles <- list.files(pattern = "DiffExpGenesAVKvs*")
AVLDiffExpGenefiles <- list.files(pattern = "DiffExpGenesAVLvs*")
AVMDiffExpGenefiles <- list.files(pattern = "DiffExpGenesAVMvs*")
AWADiffExpGenefiles <- list.files(pattern = "DiffExpGenesAWAvs*")
AWBDiffExpGenefiles <- list.files(pattern = "DiffExpGenesAWBvs*")
AWCDiffExpGenefiles <- list.files(pattern = "DiffExpGenesAWCvs*")
BAGDiffExpGenefiles <- list.files(pattern = "DiffExpGenesBAGvs*")
CANDiffExpGenefiles <- list.files(pattern = "DiffExpGenesCANvs*")
DADiffExpGenefiles <- list.files(pattern = "DiffExpGenesDAvs*")
DDDiffExpGenefiles <- list.files(pattern = "DiffExpGenesDDvs*")
DVCDiffExpGenefiles <- list.files(pattern = "DiffExpGenesDVCvs*")
I5DiffExpGenefiles <- list.files(pattern = "DiffExpGenesI5vs*")
IL1DiffExpGenefiles <- list.files(pattern = "DiffExpGenesIL1vs*")
IL2DiffExpGenefiles <- list.files(pattern = "DiffExpGenesIL2vs*")
LUADiffExpGenefiles <- list.files(pattern = "DiffExpGenesLUAvs*")
NSMDiffExpGenefiles <- list.files(pattern = "DiffExpGenesNSMvs*")
OLLDiffExpGenefiles <- list.files(pattern = "DiffExpGenesOLLvs*")
OLQDiffExpGenefiles <- list.files(pattern = "DiffExpGenesOLQvs*")
PHADiffExpGenefiles <- list.files(pattern = "DiffExpGenesPHAvs*")
PVCDiffExpGenefiles <- list.files(pattern = "DiffExpGenesPVCvs*")
PVDDiffExpGenefiles <- list.files(pattern = "DiffExpGenesPVDvs*")
PVMDiffExpGenefiles <- list.files(pattern = "DiffExpGenesPVMvs*")
RIADiffExpGenefiles <- list.files(pattern = "DiffExpGenesRIAvs*")
RICDiffExpGenefiles <- list.files(pattern = "DiffExpGenesRICvs*")
RIMDiffExpGenefiles <- list.files(pattern = "DiffExpGenesRIMvs*")
RISDiffExpGenefiles <- list.files(pattern = "DiffExpGenesRISvs*")
RMDDiffExpGenefiles <- list.files(pattern = "DiffExpGenesRMDvs*")
SMBDiffExpGenefiles <- list.files(pattern = "DiffExpGenesSMBvs*")
SMDDiffExpGenefiles <- list.files(pattern = "DiffExpGenesSMDvs*")
VBDiffExpGenefiles <- list.files(pattern = "DiffExpGenesVBvs*")
VCDiffExpGenefiles <- list.files(pattern = "DiffExpGenesVCvs*")
VDDiffExpGenefiles <- list.files(pattern = "DiffExpGenesVDvs*")
```

# Make a list of all DiffExpGenefiles for each cell type:
```{r}
DiffExpGenefileList <- rep(paste0(singletypes, "DiffExpGenefiles"), len = 46)
```

# Log2FoldChange values for unique Gene names in individual "cell type vs" files:
```{r}
setwd("D:/Zach Wolfe's DESeq analysis")

LargeGenedf <- ADLDiffExpGenefiles %>% map_dfr(read.csv) # change cell type here
LargeGenedf <- LargeGenedf[,-c(1)]
uniqueGenes <- unique(LargeGenedf$V1)

logSums <- matrix(nrow = length(uniqueGenes), ncol = 1, dimnames = list(uniqueGenes, "log2foldchange"))
logSums[is.na(logSums)] <- 0

logMeans <- matrix(nrow = length(uniqueGenes), ncol = 1, dimnames = list(uniqueGenes, "log2foldchange"))
logMeans[is.na(logMeans)] <- 0

for(i in uniqueGenes){
logSums[i,1] <- sum(LargeGenedf[LargeGenedf$V1 == i, 2])
logMeans[i,1] <- mean(LargeGenedf[LargeGenedf$V1 == i, 2])
}
```

# Save/change cell type here:
```{r}
write.csv(logSums, file = "logSumsADL_all_genes.csv")
write.csv(logMeans, file = "logMeansADL_all_genes.csv") # change cell type here
#write.csv(LargeGenedf, file = "Most interesting Genes in cell type AIY.csv")
```

# Make a heatmap for the logSums value in each gene for a given cell type:
```{r}
setwd("D:/Zach Wolfe's DESeq analysis")

logSums_gene_table <- matrix(nrow = 46934, ncol = 46, dimnames = list(genenames$Gene,singletypes))
logSums_gene_table[is.na(logSums_gene_table)] <- 0

for (i in singletypes){
logSums_for_graph <- fread(file = paste0("logSums", i, "_all_genes.csv", sep = "")) %>% as.data.frame()
rownames(logSums_for_graph) <- logSums_for_graph[,1]
  for (j in genenames){
logSums_gene_table[j,i] <- print(logSums_for_graph[j,2])
logSums_gene_table[is.na(logSums_gene_table)] <- 0
  }
}

colors <- colorRampPalette(brewer.pal(6, "PuRd"))(250)
pheatmap(logSums_gene_table,
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
         main = "Sum of All log2FoldChange Values for Unique Genes in One Cell Type vs All Other Cell Types")

logSums_gene_table <- logSums_gene_table[rowSums(logSums_gene_table) != 0,] # if you want to keep rows with 0 data, "#" this row out
maxDEcelltypes <- as.data.frame(table(colnames(logSums_gene_table))[apply(logSums_gene_table, 1, which.max)])
maxDEcelltypes <- maxDEcelltypes[order(maxDEcelltypes$Freq, decreasing = TRUE),]
minDEcelltypes <- as.data.frame(table(colnames(logSums_gene_table))[apply(logSums_gene_table, 1, which.min)])
minDEcelltypes <- minDEcelltypes[order(minDEcelltypes$Freq, decreasing = TRUE),]
maxDEGenes <- as.data.frame(table(colnames(t(logSums_gene_table)))[apply(t(logSums_gene_table), 1, which.max)])
maxDEGenes <- count(maxDEGenes$Var1)
maxDEGenes <- maxDEGenes[order(maxDEGenes$freq, decreasing = TRUE),]
minDEGenes <- as.data.frame(table(colnames(t(logSums_gene_table)))[apply(t(logSums_gene_table), 1, which.min)])
minDEGenes <- count(minDEGenes$Var1)
minDEGenes <- minDEGenes[order(minDEGenes$freq, decreasing = TRUE),]
```

# Make a log2foldchange list just for each cell type here (useful for gene ontology tools):
```{r}
for (i in 1:length(singletypes)){

j <- i + 2
gene_log2foldchange_list <- logSums_gene_table
gene_log2foldchange_list <- rownames_to_column(as.data.table(gene_log2foldchange_list, "gene_names"))
gene_log2foldchange_list <- gene_log2foldchange_list[,2:j] # column 1 contains the row_names, column j contains my cell type's log2foldchange values
gene_log2foldchange_list_cut <- subset(gene_log2foldchange_list, select = c(1,ncol(gene_log2foldchange_list)))
gene_log2foldchange_list_cut <- gene_log2foldchange_list_cut[rev(order(gene_log2foldchange_list_cut[,2])),] # reorder by cell type's log2foldchange value column
gene_log2foldchange_list_cut <- as.data.frame(gene_log2foldchange_list_cut[which(gene_log2foldchange_list_cut[,2] > 200),]) # Change log2foldchange cutoffs if desired

write.table(gene_log2foldchange_list_cut[,1], file = paste("gene log2foldchange list", singletypes[i], ".tsv"), row.names = FALSE, col.names = FALSE, quote = FALSE) # this is the list we will input into a gene enrichment tool (FuncAssociate, g:Profiler, etc.)
}
```

```{r}
write.csv(logSums_gene_table, file = "logSums_gene_table_all_genes.csv")
write.csv(gene_log2foldchange_list_cut[,1], file = "gene log2foldchange list (AFD).csv", row.names = FALSE) # this is the list we will input into a gene enrichment tool (WormCat, g:Profiler, etc.)
write.csv(maxDEcelltypes, file = "maxDEcelltypes_all_genes.csv")
write.csv(minDEcelltypes, file = "minDEcelltypes_all_genes.csv")
write.csv(maxDEGenes, file = "maxDEGenes.csv")
write.csv(minDEGenes, file = "minDEGenes.csv")
```

# Backup:
```{r}
logSums_gene_table <- read.csv("logSums_gene_table_all_genes.csv", row.names = 1)
```

# max Diff. Exp. Gene heatmap:
```{r}
maxs <- rowMaxs(as.matrix(logSums_gene_table))
logSums_gene_table <- cbind(logSums_gene_table, maxs)
logSums_gene_table <- logSums_gene_table[order(unlist(logSums_gene_table[,47]), decreasing = TRUE),]
logSums_gene_table <- logSums_gene_table[,-47]

write.csv(logSums_gene_table, file = "logSums_gene_table_all_genes.csv")

# Now may be a good time to change gene names from WBGeneXXXXXXXX to common gene names using WormBase's simple mine tool:
gene_common_names <- read.table("simplemine_results_all_genes.txt", sep = "\t", header = T) # see "gene_common_names.csv" for a list of the genes used in this portion of the analysis

# We also have some duplicate gene names - they are all called "not found". It is easier to call these by their WBGeneID, so we will correct this using this loop:
duplicates <- gene_common_names$Public.Name[duplicated(gene_common_names$Public.Name)]

for (i in duplicates) {
  indices <- which(gene_common_names$Public.Name == i)
  gene_common_names$Public.Name[indices] <- gene_common_names$Your.Input[indices]
}

rownames(logSums_gene_table) <- gene_common_names$Public.Name

col_pallete <- colorRampPalette(c("white", "white", "white", "white", "gold", "darkgoldenrod1", "darkorange1"))
colors <- col_pallete(250)
pheatmap(logSums_gene_table[1:30,],
         border_color = "black",
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         number_format = "%.2f",
         number_color = "black",
         fontsize = 7,
         fontsize_number = 6,
         angle_col = 45,
         col = colors,
         main = "Sum of All log2FoldChange Values for Most Dysregulated Genes in One Cell Type vs All Other Cell Types")
```

# Heatmap for ASK:
```{r}
# temporarily re-order logSums_gene_table by ASK:
logSums_gene_table <- logSums_gene_table[order(unlist(logSums_gene_table[,10]), decreasing = TRUE),]

# Define the specific range for the gene matrix
min_matrix <- min(logSums_gene_table[1:30,])
max_matrix <- max(logSums_gene_table[1:30,])

# Define the breaks for the gene matrix
breaks_matrix <- c(seq(min_matrix, 0, by = abs(min(logSums_gene_table[1:30,])/30)), seq(0, max_matrix, by = max(logSums_gene_table[1:30,])/30), 0.001)
breaks_matrix[breaks_matrix == 0] <- 0.001  # Set 0 as a small positive value to represent white
breaks_matrix <- unique(breaks_matrix)

# Calculate the number of breaks for negative and positive values separately
n_breaks_neg <- sum(breaks_matrix < 0)
n_breaks_pos <- length(breaks_matrix) - n_breaks_neg

# Generate the color palette for negative values in the matrix using colorRampPalette
colors_neg <- colorRampPalette(c("white"))(n = n_breaks_neg)

# Generate the color palette for positive values in the matrix using colorRampPalette
colors_pos <- colorRampPalette(c("white", "white", "lightyellow","lightyellow", "khaki1", "khaki1", "gold", "gold2", "darkgoldenrod1", "darkorange1"))(n = n_breaks_pos)

# Combine the color palettes for negative and positive values
colors_matrix <- c(colors_neg, colors_pos)

# Define breaks
breaks <- c(seq(min_matrix, 0, length.out = 25), seq(0, max_matrix, length.out = 25)[-1])
breaks <- unique(breaks)

# Define colors
colors_neg <- colorRampPalette(c("darkorchid4", "deepskyblue", "cyan", "darkslategray1","white"))(25)
colors_pos <- colorRampPalette(c("white", "lightyellow", "khaki1", "khaki1", "gold", "gold", "darkgoldenrod1", "darkorange1"))(25)

# Adjust color for zero to be white
colors <- c(colors_neg, "#FFFFFF", colors_pos)

# Create the heatmap
p <- pheatmap(logSums_gene_table[1:30,],
         border_color = "black",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = FALSE,
         number_format = "%.2f",
         number_color = "black",
         fontsize = 7,
         fontsize_row = 9.5,
         fontsize_col = 9.5,
         fontsize_number = 6,
         angle_col = 90,
         col = colors,
         breaks = breaks,
         main = "Uniqueness Index for Most Dysregulated Genes in One Cell Type vs All Other Cell Types (Sorted by ASK)")

print(p)
```

# invert for a gene-centric uniqueness index:
```{r}
# temporarily re-order logSums_gene_table by mec-2:
gene_centric_logSums_gene_table <- t(logSums_gene_table)
gene_centric_logSums_gene_table <- gene_centric_logSums_gene_table[order(unlist(gene_centric_logSums_gene_table[,2029]), decreasing = TRUE),]

# Define the specific range for the gene matrix
min_matrix <- min(gene_centric_logSums_gene_table[,c(772,2029,2035,2695,4703,5075,5085,5105,5626,6567,7128,16695)])
max_matrix <- max(gene_centric_logSums_gene_table[,c(772,2029,2035,2695,4703,5075,5085,5105,5626,6567,7128,16695)])

# Define the breaks for the gene matrix
breaks_matrix <- c(seq(min_matrix, 0, by = abs(min(gene_centric_logSums_gene_table[,c(772,2029,2035,2695,4703,5075,5085,5105,5626,6567,7128,16695)])/12)), seq(0, max_matrix, by = max(gene_centric_logSums_gene_table[,c(772,2029,2035,2695,4703,5075,5085,5105,5626,6567,7128,16695)])/12), 0.001)
breaks_matrix[breaks_matrix == 0] <- 0.001  # Set 0 as a small positive value to represent white
breaks_matrix <- unique(breaks_matrix)

# Calculate the number of breaks for negative and positive values separately
n_breaks_neg <- sum(breaks_matrix < 0)
n_breaks_pos <- length(breaks_matrix) - n_breaks_neg

# Generate the color palette for negative values in the matrix using colorRampPalette
colors_neg <- colorRampPalette(c("white"))(n = n_breaks_neg)

# Generate the color palette for positive values in the matrix using colorRampPalette
colors_pos <- colorRampPalette(c("white", "white", "lightyellow","lightyellow", "khaki1", "khaki1", "gold", "gold2", "darkgoldenrod1", "darkorange1"))(n = n_breaks_pos)

# Combine the color palettes for negative and positive values
colors_matrix <- c(colors_neg, colors_pos)

# Define breaks
breaks <- c(seq(min_matrix, 0, length.out = 25), seq(0, max_matrix, length.out = 25)[-1])
breaks <- unique(breaks)

# Define colors
colors_neg <- colorRampPalette(c("darkorchid4", "deepskyblue", "cyan", "darkslategray1","white"))(25)
colors_pos <- colorRampPalette(c("white", "lightyellow", "khaki1", "khaki1", "gold", "gold", "darkgoldenrod1", "darkorange1"))(25)

# Adjust color for zero to be white
colors <- c(colors_neg, "#FFFFFF", colors_pos)

# Create the heatmap
p <- pheatmap(gene_centric_logSums_gene_table[,c(772,2029,2035,2695,4703,5075,5085,5105,5626,6567,7128,16695)],
         border_color = "black",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = FALSE,
         number_format = "%.2f",
         number_color = "black",
         fontsize = 7,
         fontsize_row = 7,
         fontsize_col = 9.5,
         fontsize_number = 6,
         angle_col = 90,
         col = colors,
         breaks = breaks,
         main = "Uniqueness Index for Most Representative Cell Types of One Gene (mec-2) vs Other Selected Genes")
print(p)
```

# do the same for all genes in a loop:
```{r}
for (i in c(772,1562,2029,2035,2695,4703,5075,5085,5105,5626,6567,7128,16695)){
gene_centric_logSums_gene_table <- gene_centric_logSums_gene_table[order(unlist(gene_centric_logSums_gene_table[,i]), decreasing = TRUE),]

# Define the specific range for the gene matrix
min_matrix <- min(gene_centric_logSums_gene_table[,i])
max_matrix <- max(gene_centric_logSums_gene_table[,i])

# Define the breaks for the gene matrix
breaks_matrix <- c(seq(min_matrix, 0, by = abs(min(gene_centric_logSums_gene_table[,i]))), seq(0, max_matrix, by = max(gene_centric_logSums_gene_table[,i])), 0.001)
breaks_matrix[breaks_matrix == 0] <- 0.001  # Set 0 as a small positive value to represent white
breaks_matrix <- unique(breaks_matrix)

# Calculate the number of breaks for negative and positive values separately
n_breaks_neg <- sum(breaks_matrix < 0)
n_breaks_pos <- length(breaks_matrix) - n_breaks_neg

# Generate the color palette for negative values in the matrix using colorRampPalette
colors_neg <- colorRampPalette(c("white"))(n = n_breaks_neg)

# Generate the color palette for positive values in the matrix using colorRampPalette
colors_pos <- colorRampPalette(c("white", "white", "lightyellow","lightyellow", "khaki1", "khaki1", "gold", "gold2", "darkgoldenrod1", "darkorange1"))(n = n_breaks_pos)

# Combine the color palettes for negative and positive values
colors_matrix <- c(colors_neg, colors_pos)
# Define breaks
breaks <- c(seq(min_matrix, 0, length.out = 25), seq(0, max_matrix, length.out = 25)[-1])
breaks <- unique(breaks)

# Define colors
colors_neg <- colorRampPalette(c("darkorchid4", "deepskyblue", "cyan", "darkslategray1","white"))(25)
colors_pos <- colorRampPalette(c("white", "lightyellow", "khaki1", "khaki1", "gold", "gold", "darkgoldenrod1", "darkorange1"))(25)

# Adjust color for zero to be white
colors <- c(colors_neg, "#FFFFFF", colors_pos)

# force a matrix of the single column to plot
single_col_matrix <- as.matrix(gene_centric_logSums_gene_table[, i])
rownames(single_col_matrix) <- rownames(gene_centric_logSums_gene_table)

p <- pheatmap(single_col_matrix,
         border_color = "black",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = FALSE,
         number_format = "%.2f",
         number_color = "black",
         show_rownames = T,
         fontsize = 7,
         fontsize_row = 7,
         fontsize_col = 9.5,
         fontsize_number = 6,
         angle_col = 90,
         col = colors,
         breaks = breaks,
         cellwidth = 45,
         main = paste0("Uniqueness Index for ", colnames(gene_centric_logSums_gene_table)[i]))
print(p)
}
```

# save the gene-centric uniqueness index as a .csv:
```{r}
write.csv(gene_centric_logSums_gene_table, file = "gene_centric_uniqueness_index.csv", row.names = T)
```

# Check relative expression of AFD compared to other cell types:
```{r}
colSums(logSums_gene_table)["AFD"]/nrow(logSums_gene_table)
expression_value_gene <- matrix(nrow = ncol(logSums_gene_table), ncol = 2)
colnames(expression_value_gene) <- c("Cell_type", "Normalized_logSums_value")
expression_value_gene[,1] <- colnames(logSums_gene_table)

for (i in 1:ncol(logSums_gene_table)){
  expression_value_gene[i,2] <- format(colSums(logSums_gene_table)[i]/nrow(logSums_gene_table), scientific = F)
}

ggplot(data = data.frame(expression_value_gene), aes(x = expression_value_gene[,1], y = as.numeric(expression_value_gene[,2]))) +
  geom_bar(stat = "identity", fill = "lightgreen") +
  labs(title = "Normalized logSums value of genes per cell type", x = "Cell type", y = "Normalized logSums value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

# Which Genes have the highest value overall in a single cell?
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

logSums_gene_table_new <- logSums_gene_table[order(logSums_gene_table[,i], decreasing = TRUE),]
log2foldchange_summary_df[new_row, 1] <- print(rownames(logSums_gene_table_new)[1])
log2foldchange_summary_df[new_row2, 1] <- print(rownames(logSums_gene_table_new)[2])
log2foldchange_summary_df[new_row3, 1] <- print(rownames(logSums_gene_table_new)[3])
log2foldchange_summary_df[new_row4, 1] <- print(rownames(logSums_gene_table_new)[4])
log2foldchange_summary_df[new_row5, 1] <- print(rownames(logSums_gene_table_new)[5])
log2foldchange_summary_df[new_row6, 1] <- print(rownames(logSums_gene_table_new)[6])
log2foldchange_summary_df[new_row7, 1] <- print(rownames(logSums_gene_table_new)[7])
log2foldchange_summary_df[new_row8, 1] <- print(rownames(logSums_gene_table_new)[8])
log2foldchange_summary_df[new_row9, 1] <- print(rownames(logSums_gene_table_new)[9])
log2foldchange_summary_df[new_row10, 1] <- print(rownames(logSums_gene_table_new)[10])

log2foldchange_summary_df[new_row, 2] <- print(logSums_gene_table_new[1,i])
log2foldchange_summary_df[new_row2, 2] <- print(logSums_gene_table_new[2,i])
log2foldchange_summary_df[new_row3, 2] <- print(logSums_gene_table_new[3,i])
log2foldchange_summary_df[new_row4, 2] <- print(logSums_gene_table_new[4,i])
log2foldchange_summary_df[new_row5, 2] <- print(logSums_gene_table_new[5,i])
log2foldchange_summary_df[new_row6, 2] <- print(logSums_gene_table_new[6,i])
log2foldchange_summary_df[new_row7, 2] <- print(logSums_gene_table_new[7,i])
log2foldchange_summary_df[new_row8, 2] <- print(logSums_gene_table_new[8,i])
log2foldchange_summary_df[new_row9, 2] <- print(logSums_gene_table_new[9,i])
log2foldchange_summary_df[new_row10, 2] <- print(logSums_gene_table_new[10,i])
}

colnames(log2foldchange_summary_df) = c("Gene", "sum of log2foldchange")
```

```{r}
write.csv(log2foldchange_summary_df, file = "which_Gene_max_per_cell_type.csv")
```

# Now let's find which Genes are the most differentially expressed accross ALL cell types, not just in individual cells or comparisons:
```{r}
log2foldchange_per_gene_colsums <- colSums(log2foldchange_per_gene)
log2foldchange_per_gene_colsums <- log2foldchange_per_gene_colsums[order(log2foldchange_per_gene_colsums, decreasing = TRUE)]
log2foldchange_per_gene_colsums <- log2foldchange_per_gene_colsums[order(log2foldchange_per_gene_colsums, decreasing = FALSE)]
```

```{r}
write.csv(log2foldchange_per_gene_colsums, file = "log2foldchange of all Genes.csv")
```

# Looks like WBGene00010046 and WBGene00016260 are the "most intersting" for now. Let's examine them each more closely:
```{r}
colors <- colorRampPalette(brewer.pal(9, "PuRd"))(250)
interesting_Gene_table <- logSums_gene_table[c("F54D1.1", "rege-1"),]

pheatmap(interesting_Gene_table,
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
         main = "Sum of All log2FoldChange Values for Interesting Genes in One Cell Type vs All Other Cell Types")

interesting_Gene_table <- logSums_gene_table[c("F54D1.1", "rege-1", "ZK616.1", "mbl-1", "sup-12", "Y73B6BL.29"),]

pheatmap(interesting_Gene_table,
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
         main = "Sum of All log2FoldChange Values for Interesting Genes in One Cell Type vs All Other Cell Types")

pheatmap(logSums_gene_table[,c("AFD", "ASEL", "ASER")],
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
