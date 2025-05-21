## 1st ##

# See 'Type vs type DESeq comparisons.R' for statistical summaries and outputs

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# DESeq setup
```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
 
BiocManager::install("DESeq2")
library(DESeq2)
library(ggplot2)
library(stringr)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(broom)
library(ggfortify)
```

```{r}
save.image("C:/Your/Directory/.RData")
load("C:/Your/Directory/.RData")
```

# Create a count matrix using HTSeq results:
```{r}
# Identify the directory containing the .txt files
directory_path <- "C:/Your/Directory"
file_paths <- list.files(directory_path, pattern = "\\_counts.txt$", full.names = TRUE)

# Initialize variables
count_matrix <- NULL
sample_names <- NULL

# Iterate through the files and extract the counts
for (i in file_paths) {
  sample_name <- tools::file_path_sans_ext(basename(i))  # Extract sample name
  sample_names <- c(sample_names, sample_name)
  
  counts <- read_delim(i, delim = "\t", col_names = FALSE, col_types = cols_only(col_character(), col_integer()))
  count_matrix <- cbind(count_matrix, counts[[2]])
}

# Convert count matrix to matrix format, then add replicate names and gene names:
count_matrix <- as.matrix(count_matrix)
count_matrix <- count_matrix[1:(nrow(count_matrix)-5),]
colnames(count_matrix) <- str_remove(sample_names, "_counts")
genenames <- read.csv("genenames.csv")
rownames(count_matrix) <- genenames[,1]

# Save count matrix as a .csv file
write.csv(count_matrix, file = "count_matrix.csv", row.names = TRUE)

# Print the count matrix and sample names
print(count_matrix)
print(sample_names)
```

# Transform your bulk RNA-seq data:
```{r}
t_PCA <- prcomp(t(count_matrix))
```

# Unweighted for read depth or STDev (will do in DESeq later):
```{r}
print(t_PCA)

plot(t_PCA$x, type = "p", cex = 1, pch = 16)

title(main = "PCA Plot of Transposed Gene Data PCs 1-2")
```

# Add your gene names to the PCA:
```{r}
rownames(t_PCA$rotation) <- make.names(1:nrow(genenames))

list(sort(t_PCA$rotation[,1], decreasing = TRUE)[1:10])
list(sort(t_PCA$rotation[,1], decreasing = FALSE)[1:10])

list(sort(t_PCA$rotation[,2], decreasing = TRUE)[1:10])
list(sort(t_PCA$rotation[,2], decreasing = FALSE)[1:10])

list(sort(t_PCA$rotation[,3], decreasing = TRUE)[1:10])
list(sort(t_PCA$rotation[,3], decreasing = FALSE)[1:10])
```


# Make your cts object and run DESeq:
```{r}
cts <- count_matrix
# make a condition table:
condition = read.table("experimental_conditions_46_neurons.csv", header = TRUE, sep=",", row.names=1)

dds <- DESeqDataSetFromMatrix(countData = cts, colData = condition, design= ~ type)
dds <- DESeq(dds)
resultsNames(dds)

dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
#normalized_countsdata <- write.csv(normalized_counts, file="normalized_counts.csv", sep="\t", quote = F, col.names = NA) # Use for UMAP later
```

```{r}
write.csv(normalized_counts, file = "normalized_counts.csv")
```


# Re-order results
```{r}
res <- results(dds)
res[order(res$padj),]
```


# Insert desired res object below:
```{r}
with(res, plot(log2FoldChange, -log10(pvalue), pch = 18, main="", xlim=c(-4,4), ylim=c(2,35)))

with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=18, col="black"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=18, col="red"))
```

```{r}
with(subset(res, padj<.01 ), print(log2FoldChange))
with(subset(res, padj<.01 ), print(res$rownames)) # rownames are NULL
with(res, print(pvalue)) # Already re-ordered from smallest to largest (2 chunks above)
```

# Now, we are only interested in significant results:
```{r}
sig_res <- res[(res$padj<.05 & abs(res$log2FoldChange)>2 & !is.na(res$padj) & !is.na(row.names(res))) ,]
write.csv(sig_res, file="sig_results")

plotMA(res, ylim=c(-10,10))
plotMA(sig_res, ylim=c(-10,10))
```

# PCA plots:
```{r}
vsd <- vst(dds, blind=FALSE)

plotPCA(vsd, intgroup="condition") 
plotPCA(vsd, intgroup="type")

normalized_counts_assay <- assay(vsd)
```

# PCA plot using ggplot (this is the same format - ggplot - as the PCAs of AS events we will make in the future):
```{r}
pca_data <- prcomp(t(assay(vsd)))
pca_var_perc <- pca_data$sdev^2 / sum(pca_data$sdev^2) * 100
pca_df <- data.frame(
  PC1 = pca_data$x[,1],
  PC2 = pca_data$x[,2],
  CellType = rownames(colData(vsd)),
  condition = colData(vsd)$condition
)

ggplot(pca_df, aes(x = PC1, y = PC2, label = CellType, color = condition)) +
  geom_point(size = 2) +
  geom_text(hjust = 1.25, vjust = 1.25, size = 3) +
  labs(
    title = "PCA by Neuron Type (DESeq results)",
    x = paste0("Principal Component 1 (", round(pca_var_perc[1], 2), "%)"),
    y = paste0("Principal Component 2 (", round(pca_var_perc[2], 2), "%)")) +
  theme_minimal() +
  scale_color_manual(values = c("Sensory" = "dodgerblue", "Motor" = "gold3", "Interneuron" = "salmon", "Polymodal" = "darkgreen")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
```

# Heatmaps of DESeq data:
```{r}
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:25]

df <- as.data.frame(colData(dds)[,c("condition","type")])
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, fontsize = 7.5)

df <- as.data.frame(colData(dds)[1:40,c("condition","type")])
pheatmap(assay(vsd)[select,1:40], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, fontsize = 7.5)

df <- as.data.frame(colData(dds)[41:80,c("condition","type")])
pheatmap(assay(vsd)[select,41:80], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, fontsize = 7.5)

df <- as.data.frame(colData(dds)[81:120,c("condition","type")])
pheatmap(assay(vsd)[select,81:120], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, fontsize = 7.5)

df <- as.data.frame(colData(dds)[121:160,c("condition","type")])
pheatmap(assay(vsd)[select,121:160], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, fontsize = 7.5)
```

# Get sample distances:
```{r}
sampleDists <- dist(t(assay(vsd)))
```

# Heatmaps of Individual Replicates:
```{r}
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$type, sep=":")
colnames(sampleDistMatrix) <- paste(vsd$type, sep=":")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(250)
pheatmap(sampleDistMatrix[1:40,1:40], # Change rows and columns for clearer view
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         fontsize = 7)

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$type, sep=":")
colnames(sampleDistMatrix) <- paste(vsd$type, sep=":")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(250)
pheatmap(sampleDistMatrix[41:80,41:80], # Change rows and columns for clearer view
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         fontsize = 7,
         col=colors)

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$type, sep=":")
colnames(sampleDistMatrix) <- paste(vsd$type, sep=":")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(250)
pheatmap(sampleDistMatrix[81:160,81:160], # Change rows and columns for clearer view
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         fontsize = 6,
         col=colors)

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$type, sep=":")
colnames(sampleDistMatrix) <- paste(vsd$type, sep=":")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(250)
pheatmap(sampleDistMatrix, # Change rows and columns for clearer view
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         fontsize = 7,
         col=colors)
```

Copyright 2024 The Regents of the University of California

All Rights Reserved

Created by Zachery Wolfe

Department of Biochemistry

This file is part of Differential Expression in C. elegans. \
Differential Expression in C. elegans is a free pipeline: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. \
Differential Expression in C. elegans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. \
You should have received a copy of the GNU General Public License along with Differential Expression in C. elegans. If not, see <https://www.gnu.org/licenses/>.
