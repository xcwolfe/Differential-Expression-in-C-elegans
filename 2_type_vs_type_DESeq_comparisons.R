

#SECOND

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


```{r}
library(stringr)
library(magrittr)
library(dplyr)
library(pheatmap)
library(DESeq2)
library(RColorBrewer)
```

# In this script, we will be looping a data frame that contains DESeq results. It will look something like the chunk below, but at a much larger scale:
## Use this chunk to practice before committing to all cell types vs all cell types:
```{r}
# By type, ADL vs ASG:

contrast <- c("type","ADL","ASG")
res <- results(dds, contrast = contrast)

res <- res[order(res$padj),]
sum <- summary(res)

sumwrite <- write.csv(sum, file="ADLvsASGsummary.csv") #No data in excel file. See sink() below
reswrite <- write.csv(res, file="ADLvsASG.csv")

sumsink <- sink("sumADL_ASG.txt")
resprint <- print(summary(res))
sink()
```

## Basically, we need to take our DESeq results for all cell types and compare them against one another, then print the number of differentially expressed genes for each cell x cell comparison.

# Start by saving your cell types below:
```{r}
singletypes <- c("ADL","AFD","AIM","AIN","AIY","ASEL","ASER","ASG","ASI","ASK","AVA","AVE","AVG","AVH","AVK","AVL","AVM","AWA","AWB","AWC","BAG","CAN","DA","DD","DVC","I5","IL1","IL2","LUA","NSM","OLL","OLQ","PHA","PVC","PVD","PVM","RIA","RIC","RIM","RIS","RMD","SMB","SMD","VB","VC","VD")
```


# Loop through all the DESeq data:
## This took me ~48 hours to run, nearly 3,300 files in total. Use an external hard drive and make sure to keep your compurter on and connected.
```{r}
setwd("D:/Zach Wolfe's DESeq analysis") # setwd() change to external hard drive

new_results <- matrix(nrow = 46, ncol = 46, dimnames = list(singletypes, singletypes))
new_results[is.na(new_results)] <- 0

res_list <- list()

for (i in singletypes) {
  for (j in singletypes) {
    if (i != j) {
      contrast <- c("type", i, j)
      res <- results(dds, contrast = contrast)
      
      res <- res[order(abs(res$log2FoldChange), decreasing = TRUE),]
      sum <- summary(res)
      res_list[[paste(i, "vs", j, sep = "_")]] <- res
      #sink(file = paste("sum", i, "vs", j, ".txt", sep = ""))
      print(summary(res))
      #sink()
      
      new_results[i, j] <- length(which(abs(res$log2FoldChange) > 2 & res$padj < .05))
    }
    if (i == j) next
  }
}

# Write CSV files outside the loop
for (name in names(res_list)) {
  file_name <- paste("res", name, ".csv", sep = "_")
  write.csv(res_list[[name]], file = file_name)
}
```

# write/read.csv backup:
```{r}
#write.csv(new_results, file = "new_results_all_46_cell_types.csv")
#new_results <- read.csv("new_results_all_46_cell_types.csv")
```

# Heatmap plots:
```{r}
pheatmap(new_results, cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, fontsize = 7.5)

pdf(file = "heatmap with clusters.pdf", width = 10, height = 10)

colors <- colorRampPalette(brewer.pal(6, "Reds"))(300)
pheatmap(new_results,
         cluster_rows = F,
         cluster_cols = F,
         legend = F,
         fontsize_col = 7.5,
         fontsize_row = 7.5,
         col=colors,
         display_numbers = T, 
         number_format = "%.0f",
         number_color = "black",
         fontsize_number = 5)

dev.off()

colors <- colorRampPalette(brewer.pal(6, "Reds"))(300)
pheatmap(new_results,
         cluster_rows = T,
         cluster_cols = T,
         legend = F,
         fontsize_col = 9,
         fontsize_row = 9,
         col=colors,
         display_numbers = T, 
         number_format = "%.0f",
         number_color = "black",
         fontsize_number = 5)
```

Copyright 2024 The Regents of the University of California

All Rights Reserved

Created by Zachery Wolfe

Department of Biochemistry

This file is part of Differential Expression in C. elegans. \
Differential Expression in C. elegans is a free pipeline: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. \
Differential Expression in C. elegans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. \
You should have received a copy of the GNU General Public License along with Differential Expression in C. elegans. If not, see <https://www.gnu.org/licenses/>.
