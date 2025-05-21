## 3rd ##

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("M3C")
library(umap)
library(tidyverse)
```

## UMAP, or Uniform Manifold Approximation and Projection, is a machine learning method by which we can generate visuals which may (or may not) predict the likelihood of a gene being upregulated, downregulated, or more related to another gene of interest. There are important rules and limitations to using this approach; tips can be found in McInnes, L, Healy, J, UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction, ArXiv e-prints 1802.03426, 2018

## https://umap-learn.readthedocs.io/en/latest/

# Load new_results (if you haven't kept them from the second .Rmd):
```{r}
new_results <- read.csv(file = "log2foldchange_2_p-value_05.csv", header = T, row.names = 1)
results.labels <- row.names(new_results)
```

# Set seed and UMAP parameters:
```{r}
set.seed(1111)
umap.defaults

# UMAP parameters set to default:
# n_neighbors = 15, min_dist = 0.1

# low n_neighbors values will push UMAP to focus more on local structure by constraining the number of neighboring points considered when analyzing the data in high dimensions, while high values will push UMAP towards representing the big-picture structure while losing fine detail.

#  Larger values of min_dist will make UMAP pack points together more loosely, focusing instead on preserving the broad topological structure
```

# Run this again if you don't want to use assay values:
```{r}
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE) # Use for UMAP later
normalized_countsdata <- write.csv(normalized_counts, file="normalized_counts.csv", quote=F)
normalized_counts[is.na(normalized_counts)] <- 0
```

# Make sure to transpose normalized counts before running the next chunk:
```{r}
normalized_counts <- normalized_counts %>%
  t()
```

# UMAP might look prettier if we included all replicates for each neuron type
```{r, echo=FALSE}
#Perform UMAP:
custom.settings = umap.defaults
for(i in c(4:7)){
    custom.settings$n_neighbors = i # n_neighbors value. See chunk above for notes on how to choose a good n_neighbors value
print(custom.settings)

umap_results = umap(normalized_counts, config = custom.settings, labels=as.factor(row.names(normalized_counts)), controlscale=TRUE, scale=3) # Changed new_results to normalized_counts
umap_results
umap_results$layout
plot(umap_results$layout,
          main = "A UMAP Visualization of Results (colored by type)",
          xlab = "Dimension 1",
          ylab = "Dimension 2",
          pch = 16,
          col = colors(as.factor(oldtype))[80:121]
          )

legend("right", legend = unique(as.factor(oldtype)), col = colors(as.factor(oldtype))[80:121], pch = 16)

plot(umap_results$layout,
          main = "A UMAP Visualization of Results (colored by condition)",
          xlab = "Dimension 1",
          ylab = "Dimension 2",
          pch = 16,
          col = colors(as.factor(oldcondition))[5:9]
          )

legend("right", inset = c(-0.1, -0.1), legend = unique(as.factor(oldcondition)), col = colors(as.factor(oldcondition))[5:9], pch = 16)
}
```

# Every time you run the assay, you will get slightly different results for UMAP

# New attempt with Alex's Lemonade:
```{r}
# Create the data folder if it doesn't exist
if (!dir.exists("data")) {
  dir.create("data")
}

# Define the file path to the plots directory
plots_dir <- "plots"

# Create the plots folder if it doesn't exist
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Define the file path to the results directory
results_dir <- "results"

# Create the results folder if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}
```


```{r}
# Make into data frame for plotting with `ggplot2`
# The UMAP values for plotting are stored in the `layout` element
umap_plot_df <- data.frame(umap_results$layout) %>%
  # Turn sample IDs stored as row names into a column
  tibble::rownames_to_column("type") %>%
  # Add the metadata into this data frame; match by sample IDs
  dplyr::inner_join(condition, by = "type")

umap_plot_df
```

```{r}
ggplot(
  umap_plot_df,
  aes(
    x = X1,
    y = X2,
    color = type,
    label = type
  )
) +
  geom_point() +
  geom_text(size = 2, hjust = -0.5, vjust = 1)

ggplot(
  umap_plot_df,
  aes(
    x = X1,
    y = X2,
    color = condition,
    label = type
  )
) +
  geom_point() +
  geom_text(size = 2, hjust = -0.5, vjust = 1)
```

```{r}
final_annotated_umap_plot <- ggplot(
  umap_plot_df, 
  aes(
    x = X1,
    y = X2,
    color = type,
    shape = condition
  )
) +
  geom_point(size = 4) +
  scale_shape_manual(values=c(16, 17, 8, 15, 7))

final_annotated_umap_plot
```

## UMAP notes:

## The coordinates of UMAP output for any given cell can change dramatically depending on parameters, and even run to run with the same parameters (Why setting the seed is important). This means that you should not rely too heavily on the exact values of UMAP’s output.

## One particular limitation of UMAP is that while observed clusters have some meaning, the distance and density between clusters usually does not. The fact that two clusters are near each other should NOT be interpreted to mean that they are more related to each other than to more distant clusters.

Copyright 2024 The Regents of the University of California

All Rights Reserved

Created by Zachery Wolfe

Department of Biochemistry

This file is part of Differential Expression in C. elegans. \
Differential Expression in C. elegans is a free pipeline: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. \
Differential Expression in C. elegans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. \
You should have received a copy of the GNU General Public License along with Differential Expression in C. elegans. If not, see <https://www.gnu.org/licenses/>.
