## 8th ##

```{r}
library("ggplot2")
```

### Example/training:
# Create a dataframe that contains the count of AS events that pass statistical thresholds in XXX vs YYY:
# (used for linear regression later)
```{r}
# change cell type/comparisons here:
singletypes <- c("ADL","AFD","AIM","AIN","AIY","ASEL","ASER","ASG","ASI","ASK","AVA","AVE","AVG","AVH","AVK","AVL","AVM","AWA","AWB","AWC","BAG","CAN","DA","DD","DVC","I5","IL1","IL2","LUA","NSM","OLL","OLQ","PHA","PVC","PVD","PVM","RIA","RIC","RIM","RIS","RMD","SMB","SMD","VB","VC","VD")

lin_reg_df_cas_vs_int <- data.frame(comparison = character(length(singletypes)), cassette_events = numeric(length(singletypes)), intron_events = numeric(length(singletypes))) # change cell type/comparisons here

# Loop to compile cassette exon events:
for (i in 1:length(singletypes)) {
  lin_reg_df_cas_vs_int[i, 1] <- paste("ADL vs", singletypes[i])
  lin_reg_df_cas_vs_int[i, 2] <- JUM_table_cassette[1, i] # row 1 corresponds to ADL, column 2 corresponds to AFD, etc.
  lin_reg_df_cas_vs_int[i, 3] <- JUM_table_intron[1, i] # row 1 corresponds to ADL, column 2 corresponds to AFD, etc.
}
```

# Linear regression:
```{r}
summary(lm(formula = cassette_events ~ intron_events, data = lin_reg_df_cas_vs_int))

ggplot(lin_reg_df_cas_vs_int, aes(intron_events, cassette_events, label = comparison)) +
  geom_point() +
  geom_text(hjust = 0.5, vjust = -0.5, color = "blue", size = 2.5) +
  geom_smooth(method = 'lm', color='red') +
  labs(title = "Linear Regression Model of Cassette Exon Events vs Intron Retention Events in ADL vs Individual Cell Types", x = "Intron Retention Events", y = "Cassette Exon Events") +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = 'bold'), axis.title.x = element_text(size = 9, face = 'bold'), axis.title.y = element_text(size = 9, face = 'bold'))
```


# Let's make a list of all cell types vs all cell types so that we can easily create our first column in future linear regression data frames:
```{r}
comparison_list <- list()

for (i in 1:length(singletypes)) {
  for (j in 1:length(singletypes)) {
   comparison_list[[length(comparison_list) + 1]] <- paste(singletypes[i], "vs", singletypes[j])
  }
}
```

# New linear regression data frame(s) with compiled A.S. events in each column:
```{r}
lin_reg_df_AS_compilation <- data.frame(A3S_events = numeric(), A5S_events = numeric(), cassette_events = numeric(), composite_events = numeric(), intron_events = numeric(), MXE_events = numeric(), sum_of_all_A.S._events = numeric())

lin_reg_df_AS_compilation <- lin_reg_df_AS_compilation[1:(length(singletypes)^2),]
rownames(lin_reg_df_AS_compilation) <- comparison_list

# First, compile all A3S, A5S, cassette, composite, intron, and MXE events together:
for (i in 1:nrow(lin_reg_df_AS_compilation)){
lin_reg_df_AS_compilation[i,1] <- JUM_table_A3S[i]
lin_reg_df_AS_compilation[i,2] <- JUM_table_A5S[i]
lin_reg_df_AS_compilation[i,3] <- JUM_table_cassette[i]
lin_reg_df_AS_compilation[i,4] <- JUM_table_composite[i]
lin_reg_df_AS_compilation[i,5] <- JUM_table_intron[i]
lin_reg_df_AS_compilation[i,6] <- JUM_table_MXE[i]
# Let's also make a column that includes the sum of all A.S. events for each comparison:
lin_reg_df_AS_compilation[i,7] <- sum(lin_reg_df_AS_compilation[i,2:6])
}

# We now have to remove all rows which are in alphabetical order because JUM doesn't calculate these properly:
rows_to_remove <- c()

for (i in 1:nrow(lin_reg_df_AS_compilation)) {
  rowname <- rownames(lin_reg_df_AS_compilation)[i]
  row_split <- strsplit(rowname, " vs ")[[1]]
  if (row_split[1] >= row_split[2]) {
    rows_to_remove <- c(rows_to_remove, i)
  }
}

lin_reg_df_AS_compilation <- lin_reg_df_AS_compilation[-rows_to_remove,]
```

# New linear regression data frame(s) in relation to differential gene and RBP expression:
```{r}
lin_reg_df_AS_vs_DE_genes_RBPs <- data.frame(`A.S._Events` = numeric(), `Differentially_Expressed_Genes` = numeric(), `Differentially_Expressed_RBPs` = numeric())

lin_reg_df_AS_vs_DE_genes_RBPs <- lin_reg_df_AS_vs_DE_genes_RBPs[1:(length(singletypes)^2),]
rownames(lin_reg_df_AS_vs_DE_genes_RBPs) <- comparison_list

# Add all types of significantly dysregulated A.S. events together (cassettes + introns + 5' + 3' + mutually exclusive + composite)

# First, compile all A.S. events together:
row_indices <- c(1:46)  # Row indices for JUM tables
column_indices <- c(1:46)  # Column indices for JUM tables

row_counter <- 1  # Counter for lin_reg_df_AS_vs_DE_genes_RBPs row index

for (row_index in row_indices) {
  for (col_index in column_indices) {
    lin_reg_df_AS_vs_DE_genes_RBPs[row_counter, 1] <- sum(JUM_table_A3S[row_index, col_index] +
                                                         JUM_table_A5S[row_index, col_index] +
                                                         JUM_table_cassette[row_index, col_index] +
                                                         JUM_table_composite[row_index, col_index] +
                                                         JUM_table_intron[row_index, col_index] +
                                                         JUM_table_MXE[row_index, col_index])
    lin_reg_df_AS_vs_DE_genes_RBPs[row_counter, 2] <- new_results[row_index, col_index]   # new_results was created in Type vs Type DESeq comparison.Rmd
    lin_reg_df_AS_vs_DE_genes_RBPs[row_counter, 3] <- RBP_cutoff_results[row_index, col_index]   # RBP_cutoff_results was created in RNA binding protein search.Rmd
    row_counter <- row_counter + 1
  }
}


# We now have to remove all rows which are in reverse alphabetical order because JUM doesn't calculate these properly:
lin_reg_df_AS_vs_DE_genes_RBPs <- lin_reg_df_AS_vs_DE_genes_RBPs[-rows_to_remove,]
```

# Linear regression:
```{r}
lm <- summary(lm(formula = Differentially_Expressed_Genes ~ A.S._Events, data = lin_reg_df_AS_vs_DE_genes_RBPs))

ggplot(lin_reg_df_AS_vs_DE_genes_RBPs, aes(A.S._Events, Differentially_Expressed_Genes, label = rownames(lin_reg_df_AS_vs_DE_genes_RBPs))) +
  geom_point() +
  geom_text(hjust = 0.5, vjust = -0.5, color = "blue", size = 2.5) +
  annotate("text", x = max(lin_reg_df_AS_vs_DE_genes_RBPs$A.S._Events), y = max(lin_reg_df_AS_vs_DE_genes_RBPs$Differentially_Expressed_Genes), label = paste0("Adj. R-squared = ", round(lm$adj.r.squared, 3)), hjust = 1, vjust = 1, color = "black", size = 4, fontface = "bold") +
  geom_smooth(method = 'lm', color='red') +
  labs(title = "Linear Regression Model of Alternative Splicing Events vs D.E. Genes by Cell Type Comparison", x = "Alternative Splicing Events", y = "Differentially Expressed Genes") +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = 'bold'), axis.title.x = element_text(size = 9, face = 'bold'), axis.title.y = element_text(size = 9, face = 'bold'))
```

# Linear regression RBPs vs A.S. events:
```{r}
lm <- summary(lm(formula = Differentially_Expressed_RBPs ~ A.S._Events, data = lin_reg_df_AS_vs_DE_genes_RBPs))

ggplot(lin_reg_df_AS_vs_DE_genes_RBPs, aes(A.S._Events, Differentially_Expressed_RBPs, label = rownames(lin_reg_df_AS_vs_DE_genes_RBPs))) +
  geom_point() +
  geom_text(hjust = 0.5, vjust = -0.5, color = "blue", size = 2.5) +
  annotate("text", x = max(lin_reg_df_AS_vs_DE_genes_RBPs$A.S._Events), y = max(lin_reg_df_AS_vs_DE_genes_RBPs$Differentially_Expressed_RBPs), label = paste0("Adj. R-squared = ", round(lm$adj.r.squared, 3)), hjust = 1, vjust = 1, color = "black", size = 4, fontface = "bold") +
  geom_smooth(method = 'lm', color='red') +
  labs(title = "Linear Regression Model of Alternative Splicing Events vs D.E. RBPs by Cell Type Comparison", x = "Alternative Splicing Events", y = "Differentially Expressed RNA Binding Proteins") +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = 'bold'), axis.title.x = element_text(size = 9, face = 'bold'), axis.title.y = element_text(size = 9, face = 'bold'))
```

# It would be much easier for us to do linear regressions if all of the data in the chunk(s) above were in the same data frame. Let's add the D.E. gene and RBP results into the AS_compilation lin_reg data frame. But be careful - we have to order everything correctly!
```{r}
lin_reg_df_AS_compilation$Diff_Expressed_Genes = 1:nrow(lin_reg_df_AS_compilation)
lin_reg_df_AS_compilation$Diff_Expressed_RBPs = 1:nrow(lin_reg_df_AS_compilation)

lin_reg_df_AS_compilation[,8] = lin_reg_df_AS_vs_DE_genes_RBPs[,2]
lin_reg_df_AS_compilation[,9] = lin_reg_df_AS_vs_DE_genes_RBPs[,3]
```

# Now we need to do linear regression for every column vs every column:
```{r}
lin_reg_df_heatmap <- matrix(nrow = 9, ncol = 9, dimnames = list(colnames(lin_reg_df_AS_compilation),colnames(lin_reg_df_AS_compilation)))
lin_reg_df_heatmap[is.na(lin_reg_df_heatmap)] <- 0

for (i in 1:9){      # There are 8 columns in lin_reg_df_AS_compilation
  for (j in 1:9){
lin_reg_df_heatmap[i,j] <- summary(lm(formula = lin_reg_df_AS_compilation[,i] ~ lin_reg_df_AS_compilation[,j], data = lin_reg_df_AS_compilation))$adj.r.squared
  }
}
```

# Heatmap:
```{r}
colors <- colorRampPalette(brewer.pal(8, "Blues"))(50)
pheatmap(lin_reg_df_heatmap,
         border_color = "black",
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = T,
         number_format = "%.4f",
         number_color = "black",
         fontsize = 9,
         fontsize_number = 10,
         angle_col = 90,
         col = colors,
         main = "Adj. R-squared Values of Linear Regression between A.S. Events")
```

# Linear regression cassette vs composite:
```{r}
lm <- summary(lm(formula = composite_events ~ cassette_events, data = lin_reg_df_AS_compilation))

ggplot(lin_reg_df_AS_compilation, aes(cassette_events, composite_events, label = row.names(lin_reg_df_AS_compilation))) +
  geom_point() +
  geom_text(hjust = 0.5, vjust = -0.5, color = "darkgreen", size = 2.5) +
  annotate("text", x = 50, y = max(lin_reg_df_AS_compilation$composite_events), label = paste0("Adj. R-squared = ", round(lm$adj.r.squared, 3)), hjust = 1, vjust = 1, color = "black", size = 4, fontface = "bold") +
  geom_smooth(method = 'lm', color='pink') +
  labs(title = "Linear Regression Model of Composite vs Cassette Exon Events by Cell Type Comparison", x = "Cassette Events", y = "Composite Events") +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = 'bold'), axis.title.x = element_text(size = 9, face = 'bold'), axis.title.y = element_text(size = 9, face = 'bold'))
```

# Linear regression D.E. genes and RBPs:
```{r}
lm <- summary(lm(formula = Diff_Expressed_Genes ~ sum_of_all_A.S._events, data = lin_reg_df_AS_compilation))

ggplot(lin_reg_df_AS_compilation, aes(sum_of_all_A.S._events, Diff_Expressed_Genes, label = row.names(lin_reg_df_AS_compilation))) +
  geom_point() +
  geom_text(hjust = 0.5, vjust = -0.5, color = "blue", size = 2.5) +
  annotate("text", x = 1500, y = 5000, label = paste0("Adj. R-squared = ", round(lm$adj.r.squared, 3)), hjust = 1, vjust = 1, color = "black", size = 4, fontface = "bold") +
  geom_smooth(method = 'lm', color='red') +
  labs(title = "Linear Regression Model of Alternative Splicing Events vs D.E. Genes by Cell Type Comparison", x = "Alternative Splicing Events", y = "Differentially Expressed Genes") +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = 'bold'), axis.title.x = element_text(size = 9, face = 'bold'), axis.title.y = element_text(size = 9, face = 'bold'))

lm <- summary(lm(formula = Diff_Expressed_RBPs ~ sum_of_all_A.S._events, data = lin_reg_df_AS_compilation))

ggplot(lin_reg_df_AS_compilation, aes(sum_of_all_A.S._events, Diff_Expressed_RBPs, label = row.names(lin_reg_df_AS_compilation))) +
  geom_point() +
  geom_text(hjust = 0.5, vjust = -0.5, color = "blue", size = 2.5) +
  annotate("text", x = 1400, y = 120, label = paste0("Adj. R-squared = ", round(lm$adj.r.squared, 3)), hjust = 1, vjust = 1, color = "black", size = 4, fontface = "bold") +
  geom_smooth(method = 'lm', color='red') +
  labs(title = "Linear Regression Model of Alternative Splicing Events vs D.E. RNABPs by Cell Type Comparison", x = "Alternative Splicing Events", y = "Differentially Expressed RNA-Binding Proteins") +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = 'bold'), axis.title.x = element_text(size = 9, face = 'bold'), axis.title.y = element_text(size = 9, face = 'bold'))
```

# Linear regression A5S vs A3S:
```{r}
lm <- summary(lm(formula = A5S_events ~ A3S_events, data = lin_reg_df_AS_compilation))

ggplot(lin_reg_df_AS_compilation, aes(A3S_events, A5S_events, label = row.names(lin_reg_df_AS_compilation))) +
  geom_point() +
  geom_text(hjust = 0.5, vjust = -0.5, color = "chocolate3", size = 2.5) +
  annotate("text", x = 100, y = 270, label = paste0("Adj. R-squared = ", round(lm$adj.r.squared, 3)), hjust = 1, vjust = 1, color = "black", size = 4, fontface = "bold") +
  geom_smooth(method = 'lm', color='magenta') +
  labs(title = "Linear Regression Model of A5S vs A3S Events by Cell Type Comparison", x = "A3S Events", y = "A5S Events") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'), axis.title.x = element_text(size = 12, face = 'bold'), axis.title.y = element_text(size = 12, face = 'bold'))

lm <- summary(lm(formula = A5S_events ~ cassette_events, data = lin_reg_df_AS_compilation))

ggplot(lin_reg_df_AS_compilation, aes(cassette_events, A5S_events, label = row.names(lin_reg_df_AS_compilation))) +
  geom_point() +
  geom_text(hjust = 0.5, vjust = -0.5, color = "chocolate3", size = 2.5) +
  annotate("text", x = 100, y = 270, label = paste0("Adj. R-squared = ", round(lm$adj.r.squared, 3)), hjust = 1, vjust = 1, color = "black", size = 4, fontface = "bold") +
  geom_smooth(method = 'lm', color='magenta') +
  labs(title = "Linear Regression Model of A5S vs Cassette Events by Cell Type Comparison", x = "cassette Events", y = "A5S Events") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'), axis.title.x = element_text(size = 12, face = 'bold'), axis.title.y = element_text(size = 12, face = 'bold'))
```

Copyright 2024 The Regents of the University of California

All Rights Reserved

Created by Zachery Wolfe

Department of Biochemistry

This file is part of Differential Expression in C. elegans. \
Differential Expression in C. elegans is a free pipeline: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. \
Differential Expression in C. elegans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. \
You should have received a copy of the GNU General Public License along with Differential Expression in C. elegans. If not, see <https://www.gnu.org/licenses/>.
