## 7th ##

## Note from Zach: You'll notice that many of these heatmap scripts are slightly varied forms of one another. This is because the .txt final output of JUM varies by splicing event - this is particularly true for columns containing deltaPSI values.

### If running many conditions, JUM will output incorrect data for cell types that are not in alphabetical order (i.e. ASG vs AVE will output correct data, but AVE vs ASG would not). You may need to hide or move spreadsheets (or boxes in the case of heatmaps) in order to return the correct 50% of data.

```{r}
library("pheatmap")
library("RColorBrewer")
library("stringr")
library("stringi")
library("data.table")
library("tidyverse")
library("tidyr")
library("splitstackshape")
library("MatrixGenerics")
library("purrr")
library("rlang")
```

```{r}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("IsoformSwitchAnalyzeR")
BiocManager::install("sparseMatrixStats")
BiocManager::install("DelayedMatrixStats")

if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}
```

```{r}
save.image("C:/Your/directory/.RData")
load("C:/Your/directory/.RData")
```

# List your cell types/conditions:
```{r}
singletypes <- c("ADL","AFD","AIM","AIN","AIY","ASEL","ASER","ASG","ASI","ASK","AVA","AVE","AVG","AVH","AVK","AVL","AVM","AWA","AWB","AWC","BAG","CAN","DA","DD","DVC","I5","IL1","IL2","LUA","NSM","OLL","OLQ","PHA","PVC","PVD","PVM","RIA","RIC","RIM","RIS","RMD","SMB","SMD","VB","VC","VD")

singletypes_withreps <- c("ADL1", "ADL2", "ADL3", "ADL4", "AFD1", "AFD2", "AFD3", "AFD4", "AFD5", "AIM1", "AIM2", "AIM3", "AIM4", "AIN1", "AIN2", "AIN3", "AIN4", "AIN5", "AIN6", "AIY1", "AIY2", "AIY3", "ASEL1", "ASEL2", "ASEL3", "ASER1", "ASER2", "ASER3", "ASER4", "ASG1", "ASG2", "ASG3", "ASG4", "ASI1", "ASI2", "ASI3", "ASI4", "ASK1", "ASK2", "ASK3", "ASK4", "AVA1", "AVA2", "AVA3", "AVA4", "AVA5", "AVA6", "AVE1", "AVE2", "AVE3", "AVG1", "AVG2", "AVG3", "AVH1", "AVH2", "AVH3", "AVH4", "AVK1", "AVK2", "AVK3", "AVK4", "AVL1", "AVL2", "AVL3", "AVM1", "AVM2", "AVM3", "AWA1", "AWA2", "AWA3", "AWA4", "AWB1", "AWB2", "AWB3", "AWB4", "AWB5", "AWC1", "AWC2", "AWC3", "AWC4", "BAG1", "BAG2", "BAG3", "BAG4", "CAN1", "CAN2", "CAN3", "DA1", "DA2", "DA3", "DA4", "DD1", "DD2", "DD3", "DD4", "DVC1", "DVC2", "DVC3", "DVC4", "I51", "I52", "I53", "I54", "IL11", "IL12", "IL13", "IL21", "IL22", "IL23", "IL24", "LUA1", "LUA2", "LUA3", "LUA4", "NSM1", "NSM2", "NSM3", "OLL1", "OLL2", "OLQ1", "OLQ2", "OLQ3", "PHA1", "PHA2", "PHA3", "PHA4", "PVC1", "PVC2", "PVC3", "PVC4", "PVC5", "PVD1", "PVD2", "PVM1", "PVM2", "RIA1", "RIA2", "RIA3", "RIA4", "RIA5", "RIA6", "RIC1", "RIC2", "RIC3", "RIC4", "RIM1", "RIM2", "RIM3", "RIM4", "RIS1", "RIS2", "RIS3", "RMD1", "RMD2", "RMD3", "RMD4", "RMD5", "SMB1", "SMB2", "SMB3", "SMB4", "SMB5", "SMD1", "SMD2", "SMD3", "SMD4", "VB1", "VB2", "VB3", "VB4", "VC1", "VC2", "VC3", "VC4", "VC5", "VC6", "VD1", "VD2", "VD3", "VD4")
```


##You can change any of these cutoffs as desired:

# A3S events:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

JUM_table_A3S <- matrix(nrow = 46, ncol = 46, dimnames = list(singletypes,singletypes))
JUM_table_A3S[is.na(JUM_table_A3S)] <- 0

for (i in singletypes){
 for (j in singletypes){
   if (i != j){
A3S <- fread(paste0("FINAL_JUM_OUTPUT_pvalue_1", i, "vs", j, "/AS_differential_JUM_output_A3SS_events_pvalue_1_final_simplified.txt"))

A3S_split <- cSplit(A3S, 9, sep = ";", stripWhite = TRUE, type.convert = FALSE)
A3S_split[!is.finite(as.numeric(unlist(A3S_split[,9]))),] <- 0
A3S_split[is.nan(as.character(unlist(A3S_split[,9]))),] <- 0
A3S_split <- A3S_split[A3S_split$Gene != 0]
A3S_p <- A3S_split[A3S_split$pvalue < 0.05]
A3S_q <- A3S_p[A3S_p$qvalue < 0.05]

A3S_maxPSI <- A3S_q[,9:ncol(A3S_q)]
A3S_maxPSI[is.na(A3S_maxPSI[,1:ncol(A3S_maxPSI)])] <- 0
A3S_colcount <- as.numeric(ncol(A3S_maxPSI))
A3S_rowcount <- as.numeric(nrow(A3S_maxPSI))
A3S_maxPSI[A3S_maxPSI == -Inf] <- 0
A3S_maxPSI[A3S_maxPSI == Inf] <- 0
A3S_maxPSI$max <- rowMaxs(abs(as.numeric(unlist(A3S_maxPSI))), dim. = as.integer(c(A3S_rowcount, A3S_colcount)))
A3S_sig <- A3S_maxPSI[A3S_maxPSI$max > 0.1,]

JUM_table_A3S[i,j] <- print(nrow(A3S_sig))
write.csv(A3S_q, file=paste("JUM_A3S", i, "vs", j, ".csv", sep=""), row.names = FALSE, quote = FALSE) # We save as A3S_q for now because it is easier to save in a spreadsheet. We will re-eliminate deltaPSI values below 0.1 in the "combination" chunk(s) later
   }
   if (i == j) next
 }
}
```

# A5S events:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

JUM_table_A5S <- matrix(nrow = 46, ncol = 46, dimnames = list(singletypes,singletypes))
JUM_table_A5S[is.na(JUM_table_A5S)] <- 0

for (i in singletypes){
 for (j in singletypes){
   if (i != j){
A5S <- fread(paste0("FINAL_JUM_OUTPUT_pvalue_1", i, "vs", j, "/AS_differential_JUM_output_A5SS_events_pvalue_1_final_simplified.txt"))

A5S_split <- cSplit(A5S, 9, sep = ";", stripWhite = TRUE, type.convert = FALSE)
A5S_split[!is.finite(as.numeric(unlist(A5S_split[,9]))),] <- 0
A5S_split[is.nan(as.character(unlist(A5S_split[,9]))),] <- 0
A5S_split <- A5S_split[A5S_split$Gene != 0]
A5S_p <- A5S_split[A5S_split$pvalue < 0.05]
A5S_q <- A5S_p[A5S_p$qvalue < 0.05]

A5S_maxPSI <- A5S_q[,9:ncol(A5S_q)]
A5S_maxPSI[is.na(A5S_maxPSI[,1:ncol(A5S_maxPSI)])] <- 0
A5S_colcount <- as.numeric(ncol(A5S_maxPSI))
A5S_rowcount <- as.numeric(nrow(A5S_maxPSI))
A5S_maxPSI[A5S_maxPSI == -Inf] <- 0
A5S_maxPSI[A5S_maxPSI == Inf] <- 0
A5S_maxPSI$max <- rowMaxs(abs(as.numeric(unlist(A5S_maxPSI))), dim. = as.integer(c(A5S_rowcount, A5S_colcount)))
A5S_sig <- A5S_maxPSI[A5S_maxPSI$max > 0.1,]

JUM_table_A5S[i,j] <- print(nrow(A5S_sig))
write.csv(A5S_q, file=paste("JUM_A5S", i, "vs", j, ".csv", sep=""), row.names = FALSE, quote = FALSE)  # We save as A5S_q for now because it is easier to save in a spreadsheet. We will re-eliminate deltaPSI values below 0.1 in the "combination" chunk(s) later
   }
   if (i == j) next
 }
}
```

# Cassette exon events:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

JUM_table_cassette <- matrix(nrow = 46, ncol = 46, dimnames = list(singletypes,singletypes))
JUM_table_cassette[is.na(JUM_table_cassette)] <- 0

for (i in singletypes){
 for (j in singletypes){
   if (i != j){
cassette <- fread(paste0("FINAL_JUM_OUTPUT_pvalue_1", i, "vs", j, "/AS_differential_JUM_output_cassette_exon_events_pvalue_1_final_simplified.txt"))
 
cassette[!is.finite(as.numeric(unlist(cassette[,11]))),] <- 0
cassette[is.nan(as.character(unlist(cassette[,11]))),] <- 0
cassette <- cassette[cassette$Gene != 0]
cassette_p <- cassette[cassette$pvalue < 0.05,]
cassette_q <- cassette_p[cassette_p$qvalue < 0.05,]
cassette_sig <- cassette_q[abs(as.numeric(unlist(cassette_q[,11]))) > 0.1,]

JUM_table_cassette[i,j] <- print(nrow(cassette_sig))
write.csv(cassette_q, file=paste("JUM_cassette", i, "vs", j, ".csv", sep=""), row.names = FALSE, quote = FALSE)  # We save as cassette_q for now because it is easier to save in a spreadsheet. We will re-eliminate deltaPSI values below 0.1 in the "combination" chunk(s) later
   }
   if (i == j) next
 }
}
```

# Composite events:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

JUM_table_composite <- matrix(nrow = 46, ncol = 46, dimnames = list(singletypes,singletypes))
JUM_table_composite[is.na(JUM_table_composite)] <- 0

for (i in singletypes){
 for (j in singletypes){
   if (i != j){
composite <- fread(paste0("FINAL_JUM_OUTPUT_pvalue_1", i, "vs", j, "/AS_differential_JUM_output_composite_events_pvalue_1_final_simplified.txt"))

composite_split <- cSplit(composite, 8, sep = ";", stripWhite = TRUE, type.convert = FALSE)
composite_split_colcount <- as.numeric(ncol(composite_split))
composite_split[!is.finite(as.numeric(unlist(composite_split[,8]))),] <- 0
composite_split[is.nan(as.character(unlist(composite_split[,8]))),] <- 0
composite_split <- composite_split[composite_split$Gene != 0]
composite_p <- composite_split[composite_split$pvalue < 0.05]
composite_q <- composite_p[composite_p$qvalue < 0.05]
composite_q[composite_q == -Inf] <- 0
composite_q[composite_q == Inf] <- 0

composite_maxPSI <- composite_q[,8:ncol(composite_q)]
composite_maxPSI[is.na(composite_maxPSI[,1:ncol(composite_maxPSI)])] <- 0
composite_colcount <- as.numeric(ncol(composite_maxPSI))
composite_rowcount <- as.numeric(nrow(composite_maxPSI))
composite_maxPSI[composite_maxPSI == -Inf] <- 0
composite_maxPSI[composite_maxPSI == Inf] <- 0
composite_maxPSI$max <- rowMaxs(abs(as.numeric(unlist(composite_maxPSI))), dim. = as.integer(c(composite_rowcount, composite_colcount)))
composite_sig <- composite_maxPSI[composite_maxPSI$max > 0.1,]

JUM_table_composite[i,j] <- print(nrow(composite_sig))
write.csv(composite_q, file=paste("JUM_composite", i, "vs", j, ".csv", sep=""), row.names = FALSE, quote = FALSE)  # We save as composite_q for now because it is easier to save in a spreadsheet. We will re-eliminate deltaPSI values below 0.1 in the "combination" chunk(s) later
   }
   if (i == j) next
 }
}
```

# Intron retention:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

JUM_table_intron <- matrix(nrow = 46, ncol = 46, dimnames = list(singletypes,singletypes))
JUM_table_intron[is.na(JUM_table_intron)] <- 0

for (i in singletypes){
 for (j in singletypes){
   if (i != j){
intron <- fread(paste0("FINAL_JUM_OUTPUT_pvalue_1", i, "vs", j, "/AS_differential_JUM_output_intron_retention_pvalue_1_final_simplified.txt"))

intron[intron == "-INF"] <- 0
intron[intron == "INF"] <- 0
intron[is.nan(as.character(unlist(intron[,9]))),] <- 0
intron <- intron[intron$Gene != 0]
intron_p <- intron[as.numeric(intron$pvalue) < 0.05,]
intron_q <- intron_p[as.numeric(intron_p$qvalue) < 0.05,]
intron_sig <- intron_q[abs(as.numeric(unlist(intron_q[,9]))) > 0.1,]

JUM_table_intron[i,j] <- print(nrow(intron_sig))
write.csv(intron_q, file=paste("JUM_intron", i, "vs", j, ".csv", sep=""), row.names = FALSE, quote = FALSE) # We save as intron_q for now because it is easier to save in a spreadsheet. We will re-eliminate deltaPSI values below 0.1 in the "combination" chunk(s) later
   }
   if (i == j) next
 }
}
```

# MXE events:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

JUM_table_MXE <- matrix(nrow = 46, ncol = 46, dimnames = list(singletypes,singletypes))
JUM_table_MXE[is.na(JUM_table_MXE)] <- 0

for (i in singletypes){
 for (j in singletypes){
   if (i != j){
MXE <- fread(paste0("FINAL_JUM_OUTPUT_pvalue_1", i, "vs", j, "/AS_differential_JUM_output_MXE_events_pvalue_1_final_simplified.txt"), fill = TRUE)

MXE_split <- cSplit(MXE, 10, sep = ";", stripWhite = TRUE, type.convert = FALSE)
MXE_split[!is.finite(as.numeric(unlist(MXE_split[,10]))),] <- 0
MXE_split[is.nan(as.character(unlist(MXE_split[,10]))),] <- 0
MXE_split <- MXE_split[MXE_split$Gene != 0]
MXE_p <- MXE_split[MXE_split$pvalue < 0.05]
MXE_q <- MXE_p[MXE_p$qvalue < 0.05]

MXE_maxPSI <- MXE_q[,10:ncol(MXE_q)]
MXE_maxPSI[is.na(MXE_maxPSI[,1:ncol(MXE_maxPSI)])] <- 0
MXE_colcount <- as.numeric(ncol(MXE_maxPSI))
MXE_rowcount <- as.numeric(nrow(MXE_maxPSI))
MXE_maxPSI[MXE_maxPSI == -Inf] <- 0
MXE_maxPSI[MXE_maxPSI == Inf] <- 0
MXE_maxPSI$max <- rowMaxs(abs(as.numeric(unlist(MXE_maxPSI))), dim. = as.integer(c(MXE_rowcount, MXE_colcount)))
MXE_sig <- MXE_maxPSI[MXE_maxPSI$max > 0.1,]

JUM_table_MXE[i,j] <- print(nrow(MXE_sig))
write.csv(MXE_q, file=paste("JUM_MXE", i, "vs", j, ".csv", sep=""), row.names = FALSE, quote = FALSE) # We save as MXE_q for now because it is easier to save in a spreadsheet. We will re-eliminate deltaPSI values below 0.1 in the "combination" chunk(s) later
   }
   if (i == j) next
 }
}
```

# Everything that is not in alphabetical order will be marked as "NA":
```{r}
JUM_table_A3S_as_NA <- JUM_table_A3S
# Mark values as "NA" based on alphabetical order
for (i in 1:nrow(JUM_table_A3S)) {
  for (j in 1:ncol(JUM_table_A3S)) {
    if (rownames(JUM_table_A3S)[i] > colnames(JUM_table_A3S)[j]) {
      JUM_table_A3S_as_NA[i, j] <- NA
    }
  }
}

JUM_table_A5S_as_NA <- JUM_table_A5S
# Mark values as "NA" based on alphabetical order
for (i in 1:nrow(JUM_table_A5S)) {
  for (j in 1:ncol(JUM_table_A5S)) {
    if (rownames(JUM_table_A5S)[i] > colnames(JUM_table_A5S)[j]) {
      JUM_table_A5S_as_NA[i, j] <- NA
    }
  }
}

JUM_table_cassette_as_NA <- JUM_table_cassette
# Mark values as "NA" based on alphabetical order
for (i in 1:nrow(JUM_table_cassette)) {
  for (j in 1:ncol(JUM_table_cassette)) {
    if (rownames(JUM_table_cassette)[i] > colnames(JUM_table_cassette)[j]) {
      JUM_table_cassette_as_NA[i, j] <- NA
    }
  }
}

JUM_table_composite_as_NA <- JUM_table_composite
# Mark values as "NA" based on alphabetical order
for (i in 1:nrow(JUM_table_composite)) {
  for (j in 1:ncol(JUM_table_composite)) {
    if (rownames(JUM_table_composite)[i] > colnames(JUM_table_composite)[j]) {
      JUM_table_composite_as_NA[i, j] <- NA
    }
  }
}

JUM_table_intron_as_NA <- JUM_table_intron
# Mark values as "NA" based on alphabetical order
for (i in 1:nrow(JUM_table_intron)) {
  for (j in 1:ncol(JUM_table_intron)) {
    if (rownames(JUM_table_intron)[i] > colnames(JUM_table_intron)[j]) {
      JUM_table_intron_as_NA[i, j] <- NA
    }
  }
}

JUM_table_MXE_as_NA <- JUM_table_MXE
# Mark values as "NA" based on alphabetical order
for (i in 1:nrow(JUM_table_MXE)) {
  for (j in 1:ncol(JUM_table_MXE)) {
    if (rownames(JUM_table_MXE)[i] > colnames(JUM_table_MXE)[j]) {
      JUM_table_MXE_as_NA[i, j] <- NA
    }
  }
}
```


# Heatmaps:

### Note: when running many conditions, JUM may output incorrect data for cell types not in alphabetical order (i.e. ASG vs AVE will output correct data, but AVE vs ASG would not). You may need to hide boxes in heatmaps in order to plot the correct 50% (upper-right half in the heatmap) of data.
```{r}
colors <- colorRampPalette(brewer.pal(6, "Blues"))(250)
pheatmap(JUM_table_A3S_as_NA,
         #na_col = "#FFFFFF",   # Make NA values white (if desired)
         border_color = "black",
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         number_format = "%.0f",
         number_color = "black",
         fontsize = 7.5,
         fontsize_number = 6.5,
         angle_col = 90,
         col = colors,
         main = "Number of A3S Events between Cell Types")

colors <- colorRampPalette(brewer.pal(6, "Reds"))(250)
pheatmap(JUM_table_A5S_as_NA,
         border_color = "black",
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         number_format = "%.0f",
         number_color = "black",
         fontsize = 7.5,
         fontsize_number = 6.5,
         angle_col = 90,
         col = colors,
         main = "Number of A5S Events between Cell Types")

colors <- colorRampPalette(brewer.pal(6, "Greens"))(250)
pheatmap(JUM_table_cassette_as_NA,
         border_color = "black",
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         number_format = "%.0f",
         number_color = "black",
         fontsize = 7.5,
         fontsize_number = 6.5,
         angle_col = 90,
         col = colors,
         main = "Number of Cassette Exon Events between Cell Types")

colors <- colorRampPalette(brewer.pal(6, "BuPu"))(250)
pheatmap(JUM_table_composite_as_NA,
         border_color = "black",
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         number_format = "%.0f",
         number_color = "black",
         fontsize = 7.5,
         fontsize_number = 6.5,
         angle_col = 90,
         col = colors,
         main = "Number of Composite Events between Cell Types")

colors <- colorRampPalette(brewer.pal(6, "YlOrBr"))(250)
pheatmap(JUM_table_intron_as_NA,
         border_color = "black",
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         number_format = "%.0f",
         number_color = "black",
         fontsize = 7.5,
         fontsize_number = 6.5,
         angle_col = 90,
         col = colors,
         main = "Number of Intron Retention Events between Cell Types")

colors <- colorRampPalette(brewer.pal(6, "PuRd"))(250)
pheatmap(JUM_table_MXE_as_NA,
         border_color = "black",
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         number_format = "%.0f",
         number_color = "black",
         fontsize = 7.5,
         fontsize_number = 6.5,
         angle_col = 90,
         col = colors,
         main = "Number of MXE Splicing Events between Cell Types")
```
# Combine every excel file for each cell type:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

  A3S_files <- list.files(pattern = "JUM_A3S*")
  A5S_files <- list.files(pattern = "JUM_A5S*")
  cassette_files <- list.files(pattern = "JUM_cassette*")
  composite_files <- list.files(pattern = "JUM_composite*")
  intron_files <- list.files(pattern = "JUM_intron*")   
  MXE_files <- list.files(pattern = "JUM_MXE*")
  
  
# Create empty lists to store the filtered files
filtered_files_A3S <- vector("list", length(singletypes))
names(filtered_files_A3S) <- singletypes

filtered_files_A5S <- vector("list", length(singletypes))
names(filtered_files_A5S) <- singletypes

filtered_files_cassette <- vector("list", length(singletypes))
names(filtered_files_cassette) <- singletypes

filtered_files_composite <- vector("list", length(singletypes))
names(filtered_files_composite) <- singletypes

filtered_files_intron <- vector("list", length(singletypes))
names(filtered_files_intron) <- singletypes

filtered_files_MXE <- vector("list", length(singletypes))
names(filtered_files_MXE) <- singletypes

# Loop over the patterns
for (i in seq_along(singletypes)) {
  # Filter elements containing the current pattern
  filtered_files_A3S[[i]] <- A3S_files[grepl(singletypes[i], A3S_files)]
  filtered_files_A5S[[i]] <- A5S_files[grepl(singletypes[i], A5S_files)]
  filtered_files_cassette[[i]] <- cassette_files[grepl(singletypes[i], cassette_files)]
  filtered_files_composite[[i]] <- composite_files[grepl(singletypes[i], composite_files)]
  filtered_files_intron[[i]] <- intron_files[grepl(singletypes[i], intron_files)]
  filtered_files_MXE[[i]] <- MXE_files[grepl(singletypes[i], MXE_files)]
}

# Print the filtered files for each pattern
for (i in seq_along(singletypes)) {
  singletype <- singletypes[i]
  print(singletype)
  print(filtered_files_A3S[[i]])
  print(filtered_files_A5S[[i]])
  print(filtered_files_cassette[[i]])
  print(filtered_files_composite[[i]])
  print(filtered_files_intron[[i]])
  print(filtered_files_MXE[[i]])
}

# Remove files containing "PVD" from filtered_files$VD
filtered_files_A3S[["VD"]] <- filtered_files_A3S[["VD"]][!grepl("PVD", filtered_files_A3S[["VD"]])]
filtered_files_A3S[["VD"]] <- c(filtered_files_A3S$VD, "JUM_A3SPVDvsVD.csv", "JUM_A3SVDvsPVD.csv")

# Remove files containing "DVC" and "PVC" from filtered_files$VC
filtered_files_A3S[["VC"]] <- filtered_files_A3S[["VC"]][!grepl("PVC", filtered_files_A3S[["VC"]])]
filtered_files_A3S[["VC"]] <- filtered_files_A3S[["VC"]][!grepl("DVC", filtered_files_A3S[["VC"]])]
filtered_files_A3S[["VC"]] <- c(filtered_files_A3S$VC, "JUM_A3SDVCvsVC.csv", "JUM_A3SVCvsDVC.csv", "JUM_A3SPVCvsVC.csv", "JUM_A3SVCvsPVC.csv")

# Remove files containing "PVD" from filtered_files$VD
filtered_files_A5S[["VD"]] <- filtered_files_A5S[["VD"]][!grepl("PVD", filtered_files_A5S[["VD"]])]
filtered_files_A5S[["VD"]] <- c(filtered_files_A5S$VD, "JUM_A5SPVDvsVD.csv", "JUM_A5SVDvsPVD.csv")

# Remove files containing "DVC" and "PVC" from filtered_files$VC
filtered_files_A5S[["VC"]] <- filtered_files_A5S[["VC"]][!grepl("PVC", filtered_files_A5S[["VC"]])]
filtered_files_A5S[["VC"]] <- filtered_files_A5S[["VC"]][!grepl("DVC", filtered_files_A5S[["VC"]])]
filtered_files_A5S[["VC"]] <- c(filtered_files_A5S$VC, "JUM_A5SDVCvsVC.csv", "JUM_A5SVCvsDVC.csv", "JUM_A5SPVCvsVC.csv", "JUM_A5SVCvsPVC.csv")

# Remove files containing "PVD" from filtered_files$VD
filtered_files_cassette[["VD"]] <- filtered_files_cassette[["VD"]][!grepl("PVD", filtered_files_cassette[["VD"]])]
filtered_files_cassette[["VD"]] <- c(filtered_files_cassette$VD, "JUM_cassettePVDvsVD.csv", "JUM_cassetteVDvsPVD.csv")

# Remove files containing "DVC" and "PVC" from filtered_files$VC
filtered_files_cassette[["VC"]] <- filtered_files_cassette[["VC"]][!grepl("PVC", filtered_files_cassette[["VC"]])]
filtered_files_cassette[["VC"]] <- filtered_files_cassette[["VC"]][!grepl("DVC", filtered_files_cassette[["VC"]])]
filtered_files_cassette[["VC"]] <- c(filtered_files_cassette$VC, "JUM_cassetteDVCvsVC.csv", "JUM_cassetteVCvsDVC.csv", "JUM_cassettePVCvsVC.csv", "JUM_cassetteVCvsPVC.csv")

# Remove files containing "PVD" from filtered_files$VD
filtered_files_composite[["VD"]] <- filtered_files_composite[["VD"]][!grepl("PVD", filtered_files_composite[["VD"]])]
filtered_files_composite[["VD"]] <- c(filtered_files_composite$VD, "JUM_compositePVDvsVD.csv", "JUM_compositeVDvsPVD.csv")

# Remove files containing "DVC" and "PVC" from filtered_files$VC
filtered_files_composite[["VC"]] <- filtered_files_composite[["VC"]][!grepl("PVC", filtered_files_composite[["VC"]])]
filtered_files_composite[["VC"]] <- filtered_files_composite[["VC"]][!grepl("DVC", filtered_files_composite[["VC"]])]
filtered_files_composite[["VC"]] <- c(filtered_files_composite$VC, "JUM_compositeDVCvsVC.csv", "JUM_compositeVCvsDVC.csv", "JUM_compositePVCvsVC.csv", "JUM_compositeVCvsPVC.csv")

# Remove files containing "PVD" from filtered_files$VD
filtered_files_intron[["VD"]] <- filtered_files_intron[["VD"]][!grepl("PVD", filtered_files_intron[["VD"]])]
filtered_files_intron[["VD"]] <- c(filtered_files_intron$VD, "JUM_intronPVDvsVD.csv", "JUM_intronVDvsPVD.csv")

# Remove files containing "DVC" and "PVC" from filtered_files$VC
filtered_files_intron[["VC"]] <- filtered_files_intron[["VC"]][!grepl("PVC", filtered_files_intron[["VC"]])]
filtered_files_intron[["VC"]] <- filtered_files_intron[["VC"]][!grepl("DVC", filtered_files_intron[["VC"]])]
filtered_files_intron[["VC"]] <- c(filtered_files_intron$VC, "JUM_intronDVCvsVC.csv", "JUM_intronVCvsDVC.csv", "JUM_intronPVCvsVC.csv", "JUM_intronVCvsPVC.csv")

# Remove files containing "PVD" from filtered_files$VD
filtered_files_MXE[["VD"]] <- filtered_files_MXE[["VD"]][!grepl("PVD", filtered_files_MXE[["VD"]])]
filtered_files_MXE[["VD"]] <- c(filtered_files_MXE$VD, "JUM_MXEPVDvsVD.csv", "JUM_MXEVDvsPVD.csv")

# Remove files containing "DVC" and "PVC" from filtered_files$VC
filtered_files_MXE[["VC"]] <- filtered_files_MXE[["VC"]][!grepl("PVC", filtered_files_MXE[["VC"]])]
filtered_files_MXE[["VC"]] <- filtered_files_MXE[["VC"]][!grepl("DVC", filtered_files_MXE[["VC"]])]
filtered_files_MXE[["VC"]] <- c(filtered_files_MXE$VC, "JUM_MXEDVCvsVC.csv", "JUM_MXEVCvsDVC.csv", "JUM_MXEPVCvsVC.csv", "JUM_MXEVCvsPVC.csv")
```

# Write a function that returns the numeric position of an element in a vector:
```{r}
get_position <- function(vector, element) {
  position <- match(element, vector)
  if (is.na(position)) {
    return("Element not found in the vector.")
  } else {
    return(as.numeric(position))
  }
}
```

# Write a set of functions that return the elements in filtered_files in the correct alphabetical order:
```{r}
for (i in singletypes){
  if (i == "ADL"){    # ADL is first in our list of singletypes and does not require complex reordering
    ordered_files <- filtered_files_A3S[[i]][1:45]
    ordered_files <- ordered_files[order(ordered_files)]
  }
  if (i == "VC"){  # We have to arrange VC separately because we previously removed and re-added DVCs and PVCs in a unique fashion
    ordered_files <- filtered_files_A3S[[i]][c(1:42, 85, 87, 89)]
    ordered_files <- ordered_files[order(ordered_files)]
  }
  if (i == "VD"){  # We have to arrange VD separately because we previously removed and re-added PVDs in a unique fashion
    ordered_files <- filtered_files_A3S[[i]][c(1:44, 89)]
    ordered_files <- ordered_files[order(ordered_files)]
  }
  if (i != "ADL" & i != "VC" & i != "VD") {
  ordered_files <- filtered_files_A3S[[i]][order(filtered_files_A3S[[i]])]
  first_bad_file <- get_position(singletypes, i)
  last_bad_file <- get_position(singletypes, i) + get_position(singletypes, i) - 2
  ordered_files <- ordered_files[-c(first_bad_file:last_bad_file)]
  ordered_files <- ordered_files[-c(46:length(ordered_files))]
  }
  filtered_files_A3S[i] <- list(ordered_files)
}

for (i in singletypes){
  if (i == "ADL"){    # ADL is first in our list of singletypes and does not require complex reordering
    ordered_files <- filtered_files_A5S[[i]][1:45]
    ordered_files <- ordered_files[order(ordered_files)]
  }
  if (i == "VC"){  # We have to arrange VC separately because we previously removed and re-added DVCs and PVCs in a unique fashion
    ordered_files <- filtered_files_A5S[[i]][c(1:42, 85, 87, 89)]
    ordered_files <- ordered_files[order(ordered_files)]
  }
  if (i == "VD"){  # We have to arrange VD separately because we previously removed and re-added PVDs in a unique fashion
    ordered_files <- filtered_files_A5S[[i]][c(1:44, 89)]
    ordered_files <- ordered_files[order(ordered_files)]
  }
  if (i != "ADL" & i != "VC" & i != "VD") {
  ordered_files <- filtered_files_A5S[[i]][order(filtered_files_A5S[[i]])]
  first_bad_file <- get_position(singletypes, i)
  last_bad_file <- get_position(singletypes, i) + get_position(singletypes, i) - 2
  ordered_files <- ordered_files[-c(first_bad_file:last_bad_file)]
  ordered_files <- ordered_files[-c(46:length(ordered_files))]
  }
  filtered_files_A5S[i] <- list(ordered_files)
}

for (i in singletypes){
  if (i == "ADL"){    # ADL is first in our list of singletypes and does not require complex reordering
    ordered_files <- filtered_files_cassette[[i]][1:45]
    ordered_files <- ordered_files[order(ordered_files)]
  }
  if (i == "VC"){  # We have to arrange VC separately because we previously removed and re-added DVCs and PVCs in a unique fashion
    ordered_files <- filtered_files_cassette[[i]][c(1:42, 85, 87, 89)]
    ordered_files <- ordered_files[order(ordered_files)]
  }
  if (i == "VD"){  # We have to arrange VD separately because we previously removed and re-added PVDs in a unique fashion
    ordered_files <- filtered_files_cassette[[i]][c(1:44, 89)]
    ordered_files <- ordered_files[order(ordered_files)]
  }
  if (i != "ADL" & i != "VC" & i != "VD") {
  ordered_files <- filtered_files_cassette[[i]][order(filtered_files_cassette[[i]])]
  first_bad_file <- get_position(singletypes, i)
  last_bad_file <- get_position(singletypes, i) + get_position(singletypes, i) - 2
  ordered_files <- ordered_files[-c(first_bad_file:last_bad_file)]
  ordered_files <- ordered_files[-c(46:length(ordered_files))]
  }
  filtered_files_cassette[i] <- list(ordered_files)
}

for (i in singletypes){
  if (i == "ADL"){    # ADL is first in our list of singletypes and does not require complex reordering
    ordered_files <- filtered_files_composite[[i]][1:45]
    ordered_files <- ordered_files[order(ordered_files)]
  }
  if (i == "VC"){  # We have to arrange VC separately because we previously removed and re-added DVCs and PVCs in a unique fashion
    ordered_files <- filtered_files_composite[[i]][c(1:42, 85, 87, 89)]
    ordered_files <- ordered_files[order(ordered_files)]
  }
  if (i == "VD"){  # We have to arrange VD separately because we previously removed and re-added PVDs in a unique fashion
    ordered_files <- filtered_files_composite[[i]][c(1:44, 89)]
    ordered_files <- ordered_files[order(ordered_files)]
  }
  if (i != "ADL" & i != "VC" & i != "VD") {
  ordered_files <- filtered_files_composite[[i]][order(filtered_files_composite[[i]])]
  first_bad_file <- get_position(singletypes, i)
  last_bad_file <- get_position(singletypes, i) + get_position(singletypes, i) - 2
  ordered_files <- ordered_files[-c(first_bad_file:last_bad_file)]
  ordered_files <- ordered_files[-c(46:length(ordered_files))]
  }
  filtered_files_composite[i] <- list(ordered_files)
}

for (i in singletypes){
  if (i == "ADL"){    # ADL is first in our list of singletypes and does not require complex reordering
    ordered_files <- filtered_files_intron[[i]][1:45]
    ordered_files <- ordered_files[order(ordered_files)]
  }
  if (i == "VC"){  # We have to arrange VC separately because we previously removed and re-added DVCs and PVCs in a unique fashion
    ordered_files <- filtered_files_intron[[i]][c(1:42, 85, 87, 89)]
    ordered_files <- ordered_files[order(ordered_files)]
  }
  if (i == "VD"){  # We have to arrange VD separately because we previously removed and re-added PVDs in a unique fashion
    ordered_files <- filtered_files_intron[[i]][c(1:44, 89)]
    ordered_files <- ordered_files[order(ordered_files)]
  }
  if (i != "ADL" & i != "VC" & i != "VD") {
  ordered_files <- filtered_files_intron[[i]][order(filtered_files_intron[[i]])]
  first_bad_file <- get_position(singletypes, i)
  last_bad_file <- get_position(singletypes, i) + get_position(singletypes, i) - 2
  ordered_files <- ordered_files[-c(first_bad_file:last_bad_file)]
  ordered_files <- ordered_files[-c(46:length(ordered_files))]
  }
  filtered_files_intron[i] <- list(ordered_files)
}

for (i in singletypes){
  if (i == "ADL"){    # ADL is first in our list of singletypes and does not require complex reordering
    ordered_files <- filtered_files_MXE[[i]][1:45]
    ordered_files <- ordered_files[order(ordered_files)]
  }
  if (i == "VC"){  # We have to arrange VC separately because we previously removed and re-added DVCs and PVCs in a unique fashion
    ordered_files <- filtered_files_MXE[[i]][c(1:42, 85, 87, 89)]
    ordered_files <- ordered_files[order(ordered_files)]
  }
  if (i == "VD"){  # We have to arrange VD separately because we previously removed and re-added PVDs in a unique fashion
    ordered_files <- filtered_files_MXE[[i]][c(1:44, 89)]
    ordered_files <- ordered_files[order(ordered_files)]
  }
  if (i != "ADL" & i != "VC" & i != "VD") {
  ordered_files <- filtered_files_MXE[[i]][order(filtered_files_MXE[[i]])]
  first_bad_file <- get_position(singletypes, i)
  last_bad_file <- get_position(singletypes, i) + get_position(singletypes, i) - 2
  ordered_files <- ordered_files[-c(first_bad_file:last_bad_file)]
  ordered_files <- ordered_files[-c(46:length(ordered_files))]
  }
  filtered_files_MXE[i] <- list(ordered_files)
}
```

# Make a function to reverse the signs of specified columns
```{r}
# Here, num_columns represents the number of columns that WILL NOT BE reversed, starting at the beginning of the data frame:
## For example, to reverse all columns EXCEPT the first two columns, num_columns should = 2
reverse_columns <- function(data, num_columns) {
  if (nrow(data) != 0) {
    reversed_columns <- colnames(data)[-c(1:num_columns)]
    data[, (reversed_columns) := lapply(.SD, function(x) if (is.numeric(x)) -x else x), .SDcols = reversed_columns]
  }
  data
}
```

# Now we need to retain the files that are in the correct alphabetical order AND reverse the deltaPSIs for the other half of the files that are not in alphabetical order:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

# A3S:
for (i in singletypes){
  if (i != "ADL"){   # Since ADL is first in singletypes, we won't have to reverse any files for ADL
    if (i != "VD"){
  n <- get_position(singletypes, i)
  files_to_change <- filtered_files_A3S[[i]][0-n:length(filtered_files_A3S[[i]])]
    }
      if (i == "VD"){
      files_to_change <- filtered_files_A3S[[i]]
      }
  for (j in 1:length(files_to_change)){
    file <- fread(files_to_change[j])
    if (nrow(file) != 0){
    # Organize file name:
     string <- files_to_change[j]
     result <- strsplit(string, "\\.")
     before_dot <- result[[1]][1]
     after_dot <- result[[1]][2]

    file_rev <- reverse_columns(file, 8)
    write.csv(file_rev, paste0(before_dot, "_rev.", after_dot), row.names = F)
    files_to_change[j] <- basename(paste0(before_dot, "_rev.", after_dot))
    }
    if (nrow(file) == 0) next
    }
  filtered_files_A3S[[i]] <- c(files_to_change, filtered_files_A3S[[i]][n:length(filtered_files_A3S[[i]])])
  if (i == "ADL") next
    }
}

# A5S:
for (i in singletypes){
  if (i != "ADL"){   # Since ADL is first in singletypes, we won't have to reverse any files for ADL
    if (i != "VD"){
  n <- get_position(singletypes, i)
  files_to_change <- filtered_files_A5S[[i]][0-n:length(filtered_files_A5S[[i]])]
    }
      if (i == "VD"){
      files_to_change <- filtered_files_A5S[[i]]
      }
  for (j in 1:length(files_to_change)){
    file <- fread(files_to_change[j])
    if (nrow(file) != 0){
    # Organize file name:
     string <- files_to_change[j]
     result <- strsplit(string, "\\.")
     before_dot <- result[[1]][1]
     after_dot <- result[[1]][2]

    file_rev <- reverse_columns(file, 8)
    write.csv(file_rev, paste0(before_dot, "_rev.", after_dot), row.names = F, quote = F)
    files_to_change[j] <- basename(paste0(before_dot, "_rev.", after_dot))
    }
    if (nrow(file) == 0) next
    }
  filtered_files_A5S[[i]] <- c(files_to_change, filtered_files_A5S[[i]][n:length(filtered_files_A5S[[i]])])
  if (i == "ADL") next
    }
}

# cassette exons:
for (i in singletypes){
  if (i != "ADL"){   # Since ADL is first in singletypes, we won't have to reverse any files for ADL
    if (i != "VD"){
  n <- get_position(singletypes, i)
  files_to_change <- filtered_files_cassette[[i]][0-n:length(filtered_files_cassette[[i]])]
    }
      if (i == "VD"){
      files_to_change <- filtered_files_cassette[[i]]
      }
  for (j in 1:length(files_to_change)){
    file <- fread(files_to_change[j])
    if (nrow(file) != 0){
    # Organize file name:
     string <- files_to_change[j]
     result <- strsplit(string, "\\.")
     before_dot <- result[[1]][1]
     after_dot <- result[[1]][2]

    file_rev <- reverse_columns(file, 10)
    write.csv(file_rev, paste0(before_dot, "_rev.", after_dot), row.names = F)
    files_to_change[j] <- basename(paste0(before_dot, "_rev.", after_dot))
    }
    if (nrow(file) == 0) next
    }
  filtered_files_cassette[[i]] <- c(files_to_change, filtered_files_cassette[[i]][n:length(filtered_files_cassette[[i]])])
  if (i == "ADL") next
    }
}

# composite:
for (i in singletypes){
  if (i != "ADL"){   # Since ADL is first in singletypes, we won't have to reverse any files for ADL
    if (i != "VD"){
  n <- get_position(singletypes, i)
  files_to_change <- filtered_files_composite[[i]][0-n:length(filtered_files_composite[[i]])]
    }
      if (i == "VD"){
      files_to_change <- filtered_files_composite[[i]]
      }
  for (j in 1:length(files_to_change)){
    file <- fread(files_to_change[j])
    if (nrow(file) != 0){
    # Organize file name:
     string <- files_to_change[j]
     result <- strsplit(string, "\\.")
     before_dot <- result[[1]][1]
     after_dot <- result[[1]][2]

    file_rev <- reverse_columns(file, 7)
    write.csv(file_rev, paste0(before_dot, "_rev.", after_dot), row.names = F)
    files_to_change[j] <- basename(paste0(before_dot, "_rev.", after_dot))
    }
    if (nrow(file) == 0) next
    }
  filtered_files_composite[[i]] <- c(files_to_change, filtered_files_composite[[i]][n:length(filtered_files_composite[[i]])])
  if (i == "ADL") next
    }
}

# intron retention:
for (i in singletypes){
  if (i != "ADL"){   # Since ADL is first in singletypes, we won't have to reverse any files for ADL
    if (i != "VD"){
  n <- get_position(singletypes, i)
  files_to_change <- filtered_files_intron[[i]][0-n:length(filtered_files_intron[[i]])]
    }
      if (i == "VD"){
      files_to_change <- filtered_files_intron[[i]]
      }
  for (j in 1:length(files_to_change)){
    file <- fread(files_to_change[j])
    if (nrow(file) != 0){
    # Organize file name:
     string <- files_to_change[j]
     result <- strsplit(string, "\\.")
     before_dot <- result[[1]][1]
     after_dot <- result[[1]][2]

    file_rev <- reverse_columns(file, 8)
    write.csv(file_rev, paste0(before_dot, "_rev.", after_dot), row.names = F)
    files_to_change[j] <- basename(paste0(before_dot, "_rev.", after_dot))
    }
    if (nrow(file) == 0) next
    }
  filtered_files_intron[[i]] <- c(files_to_change, filtered_files_intron[[i]][n:length(filtered_files_intron[[i]])])
  if (i == "ADL") next
    }
}

# MXE:
for (i in singletypes){
  if (i != "ADL"){   # Since ADL is first in singletypes, we won't have to reverse any files for ADL
    if (i != "VD"){
  n <- get_position(singletypes, i)
  files_to_change <- filtered_files_MXE[[i]][0-n:length(filtered_files_MXE[[i]])]
    }
      if (i == "VD"){
      files_to_change <- filtered_files_MXE[[i]]
      }
  for (j in 1:length(files_to_change)){
    file <- fread(files_to_change[j])
    if (nrow(file) != 0){
    # Organize file name:
     string <- files_to_change[j]
     result <- strsplit(string, "\\.")
     before_dot <- result[[1]][1]
     after_dot <- result[[1]][2]

    file_rev <- reverse_columns(file, 9)
    write.csv(file_rev, paste0(before_dot, "_rev.", after_dot), row.names = F)
    files_to_change[j] <- basename(paste0(before_dot, "_rev.", after_dot))
    }
    if (nrow(file) == 0) next
    }
  filtered_files_MXE[[i]] <- c(files_to_change, filtered_files_MXE[[i]][n:length(filtered_files_MXE[[i]])])
  if (i == "ADL") next
    }
}

# Lastly, get rid of the final element in all the "VD" lists:
filtered_files_A3S$VD <- filtered_files_A3S$VD[-c(46)]
filtered_files_A5S$VD <- filtered_files_A5S$VD[-c(46)]
filtered_files_cassette$VD <- filtered_files_cassette$VD[-c(46)]
filtered_files_composite$VD <- filtered_files_composite$VD[-c(46)]
filtered_files_intron$VD <- filtered_files_intron$VD[-c(46)]
filtered_files_MXE$VD <- filtered_files_MXE$VD[-c(46)]
```


# Note: If there are 0 AS events in a comparison, dplyr will assume that the rows in the .csv file are logical (we want them to be character/numeric). This will prevent us from merging rows in different .csv files together in the upcoming chunks. Run this chunk to eliminate .csv files from our file list if there are 0 rows:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

# Loop over cell types and comparisons
for (cell_type in singletypes) {
   for (comparison in comparisons) {
    # Create the variable name
    list_name <- paste0("filtered_files_", comparison, "$", cell_type)
    
    # Get the list of file paths
    file_list <- eval(parse(text = list_name))
    
    # Filter out files with zero rows
    file_list_with_rows <- file_list[sapply(file_list, function(file) {
      nrows <- tryCatch(nrow(read.csv(file)), error = function(e) 0)
      nrows > 0
    })]
    
    # Assign the updated file list back to the list_name
    assign(list_name, file_list_with_rows)
  }
}
```

# deltaPSI values for unique AS events in individual "cell type vs" files:

# A3S events:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

# Loop over cell types
for (cell_type in singletypes) {
  # Create file names based on cell type and comparison
  file_name <- paste0("LargeA3Sdf_", cell_type)
  count_file_name <- paste0("A3S_event_ID_", cell_type, "_count")
  
  # Read the CSV file
  df <- get(paste0("filtered_files_A3S$", cell_type)) %>% map_dfr(read.csv)
  
  # Rename columns
  colnames_df <- colnames(df)
  deltaPSI_names <- rep(paste0("deltaPSI_", 1:(ncol(df) - 8)))
  df <- df %>% setnames(old = colnames_df[9:ncol(df)], new = deltaPSI_names)
  
  # Sort by AS_event_ID
  df <- df[order(df$AS_event_ID), ]
  
  # Calculate count
  count_df <- as.data.frame(table(df$AS_event_ID))
  count_df <- count_df[rev(order(count_df$Freq)), ]
  names(count_df) <- c("A3S_event_ID", "Freq")
  
  df[is.na(df)] <- 0
  
  df$MaxPSIs <- rowMaxs(as.numeric(unlist(df[, 9:ncol(df)])), dim. = dim(df[, 9:ncol(df)]))
  df$MinPSIs <- rowMins(as.numeric(unlist(df[, 9:ncol(df)])), dim. = dim(df[, 9:ncol(df)]))
  
  df <- df[df$MaxPSIs > 0.1 | abs(df$MinPSIs) > 0.1, ]
  
  df$first_PSI <- apply(df[, 9:ncol(df)], 1, function(x) x[x != 0][1])
  
  df <- df[which(df$MaxPSIs == (-1) * df$MinPSIs), ]
  
  unique_A3S_events <- unique(df$AS_event_ID)
  unique_A3S_genes <- unique(df$Gene)
  
  deltaPSI_sums <- df[,c("AS_event_ID", "Gene", "first_PSI", "MaxPSIs", "MinPSIs")]
  names(deltaPSI_sums) <- c("A3S_event_ID", "Gene", "first deltaPSI", "max deltaPSI", "min deltaPSI")
  deltaPSI_sums$`max delta PSI plus min delta PSI` <- deltaPSI_sums[, 4] + deltaPSI_sums[, 5]
  
  first_PSI_sums <- data.frame()
  for (i in unique_A3S_events) {
    first_PSI_sums[i, 1] <- sum(deltaPSI_sums[deltaPSI_sums$A3S_event_ID == i, 3])
  }
  
  for (i in rownames(first_PSI_sums)) {
    first_PSI_sums[i, 2] <- first(deltaPSI_sums$Gene[deltaPSI_sums$A3S_event_ID == i])
  }
  names(first_PSI_sums) <- c("sum of all first deltaPSIs for AS_event_ID", "Gene")
  
  count_df <- as.data.frame(table(df$AS_event_ID))
  rownames(count_df) <- count_df[, 1]
  for (i in rownames(count_df)) {
    count_df[i, 3] <- first_PSI_sums$Gene[count_df$Var1 == i]
    count_df[i, 4] <- first_PSI_sums$`sum of all first deltaPSIs for AS_event_ID`[count_df$Var1 == i]
  }
  names(count_df) <- c("A3S_event_ID", "Freq", "Gene", "sum of all first deltaPSIs for AS_event_ID")
  
  count_df <- count_df[rev(order(count_df$Freq)), ]
  
  # Save the dataframes to CSV files
  write.csv(df, file = paste0(file_name, ".csv"), row.names = FALSE, quote = FALSE)
  write.csv(count_df, file = paste0(count_file_name, ".csv"), row.names = FALSE, quote = FALSE)
  write.csv(first_PSI_sums, file = paste0("A3S_", cell_type, "_first_PSI_sums.csv"))
  write.csv(deltaPSI_sums, file = paste0("A3S_", cell_type, "_deltaPSI_sums.csv"), row.names = FALSE, quote = FALSE)  
}
```

# A5S events:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

# Loop over cell types
for (cell_type in singletypes) {
  # Create file names based on cell type and comparison
  file_name <- paste0("LargeA5Sdf_", cell_type)
  count_file_name <- paste0("A5S_event_ID_", cell_type, "_count")
  
  # Read the CSV file
  df <- get(paste0("filtered_files_A5S$", cell_type)) %>% map_dfr(read.csv)
  
  # Rename columns
  colnames_df <- colnames(df)
  deltaPSI_names <- rep(paste0("deltaPSI_", 1:(ncol(df) - 8)))
  df <- df %>% setnames(old = colnames_df[9:ncol(df)], new = deltaPSI_names)
  
  # Sort by AS_event_ID
  df <- df[order(df$AS_event_ID), ]
  
  # Calculate count
  count_df <- as.data.frame(table(df$AS_event_ID))
  count_df <- count_df[rev(order(count_df$Freq)), ]
  names(count_df) <- c("A5S_event_ID", "Freq")
  
  df[is.na(df)] <- 0
  
  df$MaxPSIs <- rowMaxs(as.numeric(unlist(df[, 9:ncol(df)])), dim. = dim(df[, 9:ncol(df)]))
  df$MinPSIs <- rowMins(as.numeric(unlist(df[, 9:ncol(df)])), dim. = dim(df[, 9:ncol(df)]))
  
  df <- df[df$MaxPSIs > 0.1 | abs(df$MinPSIs) > 0.1, ]
  
  df$first_PSI <- apply(df[, 9:ncol(df)], 1, function(x) x[x != 0][1])
  
  df <- df[which(df$MaxPSIs == (-1) * df$MinPSIs), ]
  
  unique_A5S_events <- unique(df$AS_event_ID)
  unique_A5S_genes <- unique(df$Gene)
  
  deltaPSI_sums <- df[,c("AS_event_ID", "Gene", "first_PSI", "MaxPSIs", "MinPSIs")]
  names(deltaPSI_sums) <- c("A5S_event_ID", "Gene", "first deltaPSI", "max deltaPSI", "min deltaPSI")
  deltaPSI_sums$`max delta PSI plus min delta PSI` <- deltaPSI_sums[, 4] + deltaPSI_sums[, 5]
  
  first_PSI_sums <- data.frame()
  for (i in unique_A5S_events) {
    first_PSI_sums[i, 1] <- sum(deltaPSI_sums[deltaPSI_sums$A5S_event_ID == i, 3])
  }
  
  for (i in rownames(first_PSI_sums)) {
    first_PSI_sums[i, 2] <- first(deltaPSI_sums$Gene[deltaPSI_sums$A5S_event_ID == i])
  }
  names(first_PSI_sums) <- c("sum of all first deltaPSIs for AS_event_ID", "Gene")
  
  count_df <- as.data.frame(table(df$AS_event_ID))
  rownames(count_df) <- count_df[, 1]
  for (i in rownames(count_df)) {
    count_df[i, 3] <- first_PSI_sums$Gene[count_df$Var1 == i]
    count_df[i, 4] <- first_PSI_sums$`sum of all first deltaPSIs for AS_event_ID`[count_df$Var1 == i]
  }
  names(count_df) <- c("A5S_event_ID", "Freq", "Gene", "sum of all first deltaPSIs for AS_event_ID")
  
  count_df <- count_df[rev(order(count_df$Freq)), ]
  
  # Save the dataframes to CSV files
  write.csv(df, file = paste0(file_name, ".csv"), row.names = FALSE, quote = FALSE)
  write.csv(count_df, file = paste0(count_file_name, ".csv"), row.names = FALSE, quote = FALSE)
  write.csv(first_PSI_sums, file = paste0("A5S_", cell_type, "_first_PSI_sums.csv"))
  write.csv(deltaPSI_sums, file = paste0("A5S_", cell_type, "_deltaPSI_sums.csv"), row.names = FALSE, quote = FALSE)  
}
```

# Cassette exons:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

# Loop over cell types
for (cell_type in singletypes) {
  # Create file names based on cell type and comparison
  file_name <- paste0("Largecassettedf_", cell_type)
  count_file_name <- paste0("cassette_event_ID_", cell_type, "_count")
  
  # Read the CSV file
  df <- get(paste0("filtered_files_cassette$", cell_type)) %>% map_dfr(read.csv)
  
  # Rename columns
  colnames_df <- colnames(df)
  deltaPSI_names <- rep(paste0("deltaPSI_", 1:(ncol(df) - 10)))
  df <- df %>% setnames(old = colnames_df[11:ncol(df)], new = deltaPSI_names)
  
  # Sort by AS_event_ID
  df <- df[order(df$AS_event_ID), ]
  
  # Calculate count
  count_df <- as.data.frame(table(df$AS_event_ID))
  count_df <- count_df[rev(order(count_df$Freq)), ]
  names(count_df) <- c("cassette_event_ID", "Freq")
  
  df[is.na(df)] <- 0
     
  df$MaxPSIs <- rowMaxs(as.numeric(unlist(df[, 11:ncol(df)])), dim. = dim(df[, 11:ncol(df)]))
  df$MinPSIs <- rowMins(as.numeric(unlist(df[, 11:ncol(df)])), dim. = dim(df[, 11:ncol(df)]))
  
  df <- df[df$MaxPSIs > 0.1 | abs(df$MinPSIs) > 0.1, ]

  df$first_PSI <- apply(df[, 11:ncol(df)], 1, function(x) x[x != 0][1])
    
  unique_cassette_events <- unique(df$AS_event_ID)
  unique_cassette_events <- na.omit(unique_cassette_events)
  unique_cassette_genes <- unique(df$Gene)
  
  deltaPSI_sums <-  df[,c("AS_event_ID", "Gene", "first_PSI", "MaxPSIs", "MinPSIs")]
  names(deltaPSI_sums) <- c("cassette_event_ID", "Gene", "first deltaPSI", "max deltaPSI", "min deltaPSI")
  deltaPSI_sums$`max delta PSI plus min delta PSI` <- deltaPSI_sums[, 4] + deltaPSI_sums[, 5]
  
  first_PSI_sums <- data.frame()
  for (i in unique_cassette_events) {
    first_PSI_sums[i, 1] <- sum(deltaPSI_sums[deltaPSI_sums$cassette_event_ID == i, 3])
  }
  
  for (i in rownames(first_PSI_sums)) {
    first_PSI_sums[i, 2] <- first(deltaPSI_sums$Gene[deltaPSI_sums$cassette_event_ID == i])
  }
  names(first_PSI_sums) <- c("sum of all first deltaPSIs for AS_event_ID", "Gene")
  
  count_df <- as.data.frame(table(df$AS_event_ID))
  rownames(count_df) <- count_df[, 1]
  for (i in rownames(count_df)) {
    count_df[i, 3] <- first_PSI_sums$Gene[count_df$Var1 == i]
    count_df[i, 4] <- first_PSI_sums$`sum of all first deltaPSIs for AS_event_ID`[count_df$Var1 == i]
  }
  names(count_df) <- c("cassette_event_ID", "Freq", "Gene", "sum of all first deltaPSIs for AS_event_ID")
  
  count_df <- count_df[rev(order(count_df$Freq)), ]
  
  # Save the dataframes to CSV files
  write.csv(df, file = paste0(file_name, ".csv"), row.names = FALSE, quote = FALSE)
  write.csv(count_df, file = paste0(count_file_name, ".csv"), row.names = FALSE, quote = FALSE)
  write.csv(first_PSI_sums, file = paste0("cassette_", cell_type, "_first_PSI_sums.csv"))
  write.csv(deltaPSI_sums, file = paste0("cassette_", cell_type, "_deltaPSI_sums.csv"), row.names = FALSE, quote = FALSE)  
}
```

# Composite events:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

# Loop over cell types
for (cell_type in singletypes) {
  # Create file names based on cell type and comparison
  file_name <- paste0("Largecompositedf_", cell_type)
  count_file_name <- paste0("composite_event_ID_", cell_type, "_count")
  
  # Read the CSV file
  df <- get(paste0("filtered_files_composite$", cell_type)) %>% map_dfr(read.csv)
  
  # Rename columns
  colnames_df <- colnames(df)
  deltaPSI_names <- rep(paste0("deltaPSI_", 1:(ncol(df) - 7)))
  df <- df %>% setnames(old = colnames_df[8:ncol(df)], new = deltaPSI_names)
  
  # Sort by AS_event_ID
  df <- df[order(df$AS_event_ID), ]
  
  # Calculate count
  count_df <- as.data.frame(table(df$AS_event_ID))
  count_df <- count_df[rev(order(count_df$Freq)), ]
  names(count_df) <- c("composite_event_ID", "Freq")
  
  df[is.na(df)] <- 0
  
  df$MaxPSIs <- rowMaxs(as.numeric(unlist(df[, 8:ncol(df)])), dim. = dim(df[, 8:ncol(df)]))
  df$MinPSIs <- rowMins(as.numeric(unlist(df[, 8:ncol(df)])), dim. = dim(df[, 8:ncol(df)]))
  
  df <- df[df$MaxPSIs > 0.1 | abs(df$MinPSIs) > 0.1, ]
  
  df$first_PSI <- apply(df[, 8:ncol(df)], 1, function(x) x[x != 0][1])
  
  df <- df[which(df$MaxPSIs == (-1) * df$MinPSIs), ]
  
  unique_composite_events <- unique(df$AS_event_ID)
  unique_composite_genes <- unique(df$Gene)
  
  deltaPSI_sums <- df[,c("AS_event_ID", "Gene", "first_PSI", "MaxPSIs", "MinPSIs")]
  names(deltaPSI_sums) <- c("composite_event_ID", "Gene", "first deltaPSI", "max deltaPSI", "min deltaPSI")
  deltaPSI_sums$`max delta PSI plus min delta PSI` <- deltaPSI_sums[, 4] + deltaPSI_sums[, 5]
  
  first_PSI_sums <- data.frame()
  for (i in unique_composite_events) {
    first_PSI_sums[i, 1] <- sum(deltaPSI_sums[deltaPSI_sums$composite_event_ID == i, 3])
  }
  
  for (i in rownames(first_PSI_sums)) {
    first_PSI_sums[i, 2] <- first(deltaPSI_sums$Gene[deltaPSI_sums$composite_event_ID == i])
  }
  names(first_PSI_sums) <- c("sum of all first deltaPSIs for AS_event_ID", "Gene")
  
  count_df <- as.data.frame(table(df$AS_event_ID))
  rownames(count_df) <- count_df[, 1]
  for (i in rownames(count_df)) {
    count_df[i, 3] <- first_PSI_sums$Gene[count_df$Var1 == i]
    count_df[i, 4] <- first_PSI_sums$`sum of all first deltaPSIs for AS_event_ID`[count_df$Var1 == i]
  }
  names(count_df) <- c("composite_event_ID", "Freq", "Gene", "sum of all first deltaPSIs for AS_event_ID")
  
  count_df <- count_df[rev(order(count_df$Freq)), ]
  
  # Save the dataframes to CSV files
  write.csv(df, file = paste0(file_name, ".csv"), row.names = FALSE, quote = FALSE)
  write.csv(count_df, file = paste0(count_file_name, ".csv"), row.names = FALSE, quote = FALSE)
  write.csv(first_PSI_sums, file = paste0("composite_", cell_type, "_first_PSI_sums.csv"))
  write.csv(deltaPSI_sums, file = paste0("composite_", cell_type, "_deltaPSI_sums.csv"), row.names = FALSE, quote = FALSE)  
}
```

# Intron Retention events:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

# Loop over cell types
for (cell_type in singletypes) {
  # Create file names based on cell type and comparison
  file_name <- paste0("Largeintrondf_", cell_type)
  count_file_name <- paste0("intron_event_ID_", cell_type, "_count")
  
  # Read the CSV file
  df <- get(paste0("filtered_files_intron$", cell_type)) %>% map_dfr(read.csv)
  
  # Rename columns
  colnames_df <- colnames(df)
  deltaPSI_names <- rep(paste0("deltaPSI_", 1:(ncol(df) - 8)))
  df <- df %>% setnames(old = colnames_df[9:ncol(df)], new = deltaPSI_names)
  
  # Sort by AS_event_ID
  df <- df[order(df$AS_event_ID), ]
  
  # Calculate count
  count_df <- as.data.frame(table(df$AS_event_ID))
  count_df <- count_df[rev(order(count_df$Freq)), ]
  names(count_df) <- c("intron_event_ID", "Freq")
  
  df[is.na(df)] <- 0

  df$MaxPSIs <- rowMaxs(as.numeric(unlist(df[, 9:ncol(df)])), dim. = dim(df[, 9:ncol(df)]))
  df$MinPSIs <- rowMins(as.numeric(unlist(df[, 9:ncol(df)])), dim. = dim(df[, 9:ncol(df)]))
  
  df <- df[df$MaxPSIs > 0.1 | abs(df$MinPSIs) > 0.1, ]

  df$first_PSI <- apply(df[, 9:ncol(df)], 1, function(x) x[x != 0][1])
    
  unique_intron_events <- unique(df$AS_event_ID)
  unique_intron_events <- na.omit(unique_intron_events)
  unique_intron_genes <- unique(df$Gene)
  
  deltaPSI_sums <- df[,c("AS_event_ID", "Gene", "first_PSI", "MaxPSIs", "MinPSIs")]
  names(deltaPSI_sums) <- c("intron_event_ID", "Gene", "first deltaPSI", "max deltaPSI", "min deltaPSI")
  deltaPSI_sums$`max delta PSI plus min delta PSI` <- deltaPSI_sums[, 4] + deltaPSI_sums[, 5]
  
  first_PSI_sums <- data.frame()
  for (i in unique_intron_events) {
    first_PSI_sums[i, 1] <- sum(deltaPSI_sums[deltaPSI_sums$intron_event_ID == i, 3])
  }
  
  for (i in rownames(first_PSI_sums)) {
    first_PSI_sums[i, 2] <- first(deltaPSI_sums$Gene[deltaPSI_sums$intron_event_ID == i])
  }
  names(first_PSI_sums) <- c("sum of all first deltaPSIs for AS_event_ID", "Gene")
  
  count_df <- as.data.frame(table(df$AS_event_ID))
  rownames(count_df) <- count_df[, 1]
  for (i in rownames(count_df)) {
    count_df[i, 3] <- first_PSI_sums$Gene[count_df$Var1 == i]
    count_df[i, 4] <- first_PSI_sums$`sum of all first deltaPSIs for AS_event_ID`[count_df$Var1 == i]
  }
  names(count_df) <- c("intron_event_ID", "Freq", "Gene", "sum of all first deltaPSIs for AS_event_ID")
  
  count_df <- count_df[rev(order(count_df$Freq)), ]
  
  # Save the dataframes to CSV files
  write.csv(df, file = paste0(file_name, ".csv"), row.names = FALSE, quote = FALSE)
  write.csv(count_df, file = paste0(count_file_name, ".csv"), row.names = FALSE, quote = FALSE)
  write.csv(first_PSI_sums, file = paste0("intron_", cell_type, "_first_PSI_sums.csv"))
  write.csv(deltaPSI_sums, file = paste0("intron_", cell_type, "_deltaPSI_sums.csv"), row.names = FALSE, quote = FALSE)  
}
```

# Mutually exclusive exons:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

# Loop over cell types
for (cell_type in singletypes) {
  # Create file names based on cell type and comparison
  file_name <- paste0("LargeMXEdf_", cell_type)
  count_file_name <- paste0("MXE_event_ID_", cell_type, "_count")
  
  # Read the CSV file
  df <- get(paste0("filtered_files_MXE$", cell_type)) %>% map_dfr(read.csv)
  
  # Rename columns
  colnames_df <- colnames(df)
  deltaPSI_names <- rep(paste0("deltaPSI_", 1:(ncol(df) - 9)))
  df <- df %>% setnames(old = colnames_df[10:ncol(df)], new = deltaPSI_names)
  
  # Sort by AS_event_ID
  df <- df[order(df$AS_event_ID), ]
  
  # Calculate count
  count_df <- as.data.frame(table(df$AS_event_ID))
  count_df <- count_df[rev(order(count_df$Freq)), ]
  names(count_df) <- c("MXE_event_ID", "Freq")
  
  df[is.na(df)] <- 0
  
  df$MaxPSIs <- rowMaxs(as.numeric(unlist(df[, 10:ncol(df)])), dim. = dim(df[, 10:ncol(df)]))
  df$MinPSIs <- rowMins(as.numeric(unlist(df[, 10:ncol(df)])), dim. = dim(df[, 10:ncol(df)]))
  
  df <- df[df$MaxPSIs > 0.1 | abs(df$MinPSIs) > 0.1, ]
  
  df$first_PSI <- apply(df[, 10:ncol(df)], 1, function(x) x[x != 0][1])
  
  df <- df[which(df$MaxPSIs == (-1) * df$MinPSIs), ]
  
  unique_MXE_events <- unique(df$AS_event_ID)
  unique_MXE_genes <- unique(df$Gene)
  
  deltaPSI_sums <- df[,c("AS_event_ID", "Gene", "first_PSI", "MaxPSIs", "MinPSIs")]
  names(deltaPSI_sums) <- c("MXE_event_ID", "Gene", "first deltaPSI", "max deltaPSI", "min deltaPSI")
  deltaPSI_sums$`max delta PSI plus min delta PSI` <- deltaPSI_sums[, 4] + deltaPSI_sums[, 5]
  
  first_PSI_sums <- data.frame()
  for (i in unique_MXE_events) {
    first_PSI_sums[i, 1] <- sum(deltaPSI_sums[deltaPSI_sums$MXE_event_ID == i, 3])
  }
  
  for (i in rownames(first_PSI_sums)) {
    first_PSI_sums[i, 2] <- first(deltaPSI_sums$Gene[deltaPSI_sums$MXE_event_ID == i])
  }
  names(first_PSI_sums) <- c("sum of all first deltaPSIs for AS_event_ID", "Gene")
  
  count_df <- as.data.frame(table(df$AS_event_ID))
  rownames(count_df) <- count_df[, 1]
  for (i in rownames(count_df)) {
    count_df[i, 3] <- first_PSI_sums$Gene[count_df$Var1 == i]
    count_df[i, 4] <- first_PSI_sums$`sum of all first deltaPSIs for AS_event_ID`[count_df$Var1 == i]
  }
  names(count_df) <- c("MXE_event_ID", "Freq", "Gene", "sum of all first deltaPSIs for AS_event_ID")
  
  count_df <- count_df[rev(order(count_df$Freq)), ]
  
  # Save the dataframes to CSV files
  write.csv(df, file = paste0(file_name, ".csv"), row.names = FALSE, quote = FALSE)
  write.csv(count_df, file = paste0(count_file_name, ".csv"), row.names = FALSE, quote = FALSE)
  write.csv(first_PSI_sums, file = paste0("MXE_", cell_type, "_first_PSI_sums.csv"))
  write.csv(deltaPSI_sums, file = paste0("MXE_", cell_type, "_deltaPSI_sums.csv"), row.names = FALSE, quote = FALSE)  
}
```

## Make a heatmap for the deltaPSIsums value in each A.S. event type for a given cell type:
# First, create a list to store the data frames
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

deltaPSI_sums_list <- list()

# read in the .csvs generated in the previous loop:

for (i in c("A3S", "A5S", "cassette", "composite", "intron", "MXE")) {
  for (j in singletypes) {
    # Generate the file name
    file_name <- paste0(i, "_event_ID_", j, "_count.csv")
    
    # Read the file using fread() and store it in the data_frames list
    deltaPSI_sums_list[[file_name]] <- fread(file_name)
  }
}
```

# Now make heatmap dataframes for each A.S. event type:
```{r}
# A3S:
A3S_events_list <- character()
for (i in singletypes) {
  A3S_events <- deltaPSI_sums_list[[paste0("A3S_event_ID_", i, "_count.csv")]]$A3S_event_ID
  A3S_events_list <- c(A3S_events_list, A3S_events)
}

row_names <- unique(A3S_events_list)
A3S_sums_heatmap <- matrix(nrow = length(row_names), ncol = 46, dimnames = list(row_names, singletypes))
A3S_sums_heatmap[is.na(A3S_sums_heatmap)] <- 0

for (i in seq_along(singletypes)) {
  A3Ssums_for_graph <- deltaPSI_sums_list[[paste0("A3S_event_ID_", singletypes[i], "_count.csv")]] %>%
    as.data.frame()
  rownames(A3Ssums_for_graph) <- A3Ssums_for_graph[, 1]

  for (j in rownames(A3Ssums_for_graph)) {
    values <- as.numeric(A3Ssums_for_graph[j, 4])
    if (!is.na(values)) {
      consolidated_row_name <- sub("\\.\\d+$", "", j)
      row_idx <- which(rownames(A3S_sums_heatmap) == consolidated_row_name)
      A3S_sums_heatmap[row_idx, i] <- sum(values, na.rm = TRUE)
    }
  }
}

maxs <- rowMaxs(as.matrix(A3S_sums_heatmap))
A3S_sums_heatmap <- cbind(A3S_sums_heatmap, maxs)
A3S_sums_heatmap <- A3S_sums_heatmap[order(unlist(A3S_sums_heatmap[,47]), decreasing = TRUE),]
A3S_sums_heatmap <- A3S_sums_heatmap[,-47]

# A5S:
A5S_events_list <- character()
for (i in singletypes) {
  A5S_events <- deltaPSI_sums_list[[paste0("A5S_event_ID_", i, "_count.csv")]]$A5S_event_ID
  A5S_events_list <- c(A5S_events_list, A5S_events)
}

row_names <- unique(A5S_events_list)
A5S_sums_heatmap <- matrix(nrow = length(row_names), ncol = 46, dimnames = list(row_names, singletypes))
A5S_sums_heatmap[is.na(A5S_sums_heatmap)] <- 0

for (i in seq_along(singletypes)) {
  A5Ssums_for_graph <- deltaPSI_sums_list[[paste0("A5S_event_ID_", singletypes[i], "_count.csv")]] %>%
    as.data.frame()
  rownames(A5Ssums_for_graph) <- A5Ssums_for_graph[, 1]

  for (j in rownames(A5Ssums_for_graph)) {
    values <- as.numeric(A5Ssums_for_graph[j, 4])
    if (!is.na(values)) {
      consolidated_row_name <- sub("\\.\\d+$", "", j)
      row_idx <- which(rownames(A5S_sums_heatmap) == consolidated_row_name)
      A5S_sums_heatmap[row_idx, i] <- sum(values, na.rm = TRUE)
    }
  }
}

maxs <- rowMaxs(as.matrix(A5S_sums_heatmap))
A5S_sums_heatmap <- cbind(A5S_sums_heatmap, maxs)
A5S_sums_heatmap <- A5S_sums_heatmap[order(unlist(A5S_sums_heatmap[,47]), decreasing = TRUE),]
A5S_sums_heatmap <- A5S_sums_heatmap[,-47]

# cassette exons:
cassette_events_list <- character()
for (i in singletypes) {
  cassette_events <- deltaPSI_sums_list[[paste0("cassette_event_ID_", i, "_count.csv")]]$cassette_event_ID
  cassette_events_list <- c(cassette_events_list, cassette_events)
}

row_names <- unique(cassette_events_list)
cassette_sums_heatmap <- matrix(nrow = length(row_names), ncol = 46, dimnames = list(row_names, singletypes))
cassette_sums_heatmap[is.na(cassette_sums_heatmap)] <- 0

for (i in seq_along(singletypes)) {
  cassettesums_for_graph <- deltaPSI_sums_list[[paste0("cassette_event_ID_", singletypes[i], "_count.csv")]] %>%
    as.data.frame()
  rownames(cassettesums_for_graph) <- cassettesums_for_graph[, 1]

  for (j in rownames(cassettesums_for_graph)) {
    values <- as.numeric(cassettesums_for_graph[j, 4])
    if (!is.na(values)) {
      consolidated_row_name <- sub("\\.\\d+$", "", j)
      row_idx <- which(rownames(cassette_sums_heatmap) == consolidated_row_name)
      cassette_sums_heatmap[row_idx, i] <- sum(values, na.rm = TRUE)
    }
  }
}

maxs <- rowMaxs(as.matrix(cassette_sums_heatmap))
cassette_sums_heatmap <- cbind(cassette_sums_heatmap, maxs)
cassette_sums_heatmap <- cassette_sums_heatmap[order(unlist(cassette_sums_heatmap[,47]), decreasing = TRUE),]
cassette_sums_heatmap <- cassette_sums_heatmap[,-47]

# composite:
composite_events_list <- character()
for (i in singletypes) {
  composite_events <- deltaPSI_sums_list[[paste0("composite_event_ID_", i, "_count.csv")]]$composite_event_ID
  composite_events_list <- c(composite_events_list, composite_events)
}

row_names <- unique(composite_events_list)
composite_sums_heatmap <- matrix(nrow = length(row_names), ncol = 46, dimnames = list(row_names, singletypes))
composite_sums_heatmap[is.na(composite_sums_heatmap)] <- 0

for (i in seq_along(singletypes)) {
  compositesums_for_graph <- deltaPSI_sums_list[[paste0("composite_event_ID_", singletypes[i], "_count.csv")]] %>%
    as.data.frame()
  rownames(compositesums_for_graph) <- compositesums_for_graph[, 1]

  for (j in rownames(compositesums_for_graph)) {
    values <- as.numeric(compositesums_for_graph[j, 4])
    if (!is.na(values)) {
      consolidated_row_name <- sub("\\.\\d+$", "", j)
      row_idx <- which(rownames(composite_sums_heatmap) == consolidated_row_name)
      composite_sums_heatmap[row_idx, i] <- sum(values, na.rm = TRUE)
    }
  }
}

maxs <- rowMaxs(as.matrix(composite_sums_heatmap))
composite_sums_heatmap <- cbind(composite_sums_heatmap, maxs)
composite_sums_heatmap <- composite_sums_heatmap[order(unlist(composite_sums_heatmap[,47]), decreasing = TRUE),]
composite_sums_heatmap <- composite_sums_heatmap[,-47]

# intron retention:
intron_events_list <- character()
for (i in singletypes) {
  intron_events <- deltaPSI_sums_list[[paste0("intron_event_ID_", i, "_count.csv")]]$intron_event_ID
  intron_events_list <- c(intron_events_list, intron_events)
}

row_names <- unique(intron_events_list)
intron_sums_heatmap <- matrix(nrow = length(row_names), ncol = 46, dimnames = list(row_names, singletypes))
intron_sums_heatmap[is.na(intron_sums_heatmap)] <- 0

for (i in seq_along(singletypes)) {
  intronsums_for_graph <- deltaPSI_sums_list[[paste0("intron_event_ID_", singletypes[i], "_count.csv")]] %>%
    as.data.frame()
  rownames(intronsums_for_graph) <- intronsums_for_graph[, 1]

  for (j in rownames(intronsums_for_graph)) {
    values <- as.numeric(intronsums_for_graph[j, 4])
    if (!is.na(values)) {
      consolidated_row_name <- sub("\\.\\d+$", "", j)
      row_idx <- which(rownames(intron_sums_heatmap) == consolidated_row_name)
      intron_sums_heatmap[row_idx, i] <- sum(values, na.rm = TRUE)
    }
  }
}

maxs <- rowMaxs(as.matrix(intron_sums_heatmap))
intron_sums_heatmap <- cbind(intron_sums_heatmap, maxs)
intron_sums_heatmap <- intron_sums_heatmap[order(unlist(intron_sums_heatmap[,47]), decreasing = TRUE),]
intron_sums_heatmap <- intron_sums_heatmap[,-47]

# MXE:
MXE_events_list <- character()
for (i in singletypes) {
  MXE_events <- deltaPSI_sums_list[[paste0("MXE_event_ID_", i, "_count.csv")]]$MXE_event_ID
  MXE_events_list <- c(MXE_events_list, MXE_events)
}

row_names <- unique(MXE_events_list)
MXE_sums_heatmap <- matrix(nrow = length(row_names), ncol = 46, dimnames = list(row_names, singletypes))
MXE_sums_heatmap[is.na(MXE_sums_heatmap)] <- 0

for (i in seq_along(singletypes)) {
  MXEsums_for_graph <- deltaPSI_sums_list[[paste0("MXE_event_ID_", singletypes[i], "_count.csv")]] %>%
    as.data.frame()
  rownames(MXEsums_for_graph) <- MXEsums_for_graph[, 1]

  for (j in rownames(MXEsums_for_graph)) {
    values <- as.numeric(MXEsums_for_graph[j, 4])
    if (!is.na(values)) {
      consolidated_row_name <- sub("\\.\\d+$", "", j)
      row_idx <- which(rownames(MXE_sums_heatmap) == consolidated_row_name)
      MXE_sums_heatmap[row_idx, i] <- sum(values, na.rm = TRUE)
    }
  }
}

maxs <- rowMaxs(as.matrix(MXE_sums_heatmap))
MXE_sums_heatmap <- cbind(MXE_sums_heatmap, maxs)
MXE_sums_heatmap <- MXE_sums_heatmap[order(unlist(MXE_sums_heatmap[,47]), decreasing = TRUE),]
MXE_sums_heatmap <- MXE_sums_heatmap[,-47]
```

# Optional: add gene names to each A.S. event ID:
```{r}
for (i in 1:46) {
  A3Ssums_for_graph <- deltaPSI_sums_list[[i]] %>%
    as.data.frame()
  rownames(A3Ssums_for_graph) <- A3Ssums_for_graph[, 1]
  
  for (j in 1:nrow(A3S_sums_heatmap)) {
    rowname <- rownames(A3S_sums_heatmap)[j]
    index <- match(rowname, A3Ssums_for_graph$A3S_event_ID)
    if (!is.na(index)) {
      A3S_event_ID <- A3Ssums_for_graph$A3S_event_ID[index]
      Gene <- A3Ssums_for_graph$Gene[index]
      combined_name <- paste(Gene, A3S_event_ID)
      rownames(A3S_sums_heatmap)[j] <- combined_name
    }
  }
}

for (i in 1:46) {
  A5Ssums_for_graph <- deltaPSI_sums_list[[i+46]] %>%
    as.data.frame()
  rownames(A5Ssums_for_graph) <- A5Ssums_for_graph[, 1]
  
  for (j in 1:nrow(A5S_sums_heatmap)) {
    rowname <- rownames(A5S_sums_heatmap)[j]
    index <- match(rowname, A5Ssums_for_graph$A5S_event_ID)
    if (!is.na(index)) {
      A5S_event_ID <- A5Ssums_for_graph$A5S_event_ID[index]
      Gene <- A5Ssums_for_graph$Gene[index]
      combined_name <- paste(Gene, A5S_event_ID)
      rownames(A5S_sums_heatmap)[j] <- combined_name
    }
  }
}

for (i in 1:46) {
  cassettesums_for_graph <- deltaPSI_sums_list[[i+92]] %>%
    as.data.frame()
  rownames(cassettesums_for_graph) <- cassettesums_for_graph[, 1]
  
  for (j in 1:nrow(cassette_sums_heatmap)) {
    rowname <- rownames(cassette_sums_heatmap)[j]
    index <- match(rowname, cassettesums_for_graph$cassette_event_ID)
    if (!is.na(index)) {
      cassette_event_ID <- cassettesums_for_graph$cassette_event_ID[index]
      Gene <- cassettesums_for_graph$Gene[index]
      combined_name <- paste(Gene, cassette_event_ID)
      rownames(cassette_sums_heatmap)[j] <- combined_name
    }
  }
}

for (i in 1:46) {
  compositesums_for_graph <- deltaPSI_sums_list[[i+138]] %>%
    as.data.frame()
  rownames(compositesums_for_graph) <- compositesums_for_graph[, 1]
  
  for (j in 1:nrow(composite_sums_heatmap)) {
    rowname <- rownames(composite_sums_heatmap)[j]
    index <- match(rowname, compositesums_for_graph$composite_event_ID)
    if (!is.na(index)) {
      composite_event_ID <- compositesums_for_graph$composite_event_ID[index]
      Gene <- compositesums_for_graph$Gene[index]
      combined_name <- paste(Gene, composite_event_ID)
      rownames(composite_sums_heatmap)[j] <- combined_name
    }
  }
}

for (i in 1:46) {
  intronsums_for_graph <- deltaPSI_sums_list[[i+184]] %>%
    as.data.frame()
  rownames(intronsums_for_graph) <- intronsums_for_graph[, 1]
  
  for (j in 1:nrow(intron_sums_heatmap)) {
    rowname <- rownames(intron_sums_heatmap)[j]
    index <- match(rowname, intronsums_for_graph$intron_event_ID)
    if (!is.na(index)) {
      intron_event_ID <- intronsums_for_graph$intron_event_ID[index]
      Gene <- intronsums_for_graph$Gene[index]
      combined_name <- paste(Gene, intron_event_ID)
      rownames(intron_sums_heatmap)[j] <- combined_name
    }
  }
}

for (i in 1:46) {
  MXEsums_for_graph <- deltaPSI_sums_list[[i+230]] %>%
    as.data.frame()
  rownames(MXEsums_for_graph) <- MXEsums_for_graph[, 1]
  
  for (j in 1:nrow(MXE_sums_heatmap)) {
    rowname <- rownames(MXE_sums_heatmap)[j]
    index <- match(rowname, MXEsums_for_graph$MXE_event_ID)
    if (!is.na(index)) {
      MXE_event_ID <- MXEsums_for_graph$MXE_event_ID[index]
      Gene <- MXEsums_for_graph$Gene[index]
      combined_name <- paste(Gene, MXE_event_ID)
      rownames(MXE_sums_heatmap)[j] <- combined_name
    }
  }
}
```

# Plot the heatmaps:
```{r}
col_pallete <- colorRampPalette(c("white", "white", "white", "gold", "gold", "darkgoldenrod1", "darkorange1"))
colors <- col_pallete(250)
pheatmap(A3S_sums_heatmap[1:30,],
         border_color = "black",
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         number_format = "%.2f",
         number_color = "black",
         fontsize = 6.5,
         fontsize_number = 6,
         angle_col = 45,
         col = colors,
         main = "Sum of All First ΔPSI Values for A3S Events in Each Cell Type")

col_pallete <- colorRampPalette(c("white", "white", "white", "gold", "gold", "darkgoldenrod1", "darkorange1"))
colors <- col_pallete(250)
pheatmap(A5S_sums_heatmap[1:30,],
         border_color = "black",
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         number_format = "%.2f",
         number_color = "black",
         fontsize = 6.5,
         fontsize_number = 6,
         angle_col = 45,
         col = colors,
         main = "Sum of All First ΔPSI Values for A5S Events in Each Cell Type")

col_pallete <- colorRampPalette(c("white", "white", "white", "gold", "gold", "gold", "darkgoldenrod1", "darkorange1", "darkorange1"))
colors <- col_pallete(250)
pheatmap(cassette_sums_heatmap[1:30,],
         border_color = "black",
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         number_format = "%.2f",
         number_color = "black",
         fontsize = 6.5,
         fontsize_number = 6,
         angle_col = 45,
         col = colors,
         main = "Sum of All First ΔPSI Values for Cassette Exon Events in Each Cell Type")

col_pallete <- colorRampPalette(c("white", "white", "white", "white", "white", "gold", "gold", "darkgoldenrod1", "darkgoldenrod1", "darkorange1"))
colors <- col_pallete(250)
pheatmap(composite_sums_heatmap[1:30,],
         border_color = "black",
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = T,
         display_numbers = F,
         number_format = "%.2f",
         number_color = "black",
         fontsize = 6.5,
         fontsize_number = 6,
         angle_col = 45,
         col = colors,
         main = "Sum of All First ΔPSI Values for Composite Events in Each Cell Type")

col_pallete <- colorRampPalette(c("white", "white", "gold", "gold", "gold", "darkgoldenrod1", "darkgoldenrod1", "darkorange1"))
colors <- col_pallete(250)
pheatmap(intron_sums_heatmap[1:30,],
         border_color = "black",
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         number_format = "%.2f",
         number_color = "black",
         fontsize = 6.5,
         fontsize_number = 6,
         angle_col = 45,
         col = colors,
         main = "Sum of All First ΔPSI Values for Intron Retention Events in Each Cell Type")

col_pallete <- colorRampPalette(c("white", "white", "white", "white", "white", "white", "white",  "gold", "gold", "darkgoldenrod1", "darkorange1"))
colors <- col_pallete(250)
pheatmap(MXE_sums_heatmap[1:30,],
         border_color = "black",
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         number_format = "%.2f",
         number_color = "black",
         fontsize = 6.5,
         fontsize_number = 6,
         angle_col = 45,
         col = colors,
         main = "Sum of All First ΔPSI Values for MXE Events in Each Cell Type")
```

# Plot the heatmaps with highly specific break points:
```{r}
# Define the specific range for the A3S matrix
min_matrix <- min(A3S_sums_heatmap[1:30,])
max_matrix <- max(A3S_sums_heatmap[1:30,])

# Define the breaks for the A3S matrix
breaks_matrix <- c(seq(min_matrix, 0, by = abs(min(A3S_sums_heatmap[1:30,])/30)), seq(0, max_matrix, by = max(A3S_sums_heatmap[1:30,])/30), 0.001)
breaks_matrix[breaks_matrix == 0] <- 0.001  # Set 0 as a small positive value to represent white
breaks_matrix <- unique(breaks_matrix)

# Calculate the number of breaks for negative and positive values separately
n_breaks_neg <- sum(breaks_matrix < 0)
n_breaks_pos <- length(breaks_matrix) - n_breaks_neg

# Generate the color palette for negative values in the matrix using colorRampPalette
colors_neg <- colorRampPalette(c("darkblue", "dodgerblue3", "dodgerblue1", "lightblue", "white"))(n = n_breaks_neg)

# Generate the color palette for positive values in the matrix using colorRampPalette
colors_pos <- colorRampPalette(c("lightyellow", "khaki", "orange", "darkorange"))(n = n_breaks_pos)

# Combine the color palettes for negative and positive values
colors_matrix <- c(colors_neg, colors_pos)

# Create the heatmap
pheatmap(A3S_sums_heatmap[1:30, ],
         border_color = "black",
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = FALSE,
         number_format = "%.2f",
         number_color = "black",
         fontsize = 6.5,
         fontsize_number = 6,
         angle_col = 90,
         col = colors_matrix,
         breaks = breaks_matrix,
         main = "Sum of All First ΔPSI Values for A3S Events in Each Cell Type")


# Define the specific range for the A5S matrix
min_matrix <- min(A5S_sums_heatmap[1:30,])
max_matrix <- max(A5S_sums_heatmap[1:30,])

# Define the breaks for the A5S matrix
breaks_matrix <- c(seq(min_matrix, 0, by = abs(min(A5S_sums_heatmap[1:30,])/30)), seq(0, max_matrix, by = max(A5S_sums_heatmap[1:30,])/30), 0.001)
breaks_matrix[breaks_matrix == 0] <- 0.001  # Set 0 as a small positive value to represent white
breaks_matrix <- unique(breaks_matrix)

# Calculate the number of breaks for negative and positive values separately
n_breaks_neg <- sum(breaks_matrix < 0)
n_breaks_pos <- length(breaks_matrix) - n_breaks_neg

# Generate the color palette for negative values in the matrix using colorRampPalette
colors_neg <- colorRampPalette(c("darkblue", "dodgerblue3", "dodgerblue1", "lightblue", "white"))(n = n_breaks_neg)

# Generate the color palette for positive values in the matrix using colorRampPalette
colors_pos <- colorRampPalette(c("lightyellow", "khaki", "orange", "darkorange"))(n = n_breaks_pos)

# Combine the color palettes for negative and positive values
colors_matrix <- c(colors_neg, colors_pos)

# Create the heatmap
pheatmap(A5S_sums_heatmap[1:30, ],
         border_color = "black",
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = FALSE,
         number_format = "%.2f",
         number_color = "black",
         fontsize = 6.5,
         fontsize_number = 6,
         angle_col = 90,
         col = colors_matrix,
         breaks = breaks_matrix,
         main = "Sum of All First ΔPSI Values for A5S Events in Each Cell Type")


# Define the specific range for the cassette matrix
min_matrix <- min(cassette_sums_heatmap[1:30,])
max_matrix <- max(cassette_sums_heatmap[1:30,])

# Define the breaks for the cassette matrix
breaks_matrix <- c(seq(min_matrix, 0, by = abs(min(cassette_sums_heatmap[1:30,])/30)), seq(0, max_matrix, by = max(cassette_sums_heatmap[1:30,])/30), 0.001)
breaks_matrix[breaks_matrix == 0] <- 0.001  # Set 0 as a small positive value to represent white
breaks_matrix <- unique(breaks_matrix)

# Calculate the number of breaks for negative and positive values separately
n_breaks_neg <- sum(breaks_matrix < 0)
n_breaks_pos <- length(breaks_matrix) - n_breaks_neg

# Generate the color palette for negative values in the matrix using colorRampPalette
colors_neg <- colorRampPalette(c("darkblue", "dodgerblue3", "dodgerblue1", "lightblue", "white"))(n = n_breaks_neg)

# Generate the color palette for positive values in the matrix using colorRampPalette
colors_pos <- colorRampPalette(c("lightyellow", "khaki", "orange", "darkorange"))(n = n_breaks_pos)

# Combine the color palettes for negative and positive values
colors_matrix <- c(colors_neg, colors_pos)

# Create the heatmap
pheatmap(cassette_sums_heatmap[1:30, ],
         border_color = "black",
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = FALSE,
         number_format = "%.2f",
         number_color = "black",
         fontsize = 6.5,
         fontsize_number = 6,
         angle_col = 90,
         col = colors_matrix,
         breaks = breaks_matrix,
         main = "Sum of All First ΔPSI Values for Cassette Exon Events in Each Cell Type")


# Define the specific range for the composite matrix
min_matrix <- min(composite_sums_heatmap[1:30,])
max_matrix <- max(composite_sums_heatmap[1:30,])

# Define the breaks for the composite matrix
breaks_matrix <- c(seq(min_matrix, 0, by = abs(min(composite_sums_heatmap[1:30,])/30)), seq(0, max_matrix, by = max(composite_sums_heatmap[1:30,])/30), 0.001)
breaks_matrix[breaks_matrix == 0] <- 0.001  # Set 0 as a small positive value to represent white
breaks_matrix <- unique(breaks_matrix)

# Calculate the number of breaks for negative and positive values separately
n_breaks_neg <- sum(breaks_matrix < 0)
n_breaks_pos <- length(breaks_matrix) - n_breaks_neg

# Generate the color palette for negative values in the matrix using colorRampPalette
colors_neg <- colorRampPalette(c("darkblue", "dodgerblue3", "dodgerblue1", "lightblue", "white"))(n = n_breaks_neg)

# Generate the color palette for positive values in the matrix using colorRampPalette
colors_pos <- colorRampPalette(c("lightyellow", "khaki", "orange", "darkorange"))(n = n_breaks_pos)

# Combine the color palettes for negative and positive values
colors_matrix <- c(colors_neg, colors_pos)

# Create the heatmap
pheatmap(composite_sums_heatmap[1:30, ],
         border_color = "black",
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = FALSE,
         number_format = "%.2f",
         number_color = "black",
         fontsize = 6,
         fontsize_number = 6,
         angle_col = 90,
         col = colors_matrix,
         breaks = breaks_matrix,
         main = "Sum of All First ΔPSI Values for Composite Events in Each Cell Type")


# Define the specific range for the intron matrix
min_matrix <- min(intron_sums_heatmap[1:30,])
max_matrix <- max(intron_sums_heatmap[1:30,])

# Define the breaks for the intron matrix
breaks_matrix <- c(seq(min_matrix, 0, by = abs(min(intron_sums_heatmap[1:30,])/30)), seq(0, max_matrix, by = max(intron_sums_heatmap[1:30,])/30), 0.001)
breaks_matrix[breaks_matrix == 0] <- 0.001  # Set 0 as a small positive value to represent white
breaks_matrix <- unique(breaks_matrix)

# Calculate the number of breaks for negative and positive values separately
n_breaks_neg <- sum(breaks_matrix < 0)
n_breaks_pos <- length(breaks_matrix) - n_breaks_neg

# Generate the color palette for negative values in the matrix using colorRampPalette
colors_neg <- colorRampPalette(c("darkblue", "dodgerblue3", "dodgerblue1", "lightblue", "white"))(n = n_breaks_neg)

# Generate the color palette for positive values in the matrix using colorRampPalette
colors_pos <- colorRampPalette(c("lightyellow", "khaki", "orange", "darkorange"))(n = n_breaks_pos)

# Combine the color palettes for negative and positive values
colors_matrix <- c(colors_neg, colors_pos)

# Create the heatmap
pheatmap(intron_sums_heatmap[1:30, ],
         border_color = "black",
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = FALSE,
         number_format = "%.2f",
         number_color = "black",
         fontsize = 6.5,
         fontsize_number = 6,
         angle_col = 90,
         col = colors_matrix,
         breaks = breaks_matrix,
         main = "Sum of All First ΔPSI Values for Intron Retention Events in Each Cell Type")


# Define the specific range for the MXE matrix
min_matrix <- min(MXE_sums_heatmap[1:30,])
max_matrix <- max(MXE_sums_heatmap[1:30,])

# Define the breaks for the MXE matrix
breaks_matrix <- c(seq(min_matrix, 0, by = abs(min(MXE_sums_heatmap[1:30,])/30)), seq(0, max_matrix, by = max(MXE_sums_heatmap[1:30,])/30), 0.001)
breaks_matrix[breaks_matrix == 0] <- 0.001  # Set 0 as a small positive value to represent white
breaks_matrix <- unique(breaks_matrix)

# Calculate the number of breaks for negative and positive values separately
n_breaks_neg <- sum(breaks_matrix < 0)
n_breaks_pos <- length(breaks_matrix) - n_breaks_neg

# Generate the color palette for negative values in the matrix using colorRampPalette
colors_neg <- colorRampPalette(c("darkblue", "dodgerblue3", "dodgerblue1", "lightblue", "white"))(n = n_breaks_neg)

# Generate the color palette for positive values in the matrix using colorRampPalette
colors_pos <- colorRampPalette(c("lightyellow", "khaki", "orange", "darkorange"))(n = n_breaks_pos)

# Combine the color palettes for negative and positive values
colors_matrix <- c(colors_neg, colors_pos)

# Create the heatmap
pheatmap(MXE_sums_heatmap[1:30, ],
         border_color = "black",
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = FALSE,
         number_format = "%.2f",
         number_color = "black",
         fontsize = 6,
         fontsize_number = 6,
         angle_col = 90,
         col = colors_matrix,
         breaks = breaks_matrix,
         main = "Sum of All First ΔPSI Values for MXE Events in Each Cell Type")
```

# detect pan-neuronal genes:
```{r}
# these pan-neuronal genes were described by Stefanakis et al. 2015:
pan_genes <- c("nsf-1", "rab-3", "ric-4", "snb-1", "unc-64", "sng-1", "snt-1", "unc-10", "unc-18", "ehs-1", "unc-11", "unc-57", "snn-1", "unc-104", "syd-2", "egl-3", "egl-21", "ric-19", "unc-31", "unc-108", "maco-1", "rgef-1", "tbb-1")

# number of cells in which pan-neuronal genes are expressed in the single-cell data:
pan_sc <- single_cell_data[single_cell_data$gene_name %in% pan_genes,]

pan_sc$count <- rowSums(pan_sc[,3:ncol(pan_sc)] != 0)
median(pan_sc$count)

A3S_sums_heatmap_pan <- A3S_sums_heatmap[A3S_sums_heatmap$single_cell_expression_count >= median(pan_sc$count),]
A5S_sums_heatmap_pan <- A5S_sums_heatmap[A5S_sums_heatmap$single_cell_expression_count >= median(pan_sc$count),]
cassette_sums_heatmap_pan <- cassette_sums_heatmap[cassette_sums_heatmap$single_cell_expression_count >= median(pan_sc$count),]
composite_sums_heatmap_pan <- composite_sums_heatmap[composite_sums_heatmap$single_cell_expression_count >= median(pan_sc$count),]
intron_sums_heatmap_pan <- intron_sums_heatmap[intron_sums_heatmap$single_cell_expression_count >= median(pan_sc$count),]
MXE_sums_heatmap_pan <- MXE_sums_heatmap[MXE_sums_heatmap$single_cell_expression_count >= median(pan_sc$count),]

combined_pan_heatmap <- rbind(A3S_sums_heatmap_pan, A5S_sums_heatmap_pan, cassette_sums_heatmap_pan, composite_sums_heatmap_pan, intron_sums_heatmap_pan, MXE_sums_heatmap_pan)
unique(combined_pan_heatmap$gene)
#combined_pan_heatmap <- combined_pan_heatmap[,-48]

# sort the pan-neuronal heatmap by absolute value rowMaxs:
maxs <- rowMaxs(abs(as.matrix(combined_pan_heatmap[1:length(singletypes)])))
combined_pan_heatmap <- cbind(combined_pan_heatmap[1:length(singletypes)], maxs)
combined_pan_heatmap <- combined_pan_heatmap[order(unlist(combined_pan_heatmap[,47]), decreasing = TRUE),]
combined_pan_heatmap <- combined_pan_heatmap[,-47]
```

# plot the pan-neuronal heatmap:
```{r}
# Define the specific range for the intron matrix
min_matrix <- min(combined_pan_heatmap[1:30,])
max_matrix <- max(combined_pan_heatmap[1:30,])

# Define the breaks for the intron matrix
breaks_matrix <- c(seq(min_matrix, 0, by = abs(min(combined_pan_heatmap[1:30,])/30)), seq(0, max_matrix, by = max(combined_pan_heatmap[1:30,])/30), 0.001)
breaks_matrix[breaks_matrix == 0] <- 0.001  # Set 0 as a small positive value to represent white
breaks_matrix <- unique(breaks_matrix)

# Calculate the number of breaks for negative and positive values separately
n_breaks_neg <- sum(breaks_matrix < 0)
n_breaks_pos <- length(breaks_matrix) - n_breaks_neg

# Generate the color palette for negative values in the matrix using colorRampPalette
colors_neg <- colorRampPalette(c("darkblue", "dodgerblue3", "dodgerblue1", "lightblue", "white"))(n = n_breaks_neg)

# Generate the color palette for positive values in the matrix using colorRampPalette
colors_pos <- colorRampPalette(c("lightyellow", "khaki", "orange", "darkorange"))(n = n_breaks_pos)

# Combine the color palettes for negative and positive values
colors_matrix <- c(colors_neg, colors_pos)

# Create the heatmap
pheatmap(combined_pan_heatmap[1:30, ],
         border_color = "black",
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = FALSE,
         number_format = "%.2f",
         number_color = "black",
         fontsize = 7,
         fontsize_number = 6,
         angle_col = 90,
         col = colors_matrix,
         breaks = breaks_matrix,
         main = "Sum of All First ΔPSI Values for AS Events of Pan-Neuronal Genes in Each Cell Type")
```

# Plots for AS events:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

ggplot(data = read.csv("A3S_event_ID_ASG_count.csv"), aes(x = A3S_event_ID, y = Freq)) +
  geom_point(size = 2) +
  geom_segment(aes(x = A3S_event_ID, xend = A3S_event_ID, y = 0, yend = Freq)) +
  labs(title = "Frequency of A3S Events in ASG cells", x = "", y = "Count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3))

ggplot(data = read.csv("A3S_event_ID_ASG_count.csv"), aes(x = Freq, y = A3S_event_ID)) +
  geom_point(size = 2) +
  geom_segment(aes(x = 0, xend = Freq, y = A3S_event_ID, yend = A3S_event_ID)) +
  labs(title = "Frequency of A3S Events in ASG cells", x = "Count", y = "") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, size = 0.3))

ggplot(data = read.csv("A3S_event_ID_ASG_count.csv"), aes(x = A3S_event_ID, y = Freq, fill = A3S_event_ID)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Frequency of A3S Events in ASG cells", x = "", y = "Count") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3), legend.position = "none")

ggplot(data = read.csv("A3S_event_ID_ASG_count.csv"), aes(x = Freq, y = A3S_event_ID, fill = A3S_event_ID)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Frequency of A3S Events in ASG cells", x = "Count", y = "") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.3), legend.position = "none")
```

# Are A3S event coordinates multiples of 3? Frame shifting vs frame preservation
```{r}
A3S_sums_heatmap$AS_event_ID <- sapply(strsplit(rownames(A3S_sums_heatmap), " "), `[`, 2)
for (i in 1:nrow(A3S_sums_heatmap)){
elements <- unlist(strsplit(A3S_sums_heatmap$AS_event_ID[i], "_"))
if (elements[2] == "+"){    ### if positive-stranded AS event, extract last 2 coordinates of AS event
shift <- print(as.numeric(elements[5])-as.numeric(elements[4]))
if (shift %% 3 == 0){
  A3S_sums_heatmap$shift_or_preserve[i] <- "preserve"
}
else A3S_sums_heatmap$shift_or_preserve[i] <- "shift"
}
if (elements[2] == "-"){    ### if negative-stranded AS event, extract first 2 coordinates of AS event
shift <- print(as.numeric(elements[4])-as.numeric(elements[3]))
if (shift %% 3 == 0){
  A3S_sums_heatmap$shift_or_preserve[i] <- "preserve"
}
else A3S_sums_heatmap$shift_or_preserve[i] <- "shift"
}
}

shift_percentage <- sum(A3S_sums_heatmap$shift_or_preserve == "shift") / nrow(A3S_sums_heatmap) * 100
print(noquote(paste0(shift_percentage, "% of all alternative 3' splice sites are frame-shifting")))
preserve_percentage <- sum(A3S_sums_heatmap$shift_or_preserve == "preserve") / nrow(A3S_sums_heatmap) * 100
print(noquote(paste0(preserve_percentage, "% of all alternative 3' splice sites are frame-preserving")))

# OLL A3S investigation:
filtered_df <- subset(A3S_sums_heatmap, OLL > 0)

shift_percentage <- sum(filtered_df$shift_or_preserve == "shift") / nrow(filtered_df) * 100
print(noquote(paste0(shift_percentage, "% of alternative 3' splice sites in OLL are frame-shifting")))
preserve_percentage <- sum(filtered_df$shift_or_preserve == "preserve") / nrow(filtered_df) * 100
print(noquote(paste0(preserve_percentage, "% of alternative 3' splice sites in OLL are frame-preserving")))

elements_list <- unlist(filtered_df$AS_event_ID)

# Initialize a vector to store the nt differences
diff_values <- numeric(length(elements_list))

for (i in 1:length(elements_list)) {
  elements <- elements_list[[i]]
  elements <- unlist(strsplit(elements, "_"))
  if (elements[2] == "+") {
    result <- as.numeric(elements[5]) - as.numeric(elements[4])
  } else if (elements[2] == "-") {
    result <- as.numeric(elements[4]) - as.numeric(elements[3])
  }
  diff_values[i] <- result
}

# stats:
median_diff <- median(diff_values)
max_diff <- max(diff_values)
min_diff <- min(diff_values)

print(noquote(paste("The median size difference in significant alternative 3' splice sites in OLL neurons is:", median_diff, "nucleotides")))
print(noquote(paste("The maximum size difference in significant alternative 3' splice sites in OLL neurons is:", max_diff, "nucleotides")))
print(noquote(paste("The minimum size difference in significant alternative 3' splice sites in OLL neurons is:", min_diff, "nucleotides")))
```

# how about retained introns?
```{r}
intron_sums_heatmap$AS_event_ID <- sapply(strsplit(rownames(intron_sums_heatmap), " "), `[`, 2)
for (i in 1:nrow(intron_sums_heatmap)){
elements <- unlist(strsplit(intron_sums_heatmap$AS_event_ID[i], "_"))
shift <- print(as.numeric(elements[4])-as.numeric(elements[3]))
if (shift %% 3 == 0){
  intron_sums_heatmap$shift_or_preserve[i] <- "preserve"
}
else intron_sums_heatmap$shift_or_preserve[i] <- "shift"
}

#### this script continues in "12. Identifying stop codons and reading frames in intron retention and A3S events.Rmd"
```

# PCAs of different splice types:
```{r}
### A3S events:
mega_A3S_PCA <- sapply(mega_A3S[,5:50], function(x) as.numeric(gsub("%", "", x)))
mega_A3S_PCA <- mega_A3S_PCA[rowSums(is.na(mega_A3S_PCA)) != ncol(mega_A3S_PCA), ]

for(i in 1:ncol(mega_A3S_PCA)) {
  mega_A3S_PCA[is.na(mega_A3S_PCA[,i]), i] <- mean(mega_A3S_PCA[,i], na.rm = TRUE)
}

data_t <- t(mega_A3S_PCA)
data_scaled <- scale(data_t)
pca_res <- prcomp(data_scaled, center = TRUE, scale. = TRUE)
pca_var <- pca_res$sdev^2
pca_var_perc <- pca_var / sum(pca_var) * 100
pca_data <- data.frame(pca_res$x)
pca_data$CellType <- rownames(pca_data)

pca_data <- pca_data %>%
  mutate(CellTypeBase = ifelse(grepl("\\d$", CellType), CellType, sub(".*(\\d+)$", "\\1", CellType)))

pca_data <- left_join(pca_data, condition, by = c("CellTypeBase" = "type"))

# plot the PCA:
ggplot(pca_data, aes(x = PC1, y = PC2, label = CellType, color = condition)) +
  geom_point(size = 2) +
  geom_text(aes(label = CellType), hjust = 1.25, vjust = 1.25, size = 3) +
  labs(
    title = "PCA of A3S Events (upstreamness)",
    x = paste0("Principal Component 1 (", round(pca_var_perc[1], 2), "%)"),
    y = paste0("Principal Component 2 (", round(pca_var_perc[2], 2), "%)")
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Sensory" = "dodgerblue", "Motor" = "gold3", "Interneuron" = "salmon", "Polymodal" = "darkgreen")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

### A5S events:
mega_A5S_PCA <- sapply(mega_A5S[,5:50], function(x) as.numeric(gsub("%", "", x)))
mega_A5S_PCA <- mega_A5S_PCA[rowSums(is.na(mega_A5S_PCA)) != ncol(mega_A5S_PCA), ]

for(i in 1:ncol(mega_A5S_PCA)) {
  mega_A5S_PCA[is.na(mega_A5S_PCA[,i]), i] <- mean(mega_A5S_PCA[,i], na.rm = TRUE)
}

data_t <- t(mega_A5S_PCA)
data_scaled <- scale(data_t)
pca_res <- prcomp(data_scaled, center = TRUE, scale. = TRUE)
pca_var <- pca_res$sdev^2
pca_var_perc <- pca_var / sum(pca_var) * 100
pca_data <- data.frame(pca_res$x)
pca_data$CellType <- rownames(pca_data)

pca_data <- pca_data %>%
  mutate(CellTypeBase = ifelse(grepl("\\d$", CellType), CellType, sub(".*(\\d+)$", "\\1", CellType)))

pca_data <- left_join(pca_data, condition, by = c("CellTypeBase" = "type"))

# plot the PCA:
ggplot(pca_data, aes(x = PC1, y = PC2, label = CellType, color = condition)) +
  geom_point(size = 2) +
  geom_text(aes(label = CellType), hjust = 1.25, vjust = 1.25, size = 3) +
  labs(
    title = "PCA of A5S Events (upstreamness)",
    x = paste0("Principal Component 1 (", round(pca_var_perc[1], 2), "%)"),
    y = paste0("Principal Component 2 (", round(pca_var_perc[2], 2), "%)")
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Sensory" = "dodgerblue", "Motor" = "gold3", "Interneuron" = "salmon", "Polymodal" = "darkgreen")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

### cassette events:
mega_cassette_PCA <- sapply(mega_cassette[,5:50], function(x) as.numeric(gsub("%", "", x)))
mega_cassette_PCA <- mega_cassette_PCA[rowSums(is.na(mega_cassette_PCA)) != ncol(mega_cassette_PCA), ]

for(i in 1:ncol(mega_cassette_PCA)) {
  mega_cassette_PCA[is.na(mega_cassette_PCA[,i]), i] <- mean(mega_cassette_PCA[,i], na.rm = TRUE)
}

data_t <- t(mega_cassette_PCA)
data_scaled <- scale(data_t)
pca_res <- prcomp(data_scaled, center = TRUE, scale. = TRUE)
pca_var <- pca_res$sdev^2
pca_var_perc <- pca_var / sum(pca_var) * 100
pca_data <- data.frame(pca_res$x)
pca_data$CellType <- rownames(pca_data)

pca_data <- pca_data %>%
  mutate(CellTypeBase = ifelse(grepl("\\d$", CellType), CellType, sub(".*(\\d+)$", "\\1", CellType)))

pca_data <- left_join(pca_data, condition, by = c("CellTypeBase" = "type"))

# plot the PCA:
ggplot(pca_data, aes(x = PC1, y = PC2, label = CellType, color = condition)) +
  geom_point(size = 2) +
  geom_text(aes(label = CellType), hjust = 1.25, vjust = 1.25, size = 3) +
  labs(
    title = "PCA of Cassette Events",
    x = paste0("Principal Component 1 (", round(pca_var_perc[1], 2), "%)"),
    y = paste0("Principal Component 2 (", round(pca_var_perc[2], 2), "%)")
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Sensory" = "dodgerblue", "Motor" = "gold3", "Interneuron" = "salmon", "Polymodal" = "darkgreen")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

### composite events:
mega_composite_PCA <- sapply(mega_composite[,5:50], function(x) as.numeric(gsub("%", "", x)))
mega_composite_PCA <- mega_composite_PCA[rowSums(is.na(mega_composite_PCA)) != ncol(mega_composite_PCA), ]

for(i in 1:ncol(mega_composite_PCA)) {
  mega_composite_PCA[is.na(mega_composite_PCA[,i]), i] <- mean(mega_composite_PCA[,i], na.rm = TRUE)
}

data_t <- t(mega_composite_PCA)
data_scaled <- scale(data_t)
pca_res <- prcomp(data_scaled, center = TRUE, scale. = TRUE)
pca_var <- pca_res$sdev^2
pca_var_perc <- pca_var / sum(pca_var) * 100
pca_data <- data.frame(pca_res$x)
pca_data$CellType <- rownames(pca_data)

pca_data <- pca_data %>%
  mutate(CellTypeBase = ifelse(grepl("\\d$", CellType), CellType, sub(".*(\\d+)$", "\\1", CellType)))

pca_data <- left_join(pca_data, condition, by = c("CellTypeBase" = "type"))

# plot the PCA:
ggplot(pca_data, aes(x = PC1, y = PC2, label = CellType, color = condition)) +
  geom_point(size = 2) +
  geom_text(aes(label = CellType), hjust = 1.25, vjust = 1.25, size = 3) +
  labs(
    title = "PCA of Composite Events",
    x = paste0("Principal Component 1 (", round(pca_var_perc[1], 2), "%)"),
    y = paste0("Principal Component 2 (", round(pca_var_perc[2], 2), "%)")
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Sensory" = "dodgerblue", "Motor" = "gold3", "Interneuron" = "salmon", "Polymodal" = "darkgreen")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

### intron events:
mega_intron_PCA <- sapply(mega_intron[,5:50], function(x) as.numeric(gsub("%", "", x)))
mega_intron_PCA <- mega_intron_PCA[rowSums(is.na(mega_intron_PCA)) != ncol(mega_intron_PCA), ]

for(i in 1:ncol(mega_intron_PCA)) {
  mega_intron_PCA[is.na(mega_intron_PCA[,i]), i] <- mean(mega_intron_PCA[,i], na.rm = TRUE)
}

data_t <- t(mega_intron_PCA)
data_scaled <- scale(data_t)
pca_res <- prcomp(data_scaled, center = TRUE, scale. = TRUE)
pca_var <- pca_res$sdev^2
pca_var_perc <- pca_var / sum(pca_var) * 100
pca_data <- data.frame(pca_res$x)
pca_data$CellType <- rownames(pca_data)

pca_data <- pca_data %>%
  mutate(CellTypeBase = ifelse(grepl("\\d$", CellType), CellType, sub(".*(\\d+)$", "\\1", CellType)))

pca_data <- left_join(pca_data, condition, by = c("CellTypeBase" = "type"))

# plot the PCA:
ggplot(pca_data, aes(x = PC1, y = PC2, label = CellType, color = condition)) +
  geom_point(size = 2) +
  geom_text(aes(label = CellType), hjust = 1.25, vjust = 1.25, size = 3) +
  labs(
    title = "PCA of Intron Retention Events",
    x = paste0("Principal Component 1 (", round(pca_var_perc[1], 2), "%)"),
    y = paste0("Principal Component 2 (", round(pca_var_perc[2], 2), "%)")
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Sensory" = "dodgerblue", "Motor" = "gold3", "Interneuron" = "salmon", "Polymodal" = "darkgreen")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

### MXE events:
mega_MXE_PCA <- sapply(mega_MXE[,5:50], function(x) as.numeric(gsub("%", "", x)))
mega_MXE_PCA <- mega_MXE_PCA[rowSums(is.na(mega_MXE_PCA)) != ncol(mega_MXE_PCA), ]

for(i in 1:ncol(mega_MXE_PCA)) {
  mega_MXE_PCA[is.na(mega_MXE_PCA[,i]), i] <- mean(mega_MXE_PCA[,i], na.rm = TRUE)
}

data_t <- t(mega_MXE_PCA)
data_scaled <- scale(data_t)
pca_res <- prcomp(data_scaled, center = TRUE, scale. = TRUE)
pca_var <- pca_res$sdev^2
pca_var_perc <- pca_var / sum(pca_var) * 100
pca_data <- data.frame(pca_res$x)
pca_data$CellType <- rownames(pca_data)

pca_data <- pca_data %>%
  mutate(CellTypeBase = ifelse(grepl("\\d$", CellType), CellType, sub(".*(\\d+)$", "\\1", CellType)))

pca_data <- left_join(pca_data, condition, by = c("CellTypeBase" = "type"))

# plot the PCA:
ggplot(pca_data, aes(x = PC1, y = PC2, label = CellType, color = condition)) +
  geom_point(size = 2) +
  geom_text(aes(label = CellType), hjust = 1.25, vjust = 1.25, size = 3) +
  labs(
    title = "PCA of MXE Events (upstreamness)",
    x = paste0("Principal Component 1 (", round(pca_var_perc[1], 2), "%)"),
    y = paste0("Principal Component 2 (", round(pca_var_perc[2], 2), "%)")
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Sensory" = "dodgerblue", "Motor" = "gold3", "Interneuron" = "salmon", "Polymodal" = "darkgreen")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

### All AS events combined:
all_datasets <- list(mega_A3S, mega_A5S, mega_cassette, mega_composite, mega_intron, mega_MXE)
combined_data <- do.call(rbind, all_datasets)

# Convert necessary columns to numeric, replacing '%' and handling NA values
mega_ALL_PCA <- sapply(combined_data[, 5:50], function(x) as.numeric(gsub("%", "", x)))

mega_ALL_PCA <- mega_ALL_PCA[rowSums(is.na(mega_ALL_PCA)) != ncol(mega_ALL_PCA), ]

for(i in 1:ncol(mega_ALL_PCA)) {
  mega_ALL_PCA[is.na(mega_ALL_PCA[,i]), i] <- mean(mega_ALL_PCA[,i], na.rm = TRUE)
}

data_t <- t(mega_ALL_PCA)
data_scaled <- scale(data_t)
pca_res <- prcomp(data_scaled, center = TRUE, scale. = TRUE)
pca_var <- pca_res$sdev^2
pca_var_perc <- pca_var / sum(pca_var) * 100
pca_data <- data.frame(pca_res$x)
pca_data$CellType <- rownames(pca_data)

pca_data <- pca_data %>%
  mutate(CellTypeBase = ifelse(grepl("\\d$", CellType), CellType, sub(".*(\\d+)$", "\\1", CellType)))

pca_data <- left_join(pca_data, condition, by = c("CellTypeBase" = "type"))

# plot the PCA:
ggplot(pca_data, aes(x = PC1, y = PC2, label = CellType, color = condition)) +
  geom_point(size = 2) +
  geom_text(aes(label = CellType), hjust = 1.25, vjust = 1.25, size = 3) +
  labs(
    title = "PCA of ALL AS Events Combined",
    x = paste0("Principal Component 1 (", round(pca_var_perc[1], 2), "%)"),
    y = paste0("Principal Component 2 (", round(pca_var_perc[2], 2), "%)")
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Sensory" = "dodgerblue", "Motor" = "gold3", "Interneuron" = "salmon", "Polymodal" = "darkgreen")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

### All AS events combined (minus composite events):
all_datasets <- list(mega_A3S, mega_A5S, mega_cassette, mega_intron, mega_MXE)
combined_data <- do.call(rbind, all_datasets)

# Convert necessary columns to numeric, replacing '%' and handling NA values
mega_ALL_PCA <- sapply(combined_data[, 5:50], function(x) as.numeric(gsub("%", "", x)))

mega_ALL_PCA <- mega_ALL_PCA[rowSums(is.na(mega_ALL_PCA)) != ncol(mega_ALL_PCA), ]

for(i in 1:ncol(mega_ALL_PCA)) {
  mega_ALL_PCA[is.na(mega_ALL_PCA[,i]), i] <- mean(mega_ALL_PCA[,i], na.rm = TRUE)
}

data_t <- t(mega_ALL_PCA)
data_scaled <- scale(data_t)
pca_res <- prcomp(data_scaled, center = TRUE, scale. = TRUE)
pca_var <- pca_res$sdev^2
pca_var_perc <- pca_var / sum(pca_var) * 100
pca_data <- data.frame(pca_res$x)
pca_data$CellType <- rownames(pca_data)

pca_data <- pca_data %>%
  mutate(CellTypeBase = ifelse(grepl("\\d$", CellType), CellType, sub(".*(\\d+)$", "\\1", CellType)))

pca_data <- left_join(pca_data, condition, by = c("CellTypeBase" = "type"))

# plot the PCA:
ggplot(pca_data, aes(x = PC1, y = PC2, label = CellType, color = condition)) +
  geom_point(size = 2) +
  geom_text(aes(label = CellType), hjust = 1.25, vjust = 1.25, size = 3) +
  labs(
    title = "PCA of ALL AS Events Combined (minus composite events)",
    x = paste0("Principal Component 1 (", round(pca_var_perc[1], 2), "%)"),
    y = paste0("Principal Component 2 (", round(pca_var_perc[2], 2), "%)")
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Sensory" = "dodgerblue", "Motor" = "gold3", "Interneuron" = "salmon", "Polymodal" = "darkgreen")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
```

# Combine every cell type for each A.S. event into one data frame to get the "most interesting" A.S. events in all comparisons:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

A3S_events_dataframe <- data.frame()

# Loop over the comparisons
for (i in singletypes) {
  # Construct the file name
  file_name <- paste0("A3S_event_ID_", i, "_count.csv")
  
  # Check if the file exists
  if (file.exists(file_name)) {
    # Read the CSV file
    df <- read.csv(file_name)
    
    # Append the data frame to A3S_events_dataframe
    A3S_events_dataframe <- rbind(A3S_events_dataframe, df)
  }
}

A3S_events_dataframe <- ddply(A3S_events_dataframe, c("A3S_event_ID", "Gene"), numcolwise(sum))
A3S_events_dataframe$Freq <- A3S_events_dataframe$Freq/2 # We have to divide by 2 to get the true frequency of A.S. events because there are 2 datasets for each comparison
A3S_events_dataframe <- A3S_events_dataframe[rev(order(A3S_events_dataframe$Freq)),]

# Create an empty data frame to store the results
A5S_events_dataframe <- data.frame()

# Loop over the comparisons
for (i in singletypes) {
  # Construct the file name
  file_name <- paste0("A5S_event_ID_", i, "_count.csv")
  
  # Check if the file exists
  if (file.exists(file_name)) {
    # Read the CSV file
    df <- read.csv(file_name)
    
    # Append the data frame to A5S_events_dataframe
    A5S_events_dataframe <- rbind(A5S_events_dataframe, df)
  }
}

A5S_events_dataframe <- ddply(A5S_events_dataframe, c("A5S_event_ID", "Gene"), numcolwise(sum))
A5S_events_dataframe$Freq <- A5S_events_dataframe$Freq/2
A5S_events_dataframe <- A5S_events_dataframe[rev(order(A5S_events_dataframe$Freq)),]

# Create an empty data frame to store the results
cassette_events_dataframe <- data.frame()

# Loop over the comparisons
for (i in singletypes) {
  # Construct the file name
  file_name <- paste0("cassette_event_ID_", i, "_count.csv")
  
  # Check if the file exists
  if (file.exists(file_name)) {
    # Read the CSV file
    df <- read.csv(file_name)
    
    # Append the data frame to cassette_events_dataframe
    cassette_events_dataframe <- rbind(cassette_events_dataframe, df)
  }
}

cassette_events_dataframe <- ddply(cassette_events_dataframe, c("cassette_event_ID", "Gene"), numcolwise(sum))
cassette_events_dataframe$Freq <- cassette_events_dataframe$Freq/2
cassette_events_dataframe <- cassette_events_dataframe[rev(order(cassette_events_dataframe$Freq)),]

# Create an empty data frame to store the results
composite_events_dataframe <- data.frame()

# Loop over the comparisons
for (i in singletypes) {
  # Construct the file name
  file_name <- paste0("composite_event_ID_", i, "_count.csv")
  
  # Check if the file exists
  if (file.exists(file_name)) {
    # Read the CSV file
    df <- read.csv(file_name)
    
    # Append the data frame to composite_events_dataframe
    composite_events_dataframe <- rbind(composite_events_dataframe, df)
  }
}

composite_events_dataframe <- ddply(composite_events_dataframe, c("composite_event_ID", "Gene"), numcolwise(sum))
composite_events_dataframe$Freq <- composite_events_dataframe$Freq/2
composite_events_dataframe <- composite_events_dataframe[rev(order(composite_events_dataframe$Freq)),]

# Create an empty data frame to store the results
intron_events_dataframe <- data.frame()

# Loop over the comparisons
for (i in singletypes) {
  # Construct the file name
  file_name <- paste0("intron_event_ID_", i, "_count.csv")
  
  # Check if the file exists
  if (file.exists(file_name)) {
    # Read the CSV file
    df <- read.csv(file_name)
    
    # Append the data frame to intron_events_dataframe
    intron_events_dataframe <- rbind(intron_events_dataframe, df)
  }
}

intron_events_dataframe <- ddply(intron_events_dataframe, c("intron_event_ID", "Gene"), numcolwise(sum))
intron_events_dataframe$Freq <- intron_events_dataframe$Freq/2
intron_events_dataframe <- intron_events_dataframe[rev(order(intron_events_dataframe$Freq)),]

# Create an empty data frame to store the results
MXE_events_dataframe <- data.frame()

# Loop over the comparisons
for (i in singletypes) {
  # Construct the file name
  file_name <- paste0("MXE_event_ID_", i, "_count.csv")
  
  # Check if the file exists
  if (file.exists(file_name)) {
    # Read the CSV file
    df <- read.csv(file_name)
    
    # Append the data frame to MXE_events_dataframe
    MXE_events_dataframe <- rbind(MXE_events_dataframe, df)
  }
}

MXE_events_dataframe <- ddply(MXE_events_dataframe, c("MXE_event_ID", "Gene"), numcolwise(sum))
MXE_events_dataframe$Freq <- MXE_events_dataframe$Freq/2
MXE_events_dataframe <- MXE_events_dataframe[rev(order(MXE_events_dataframe$Freq)),]
```

# Plots of A.S. events in all comparisons:
```{r}
ggplot(data = A3S_events_dataframe, aes(x = Freq, y = A3S_event_ID, fill = A3S_event_ID)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Frequency of A3S Events in All Cells", x = "Count", y = "") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.3), legend.position = "none")

ggplot(data = A3S_events_dataframe, aes(x = Freq, y = Gene, fill = Gene)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Frequency of A3S Events in All Cells", x = "Count", y = "") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.3), legend.position = "none")

ggplot(data = A5S_events_dataframe, aes(x = Freq, y = A5S_event_ID, fill = A5S_event_ID)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Frequency of A5S Events in All Cells", x = "Count", y = "") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.3), legend.position = "none")

ggplot(data = A5S_events_dataframe, aes(x = Freq, y = Gene, fill = Gene)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Frequency of A5S Events in All Cells", x = "Count", y = "") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.3), legend.position = "none")

ggplot(data = cassette_events_dataframe, aes(x = Freq, y = cassette_event_ID, fill = cassette_event_ID)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Frequency of Cassette Events in All Cells", x = "Count", y = "") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.3), legend.position = "none")

ggplot(data = cassette_events_dataframe, aes(x = Freq, y = Gene, fill = Gene)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Frequency of Cassette Events in All Cells", x = "Count", y = "") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.3), legend.position = "none")

ggplot(data = composite_events_dataframe, aes(x = Freq, y = composite_event_ID, fill = composite_event_ID)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Frequency of Composite Events in All Cells", x = "Count", y = "") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.3), legend.position = "none")

ggplot(data = composite_events_dataframe, aes(x = Freq, y = Gene, fill = Gene)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Frequency of Composite Events in All Cells", x = "Count", y = "") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.3), legend.position = "none")

ggplot(data = intron_events_dataframe, aes(x = Freq, y = intron_event_ID, fill = intron_event_ID)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Frequency of Intron Retention Events in All Cells", x = "Count", y = "") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.3), legend.position = "none")

ggplot(data = intron_events_dataframe, aes(x = Freq, y = Gene, fill = Gene)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Frequency of Intron Retention Events in All Cells", x = "Count", y = "") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.3), legend.position = "none")

ggplot(data = MXE_events_dataframe, aes(x = Freq, y = MXE_event_ID, fill = MXE_event_ID)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Frequency of MXE Events in All Cells", x = "Count", y = "") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.3), legend.position = "none")

ggplot(data = MXE_events_dataframe, aes(x = Freq, y = Gene, fill = Gene)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Frequency of MXE Events in All Cells", x = "Count", y = "") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.3), legend.position = "none")
```


# Pie charts:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

for (i in singletypes){
A3S_event_ID_count <- fread(paste0("A3S_event_ID_", i, "_count.csv"))
A5S_event_ID_count <- fread(paste0("A5S_event_ID_", i, "_count.csv"))
cassette_event_ID_count <- fread(paste0("cassette_event_ID_", i, "_count.csv"))
composite_event_ID_count <- fread(paste0("composite_event_ID_", i, "_count.csv"))
intron_event_ID_count <- fread(paste0("intron_event_ID_", i, "_count.csv"))
MXE_event_ID_count <- fread(paste0("MXE_event_ID_", i, "_count.csv"))

values <- c(nrow(A3S_event_ID_count), nrow(A5S_event_ID_count), nrow(cassette_event_ID_count), nrow(composite_event_ID_count), nrow(intron_event_ID_count), nrow(MXE_event_ID_count))

piepercent <- round(100*values/sum(values), 2)
pie_df <- data.frame(piepercent, AStype = c("A3S events", "A5S events", "Cassette exon events", "Composite events", "Intron retention events", "MXE events"))

ggplot(pie_df, aes(x = " ", y = piepercent, fill = c("A3S events", "A5S events", "Cassette exon events", "Composite events", "Intron retention events", "MXE events"))) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0, direction = -1) +
  labs(title = paste0("Percentage of A.S. Events in Comparisons with ", i), fill = "A.S. Type") +
  geom_text(aes(label = paste0(round(piepercent, 2), "%")), position = position_stack(vjust = 0.5), size = 3.5) +
  scale_fill_manual(values = c("cadetblue1", "brown2", "darkolivegreen3", "violet", "gold1", "hotpink4")) +
  theme_void()
}

A3S_event_ID_count <- fread(paste0("JUM_A3SADLvsOLL.csv"))
A5S_event_ID_count <- fread(paste0("JUM_A5SADLvsOLL.csv"))
cassette_event_ID_count <- fread(paste0("JUM_cassetteADLvsOLL.csv"))
composite_event_ID_count <- fread(paste0("JUM_compositeADLvsOLL.csv"))
intron_event_ID_count <- fread(paste0("JUM_intronADLvsOLL.csv"))
MXE_event_ID_count <- fread(paste0("JUM_MXEADLvsOLL.csv"))

values <- c(nrow(A3S_event_ID_count), nrow(A5S_event_ID_count), nrow(cassette_event_ID_count), nrow(composite_event_ID_count), nrow(intron_event_ID_count), nrow(MXE_event_ID_count))

piepercent <- round(100*values/sum(values), 2)
pie_df <- data.frame(piepercent, AStype = c("A3S events", "A5S events", "Cassette exon events", "Composite events", "Intron retention events", "MXE events"))

ggplot(pie_df, aes(x = " ", y = piepercent, fill = c("A3S events", "A5S events", "Cassette exon events", "Composite events", "Intron retention events", "MXE events"))) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0, direction = -1) +
  labs(title = "Percentage of A.S. Event Types in ADL vs OLL", fill = "A.S. Type") +
  geom_text(aes(label = paste0(round(piepercent, 2), "%")), position = position_stack(vjust = 0.5), size = 3.5) +
  scale_fill_manual(values = c("cadetblue1", "brown2", "darkolivegreen3", "violet", "gold1", "hotpink4")) +
  theme_void()

# For all A.S. events combined:
A3S_event_ID_count <- nrow(A3S_events_dataframe[A3S_events_dataframe$Freq > 5,])
A5S_event_ID_count <- nrow(A5S_events_dataframe[A5S_events_dataframe$Freq > 5,])
cassette_event_ID_count <- nrow(cassette_events_dataframe[cassette_events_dataframe$Freq > 5,])
composite_event_ID_count <- nrow(composite_events_dataframe[composite_events_dataframe$Freq > 5,])
intron_event_ID_count <- nrow(intron_events_dataframe[intron_events_dataframe$Freq > 5,])
MXE_event_ID_count <- nrow(MXE_events_dataframe[MXE_events_dataframe$Freq > 5,])

values <- c(A3S_event_ID_count, A5S_event_ID_count, cassette_event_ID_count, composite_event_ID_count, intron_event_ID_count, MXE_event_ID_count)

piepercent <- round(100*values/sum(values), 2)
pie_df <- data.frame(piepercent, AStype = c("A3S events", "A5S events", "Cassette exon events", "Composite events", "Intron retention events", "MXE events"))

colors <- brewer.pal(n = 6, name = "Pastel1")

ggplot(pie_df, aes(x = " ", y = piepercent, fill = c("A3S events", "A5S events", "Cassette exon events", "Composite events", "Intron retention events", "MXE events"))) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0, direction = -1) +
  labs(title = bquote(bold("Percentage of A.S. Event Types (frequency > 5) in All Comparisons Combined")), fill = "A.S. Type", ) +
  geom_text(aes(label = paste0(round(piepercent, 2), "%")), position = position_stack(vjust = 0.55), size = 4.8) +
  scale_fill_manual(values = colors) +
  theme_void()
```

Copyright 2024 The Regents of the University of California

All Rights Reserved

Created by Zachery Wolfe

Department of Biochemistry

This file is part of Differential Expression in C. elegans. \
Differential Expression in C. elegans is a free pipeline: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. \
Differential Expression in C. elegans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. \
You should have received a copy of the GNU General Public License along with Differential Expression in C. elegans. If not, see <https://www.gnu.org/licenses/>.
