#ELEVENTH

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
library("pheatmap")
library("RColorBrewer")
library("tidyverse")
library("data.table")
library("purrr")
library("foreach")
library("parallel")
```

# List your cell types/conditions:
```{r}
singletypes <- c("ADL","AFD","AIM","AIN","AIY","ASEL","ASER","ASG","ASI","ASK","AVA","AVE","AVG","AVH","AVK","AVL","AVM","AWA","AWB","AWC","BAG","CAN","DA","DD","DVC","I5","IL1","IL2","LUA","NSM","OLL","OLQ","PHA","PVC","PVD","PVM","RIA","RIC","RIM","RIS","RMD","SMB","SMD","VB","VC","VD")

singletypes_withreps <- c("ADL1", "ADL2", "ADL3", "ADL4", "AFD1", "AFD2", "AFD3", "AFD4", "AFD5", "AIM1", "AIM2", "AIM3", "AIM4", "AIN1", "AIN2", "AIN3", "AIN4", "AIN5", "AIN6", "AIY1", "AIY2", "AIY3", "ASEL1", "ASEL2", "ASEL3", "ASER1", "ASER2", "ASER3", "ASER4", "ASG1", "ASG2", "ASG3", "ASG4", "ASI1", "ASI2", "ASI3", "ASI4", "ASK1", "ASK2", "ASK3", "ASK4", "AVA1", "AVA2", "AVA3", "AVA4", "AVA5", "AVA6", "AVE1", "AVE2", "AVE3", "AVG1", "AVG2", "AVG3", "AVH1", "AVH2", "AVH3", "AVH4", "AVK1", "AVK2", "AVK3", "AVK4", "AVL1", "AVL2", "AVL3", "AVM1", "AVM2", "AVM3", "AWA1", "AWA2", "AWA3", "AWA4", "AWB1", "AWB2", "AWB3", "AWB4", "AWB5", "AWC1", "AWC2", "AWC3", "AWC4", "BAG1", "BAG2", "BAG3", "BAG4", "CAN1", "CAN2", "CAN3", "DA1", "DA2", "DA3", "DA4", "DD1", "DD2", "DD3", "DD4", "DVC1", "DVC2", "DVC3", "DVC4", "I51", "I52", "I53", "I54", "IL11", "IL12", "IL13", "IL21", "IL22", "IL23", "IL24", "LUA1", "LUA2", "LUA3", "LUA4", "NSM1", "NSM2", "NSM3", "OLL1", "OLL2", "OLQ1", "OLQ2", "OLQ3", "PHA1", "PHA2", "PHA3", "PHA4", "PVC1", "PVC2", "PVC3", "PVC4", "PVC5", "PVD1", "PVD2", "PVM1", "PVM2", "RIA1", "RIA2", "RIA3", "RIA4", "RIA5", "RIA6", "RIC1", "RIC2", "RIC3", "RIC4", "RIM1", "RIM2", "RIM3", "RIM4", "RIS1", "RIS2", "RIS3", "RMD1", "RMD2", "RMD3", "RMD4", "RMD5", "SMB1", "SMB2", "SMB3", "SMB4", "SMB5", "SMD1", "SMD2", "SMD3", "SMD4", "VB1", "VB2", "VB3", "VB4", "VC1", "VC2", "VC3", "VC4", "VC5", "VC6", "VD1", "VD2", "VD3", "VD4")
```

# Read JUM excel files created in 7. Cell type vs cell type heatmap script for JUM analysis.Rmd and create .bed file which will be used for MEME:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

for (i in singletypes){
  for (j in c("A3S", "A5S", "cassette", "intron")){
  df <- deltaPSI_sums_list[[paste0(j, "_event_ID_", i, "_count.csv")]]
  df <- df[order(-df[,4]),]
  df <- df[df$`sum of all first deltaPSIs for AS_event_ID` > 5,]   ## change sum of all first deltaPSIs threshold here
  if (nrow(df) > 0){
  result <- lapply(df[[paste0(j, "_event_ID")]], function(x) unlist(strsplit(x, split = "[_]")))
  result <- lapply(result, function(x) x[x != ""])
  
# The first three required BED fields are:

# chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random)
# chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
# chromEnd - The ending position of the feature in the chromosome or scaffold.
  
# The 9 additional optional BED fields are:
 
# name - Defines the name of the BED line. The MEME Suite tools ignore the contents of this field.
# score - A score between 0 and 1000.
# strand - Defines the strand. Either "." (=no strand) or "+" or "-".
# thickStart - The MEME Suite tools ignore the contents of this field.
# thickEnd - The MEME Suite tools ignore the contents of this field.
# itemRgb - The MEME Suite tools ignore the contents of this field.
# blockCount - The MEME Suite tools ignore the contents of this field.
# blockSizes - The MEME Suite tools ignore the contents of this field.
# blockStarts - The MEME Suite tools ignore the contents of this field.
  
  for (r in 1:nrow(df)){
  df$chrom[r] <- result[[r]][1]
  df$chromStart[r] <- result[[r]][3]
  df$chromEnd[r] <- result[[r]][length(result[[1]])]
  df$strand[r] <- result[[r]][2]
    }

  df <- df[,-c(1:4)]
  df$new_col <- c("")
  df$new_col2 <- c("")
  df <- df[, c(1,2,3,6,5,4)]
  
  # MEME will not accept sequence lengths longer than 100,000 nucletotides. Let's eliminate them:
  df$chromEnd <- as.numeric(df$chromEnd)
  df$chromStart <- as.numeric(df$chromStart)
  df <- subset(df, chromEnd - chromStart <= 99999)
  # We should also eliminate MtDNA since MEME does not incorporate MtDNA into the genome(s) it uses:
  df <- subset(df, chrom != "ChrMtDNA")
  
  colnames(df) <- c("chrom", "chromStart", "chromEnd", "", "", "Strand")
  }
  if (nrow(df) == 0) next
  write.table(df, file = paste0(i, "_", j, "_MEME.tsv"), sep = "\t", quote = F, row.names = F, col.names = F)
  }
}
```

# Create BACKGROUND .bed file which will be used for MEME:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

for (i in singletypes){
  for (j in c("A3S", "A5S", "cassette", "intron")){
  df <- deltaPSI_sums_list[[paste0(j, "_event_ID_", i, "_count.csv")]]
  df <- df[order(-df[,4]),]
  df <- df[df$`sum of all first deltaPSIs for AS_event_ID` < 1,]   ## change sum of all first deltaPSIs threshold here - should be non-up and non-downregulated events for background
  df <- df[df$`sum of all first deltaPSIs for AS_event_ID` > -1,]   ## change sum of all first deltaPSIs threshold here - should be non-up and non-downregulated events for background
  if (nrow(df) > 0){
  result <- lapply(df[[paste0(j, "_event_ID")]], function(x) unlist(strsplit(x, split = "[_]")))
  result <- lapply(result, function(x) x[x != ""])
  
# The first three required BED fields are:

# chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random)
# chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
# chromEnd - The ending position of the feature in the chromosome or scaffold.
  
# The 9 additional optional BED fields are:
 
# name - Defines the name of the BED line. The MEME Suite tools ignore the contents of this field.
# score - A score between 0 and 1000.
# strand - Defines the strand. Either "." (=no strand) or "+" or "-".
# thickStart - The MEME Suite tools ignore the contents of this field.
# thickEnd - The MEME Suite tools ignore the contents of this field.
# itemRgb - The MEME Suite tools ignore the contents of this field.
# blockCount - The MEME Suite tools ignore the contents of this field.
# blockSizes - The MEME Suite tools ignore the contents of this field.
# blockStarts - The MEME Suite tools ignore the contents of this field.
  
  for (r in 1:nrow(df)){
  df$chrom[r] <- result[[r]][1]
  df$chromStart[r] <- result[[r]][3]
  df$chromEnd[r] <- result[[r]][length(result[[1]])]
  df$strand[r] <- result[[r]][2]
    }

  df <- df[,-c(1:4)]
  df$new_col <- c("")
  df$new_col2 <- c("")
  df <- df[, c(1,2,3,6,5,4)]

  
  # MEME will not accept sequence lengths longer than 100,000 nucletotides. Let's eliminate them:
  df$chromEnd <- as.numeric(df$chromEnd)
  df$chromStart <- as.numeric(df$chromStart)
  df <- subset(df, chromEnd - chromStart <= 99999)
  # We should also eliminate MtDNA since MEME does not incorporate MtDNA into the genome(s) it uses:
  df <- subset(df, chrom != "ChrMtDNA")
  
  colnames(df) <- c("chrom", "chromStart", "chromEnd", "", "", "Strand")
  }
  if (nrow(df) == 0) next
  write.table(df, file = paste0(i, "_", j, "_MEME_background.tsv"), sep = "\t", quote = F, row.names = F, col.names = F)
  }
}
```

Copyright 2024 The Regents of the University of California

All Rights Reserved

Created by Zachery Wolfe

Department of Biochemistry

This file is part of Differential Expression in C. elegans. \
Differential Expression in C. elegans is a free pipeline: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. \
Differential Expression in C. elegans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. \
You should have received a copy of the GNU General Public License along with Differential Expression in C. elegans. If not, see <https://www.gnu.org/licenses/>.
