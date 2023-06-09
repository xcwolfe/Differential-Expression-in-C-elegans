---
title: "Creating new JUM txts"
output: html_document
---

#SIXTH

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
library(tidyr)
library(dplyr)
```

# Load your cell types and save them as a character vector:
```{r}
singletypes <- c("ADL","AFD","AIM","AIN","AIY","ASEL","ASER","ASG","ASI","ASK","AVA","AVE","AVG","AVH","AVK","AVL","AVM","AWA","AWB","AWC","BAG","CAN","DA","DD","DVC","I5","IL1","IL2","LUA","NSM","OLL","OLQ","PHA","PVC","PVD","PVM","RIA","RIC","RIM","RIS","RMD","SMB","SMD","VB","VC","VD")

# Create a data frame with your cell types and number of replicates per cell type (used to print --Condition1_fileNum_threshold later):

type = as.factor(c(rep("ADL",4),rep("AFD",5),rep("AIM",4),rep("AIN",6),rep("AIY",3),rep("ASEL",3),rep("ASER",4),rep("ASG",4),rep("ASI",4),rep("ASK",4),rep("AVA",6),rep("AVE",3),rep("AVG",3),rep("AVH",4),rep("AVK",4),rep("AVL",3),rep("AVM",3),rep("AWA",4),rep("AWB",5),rep("AWC",4),rep("BAG",4),rep("CAN",3),rep("DA",4),rep("DD",4),rep("DVC",4),rep("I5",4),rep("IL1",3),rep("IL2",4),rep("LUA",4),rep("NSM",3),rep("OLL",2),rep("OLQ",3),rep("PHA",4),rep("PVC",5),rep("PVD",2),rep("PVM",2),rep("RIA",6),rep("RIC",4),rep("RIM",4),rep("RIS",3),rep("RMD",5),rep("SMB",5),rep("SMD",4),rep("VB",4),rep("VC",6),rep("VD",4)))

typescounts <- c(4,5,4,6,3,3,4,4,4,4,6,3,3,4,4,3,3,4,5,4,4,3,4,4,4,4,3,4,4,3,2,3,4,5,2,2,6,4,4,3,5,5,4,4,6,4)
combinedtype = as.data.frame(list(Type = singletypes, Count = as.numeric(typescounts-1), TrueCount = typescounts))
rownames(combinedtype) = combinedtype$Type
```

# .txt files that contain all runs for each step of JUM:
```{r}
setwd("D:/Zach Wolfe's DESeq analysis")

for(i in singletypes){
  for(j in singletypes){
    if(i != j){   
     sink(file=paste0("JUM_A_runs", i, "vs", j, ".txt", sep=""))
      cat(paste("sbatch JUM_A_", i, "vs", j, ".txt ", sep=""))
      sink()}
    if (i == j) next
  }
}

for(i in singletypes){
  for(j in singletypes){
    if(i != j){   
      
    sink(file=paste0("JUM_B_runs", i, "vs", j, ".txt", sep=""))
      cat(paste("sbatch JUM_B_", i, "vs", j, ".txt ", sep=""))
      sink()}
    if (i == j) next
  }
}

for(i in singletypes){
  for(j in singletypes){
    if(i != j){   
      
     sink(file=paste("JUM_C_runs", i, "vs", j, ".txt", sep=""))
      cat(paste("sbatch JUM_C_", i, "vs", j, ".txt ", sep=""))
      sink()}
    if (i == j) next
  }
}
```

# JUM_A:
```{r}
setwd("D:/Zach Wolfe's DESeq analysis")

for(i in singletypes){
  for(j in singletypes){
    if(i != j){   
      
sink(file=paste("JUM_A_", i, "vs", j, ".txt", sep=""))
      cat(paste("#!/bin/bash",
"",
"",
"#SBATCH -N 1",
"#SBATCH -o output.out",
"#SBATCH -e output.err",
"#SBATCH -p highmem",
"",
"module load gcc",
"module load bedtools2",
"module load samtools",
"module load perl",
"",
"", sep = "\n"),
paste0(" /scratch/group/norrislab/JUM/JUM_2.02/JUM_A.sh --Folder /scratch/group/norrislab/JUM/JUM_2.02 --JuncThreshold 5 --Condition1_fileNum_threshold ", combinedtype[i,2], " --Condition2_fileNum_threshold ", combinedtype[j,2], " --IRthreshold 5 --Readlength 100 --Thread ", combinedtype[i,3], " --Condition1SampleName ", i, "1,", i, "2,", i, "3 --Condition2SampleName ", j, "1,", j, "2,", j, "3"))
sink()
    }
    if (i == j) next
  }
}
```

# Run R script step:
```{r}
setwd("D:/Zach Wolfe's DESeq analysis")

for(i in singletypes){
  for(j in singletypes){
    if(i != j){   
x <- as.character(c(paste0(i, 1:3, "   ", i), paste0(j, 1:3, "   ", j)))

# experiment_design file:
sink(file=paste("experimental_design", i, "vs", j, ".txt", sep = ""))
     cat(paste("condition",
               x[1], x[2], x[3], x[4], x[5], x[6], sep = "\n"))
     sink()
    }
    if (i == j) next
  }
}
```

```{r}
setwd("D:/Zach Wolfe's DESeq analysis")

for(i in singletypes){
  for(j in singletypes){
    if(i != j){   

cat(file = paste0("Rscript", i, "vs", j, ".txt"),
      paste0("Rscript /scratch/group/norrislab/JUM/JUM_2.02/R_script_JUM.R experimental_design", i, "vs", j, ".txt >outputFile.Rout 2> errorFile.Rout"))
    }
    if (i == j) next}
}
```

# JUM_B:
```{r}
setwd("D:/Zach Wolfe's DESeq analysis")

for(i in singletypes){
  for(j in singletypes){
    if(i != j){   
      
sink(file=paste("JUM_B_", i, "vs", j, ".txt", sep=""))
      cat(paste("#!/bin/bash",
"",
"",
"#SBATCH -N 1",
"#SBATCH -o output.out",
"#SBATCH -e output.err",
"#SBATCH -p standard-s",
"",
"module load gcc",
"module load bedtools2",
"module load samtools",
"module load perl",
"",
"", sep = "\n"),
paste0(" /scratch/group/norrislab/JUM/JUM_2.02/JUM_B.sh --Folder /scratch/group/norrislab/JUM/JUM_2.02 --Test pvalue --Cutoff 1 --TotalFileNum ", combinedtype[i,3]+combinedtype[j,3], " --Condition1_fileNum_threshold ", combinedtype[i,2], " --Condition2_fileNum_threshold ", combinedtype[j,2], " --Condition1SampleName ", i, "1,", i, "2,", i, "3 --Condition2SampleName ", j, "1,", j, "2,", j, "3")) # Change pvalue cutoff as needed
sink()
    }
    if (i == j) next
  }
}
```

# JUM_C:
```{r}
setwd("D:/Zach Wolfe's DESeq analysis")

for(i in singletypes){
  for(j in singletypes){
    if(i != j){   
      
sink(file=paste("JUM_C_", i, "vs", j, ".txt", sep=""))
      cat(paste("#!/bin/bash",
"",
"",
"#SBATCH -N 1",
"#SBATCH -o output.out",
"#SBATCH -e output.err",
"#SBATCH -p standard-s",
"",
"module load gcc",
"module load bedtools2",
"module load samtools",
"module load perl",
"",
"", sep = "\n"),
paste(" /scratch/group/norrislab/JUM/JUM_2.02/JUM_C.sh --Folder /scratch/group/norrislab/JUM/JUM_2.02 --Test pvalue --Cutoff 1 --TotalCondition1FileNum", combinedtype[i,3], "--TotalCondition2FileNum",   combinedtype[j,3], " --REF refFlat.txt"))
sink()
    }
    if (i == j) next
  }
}
```
