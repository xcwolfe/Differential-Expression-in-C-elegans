---
title: "Creating new JUM txts"
output: html_document
date: "2022-08-29"
---

#FIFTH

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
library(tidyr)
library(dplyr)
```


```{r}
singletypesprac <-  c("ASG","AVE","AVG","AWA","AWB","DD","PVD","VD")

# Create a data frame with your cell types and number of replicates per cell type (used to print --Condition1_fileNum_threshold later):
oldtypesprac = c("ASG","AVE","AVG","AWA","AWB","DD","PVD","VD")
oldtypespraccounts = c(3,2,2,3,4,2,1,3)
combinedtype = as.data.frame(list(Type = oldtypesprac, Count = oldtypespraccounts, TrueCount = as.numeric(oldtypespraccounts+1)))
rownames(combinedtype) = combinedtype$Type
```

# Unix runs files:
```{r}
setwd("D:/Zach Wolfe's DESeq analysis")

for(i in singletypesprac){
  for(j in singletypesprac){
    if(i != j){   
     sink(file=paste0("JUM_A_runs", i, "vs", j, ".txt", sep=""))
      cat(paste("sbatch JUM_A_", i, "vs", j, ".txt ", sep=""))
      sink()}
    if (i == j) next
  }
}

for(i in singletypesprac){
  for(j in singletypesprac){
    if(i != j){   
      
    sink(file=paste0("JUM_B_runs", i, "vs", j, ".txt", sep=""))
      cat(paste("sbatch JUM_B_", i, "vs", j, ".txt ", sep=""))
      sink()}
    if (i == j) next
  }
}

for(i in singletypesprac){
  for(j in singletypesprac){
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

for(i in singletypesprac){
  for(j in singletypesprac){
    if(i != j){   
      
sink(file=paste("JUM_A_", i, "vs", j, ".txt", sep=""))
      cat(paste("#!/bin/bash",
"",
"",
"#SBATCH -N 1",
"#SBATCH -t 1200",
"#SBATCH -o output.out",
"#SBATCH -e output.err",
"#SBATCH -p high-mem-1",
"",
"module load bedtools2",
"module load samtools",
"module load gcc-6.3",
"module load perl/5.24.1-44zo6ky",
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

for(i in singletypesprac){
  for(j in singletypesprac){
    if(i != j){   
     
# z <- matrix(combinedtype$TrueCount, ncol = 2, nrow = 1)
# colnames(z) <- c(i, j)
# rownames(z) <- c("Replicates")
# imult <- z[1,i]
# jmult <- z[1,j]
# combinedmults <- as.numeric(c(imult, jmult))
x <- as.character(c(paste0(i, 1:3, "   ", i), paste0(j, 1:3, "   ", j)))
# y <- as.character(c(i,i,i,i,j,j,j,j))

# experiment_design file:
sink(file=paste("experimental_design", i, "vs", j, ".txt", sep = ""))
     # cat(cbind(x,y))   
     # cat(t(rbind(x,y))) # Fix this. Get it into a table!
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

for(i in singletypesprac){
  for(j in singletypesprac){
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

for(i in singletypesprac){
  for(j in singletypesprac){
    if(i != j){   
      
sink(file=paste("JUM_B_", i, "vs", j, ".txt", sep=""))
      cat(paste("#!/bin/bash",
"",
"",
"#SBATCH -N 1",
"#SBATCH -t 1000",
"#SBATCH -o output.out",
"#SBATCH -e output.err",
"#SBATCH -p standard-mem-s",
"",
"module load bedtools2",
"module load samtools",
"module load gcc-6.3",
"module load perl/5.24.1-44zo6ky",
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

for(i in singletypesprac){
  for(j in singletypesprac){
    if(i != j){   
      
sink(file=paste("JUM_C_", i, "vs", j, ".txt", sep=""))
      cat(paste("#!/bin/bash",
"",
"",
"#SBATCH -N 1",
"#SBATCH -t 50",
"#SBATCH -o output.out",
"#SBATCH -e output.err",
"#SBATCH -p standard-mem-s",
"",
"module load bedtools2",
"module load samtools",
"module load gcc-6.3",
"module load perl/5.24.1-44zo6ky",
"",
"", sep = "\n"),
paste(" /scratch/group/norrislab/JUM/JUM_2.02/JUM_C.sh --Folder /scratch/group/norrislab/JUM/JUM_2.02 --Test pvalue --Cutoff 1 --TotalCondition1FileNum", combinedtype[i,3], "--TotalCondition2FileNum",   combinedtype[j,3], " --REF refFlat.txt"))
sink()
    }
    if (i == j) next
  }
}
```
