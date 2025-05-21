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
setwd("D:/Zach Wolfe's JUM analysis")

for(i in singletypes){
  for(j in singletypes){
    if(i != j){   
      
    type1 <- i
    count1 <- combinedtype$Count[combinedtype$Type == i]
        
    type2 <- j
    count2 <- combinedtype$Count[combinedtype$Type == j]
        
        sample_numbers_condition1 <- paste0(type1, 1:count1)
        sample_numbers_condition2 <- paste0(type2, 1:count2)
      
sink(file=paste("JUM_A_", i, "vs", j, ".txt", sep=""))
      cat(paste("#!/bin/bash",
"",
"",
"#SBATCH -N 1",
"#SBATCH -o output.out",
"#SBATCH -e output.err",
"#SBATCH -p highmem",
"#SBATCH --mem=500G",
"",
"module load gcc",
"module load bedtools2",
"module load samtools",
"module load perl",
"",
"", sep = "\n"),
paste0(" /lustre/work/client/group/norrislab/JUM/JUM_2.02/JUM_A.sh --Folder /lustre/work/client/group/norrislab/JUM/JUM_2.02 --JuncThreshold 5 --Condition1_fileNum_threshold ", combinedtype[i,2], " --Condition2_fileNum_threshold ", combinedtype[j,2], " --IRthreshold 5 --Readlength 100 --Thread ", combinedtype[i,3], paste(" --Condition1SampleName", paste(sample_numbers_condition1, collapse = ","), "--Condition2SampleName", paste(sample_numbers_condition2, collapse = ","), sep = " ")))
sink()
    }
    if (i == j) next
  }
}
```

# Run R script step:
# This is the experiment_design .txt file:
```{r}
setwd("D:/Zach Wolfe's JUM analysis")

for (i in 1:length(singletypes)) {
  type1 <- singletypes[i]
  
  for (j in 1:length(singletypes)) {
    type2 <- singletypes[j]
    
    if (i != j) {
      count1 <- combinedtype$TrueCount[combinedtype$Type == type1]
      count2 <- combinedtype$TrueCount[combinedtype$Type == type2]
      
      table_data <- data.frame(condition = c(paste0(type1, 1:count1), paste0(type2, 1:count2)), type = rep(c(type1, type2), c(count1, count2)))
      table_name <- paste0("experiment_design", type1, "vs", type2)
      colnames(table_data) <- c("", "condition")
      
      #print(table_data)
      assign(table_name, table_data)
      write.table(table_data, file = paste(table_name, ".txt", sep = ""), quote = FALSE, row.names = FALSE)
    }
  }
}
```

# This is the actual R script .txt file:
```{r}
setwd("D:/Zach Wolfe's JUM analysis")

for (i in singletypes) {
  for (j in singletypes) {
    if (i != j) {
      script <- paste0(
        "#!/bin/bash",
        "\n",
        "\n",
        "#SBATCH -N 1",
        "\n",
        "#SBATCH -o output.out",
        "\n",
        "#SBATCH -e output.err",
        "\n",
        "#SBATCH -p standard-s",
        "\n",
        "#SBATCH --mem=50G",
        "\n",
        "\n",
        "",
        "module load gcc",
        "\n",
        "module load R",
        "\n",
        "\n",
        "",
        "",
        "Rscript /lustre/work/client/group/norrislab/JUM/JUM_2.02/R_script_JUM.R experiment_design", i, "vs", j, ".txt >outputFile.Rout 2> errorFile.Rout"
      )
      
      cat(script, file = paste0("Rscript", i, "vs", j, ".txt"))
    }
    if (i == j) next
  }
}
```

# JUM_B:
```{r}
setwd("D:/Zach Wolfe's DESeq analysis")

for (i in singletypes) {
  for (j in singletypes) {
    if (i != j) {   
      
      type1 <- i
      count1 <- combinedtype$TrueCount[combinedtype$Type == i]
      
      type2 <- j
      count2 <- combinedtype$TrueCount[combinedtype$Type == j]
      
      sample_numbers_condition1 <- paste0(type1, 1:count1)
      sample_numbers_condition2 <- paste0(type2, 1:count2)
      
sink(file = paste("JUM_B_", i, "vs", j, ".txt", sep = ""))
cat(
paste(
"#!/bin/bash",
"",
"",
"#SBATCH -N 1",
"#SBATCH -o output.out",
"#SBATCH -e output.err",
"#SBATCH -p standard-s",
"#SBATCH --mem=500G",
 "",
"module load gcc",
"module load bedtools2",
"module load samtools",
"module load perl",
"",
"", sep = "\n"),
paste0(" /lustre/work/client/group/norrislab/JUM/JUM_2.02/JUM_B.sh --Folder /lustre/work/client/group/norrislab/JUM/JUM_2.02 --Test pvalue --Cutoff 1 --TotalFileNum ", combinedtype[i, 3] + combinedtype[j, 3], " --Condition1_fileNum_threshold ", combinedtype[i, 2], " --Condition2_fileNum_threshold ", combinedtype[j, 2], " --Condition1SampleName ", paste(sample_numbers_condition1, collapse = ","), " --Condition2SampleName ", paste(sample_numbers_condition2, collapse = ",")))
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
"#SBATCH --mem=50G",
"",
"module load gcc",
"module load bedtools2",
"module load samtools",
"module load perl",
"",
"", sep = "\n"),
paste(" /lustre/work/client/group/norrislab/JUM/JUM_2.02/JUM_C.sh --Folder /lustre/work/client/group/norrislab/JUM/JUM_2.02 --Test pvalue --Cutoff 1 --TotalCondition1FileNum", combinedtype[i,3], "--TotalCondition2FileNum",   combinedtype[j,3], " --REF refFlat.txt"))
sink()
    }
    if (i == j) next
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
