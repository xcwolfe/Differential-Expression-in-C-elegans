# STAR and JUM workflow
This pipeline illustrates how to run STAR and JUM using a large cluster computer (with hopefully lots of memory). We will be utilizing Unix (C language) and later R to analyze our JUM data

Navigate to your working directory:

```
$ cd /work/users/zwolfe
```

Make a folder to put all your fastq files in:

```
$ mkdir srrfiles
$ cd srrfiles/
```

Now we need to import our SRA toolkit:

```
$ export MODULEPATH="${MODULEPATH}:/scratch/group/norrislab/modules/sratoolkit"
$ module load sratoolkit
```

Setup the parameters for the cluster computer:

`$ vdb-config –interactive`

  a.	Enable remote access [X]
  
  b.	Cache tab
  
    i. Enable local file-sharing [X]
    
    ii. choose user-repository: /work/users/zwolfe/srrfiles
    
    iii. Save
    
  c.	AWS tab
  
    i. report cloud instance identity [X]
    
    ii. Save and exit
    
 ```
$ fastq-dump -I --split-files SRRXXXXXX SRRXXXXXX SRRXXXXXX…
$ du -sh *
```

  a. This is to verify that all file sizes are correct. Some may download incorrectly or partially due to timeouts on HPC.
      
Create a new directory that includes your genome and initialize STAR:

Go to http://genome.ucsc.edu/cgi-bin/hgTables and download the latest *C. elegans* genome as a .gtf file:

```
$ mkdir genome
$ cd genome/
$ wget ftp://ftp.ensembl.org/pub/release-107/fasta/caenorhabditis_elegans/dna/*.fa.gz
$ gunzip *.gz
```

Now we can initialize STAR to generate Genome and SA output files:
    
```
$ dos2unix initialize_star.txt
$ sbatch initialize_star.txt
```

Now, let's actually run STAR on each file/cell type:

```
$ export MODULEPATH="${MODULEPATH}:/scratch/group/norrislab/modules/sratoolkit"
$ module load star
$ mkdir XXX
$ mkdir rep1 rep2 rep3
$ cp /scratch/group/norrislab/Zach/run_STAR_chimeras.txt
```

run_STAR_chimeras.txt:
```
#!/bin/bash


#SBATCH -N 1
#SBATCH -t 1400
#SBATCH -o output.out
#SBATCH -e output.err
#SBATCH -p standard-mem-s

module load star

STAR --genomeDir /scratch/group/norrislab/genome --readFilesIn PAN_SRR13995337*.fastq --outSAMstrandField intronMotif --outReadsUnmapped Fastx --runThreadN 20 --chimSegmentMin 15
```

  a. Add parameters --chimOutType Junctions SeparateSAMold
  
  b. Make sure to change this file every time to reflect the names of the neuron type (ASG, etc.)
  
```      
$ cd rep1/
$ cp /work/users/zwolfe/srrfiles/SRRXXXXXXXX*.fastq
$ mv SRR13995310_1.fastq XXX_SRR13995310_1.fastq
```

  a. XXX = neuron type
  
```
$ dos2unix ../run_STAR_chiemras.txt
$ sbatch ../run_STAR_chimeras.txt (once you are in the rep/ folder)
```

Run samtools:
     
```     
$ cd /work/users/zwolfe/srrfiles/XXX/rep1
```

  a. add run_samtools.txt (from Olivia) to each folder:
  
```
#!/bin/bash

#SBATCH -N 2
#SBATCH -t 180
#SBATCH -o output.out
#SBATCH -e output.err
#SBATCH -p standard-mem-s

module load samtools

samtools index Aligned.sortedByCoord.out.bam
```

```
$ dos2unix run_samtools.txt
$ sbatch run_samtools.txt 
```

(Do this in each folder)

RENAME Aligned.out.sam, Aligned.out_sorted.bam, and SJ.out.tab files:

  a. ASG1Aligned.out.sam
  
Put these renamed files, ALONG WITH all JUM and experimental_design files, in their own JUM_analysis folder 

  a. JUM_analysis_ASGvsAVE
  

Change run_JUM_A to reflect appropriate neuron type:

  a. ASG1, AVE1, etc.
  
    i. Make sure reps are separated by a comma and do not contain spaces.
    
  b. repeat for EVERY comparison
  
```
$ dos2unix run_JUM_A.txt
$ sbatch run_JUM_A.txt
$ cd /JUM_diff
```

Create experimental_design.txt:
`  condition`

`ctrl1 control`

`ctrl2 control`

`ctrl3 control`

`treat1 treatment`

`treat2 treatment`

`treat3 treatment`

Load r module:

```
$ module load rstudio`
$ Rscript /scratch/group/norrislab/JUM/JUM_2.02/R_script_JUM.R experimental_design.txt >outputFile.Rout 2> errorFile.Rout
```

  a.	Make sure the R script outputs a file called AS_differential.txt
  
Change run_JUM_B to reflect appropriate neuron type:

  a. ASG1, AVE1, etc.
  
  b. -- TotalFileNum (depends on the total number of files you have between your two cell types)
  
  c. -- Condition1_FileNum_threshold/TotalCondition2_FileNum_threshold (depends on the total number of files you have between your two cell types; should be N - 1)
  
  d. -- Cutoff to your desired p-value (default is 1 to see all splicing events)
  
  e. repeat for EVERY comparison
  
  f. NOTE: If you want to do multiple p-value comparisons, make sure to delete the previous JUM_B runs or you will create a downward spiral of directories
  
```
$ dos2unix run_JUM_B.txt
$ sbatch run_JUM_B.txt
```

Download and moved refFlat.txt and run_JUM_C.txt to FINAL_OUTPUT_pvalue_1 directory 

```
$ cd FINAL_OUTPUT_pvalue_1
```

Change run_JUM_C to reflect appropriate neuron type:

  a. ASG1, AVE1, etc.
  
  b. -- TotalFileNum (depends on the total number of files you have between your two cell types)
  
  c. -- TotalCondition1_FileNum/TotalCondition2_FileNum (depends on the total number of files you have between your two cell types)
  
  d. repeat for EVERY comparison
  
```
$ dos2unix run_JUM_C.txt
$ sbatch run_JUM_C.txt
$ cd /work/users/zwolfe/srrfiles/JUM_analysis/JUM_diff/FINAL_JUM_OUTPUT_pvalue_1/
```

  a. This is where the detailed and simplified summaries of each AS event are
  
When complete with the previous step, it is wise to save and export this unique FINAL_JUM_OUTPUT folder because it will be replaced when you run JUM through again.
I have simplified this process by generating the necessary .txt files in R. See JUM file automation.Rmd

  a. I put the entire output (all .txt files) in each of the following folders for smooth runs: Jum_analysis, JUM_diff, temp_JUM_C_run


Download JUM results from the cluster:

  a. WinSCP to D:/Zach Wolfe’s JUM analysis
  
  b. Make a NEW DIRECTORY for each JUM comparison. For example: D:\Zach Wolfe's JUM analysis\FINAL_JUM_OUTPUT_pvalue_1ASGvsAVE

Open the downloaded results in R/Rstudio. See “Heatmap Script for JUM analysis.Rmd”

