Before you begin… 

How many samples do you have? If you have more than 2 experimental samples, you will save yourself quite a bit of time by batch running many looped JUM and DCC scripts. See “JUM and DCC megaloops” for more information. 

 

STAR on all neurons using all 180 SRA files: 

Navigate to your work/group/norrislab directory which you plan to download all SRA files into (we need scratch space as this will take lots of memory): 

```
cd /work/group/norrislab/Zach/ 
```
 

Create a file called sra_accessions.txt which is a list of all your SRR numbers for your designated project and move/drag it to your scratch space: 
```
SRR24086884 
SRR24086885 
SRR24086886 
SRR24086887 
…
```

 

Convert this file from dos format to Unix format. NOTE: You will have to do this whenever you create a file on your personal computer that you plan to run on ManeFrame: 

```
module load gcc dos2unix 
dos2unix sra_accessions.txt
```

Load necessary modules for download: 
```
module load gcc 
module load spack 
source <(spack module tcl loads --dependencies entrezdirect@10.7.20190114) 
source <(spack module tcl loads --dependencies sratoolkiit@2.10.9-gcc-9.2.0-qiw3o6v) 
source <(spack module tcl loads --dependencies parallel@20210922) 
module load sratoolkit-2.10.9-gcc-9.2.0-qiw3o6v
```

Now use a loop of the fastq-dump command to download the .fastq files listed in your project (dos2unix and sbatch this script): 

```
#!/bin/bash 
#SBATCH -N 1 
#SBATCH -o output.out 
#SBATCH -e output.err 
#SBATCH -p highmem  

module load spack 
module load sratoolkit/3.0.0-uicqmaw 

while read -r line; do 

  fastq-dump --split-files "$line" 

done < sra_accessions.txt 
```
 

Perform STAR on the 2 biggest .fastq files for each SRR number (dos2unix and sbatch this script): 

WARNING: Sometimes we have more than 2 .fastq files present on SRA. Only the two biggest .fastq files are important to us when we run STAR. Make sure you modify this script to select only the 2 biggest .fastq files for each replicate from your .fastq dump: 

 
```
#!/bin/bash 
#SBATCH -N 1 
#SBATCH -o output.out 
#SBATCH -e output.err 
#SBATCH -p highmem 

module load gcc star 

STAR --runThreadN 10 --genomeDir /work/users/zwolfe/genome_109/ --outSAMtype BAM SortedByCoordinate --readFilesIn *_1.fastq *_2.fastq --outSAMstrandField intronMotif --readFilesCommand cat --outReadsUnmapped Fastx --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMatchNmin 1 --outFilterMismatchNmax 2 --chimSegmentMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15 --alignTranscriptsPerReadNmax 200000 --alignTranscriptsPerWindowNmax 20000 --outFileNamePrefix ADL1 
```
 

Once STAR is done running for one replicate, make sure you rename the output files (chimeric.out.sam, chimeric.out.junction, etc) to reflect which replicate they correspond to (ADL1chimeric.out.sam, ADL1chimeric.out.junction, etc). If you don’t do this but stay in the same working directory, STAR will overwrite all of your files! 

- OR use the --outFileNamePrefix flag above 


Run samtools on the .bam output files (dos2unix and sbatch this script): 
 
```
#!/bin/bash 
#SBATCH -N 2 
#SBATCH -o output.out 
#SBATCH -e output.err 
#SBATCH -p standard-s

module load gcc 
module load samtools 

samtools index NSM2Aligned.sortedByCoord.out.bam  
```

JUM will require us to create an Aligned.out_sorted.bam file whereas DCC prefers Aligned.sortedByCoord.out.bam, so we will do both (dos2unix and sbatch this script – could also be combined with the previous samtools script): 

```
#!/bin/bash 
#SBATCH -N 1 
#SBATCH -o output.out 
#SBATCH -e output.err 
#SBATCH -p standard-s 

module load gcc samtools  

# List of unique prefixes for the BAM files 
prefixes=("ADL1" "ADL2" …"VD4") 

# Make files with the given prefixes 
for prefix in "${prefixes[@]}"; do 
files=("${prefix}Aligned.sortedByCoord.out.bam") 
if [[ ${#files[@]} -gt 0 ]]; then 
samtools view -b "${files[0]}" | sed "s/sortedByCoord.out/${prefix}_out_sorted/" | samtools reheader - "${files[0]}" > "${prefix}Aligned.out_sorted.bam" 
fi 
done 
```
 
We also need an aligned.out.sam file for each replicate. I am going to create a loop to make an aligned.out.sam file for each Aligned.sortedByCoord.out.bam file corresponding to each replicate (dos2unix and sbatch this script – could also be combined with the previous samtools script): 

```
#!/bin/bash 
#SBATCH -N 2 
#SBATCH -o output.out 
#SBATCH -e output.err 
#SBATCH -p standard-s 

module load gcc 
module load samtools 

# List of unique prefixes for the BAM files 
prefixes=("ADL1" "ADL2" …"VD4") 

# Loop through the prefixes and convert BAM to SAM 
for prefix in "${prefixes[@]}" 
do 
    bam_file="${prefix}Aligned.sortedByCoord.out.bam" 
    sam_file="${prefix}Aligned.out.sam" 
    samtools view -h -o "$sam_file" "$bam_file" 
done 
```

We now have all the files we need to perform DESeq and JUM. 

Once all of these files have been loaded into your scratch space, move each STAR output file (corresponding to one replicate) into its own folder; I have decided to execute this using a loop (dos2unix and sbatch this script): 

```
#!/bin/bash 
#SBATCH -N 1 
#SBATCH -o output.out 
#SBATCH -e output.err 
#SBATCH -p standard-s 

source_dir="/lustre/work/client/group/norrislab/Zach" 
destination_dir="/lustre/work/client/group/norrislab/Zach/STAR_results" 
prefixes=("ADL" "AFD" "AIM" "AIN" "AIY" "ASEL" "ASER" "ASG" "ASI" "ASK" "AVA" "AVE" "AVG" "AVH" "AVK" "AVL" "AVM" "AWA" "AWB" "AWC" "BAG" "CAN" "DA" "DD" "DVC" "I5" "IL1" "IL2" "LUA" "NSM" "OLL" "OLQ" "PHA" "PVC" "PVD" "PVM" "RIA" "RIC" "RIM" "RIS" "RMD" "SMB" "SMD" "VB" "VC" "VD") 

# Move files with the given prefixes 
for prefix1 in "${prefixes[@]}"; do 
  for prefix2 in "${prefixes[@]}"; do 
    if [[ "$prefix1" != "$prefix2" ]]; then 
      sub_dir="${prefix1}vs${prefix2}" 
      mkdir -p "$destination_dir/$sub_dir" 
      cp "$source_dir/$prefix1"* "$destination_dir/$sub_dir/" 
      cp "$source_dir/$prefix2"* "$destination_dir/$sub_dir/" 
    fi 
  done 
done  
```

We also need to run HTSeq to generate counts from our STAR outputs. We can start by setting up a conda environment and pip installing HTSeq (if HTSeq is already installed, you can ignore this section): 

```
module load conda 
conda create -n htseq_env -c conda-forge python=3.9 
pip install HTSeq 
```

Here is my command for running HTSeq (dos2unix and sbatch this script): 

```
#!/bin/bash 
#SBATCH -N 1 
#SBATCH -o A.HTSeq.output.out 
#SBATCH -e output.err 
#SBATCH -p standard-s 
 
module load gcc 
module load python/3.10.8-dkpz5k5 
 
htseq-count --stranded=no -f bam /lustre/work/client/group/norrislab/Zach/ADL1Aligned.sortedByCoord.out.bam /lustre/work/client/users/zwolfe/genome/Caenorhabditis_elegans.WBcel235.107.gtf > ADL1_counts.txt
```

We now need to combine our HTSeq results to create a count matrix. We can do this in R; see CenGen bulk RNAseq data PCA and DESeq.Rmd. Make sure all of your count files are in the same directory before running the R script(s). 


JUM on all neurons using all STAR and samtools files: 

Once your STAR output files have all been moved to the appropriate destination folder, we will be ready to run JUM. Let’s start with JUM_A. These JUM.txt scripts have all been created using JUM file automation.Rmd (drag each JUM_A.txt into the appropriate comparison folder and dos2unix and sbatch this script): 

```
#!/bin/bash 
#SBATCH -N 1 
#SBATCH -o output.out 
#SBATCH -e output.err 
#SBATCH -p highmem 

module load gcc 
module load bedtools2 
module load samtools 
module load perl 

  /lustre/work/client/group/norrislab/JUM/JUM_2.02/JUM_A.sh --Folder /lustre/work/client/group/norrislab/JUM/JUM_2.02 --JuncThreshold 5 --Condition1_fileNum_threshold 3 --Condition2_fileNum_threshold 4 --IRthreshold 5 --Readlength 100 --Thread 4 --Condition1SampleName ADL1,ADL2,ADL3,ADL4 --Condition2SampleName AFD1,AFD2,AFD3,AFD4,AFD5
```
 
Before running the Rscript step, you must download the latest version of BiocManager that is compatible with R. I have changed the Rscript step on my own and moved it to /lustre/work/client/group/norrislab/Zach like so (dos2unix this script): 

```
args = commandArgs(trailingOnly=TRUE) 

# test if experiment design input file is provided; if not, return an error 

if (length(args)==0) { 
  stop("An experiment design file must be supplied (input file).n", call.=FALSE) 
} 


#dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE) 

if (!require(DEXSeq)) biocManager("DEXSeq", dependencies=TRUE) 

#if (!require(DEXSeq)) biocManager("DEXSeq", suppressUpdates = TRUE, lib = Sys.getenv("R_LIBS_USER"), dependencies=TRUE) 

library(DEXSeq) 

annotationfile = file.path("combined_AS_JUM.gff") 
sampledesign = read.table(args[1], header=TRUE) 
sampledesign 
fullFilenames <- dir(pattern ="combined_count.txt") 
fullFilenames  

JUM = DEXSeqDataSetFromHTSeq( 
fullFilenames, 
sampleData=sampledesign, 
design= ~ sample + exon + condition:exon, 
flattenedfile=annotationfile) 

JUM = estimateSizeFactors(JUM) 
JUM = estimateDispersions(JUM) 
JUM = testForDEU(JUM) 
JUM = estimateExonFoldChanges(JUM, fitExpToVar="condition") 
dxr1 = DEXSeqResults(JUM) 
dxr1_sub <- dxr1[,1:12] 
transcripts <- as.character(dxr1$transcripts) 
write.table(cbind(dxr1_sub, transcripts), "AS_differential.txt", sep="\t", quote=F) 
q() 
```

Change your directory to JUM_diff (cd JUM_diff). Drag each experiment_design.txt and Rscript.txt into the appropriate comparison folder and dos2unix and sbatch this script:  

```
#!/bin/bash 
#SBATCH -N 1 
#SBATCH -o output.out 
#SBATCH -e output.err 
#SBATCH -p standard-s 
#SBATCH --mem=50G 

module load gcc 
module load R 

Rscript /lustre/work/client/group/norrislab/Zach/R_script_JUM.txt experiment_designADLvsAFD.txt >outputFile.Rout 2> errorFile.Rout 
```

Before moving on to JUM_B, you may have to install the Descriptive perl package. To do this simply download the perl module and run the function below: 
 
```
module load gcc perl 
cpan -i Statistics::Descriptive 
cpan -i List::MoreUtils 
```

Now that the Rscript step is done and the descriptive package is loaded, we can run JUM_B (dos2unix and sbatch this script): 

```
#!/bin/bash 
#SBATCH -N 1 
#SBATCH -o output.out 
#SBATCH -e output.err 
#SBATCH -p standard-s 
#SBATCH --mem=500G  

module load gcc 
module load bedtools2 
module load samtools 
module load perl 
 
  /lustre/work/client/group/norrislab/JUM/JUM_2.02/JUM_B.sh --Folder /lustre/work/client/group/norrislab/JUM/JUM_2.02 --Test pvalue --Cutoff 1 --TotalFileNum 9 --Condition1_fileNum_threshold 3 --Condition2_fileNum_threshold 4 --Condition1SampleName ADL1,ADL2,ADL3,ADL4 --Condition2SampleName AFD1,AFD2,AFD3,AFD4,AFD5 
```

Finally, after JUM_B is done, change directories to the FINAL_JUM_OUTPUT_pvalue_1 directory (cd FINAL_JUM_OUTPUT_pvalue_1). Drag in refFlat.txt and run JUM_C (dos2unix and sbatch this script): 

```
#!/bin/bash 
#SBATCH -N 1 
#SBATCH -o output.out 
#SBATCH -e output.err 
#SBATCH -p standard-s 
#SBATCH --mem=50G 

module load gcc 
module load bedtools2 
module load samtools 
module load perl 


  /lustre/work/client/group/norrislab/JUM/JUM_2.02/JUM_C.sh --Folder /lustre/work/client/group/norrislab/JUM/JUM_2.02 --Test pvalue --Cutoff 1 --TotalCondition1FileNum 4 --TotalCondition2FileNum 5  --REF refFlat.txt 
 ```

You may find it worthwhile to place all FINAL_JUM_OUTPUT_pvalue_1 directories in a single overarching directory. I have written this (optional) loop if you want to organize your files like me: 

```
#!/bin/bash 
#SBATCH -N 1 
#SBATCH -o output.out 
#SBATCH -e output.err 
#SBATCH -p standard-s 


base_path="/lustre/work/client/group/norrislab/Zach/STAR_results" 

prefixes=("ADL" "AFD" "AIM" "AIN" "AIY" "ASEL" "ASER" "ASG" "ASI" "ASK" "AVA" "AVE" "AVG" "AVH" "AVK" "AVL" "AVM" "AWA" "AWB" "AWC" "BAG" "CAN" "DA" "DD" "DVC" "I5" "IL1" "IL2" "LUA" "NSM" "OLL" "OLQ" "PHA" "PVC" "PVD" "PVM" "RIA" "RIC" "RIM" "RIS" "RMD" "SMB" "SMD" "VB" "VC" "VD") 

for prefix1 in "${prefixes[@]}"; do 
 for prefix2 in "${prefixes[@]}"; do 
  if [[ "$prefix1" != "$prefix2" ]]; then 
            source_directory="$base_path/${prefix1}vs${prefix2}/JUM_diff/FINAL_JUM_OUTPUT_pvalue_1" 
            destination_directory="$base_path/${prefix1}vs${prefix2}/JUM_diff/FINAL_JUM_OUTPUT_pvalue_1${prefix1}vs${prefix2}" 

# Create the destination directory 
 mkdir -p "$destination_directory" 

# Copy all files from the source directory to the destination directory 
 cp -r "$source_directory"/* "$destination_directory"   

# Copy newly created detination directory to the newly created JUM_results folder: 
 cp -r "$destination_directory" "/lustre/work/client/group/norrislab/Zach/JUM_results" 
  fi 
 done 
done 
```
 
DCC on all neurons using all STAR and samtools files: 

We can also run DCC. If you haven’t already, install DCC and its dependencies using micromamba:  

```
cd $WORK; curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba 
./bin/micromamba shell init -s bash -p ~/micromamba 
source ~/.bashrc 
micromamba create --name pysam_env -c defaults -c bioconda -c r -c conda-forge --override-channels  python=3.8 numpy pip jupyterlab matplotlib pandas pysam htseq 
git clone https://github.com/dieterich-lab/DCC.git 
cd DCC 

module load gcc 
micromamba activate pysam_env 
python setup.py install –user 
```

We can go back to /lustre/work/client/group/norrislab/Zach and create a new directory called /STAR_results_single_mates in which we will run STAR on individual mates of each replicate (dos2unix and sbatch this script): 

```
#!/bin/bash 
#SBATCH -N 1 
#SBATCH -o output.out 
#SBATCH -e output.err 
#SBATCH -p highmem 

module load gcc star 

STAR --runThreadN 10 --genomeDir /work/users/zwolfe/genome/ --outSAMtype BAM SortedByCoordinate --readFilesIn /lustre/work/client/group/norrislab/Zach/SRR24086885_1.fastq --outSAMstrandField intronMotif --readFilesCommand cat --outReadsUnmapped Fastx --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMatchNmin 1 --outFilterMismatchNmax 2 --chimSegmentMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15 --alignTranscriptsPerReadNmax 200000 --alignTranscriptsPerWindowNmax 20000 \ 
 --outFileNamePrefix OLL1_1 
```

We also need to create some .txt files that DCC reads as part of its function. These are called samplesheet.txt, mate1.txt, and mate2.txt. We need to create these for each replicate: 
samplesheet.txt will represent the file path of the chimeric.out.junction file from the combined STAR mapping (remember to dos2unix): 

```
/lustre/work/client/group/norrislab/Zach/OLQ1Chimeric.out.junction 
```

mate1.txt will represent the file path of the chimeric.out.junction file from the STAR mapping of the first mate only (remember to dos2unix): 

```
/lustre/work/client/group/norrislab/Zach/STAR_results_single_mates/OLQ1_1Chimeric.out.junction 
```

mate2.txt will represent the file path of the chimeric.out.junction file from the STAR mapping of the second mate only (remember to dos2unix): 

```
/lustre/work/client/group/norrislab/Zach/STAR_results_single_mates/OLQ1_2Chimeric.out.junction 
```

Finally, we need to create a batch script called DCC.txt which will actually run DCC on our samples (dos2unix and sbatch this script): 

```
#!/bin/bash 
#SBATCH -N 1 
#SBATCH -o output.out 
#SBATCH -e output.err 
#SBATCH -p highmem 
#SBATCH --mem=200G 

DCC @samplesheet.txt -mt1 @mate1.txt -mt2 @mate2.txt -D -an /work/users/zwolfe/genome/Caenorhabditis_elegans.WBcel235.107.gtf -B /lustre/work/client/group/norrislab/Zach/OLQ1Aligned.sortedByCoord.out.bam -Pi -N -F -M -Nr 5 1 -fg -G -A /work/users/zwolfe/genome/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa  
```

It would be less laborious to perform all of these functions in a loop at once. I have provided the loop to create all of these files and run DCC here – just make sure you are in your pysam_env before running (dos2unix and sbatch this script): 

```
#!/bin/bash 
#SBATCH -N 1 
#SBATCH -o output.out 
#SBATCH -e output.err 
#SBATCH -p highmem 
#SBATCH --mem=200G 

prefixes=(“ADL1" "ADL2" …"VD4") 

for prefix in "${prefixes[@]}"; do 
  samplesheet_file="/lustre/work/client/group/norrislab/Zach/${prefix}Chimeric.out.junction" 
  mate1_file="/lustre/work/client/group/norrislab/Zach/STAR_results_single_mates/${prefix}_1Chimeric.out.junction" 
  mate2_file="/lustre/work/client/group/norrislab/Zach/STAR_results_single_mates/${prefix}_2Chimeric.out.junction" 
  output_dir="/lustre/work/client/group/norrislab/Zach/STAR_results_single_mates/DCC/${prefix}_DCC_output"
  mkdir -p "$output_dir"  # Create the output directory if it doesn't exist 

  DCC "$samplesheet_file" -mt1 "$mate1_file" -mt2 "$mate2_file" -D -an /work/users/zwolfe/genome/Caenorhabditis_elegans.WBcel235.107.gtf -B /lustre/work/client/group/norrislab/Zach/${prefix}Aligned.sortedByCoord.out.bam -Pi -N -F -M -Nr 5 1 -fg -G -A /work/users/zwolfe/genome/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa 

  # Move the output files to the output directory 
  mv LinearCount "$output_dir" 
  mv CircSkipJunctions "$output_dir" 
  mv CircCoordinates "$output_dir" 
  mv CircRNACount "$output_dir" 
  mv DCC-*.log "$output_dir" 
done 
```

It will be worthwhile to rename these DCC output files with their replicate name as their respective prefix - (dos2unix and sbatch this script): 

WARNING: Sometimes no circular RNA is detected in a given replicate. If that is the case, no output files will be generated. Make sure to exclude these replicates from the mv loop below: 

```
#!/bin/bash 
#SBATCH -N 1 
#SBATCH -o output.out 
#SBATCH -e output.err 
#SBATCH -p standard-s 
#SBATCH --mem=20G 

prefixes=(“ADL1" "ADL2" …"VD4") 

for prefix in "${prefixes[@]}"; do 
  output_dir="/lustre/work/client/group/norrislab/Zach/STAR_results_single_mates/DCC/${prefix}_DCC_output" 

  # Rename specific output files 
  mv "$output_dir/CircSkipJunctions" "$output_dir/${prefix}CircSkipJunctions" 
  mv "$output_dir/LinearCount" "$output_dir/${prefix}LinearCount" 
  mv "$output_dir/CircRNACount" "$output_dir/${prefix}CircRNACount" 
  mv "$output_dir/CircCoordinates" "$output_dir/${prefix}CircCoordinates" 
done 
``` 

Now we can move these output files to our hard drive and run them in the R script circRNA filtering and testing.Rmd 

Copyright 2024 The Regents of the University of California

All Rights Reserved

Created by Zachery Wolfe

Department of Biochemistry

This file is part of Differential Expression in C. elegans. \
Differential Expression in C. elegans is a free pipeline: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. \
Differential Expression in C. elegans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. \
You should have received a copy of the GNU General Public License along with Differential Expression in C. elegans. If not, see <https://www.gnu.org/licenses/>.
