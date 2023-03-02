# DCC pipeline:

 ```
$ mkdir genome

$ cd genome/

$ wget ftp://ftp.ensembl.org/pub/release-107/fasta/caenorhabditis_elegans/dna/*.fa.gz

$ gunzip *.gz
 ```

Go to http://useast.ensembl.org/Caenorhabditis_elegans/Info/Index and download the latest *C. elegans* genome as a .gtf file:

`$ ftp://ftp.ensembl.org/pub/release-107/gtf/caenorhabditis_elegans/*.gz`

Initialize STAR on your genome:

  a.	create initialize_star.txt:
 
```
#!/bin/bash

#SBATCH -N 1
#SBATCH -t 1000
#SBATCH -o output.out
#SBATCH -e output.err
#SBATCH -p high-mem-1
module load star
STAR \
 --runThreadN 1 \
 --runMode genomeGenerate \
 --genomeDir /work/users/zwolfe/genome \
 --genomeFastaFiles /work/users/zwolfe/genome/*.fa \
 --sjdbGTFfile /work/users/zwolfe/genome/c_elegans.gtf \
 --sjdbOverhang 99

$ dos2unix initialize_star.txt

$ sbatch initialize_star.txt
 ```
  
Installed DCC and its dependencies:

```
$ cd $WORK; curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba

$ ./bin/micromamba shell init -s bash -p ~/micromamba

$ source ~/.bashrc

$ micromamba create --name pysam_env -c defaults -c bioconda -c r -c conda-forge --override-channels  python=3.8 numpy pip jupyterlab matplotlib pandas pysam htseq

$ git clone https://github.com/dieterich-lab/DCC.git

$ cd DCC

$ module load gcc-9.2

$ micromamba activate pysam_env

$ python setup.py install –user
 ```

Create ASG1 directory in DCC folder:

`$ mkdir ASG1`

Create mate1 and mate2 folders in ASG1 directory:

```
$ cd ASG1

$ mkdir mate1

$ mkdir mate2
 ```

Moved ASG_SRR13995310_01.fastq and ASG_SRR13995310_02.fastq files to ASG1 directory

Moved ASG_SRR13995310_01.fastq to mate1 folder and ASG_SRR13995310_02.fastq to mate2 folder

Create samplesheet.txt  and move into ASG1 directory:

```
/work/users/zwolfe/DCC/ASG1/Chimeric.out.junction
```

Create mate1.txt  and move into ASG1 directory:

```
/work/users/zwolfe/DCC/ASG1/mate1/Chimeric.out.junction
```

Create mate2.txt  and move into ASG1 directory:

```
/work/users/zwolfe/DCC/ASG1/mate2/Chimeric.out.junction
```
```
$ Module load star
$ Module load samtools
$ micromamba activate pysam_env
```

Create STAR_in_DCC.txt :

```
#!/bin/bash

#SBATCH -N 1
#SBATCH -t 10000
#SBATCH -o output.out
#SBATCH -e output.err
#SBATCH -p high-mem-1

module load star
STAR --runThreadN 10 --genomeDir /work/users/zwolfe/genome/ --outSAMtype BAM SortedByCoordinate --readFilesIn *_1.fastq *_2.fastq --readFilesCommand cat --outReadsUnmapped Fastx --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMatchNmin 1 --outFilterMismatchNmax 2 --chimSegmentMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15 --alignTranscriptsPerReadNmax 200000 --alignTranscriptsPerWindowNmax 20000 

$ dos2unix STAR_in_DCC.txt
$ sbatch STAR_in_DCC.txt
```

```
$ cd mate1
```

Create STAR_in_DCC_mate_1.txt:

```
#!/bin/bash

#SBATCH -N 1
#SBATCH -t 10000
#SBATCH -o output.out
#SBATCH -e output.err
#SBATCH -p high-mem-1

module load star
STAR --runThreadN 10 --genomeDir /work/users/zwolfe/genome/ --outSAMtype None --readFilesIn *_1.fastq --readFilesCommand cat --outReadsUnmapped Fastx --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMatchNmin 1 --outFilterMismatchNmax 2 --chimSegmentMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15 --alignTranscriptsPerReadNmax 200000 --alignTranscriptsPerWindowNmax 20000 

$ dos2unix STAR_in_DCC_mate_1.txt
$ sbatch STAR_in_DCC_mate_1.txt

$ cd .. 
```

```
$ cd mate2
```

Create STAR_in_DCC_mate_2.txt:

```
#!/bin/bash

#SBATCH -N 1
#SBATCH -t 10000
#SBATCH -o output.out
#SBATCH -e output.err
#SBATCH -p high-mem-1

module load star
STAR --runThreadN 10 --genomeDir /work/users/zwolfe/genome/ --outSAMtype None --readFilesIn *_2.fastq --readFilesCommand cat --outReadsUnmapped Fastx --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMatchNmin 1 --outFilterMismatchNmax 2 --chimSegmentMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15 --alignTranscriptsPerReadNmax 200000 --alignTranscriptsPerWindowNmax 20000 

$ dos2unix STAR_in_DCC_mate_2.txt
$ sbatch STAR_in_DCC_mate_2.txt
```

```
33.	cd .. 
```

Create run_samtools_DCC.txt:

```
#!/bin/bash

#SBATCH -N 2
#SBATCH -t 180
#SBATCH -o output.out
#SBATCH -e output.err
#SBATCH -p standard-mem-s

module load samtools

samtools index Aligned.sortedByCoord.out.bam

$ dos2unix run_samtools_DCC.txt
$ sbatch run_samtools_DCC.txt
```

```
$ dos2unix DCC_ASG.txt 
$ dos2unix samplesheet.txt
$ dos2unix mate1.txt
$ dos2unix mate2.txt
$ sbatch DCC_ASG.txt 

Once you have completed each DCC run, rename your DCC output files with the cell type and replicate number as a prefix for the output files:
a.	CircSkipJunctions = ASG1CircSkipJunctions
b.	LinearCount = ASG1LinearCount
c.	CircRNACount = ASG1CircRNACount
d.	CircCoordinates = ASG1CircCoordinates

Download each output file to your R working directory
