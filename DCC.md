# DCC pipeline:

`mkdir genome
cd genome/
wget ftp://ftp.ensembl.org/pub/release-107/fasta/caenorhabditis_elegans/dna/*.fa.gz
gunzip *.gz`
5.	go to http://useast.ensembl.org/Caenorhabditis_elegans/Info/Index and download the latest C. elegans genome as a .gtf file
6.	gunzip ftp://ftp.ensembl.org/pub/release-107/gtf/caenorhabditis_elegans/*.gz
7.	Initialize STAR on your genome:
a.	create initialize_star.txt 
b.	dos2unix initialize_star.txt
c.	sbatch initialize_star.txt
8.	Installed DCC and its dependencies:
a.	cd $WORK; curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
b.	./bin/micromamba shell init -s bash -p ~/micromamba
c.	source ~/.bashrc
d.	micromamba create --name pysam_env -c defaults -c bioconda -c r -c conda-forge --override-channels  python=3.8 numpy pip jupyterlab matplotlib pandas pysam htseq
e.	git clone https://github.com/dieterich-lab/DCC.git
f.	cd DCC
g.	module load gcc-9.2
h.	micromamba activate pysam_env
i.	python setup.py install –user

9.	Created ASG1 directory in DCC folder:
a.	mkdir ASG1
10.	Created mate1 and mate2 folders in ASG1 directory:
a.	cd ASG1
b.	mkdir mate1
c.	mkdir mate2
11.	Moved ASG_SRR13995310_01.fastq and ASG_SRR13995310_02.fastq files to ASG1 directory.
12.	Moved ASG_SRR13995310_01.fastq to mate1 folder and ASG_SRR13995310_02.fastq to mate2 folder
13.	Create samplesheet.txt  and move into ASG1 directory
14.	Create mate1.txt  and move into ASG1 directory
15.	Create mate2.txt  and move into ASG1 directory
16.	Module load star
17.	Module load samtools
18.	micromamba activate pysam_env 
19.	Create STAR_in_DCC.txt 
20.	dos2unix STAR_in_DCC.txt
21.	sbatch STAR_in_DCC.txt
22.	cd mate1
23.	Create STAR_in_DCC_mate_1.txt 
24.	dos2unix STAR_in_DCC_mate_1.txt
25.	sbatch STAR_in_DCC_mate_1.txt
26.	cd .. 
27.	cd mate2
28.	Create STAR_in_DCC_mate_2.txt 
29.	dos2unix STAR_in_DCC_mate_2.txt
30.	sbatch STAR_in_DCC_mate_2.txt
31.	cd .. 
32.	Create run_samtools_DCC.txt 
33.	dos2unix run_samtools_DCC.txt
34.	sbatch run_samtools_DCC.txt
35.	dos2unix DCC_ASG.txt 
36.	dos2unix samplesheet.txt
37.	dos2unix mate1.txt
38.	dos2unix mate2.txt
39.	sbatch DCC_ASG.txt 
40.	Once you have completed each DCC run, rename your DCC output files with the cell type and replicate number as a prefix for the output files:
a.	CircSkipJunctions = ASG1CircSkipJunctions
b.	LinearCount = ASG1LinearCount
c.	CircRNACount = ASG1CircRNACount
d.	CircCoordinates = ASG1CircCoordinates
41.	Download each output file to your R working directory