# Differential Expression in *C. elegans*
This repository contains the scripts, commands, and input files necessary to analyze differential gene, RNA-binding protein, and A.S. expression in *C. elegans* (as of September 2024).

System requirements:
R version 4.4.1 (Race for Your Life) or greater
Rstudio version 2024.09.0+375 or greater
A command-line interface (shell) capable of running the UNIX operating system

R can be downloaded and installed from the comprehensive R Archive Network (https://cran.r-project.org/)

I have ordered the R scripts in chronological order (1-12). I recommend running all scripts in this order to generate the correct files and file types.

You can run STAR and JUM.md any time before R script #6 (JUM file automation.R), although I would recommend running it as soon as possible. If you are starting this pipeline with .fastq or .fastq.gz files, you *have* to run STAR and HTSeq before the first script.

The Bulk RNA-seq data I use (as an example) comes from Barrett et al. 2022 (CeNGEN). A sample count matrix is provided for all 46 cell types (Barrett_et_al_2022_CeNGEN_bulk_RNAseq_data.csv). This repository contains an additional dataset (Barrett_et_al_2022_CeNGEN_practice_RNAseq_data.csv) which can be used as a practice dataset; it contains only 8 cell types as opposed to the 46 used in the full dataset.

I have also acquired the names of 31173 known *C. elegans* genes (genenames.csv) and 484 known *C. elegans* RNA-binding proteins (simplemine_results.csv) from WormBase

You will need to adjust the parameters in these files for the specific cell types you are interested in - I was interested in 46 *C. elegans* neurons, for instance, so that is how my pipeline has been constructed.

Zachery Wolfe, Ph.D.

Postdoctoral Researcher

University of California, Riverside

zacheryw@ucr.edu

