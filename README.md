# Differential Expression in *C. elegans*
This repository contains the scripts, commands, and input files necessary to analyze differential gene, RNA-binding protein, and A.S. expression in *C. elegans* (as of October 2024).

System requirements:
R version 4.4.1 (Race for Your Life) or greater \
Rstudio version 2024.09.0+375 or greater \
A command-line interface (shell) capable of running the UNIX operating system \
JUM V3.0.0 or greater

R can be downloaded and installed from the comprehensive R Archive Network (https://cran.r-project.org/) \
JUM and its dependent scripts can be downloaded from https://github.com/qqwang-berkeley/JUM \
All system requirements can be downloaded within a matter of minutes on a typical desktop computer

I have ordered the R scripts in chronological order (1-12). I recommend running all scripts in this order to generate the correct files and file types.

You can run STAR and JUM.md any time before R script #6 (JUM file automation.R), although I would recommend running it as soon as possible. If you are starting this pipeline with .fastq or .fastq.gz files, you *have* to run STAR and HTSeq before the first script.

The Bulk RNA-seq data I use (as an example) comes from Barrett et al. 2022 (CeNGEN). A sample count matrix is provided for all 46 cell types (Barrett_et_al_2022_CeNGEN_bulk_RNAseq_data.csv). This repository contains an additional dataset (Barrett_et_al_2022_CeNGEN_practice_RNAseq_data.csv) which can be used as a practice dataset; it contains only 8 cell types as opposed to the 46 used in the full dataset.

If using the practice dataset, you should expect an output that includes the following (from .Rmds 1, 2, 4, and 5):
1. Associated PCAs (from the underlying provided count matrix) for 8 cell types
2. Differential expression of detected genes, including RNA-binding proteins, between the 8 cell types (both in list and heatmap form)
3. Uniqueness indices and associated heatmaps for all genes, including RNA-binding proteins, representative of the 8 cell types \

These four .Rmds can be run on the practice dataset within a matter of hours; perhaps less if using a computer with a powerful processor (or in an HPC environment)

I have also provided the names of 31173 known *C. elegans* genes (genenames.csv) and 484 known *C. elegans* RNA-binding proteins (simplemine_results.csv) from WormBase

You will need to adjust the parameters in these markdowns for the specific cell types you are interested in - I was interested in 46 *C. elegans* neurons, for instance, so that is how my pipeline has been constructed.

Please email me or comment in the issues tab if you have any questions or find any errors; I will address them as soon as possible.

Zachery Wolfe, Ph.D.

Postdoctoral Researcher

University of California, Riverside

zacheryw@ucr.edu

Copyright 2024 The Regents of the University of California

All Rights Reserved

Created by Zachery Wolfe

Department of Biochemistry
