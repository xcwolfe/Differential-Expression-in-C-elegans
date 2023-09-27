# Differential Expression in *C. elegans*
This repository contains the scripts, commands, and input files necessary to analyze differential gene, RNA-binding protein, and A.S. expression in *C. elegans* (as of September 2023).

I have ordered the R scripts in chronological order (1-9). I recommend running all scripts in this order to generate the correct files and file types.

You can run STAR and JUM workflow any time before R script #6 (JUM file automation.R), although I would recommend running it as soon as possible. If you are starting this pipeline with .fastq or .fastq.gz files, you *have* to run STAR before the first script.

The Bulk RNA-seq data I use comes from Barrett et al. 2022 (CeNGEN)

I have also acquired the names of 31173 known *C. elegans* genes (genenames.csv) and 484 known *C. elegans* RNA-binding proteins (simplemine_results.csv) from WormBase

You will need to adjust the parameters in these files for the specific cell types you are interested in - I was interested in 46 *C. elegans* neurons, for instance, so that is how my pipeline has been constructed.

Zachery Wolfe, Ph.D.

Postdoctoral Researcher

Southern Methodist University

zwolfe@smu.edu

