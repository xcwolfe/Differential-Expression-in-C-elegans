# Differential-Expression-in-C-elegans
This repository contains the scripts and commands necessary to analyze differential gene, RNA-binding protein, and A.S. expression in C. elegans as of February 2023.

I have ordered the R scripts in chronological order (1-7). I recommend running all scripts in this order to generate the correct files and file types.

You can run JUM workflow any time before R script #5 (JUM file automation.R)

The Bulk RNA-seq data I use comes from Barrett et al. 2022 (CeNGEN)
I have also acquired the names of 484 known C. elegans RNA-binding proteins from WormBase (simplemine_results.csv)

You will need to adjust the parameters in these files for the specific cell types you are interested in - I was interested in 41 C. elegans neurons, for instance, so that is how my pipeline has been constructed.

Zachery Wolfe, Ph.D.

Postdoctoral Researcher

Southern Methodist University

zwolfe@smu.edu

