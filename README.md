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

I have ordered the R scripts (1-12). I recommend running all scripts in this order to generate the correct files and file formats. There is a bonus script (#13) on evolutionary conservation that ultimately was not included in the final manuscript: Wolfe, Z., Liska, D., & Norris, A. (2025). Deep transcriptomics reveals cell-specific isoforms of pan-neuronal genes. *Nature Communications, 16*(1), 1-14.

You can run STAR_and_JUM.md any time before R script #6 (JUM_file_automation.R), although I would recommend running it as soon as possible. If you are starting this pipeline with .fastq or .fastq.gz files, you *have* to run STAR and HTSeq before the first script.

The Bulk RNA-seq data I use (as an example) comes from Barrett et al. 2022 (CeNGEN) (PRJNA952691). A sample count matrix is provided for all 46 cell types (Barrett_et_al_2022_CeNGEN_bulk_RNAseq_data.csv). This repository contains an additional dataset (Barrett_et_al_2022_CeNGEN_practice_RNAseq_data.csv) which can be used as a practice dataset; it contains only 8 cell types as opposed to the 46 used in the full dataset.

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

Academic Software License: Copyright 2024 UCR (“Institution”). Academic or nonprofit researchers are permitted to use this Software (as defined below) subject to Paragraphs 1-4:

1. Institution hereby grants to you free of charge, so long as you are an academic or nonprofit researcher, a nonexclusive license under Institution’s copyright ownership interest in this software and any derivative works made by you thereof (collectively, the “Software”) to use, copy, and make derivative works of the Software solely for educational or academic research purposes, and to distribute such Software free of charge to other academic or nonprofit researchers for their educational or academic research purposes, in all cases subject to the terms of this Academic Software License. Except as granted herein, all rights are reserved by Institution, including the right to pursue patent protection of the Software.

2. Any distribution of copies of this Software -- including any derivative works made by you thereof -- must include a copy (including the copyright notice above), and be made subject to the terms, of this Academic Software License; failure by you to adhere to the requirements in Paragraphs 1 and 2 will result in immediate termination of the license granted to you pursuant to this Academic Software License effective as of the date you first used the Software.

3. IN NO EVENT WILL INSTITUTION BE LIABLE TO ANY ENTITY OR PERSON FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF INSTITUTION HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. INSTITUTION SPECIFICALLY DISCLAIMS ANY AND ALL WARRANTIES, EXPRESS AND IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE IS PROVIDED “AS IS.” INSTITUTION HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS OF THIS SOFTWARE.

4. Any academic or scholarly publication arising from the use of this Software or any derivative works thereof will include the following acknowledgment:  The Software used in this research was created by Zachery Wolfe of University of California, Riverside. Copyright 2024 UCR.

Commercial entities: please contact zacheryw@ucr.edu for licensing opportunities.

This file is part of Differential Expression in C. elegans. \
Differential Expression in C. elegans is a free pipeline: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. \
Differential Expression in C. elegans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. \
You should have received a copy of the GNU General Public License along with Differential Expression in C. elegans. If not, see <https://www.gnu.org/licenses/>.

