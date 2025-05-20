SARS/MERS Genomic Analysis with R
This repository contains an R script (lab8.R) for performing phylogenetic analysis on SARS-CoV-2 and MERS genomic sequences. The script processes FASTA files, conducts multiple sequence alignments (MSA) with ClustalW, constructs phylogenetic trees using PhyML, and visualizes results with ggplot2. This project showcases skills in R programming, bioinformatics workflows, and data visualization, as part of my bioinformatics coursework at Langara College.
Overview

Purpose: Investigate evolutionary relationships between SARS-CoV-2 and MERS sequences using phylogenetic tree construction and visualization.
Features:
Processes SARS-CoV-2 and MERS FASTA sequences for alignment and phylogenetic analysis.
Performs MSA using clustalw2 via the ape package.
Constructs phylogenetic trees with phyml using the phangorn package.
Preprocesses sequences with seqinr for efficient data handling.
Visualizes phylogenetic relationships using ggplot2 and manages data with dplyr.
Automates dependency installation with bash scripts for reproducibility.


Date: Developed in Fall 2024.

Prerequisites

R 4.0+
Required R packages:
ape (for sequence alignment)
phangorn (for phylogenetic tree construction)
seqinr (for sequence preprocessing)
ggplot2 (for visualization)
dplyr (for data wrangling)


External tools:
clustalw2 (for multiple sequence alignment)
phyml (for phylogenetic tree construction)


A working internet connection to download FASTA files (or local copies).

Installation

Clone the Repository:
git clone https://github.com/ilaypaz/Bioinformatics.git
cd Bioinformatics


Install R Packages:

Open R or RStudio and run:install.packages(c("ape", "phangorn", "seqinr", "ggplot2", "dplyr"))




Install External Tools:

Using conda (recommended):conda install -c bioconda clustalw phyml


Alternatively, download and install manually:
ClustalW
PhyML




Verify Tools:

Ensure clustalw2 and phyml are in your PATH. Test with:clustalw2 -version
phyml --version





Usage

Prepare FASTA Files:

Place SARS-CoV-2 and MERS FASTA files in a data directory within the project folder, or update the script to point to your file locations.
Example files can be sourced from public databases like NCBI GenBank.


Run the Script:

In R or RStudio, set the working directory and run:setwd("path/to/Bioinformatics")
source("lab8.R")


Alternatively, from the command line:Rscript lab8.R




Expected Output:

Alignment Files: SARS_CoV_2_msa.fasta, MERS_msa.fasta
Tree Files: SARS_CoV_2_tree.tre, MERS_tree.tre
Plots: SARS_CoV_2_dendrogram.pdf, MERS_dendrogram.pdf
Console Output: Sequence IDs, alignment summaries, and phylogenetic tree statistics.


Notes:

The script includes error handling for missing dependencies or files.
Adjust file paths in the script if your FASTA files are stored elsewhere.
Bash scripts for installing clustalw2 and phyml are included for macOS; modify for other operating systems if needed.



Outputs and Analysis

Phylogenetic Trees: Trees constructed with phyml provide insights into the evolutionary relationships between SARS-CoV-2 and MERS sequences.
Visualizations: Publication-quality dendrograms generated with ggplot2, showing sequence divergence and clustering.
Data Management: Cleaned and processed datasets using dplyr, ensuring efficient analysis workflows.

Acknowledgments

FASTA sequences were sourced from NCBI GenBank.
Utilizes R packages ape, phangorn, seqinr, ggplot2, and dplyr for analysis and visualization.
Developed as part of bioinformatics coursework at Langara College, Fall 2024.

Contact

Author: Ilay Paz
Email: pazlilay@gmail.com
GitHub: github.com/ilaypaz

Last updated: May 20, 2025, 01:53 AM PDT
