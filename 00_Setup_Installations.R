# Coastal MicroEco Lab dada2 Tutorial
# 00_Prep work/Installations

# This tutorial has assumed you have already installed R and RStudio
# See "git_token_classic.R" to setup token for git

# Bioconductor----
# Install bioconductor: https://www.bioconductor.org/install/
# Bioconductor is an open-source, R-based project that provides tools for the analysis and 
# comprehension of high-throughput genomic and molecular biology data. I
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.22")

# dada2----
# Install dada2 (using BioCManager): https://benjjneb.github.io/dada2/dada-installation.html
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.22") # change version as necessary


