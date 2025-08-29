# SPIEC-EASI Application Project

# Packages installation
This project requires a mix of CRAN and Bioconductor packages.
Below is a breakdown of which packages come from where, and how to install them.

## 📦 CRAN Packages
install with
```r{}
install.packages(c(
  "igraph", 
  "MASS", 
  "pheatmap", 
  "huge", 
  "pulsar", 
  "MLmetrics", 
  "brew", 
  "tidyverse", 
  "ggplot2", 
  "dplyr", 
  "tibble", 
  "tidyr", 
  "Matrix"
))
```
## 🧬 Bioconductor Packages
First, install Bioconductor’s manager if you don’t have it yet:
```r{}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```
Then install the required Bioconductor packages:
```r{}
BiocManager::install(c(
  "phyloseq", 
  "microbiome", 
  "SpiecEasi", 
  "NetCoMi", 
  "limma", 
  "ComplexHeatmap", 
  "circlize"
))
```


# 📜 R_Scripts/
In this folder you can find all used R Scripts for this seminar paper.
* ApplicationNet.R – runs the main SPIEC-EASI network application on microbiome data
* SimulationHUB.R – simulates two-hub-based networks for testing network inference methods
* SimulationTenHUB.R – simulates ten-hub-based network
  
# 📊 Plots/
In this folder you can find all plots used in this paper. 
Some of the plots follow a certain scheme:
* p40n200_* ->low-dimensional plots from the two-hub simulation
* p200n200_* -> high-dimensional plots from the ten-hub simulation

# 🦠 agp_obese_underweigt/
This folder contains the code for the data aggregation (in the subfolder analysis/).
Since the original dataset is too large to upload to github, the already processed data is included in the subfolder data/.
