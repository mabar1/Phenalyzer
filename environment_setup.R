# environment_setup.R
# Run this script before using Phenalyzer_v1.R to ensure all required packages are installed and loaded.

# Helper function to install CRAN packages if missing
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Helper function to install Bioconductor packages if missing
install_bioc_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(pkg)
  }
}

# Helper function to install GitHub packages if missing
install_github_if_missing <- function(repo, pkg = NULL) {
  if (is.null(pkg)) pkg <- sub(".*/", "", repo)
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
    remotes::install_github(repo)
  }
}

install_github_devtools_if_missing <- function(repo, pkg = NULL) {
  if (is.null(pkg)) pkg <- sub(".*/", "", repo)
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
    devtools::install_github(repo)
  }
}


# CRAN packages
cran_packages <- c(
  "devtools", "rstudioapi", "progress", "data.table", "dplyr", "tibble", "ggplot2", "cowplot",
  "ggridges", "scico", "DescTools", "Polychrome", "stringr", "gridExtra",
  "patchwork", "tidyr", "magrittr", "forcats", "remotes", "uwot", "ggsignif", "ggrepel",
  "raster","tiff","sf","s2","stars","exactextractr","tidyr","readxl","ggbeeswarm",
   "fastcluster", "dendsort",  "ggh4x","tidygam","tidyverse"

)

# Bioconductor packages
bioc_packages <- c(
   "ComplexHeatmap", 
   "circlize", 
   "flowCore", 
   "cytolib", 
   "FlowSOM",
   "rhdf5", 
   "HDF5Array",
   "JinmiaoChenLab/cytofkit"
)

# GitHub packages
github_packages <- list(
  Rphenoannoy = "stuchly/Rphenoannoy",
  ggpubr = "kassambara/ggpubr",
  spectre = "immunedynamics/spectre" # <- make sure to install this last since its dependencies are fussy!
  
   #Leiden: TomKellyGenetics Leiden engine (if needed)
   #leiden = "TomKellyGenetics/leiden"
)

# List of GitHub packages to install with devtools
github_devtools_packages <- list(
  ggpubr = "kassambara/ggpubr"
)


# Optional: install nortest if you want to use ad.test (mentioned in comments)
if (!requireNamespace("nortest", quietly = TRUE)) install.packages("nortest")


# Install CRAN packages
for (pkg in cran_packages) install_if_missing(pkg)

# Install Bioconductor packages
for (pkg in bioc_packages) install_bioc_if_missing(pkg)

# Install GitHub packages
for (pkg in github_packages) install_github_if_missing(pkg)

# Install GitHub packages with devtools
for (pkg in github_devtools_packages) install_github_devtools_if_missing(pkg)

# Load all packages (suppress warnings for already loaded)
all_packages <- c(
  cran_packages,
  bioc_packages,
  names(github_packages)
)
for (pkg in all_packages) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Additional: install tidyverse if not present (covers dplyr, tibble, tidyr, ggplot2, etc.)
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")



cat("All required packages are installed and loaded.\n")