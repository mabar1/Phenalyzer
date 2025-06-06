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

# CRAN packages
cran_packages <- c(
  "rstudioapi", "progress", "data.table", "dplyr", "tibble", "ggplot2", "ggpubr", "cowplot",
  "ggridges", "scico", "quasirandom", "DescTools", "Polychrome", "stringr", "gridExtra",
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
  spectre = "immunedynamics/spectre",
  Rphenoannoy = "stuchly/Rphenoannoy",
   #Leiden: TomKellyGenetics Leiden engine (if needed)
   leiden = "TomKellyGenetics/leiden"
)


# Optional: install nortest if you want to use ad.test (mentioned in comments)
if (!requireNamespace("nortest", quietly = TRUE)) install.packages("nortest")


# Install CRAN packages
for (pkg in cran_packages) install_if_missing(pkg)

# Install Bioconductor packages
for (pkg in bioc_packages) install_bioc_if_missing(pkg)

# Install GitHub packages
for (pkg in github_packages) install_github_if_missing(pkg)

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