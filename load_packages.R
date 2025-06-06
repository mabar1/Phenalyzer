### Load libraries
  
  ## load all packages of spectre
  #if(!require('remotes')) {install.packages('remotes')} # Installs the package 'remotes'
  #cytolib flowCore FlowSOM does not work
  #BiocManager::install("cytolib")
  # BiocManager::install("flowCore")
  #BiocManager::install("FlowSOM")
  #remotes::install_github(repo = "immunedynamics/spectre") # Install the Spectre package
library(Spectre)

#Spectre::package.check(type = 'spatial')
#Spectre::package.load(type = 'spatial')

        #1: In library(package, lib.loc = lib.loc, character.only = TRUE, logical.return = TRUE,  :
        #there is no package called ‘raster’
        #install.packages("raster")
        #: In library(package, lib.loc = lib.loc, character.only = TRUE, logical.return = TRUE,  :
         # there is no package called ‘tiff’
          #there is no package called ‘exactextractr’
           # there is no package called ‘sf’
            #ere is no package called ‘stars’
             # there is no package called ‘s2’
             # same for all of them, tiff, stars, s2, exactextractr, raster
      #  1: In library(package, lib.loc = lib.loc, character.only = TRUE, logical.return = TRUE,  :
  #there is no package called ‘rhdf5’
#2: In library(package, lib.loc = lib.loc, character.only = TRUE, logical.return = TRUE,  :
#  there is no package called ‘HDF5Array’


        # install.packages("qs")
        #install.packages("sf")
library(qs)
library(sf)

# we need to tweak the pre-defined functions once again:
# we gonna bring in here the adapted color.plot function from shinyIMC:
library(ggrepel) # repelling labels in make.colour.plot
library(scico) # just like brewer, another package to deliver color grads. We need that for the scatter plots of two variables, along with the stat_density_2d lines
library(Polychrome) # another color palette aiming for colors maximally separated
library(ggbeeswarm) # call quasirandom 
#library(tidyverse) # use pipes to handle data.frames
library(stringr) # use strings to process conditions. we use str_detect to pull out strings from Diet vector

# use Bioconductor to get these in R 4.2:
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
library(cytofkit) # rphenograph this masks now Rphenograph function. see https://github.com/JinmiaoChenLab/cytofkit BiocManager::install("JinmiaoChenLab/cytofkit") 
# Rphenoannoy is installed via git:
#devtools::install_github("stuchly/Rphenoannoy")
library(Rphenoannoy) # faster Phenograph engine 
#library("leiden") # install.packages("leiden") https://cran.r-project.org/web/packages/leiden/vignettes/run_leiden.html
library(ComplexHeatmap) # use complexheatmap to plot hm. install via Biocunductor BiocManager::install("ComplexHeatmap")

library(uwot) # umap


library(dendsort) # call in hclust and sort to find the two most similar clusters to collapse
library(fastcluster) # lets use a faster engine than dendsort in the engine
suppressPackageStartupMessages(library(circlize)) # we use a col grad from here to adjust blue-red according to the z-score range in the heatmap

library(progress) # for progress bars: prepare them by  pb <- progress_bar$new(total = amount.of.TMAs ) https://cran.r-project.org/web/packages/progress/progress.pdf

library(ggh4x) # use it for the background strip colors in facet_wrap2

# Fun with GAM models to batch normalize 
library(mgcv) # use additive models for the gam() functoin


#library(tidymv) # use this for plot_smooths( in ggplot
#remotes::install_github("stefanocoretta/tidygam@v1.0.0")
library(tidygam)

library(tibble)
library(tidyverse)
        
library(readxl)