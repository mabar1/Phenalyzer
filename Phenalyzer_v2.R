###################################################################################
### Spectre: spatial 2 - cellular analysis
####
# This is v19 brought over from the Felsenstein-steinbock project on 20250328 and adapted to fit into the PCF workflow via QuPath Cellpose Segmentation
####
#   
#   
#
#   Here we have to load the metadata file that is fixed in the Analysis folder. The same file is used by manual gating method
#   There is one and only one file like that. As of v5 this is not entirely true and the multisegmented line to gate into CD45+ cells is defined in this code.
#
#  WATCH OUT
#  late OTs (like OT17) have no longer CD62L_Lu175.tiff in!! you need to keep that empty _Lu175.tiff in the ROI folders or spectre list lengths of ROI-wise markers are not consistent!!
#   -> we patch that by not taking Lu175 into further transformation and scaling but leave it only in the raw signals. You could also delete CD62L_Lu175 along with x_Lu175, sure...
#  
#
#
#  changes in v1: - subset for the CD45 Gate and segmentation area penalties (permanent removal of these cells)
#                 - staging so that I can run the script in chunks and control what is handed over into next stage
#                 - add UMAP and Phenograph and return the data into cell.dat
#                 - make sure that when re-running cluster annotation, the column is deleted beforehand
#                 - adapt plotting parameters to fit our purpose
#
#  changes in v2: - include outlier columns to flag 1) CD45-gate penalty, 2) segmentation-area cells. Flagged cells remain in the cell.dat file so that:
#                 - option to bring back INS high cells back into the cell.dat file after clustering, marked with an own annotated.cluster name for spatial analysis
#                 - Quality control: check inside-gate but outside-area cells if they are truly rubbish
#                 - Quality control: have a look at the gate(s) in scatter plots to make sure the gates fit on all batches
#
#  changes in v3: - adapt the code for spectre package v1.0.0 
#                 - bring back tweaked functions from shinyIMC
#                 - start generalizing drawing scatter plots
#
#  changes in v4: - change the gating routine to hierarchical gating and clustering
#                 - change the string definition of high-pass gate (_Hgate) or low-pass gate (_Lgate) in metadata to allow automatic gating
#                 - change the dataframe handling during clustering the gated subsets
#
#  changes in v5: - use multisegmented line to gate into CD45+ cells
#                 - collapse the clusters according to dendrogram levels and check differeneces in WD vs ND
#                 - renumber the collapsed clusters so that there are no skipped numbers
#
#  changes in v6: - harmonize cluster collapse strategy and find a single number read-out for "cluster homogeneity"
#                 - clear out the comments for Laia
#
#  changes in v7: - automate cluster collapse routine
#                 - clear out old and experimental code
#
#  changes in v8: - move to complexheatmap engine for good, using z-scored data before plotting them
#                 - Use an engine that does transformation/normalization in one go: get rid of the long appended column names
#                 - move to z-scored transformed data for clustering for good
#                 - fix a bug in the collapse routine to allow no-collapse at all with level 0
#
#  changes in v9: - Follow roughly Lev's idea to include batch normalization for selected markers present on the pos. ctr. Check in STAGE 0, scale before asinh transf.
#                 - Switch from geom_density_ridges to facet_grid engine that is 30x faster in plotting
#                 - finish batchnormalization routine without piping it downstream into the gating routine (!!gating with raw, cluster with scaled data!!)
#
#  changes in v10:- Batch normalize routine creates a copy of cell.dat, scaled.cell.dat with the same folder names to pipe it downstream into gating without hickupps
#
#  changes in v11:- Added new STAGE for Proximity Maps based on EDM matrices created by imageJ
#                 - Clean up STAGES following collapse which were out-of-date
#                 - STAGE 2 can be re-run a second time to delete clusters from the subgates
#                 - STAGE 7 sources the EDM files and does proximity calcluations to these binary masks
#                 - rebuild.cell.data.anew=0 condition is now working again automatically, handling all the different files from different STAGEs
#
#  changes in v12:- clean up redundant variables
#                 - Update talk-back in the end of STAGEs
#                 - implement spillover correction via excel sheet in percent. Empty entries are filled in with 0.
#
#  changes in v13:- simplify batch normalization routine in STAGE 0 run-through 
#                 - homogenize polygon line design
#                 - adapt make.scatter.plot to deal with the case if a marker is used multiple times to gate
#
#  changes in v14:- clean up code in later STAGES, move lymphnode masking up in the code to allow clean analysis of on-ROI-compositions
#                 - STAGE 4: is now masking lymph nodes and calculating cluster-wise distances to EDM matrices
#                 - clean up all code titles
#                 - change imageJ lymphnode code and produce EDM: now we can extend the lymphnode mask by x pixel in case the mask did not fit well in certain ROIs
#  changes in v15:- beside inside.islet, we can now define also a second area, the rim around the mask to a defined pixel extension
#                 - update scatter.plot functions so that they automatically recognize if the gate flags are coming from a threshold or a multi-segmented line and set colors accordingly
#                 - fix the gating step plotting routine
#  changes in v16:- things get messy with v15.2 in BeCeFa that has partly updated v15 code in Laias Project. So as we run through M.Felsenstein Data, lets pool the both here
#  changes in v17:- clear out old code and re-start the analysis IMC-friendly way since the old paper got rejected. See notes in MB-2-67 (5th of Oct 23)
#  changes in v18:- moved over from v17_home in the Felsenstein Ilastik folder
#  changes in v19:- bring in Juliettes Code from v18
#                 - marker selection to process is no longer numbers but the marker names in spatial.dat 
#  - MOVE to PhenoCycler data and reset version -
#  changes in v0: - adapt column name handling
#                 - change how cluster column names and markercomparison list is built
#                 - QC include marker histograms 
#
#
#
###################################################################################

#rm(list=ls()) # if needed

########################## WORKFLOW GLOBAL CONTROL ##########################
#
#---- whichSTAGE are we running? ---
# 0   - read in data
#     - map metadata
#   -> Be aware that when you change metadata excel file, you will need start again at STAGE 0 to bring the new info into an unmapped "cell.dat"
#     - remove bad ROIs or even TMAs
#     - correct the raw data with a spillover correction matrix
#     - correct batch effect by normalizing a selection of markes on a positive control that is ablated on every ROI:
#                  - STAGE 0 needs to be run several times if you want batch normalization.
#                  - in the first run we just aim to connect metadata as to find the permanent column numbers needed to define "cellular.cols" in the end
#                  - the second run collects raw signal percentiles of each marker given in batch.normalization.matrix for batch normalization. 
#                    Check the outputs and decide on which percentile you want each marker to be normalized and put the number in "batch.normalization.matrix".
#                  - the third run applies the noted percentiles and normalizes the raw data accordingly 
#                    Even if you wrote percentiles into "batch.normalization.matrix" from the beginning, these numbers were not read in until third run of STAGE 0
#   -> define "cellular.cols" in the code: We need to drop markers that are not present on all ROIs, and chose all markers we want to proceed with
#
# 1   - flag segmentation area outlier according to metadata high/lowpass thresholds
#     - flag cells according to metadata high/lowpass thresholds or by supplying multi-segmented gate lines via "polygon.gates" list
#     - plot these gates for every TMA. Define which markers to plot via "markercomparision" list
#     - QC plots that will help finding the correct thresholds:
#                  - Cell segmentation area outliers per TMA as per thresholds read from metadata
#                  - All gates (metadata or polygon.gates) plotted on entire data of every TMA
#                  - All gates (metadata or polygon.gates) sequentially plotted for each hierarchical gating step
#                    These plots can be of use if the first subset blocks a clear view to set a proper threshold for a later gate
#       Adapt the numbers and re-run STAGE 0 to bring the new thresholds read from metadata!!
#     - Setup hierarchical gating with flags just written, create gated subsets
#     - Perform asinh transformation with the chosen cofactor followed by either min-max or z-score normalization
#  -> define "cluster.cols" in the code: select marker channel signals that go into clustering of every subset
#
# 2   - cluster and dim-reduction
#  -> Optionally: re-run with a non-empty "delete.clustersinsubsets" list to remove defined clusters for good. There is no way back.
#
# 3 Watch out, this stage will overwrite "gated.cell.dat" that entered STAGE2. Collapsed clusters are gone for good. There is no way back. 
#     - (sequentially) collapse clusters in each subset
#     - rename the cluster numbers into three-digit format and then:
#     - pool the subsets back into one "gated.cell.dat" object. 
#     - annotate cluster numbers where needed. Leave the other untouched and they will keep their three-digit number (yes, this is different to how spectre handles it)
# 
# 4   - Load in masks: binary to mask lymphnodes
#     - Distance Maps for spatial analysis
#
# 5 Watch out, is time consuming if executed on large data (>100 ROIs)! 
#     - spatial plotting of clusters (or subsets thereof) on all ROIs (or a subset of) 
# 
# 6   - project-specific plots
# 
# 7   - Proximity-index routine

# Start with STAGE 0 and increase one by one:
whichSTAGE <- 0









########################## STAGE 0 ##########################

if(whichSTAGE == 0){

#--------------------------------------
# Set Working Directory of the code
# and source the setting.
#--------------------------------------

dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory
source("global_settings.R") 
source("load_packages.R") # load all packages we need for the script to run
# now that we loaded spectre, lets make sure the default functions get overwritten with the adapted functions
source(paste0("utils_v",utils_version,".R")) # we need this line here when we run STAGE 0 the first time.
#--------------------------------------
# Set Output Directory Tree
#--------------------------------------
date <- gsub('-','',strsplit(x = as.character( Sys.time() ), split = ' ')[[1]][1])
setwd(PrimaryDirectory)
dir.create("output")
setwd("output")

dir.create(date)
setwd(date)
OutputDirectory <- getwd()

dir.create("Data")
setwd("Data")
OutputDataDirectory <- getwd()


  


if (!exists("how.often.ran.STAGE0")){
  #first time you run STAGE 0 create the counter
  how.often.ran.STAGE0<-0
}else{
  # after first run, increase counter
  how.often.ran.STAGE0<-how.often.ran.STAGE0+1
}


# Load Metadata file 
setwd(PrimaryDirectory)
sample.meta <- read_excel(metadatapath)

# check if tsv.file.name column ends with .tsv.
# if not, add it:
if(!any(grepl("\\.tsv$", sample.meta$tsv.file.name))){
  sample.meta$tsv.file.name <- paste0(sample.meta$tsv.file.name, ".tsv")
}

















if(rebuild.cell.data.anew == 1){


#--------------------------------------
# Set Input Directory
# if we rebuild anew, the input directory is the QuPath_tsv_inputs folder
#--------------------------------------
setwd(PrimaryDirectory)
setwd(paste0("QuPath_tsv_inputs/")) # all tsv files will load
InputDirectory <- getwd()
InputDirectory



#--------------------------------------
# load tsv files from QuPath Cellpose Segmentation
#--------------------------------------
setwd(InputDirectory)
# Merge all tsv files in all.tsv.files
all.tsv.files <- list.files(pattern = "\\.tsv$") 
# Read all files into a list of data frames
tsv_list <- lapply(seq_along(all.tsv.files), function(i) {
  df <- readr::read_tsv(all.tsv.files[i])
  df <- tibble::add_column(df, "tsv file" = all.tsv.files[i], .after = 2)
  df <- df %>% dplyr::rename_all(~ make.names(.))
  df
})

# Generalized code to compare all n list entries and print which names are different between any of the list entries
all_names <- lapply(tsv_list, names)
unique_names <- unique(unlist(all_names))
for (i in seq_along(tsv_list)) {
  missing_in_i <- setdiff(unique_names, names(tsv_list[[i]]))
  if (length(missing_in_i) > 0) {
    cat(sprintf("tsv_list[[%d]] (%s) is missing columns: %s\n", 
                i, all.tsv.files[i], paste(missing_in_i, collapse = ", ")))
  } else {
    cat(sprintf("tsv_list[[%d]] (%s) has all columns present.\n", 
                i, all.tsv.files[i]))
  }
}


# Append all data frames row-wise
PCF_data <- dplyr::bind_rows(tsv_list)
rm(tsv_list)


#--------------------------------------
# Prepare the marker names  that we process
#--------------------------------------
# before we even do anything, 
# lets make sure there is no typo in the Markers2Process vector:

# lets assemble the column names of the markers we will process.
# we do not just asseble the columns based on Markers2Process, but
# we attempt to pull these assebled columns from PCF_data...
cellular.cols <- names(
PCF_data %>%
  dplyr::select(
    starts_with(
      paste0(
        Markers2Process,
        qupath.separator,
        cellular.segment.compartment,
        qupath.separator,
        cellular.segment.readout
      )
    )
  )
)

# ...and then can check if all markers are present in the data:
if(length(Markers2Process) != length(cellular.cols)){
  cat(paste0("\n\n\n----------------------------------- \n ERROR in STAGE0 start-up   \n 
  Not all markers to process further were found in your data \n  
  Check your Markers2Process vector and make sure it matches the marker spelling \n 
  -> Seems you misspelled ",length(Markers2Process)-length(cellular.cols)," of the following entries? \n")) 
  temp.string <- names(
  PCF_data %>%
    dplyr::select(
      -starts_with(
        paste0(
          Markers2Process,
          qupath.separator,
          cellular.segment.compartment,
          qupath.separator,
          cellular.segment.readout
        )
      )
    ) %>%
    # lets remove all metadata columns we just added
    dplyr::select(-starts_with(names(sample.meta))) %>%
    # lets also remove the QuPath-specific columns
    dplyr::select(
      -starts_with(
        c(
          "Image", "Object.ID", "Object.type", "Name", "Classification",
          "TMA.core", "Parent", "ROI", "Centroid.X.µm", "Centroid.Y.µm",
          "tsv.file", "Nucleus.Cell.area.ratio"
        )
      )
    ) %>%
    dplyr::select(-ends_with(c("Std.Dev.", "Variance")))
)

  temp.string <- unique( sub("\\.\\..*", "", temp.string) )
  print(  as.matrix( temp.string  ))

  cat("-----------------------------------\n\n")
  stop()
}


#--------------------------------------
# Tidy up PCF_data before proceeding
#--------------------------------------
# the groovy script that runs CellPose in QuPath is heavily adapted
# consequently, the column names in tsv files need clean up:

# 1) due to assembling the cell cytoplasm and nucleus manually, 
# TMA.column is empty and TMA info added into "Name" colum
# fix that:
PCF_data <- PCF_data %>%
  dplyr::select(-TMA.core) %>%
  dplyr::rename(TMA.core= Name)

# 2) while assembling the cells, we also added non-nucleated cell profiles
# if you did NOT do a object classification in QuPath, 
# "Classification" column is still populated with
# this information: "nucleated profile" and "anucleated profile". 
# if you DO classify the object, these tags were overwritten by the Object Classifier of QuPath.
# but dont worry, we bring back the nucleation status by checking the intensities of Cytoplasm:
# anucleated profiles are composed of nuc and cyto being the same polygon
# so when shape analyis ran, the Cytoplasm reads became NA in these cells
# we take advantage of this fact and create a new column called "Object.cross.section":

# Check and process the Classification column in PCF_data
if ("Classification" %in% names(PCF_data)) {
  unique_vals <- unique(na.omit(PCF_data$Classification))
  allowed_vals <- c("nucleated profile", "anucleated profile")
  if (all(unique_vals %in% allowed_vals) && length(unique_vals) <= 2) {
    # No Object classification ran, profile annotations are conserved 
    # only the column has wrong name:
    PCF_data <- PCF_data %>%
      dplyr::rename(Object.cross.section = Classification)
  } else {
    # Other values present because Object Classifier ran.
    # create Object.cross.section based on ..Cytoplasm.. columns
    cytoplasm_cols <-  grep(paste0(qupath.separator,"Cytoplasm",qupath.separator)   , names(PCF_data), value = TRUE)
    PCF_data$Object.cross.section <- ifelse(
      apply(PCF_data[, cytoplasm_cols, drop = FALSE], 1, function(row) any(is.na(row))), #dataframe
      #apply(PCF_data[, ..cytoplasm_cols], 1, function(row) any(is.na(row))), #data.table
      "anucleated profile", "nucleated profile"
    )
  }
}


# 3) Something in the groovy script, or the export itself creates a bunch of ghost cells without any reads
# these cells carry a centroid and Object.type Cell, they are assigned to TMA.cores,
# I checked the CEntroid.X.µm and Centroid.Y.µm of a couple of these cells, none of them is drawn in QuPath 
# They even carry Clasisfication, but their marker values are filled with NA values.
#cells_with_all_na <- PCF_data[apply(PCF_data[, ..cellular.cols], 1, function(row) all(is.na(row))), ]
# we remove them if all markers to process are NA:
cells_with_all_na <- PCF_data[apply(PCF_data[, cellular.cols, drop = FALSE], 1, function(row) all(is.na(row))), ]

if (ncol(cells_with_all_na) > 0) {
  message(paste0(ncol(cells_with_all_na), " cells were found that have all NA values in the Markers2Process columns. Removing them..."))
  print(head(cells_with_all_na))
  # Remove cells where all Markers2Process columns are NA
#PCF_data <- PCF_data[!apply(PCF_data[, ..cellular.cols], 1, function(row) all(is.na(row)))]
PCF_data <- PCF_data[!apply(PCF_data[, cellular.cols], 1, function(row) all(is.na(row))), ]
}
rm(cells_with_all_na)


#--------------------------------------
# Mapping Metadata to the cellular data
#--------------------------------------

# check if any of the tsv files has NA in the TMA.core column:
if(any(is.na(PCF_data$TMA.core))){
  message("TMA.core column in PCF_data has NA values")
  message("-> Assuming there was no TMA present and mapping metadata by tsv file names:")

  # remove the TMA.core column from sample.meta. 
  #We got this column already in PCFdata and would not be able to map 
  #if more than one column name are shared among the two files:
  sample.meta <- sample.meta %>% dplyr::select(-TMA.core)

  #ready to map meta into df
  PCF_data <- do.add.cols(as.data.frame(PCF_data), 
    base.col = 'tsv.file', 
    add.dat = sample.meta, 
    add.by = 'tsv.file.name'
  )
  message("Mapping metadata to the cellular data complete")

}else {

  # normal case: we have TMA.core column in PCF_data
  message("-> TMA core IDs found, mapping each core:")
  # sample.meta tsv.file.name column is already set, so is PCF_data$tsv.file

  # We need a unique TMA - core ID to tell the cores apart:
  # paste TMA.ID and TMA.core into a new column called TMA.ID_core 
  # for PCF_data this is straightforward:.after
  PCF_data$TMA.ID_core <- paste(PCF_data$tsv.file,PCF_data$TMA.core, sep="_")
  # in sample.meta, lets make sure the TMA.core column comes properly formatted.
  # we need A-1, not A1 oder A.1 or whatsoever.

  if ("TMA.core" %in% names(sample.meta)) {
    # Regular expression to match a single letter, a hyphen, and another letter
    if (all(grepl("^[A-Za-z]-[0-9]$", sample.meta$TMA.core))) {
      # all good, prepare the column to map:
      sample.meta$TMA.ID_core <- paste(sample.meta$tsv.file.name, sample.meta$TMA.core, sep = "_")
      # remove the TMA.core column from sample.meta for upcoming mapping:
      sample.meta <- sample.meta %>% dplyr::select(-TMA.core)
      # ready to map meta into df
      PCF_data <- do.add.cols(as.data.frame(PCF_data), base.col = 'TMA.ID_core', add.dat = sample.meta, add.by = 'TMA.ID_core')
      message("Mapping metadata to the cellular data complete")
    } else {
      # lets attempt to fix TMA.core column, maybe the column is in A1 format:
      sample.meta$TMA.core <- gsub("([A-Za-z]+)([0-9]+)", "\\1-\\2", sample.meta$TMA.core)
      # Regular expression to match a single letter, a hyphen, and another letter
      if (all(grepl("^[A-Za-z]-[0-9]$", sample.meta$TMA.core))) {
        # fixing worked. proceed with mapping:
        # prepare the column to map:
        sample.meta$TMA.ID_core <- paste(sample.meta$tsv.file.name, sample.meta$TMA.core, sep = "_")
        sample.meta <- sample.meta %>% dplyr::select(-TMA.core)
        PCF_data <- do.add.cols(as.data.frame(PCF_data), base.col = 'TMA.ID_core', add.dat = sample.meta, add.by = 'TMA.ID_core')
        message("Fixed the TMA.core format. Mapping metadata to the cellular data complete")
      }
    }
  } else {
    message("Mapping metadata failed. Make sure your metadata contains TMA.core and tsv.file.name matches the file names in the QuPath_tsv_inputs folder.")
    stop()
  }
  # end mapping: PCF_data$TMA.core was there.

# now that we mapped batch numbers, we could clear entire cores:
#e.g. PCF_data <- PCF_data[ grep("20250228_27X_bleaching trial (control).tsv_B-2",  PCF_data$TMA.ID_core, invert = TRUE) , ]


# We will now switch to the spectre ecosystem:
# start off by renaming to cell.dat:
cell.dat <- PCF_data
setwd(OutputDataDirectory)
qsave(cell.dat, "unscaled.cell.dat.qs")
rm(PCF_data)



# all settings done. 
# initialize following stage 0 varibles :
# all_images, amount.of.TMAs, all_cores, amount.of.cores, markercomparision
list2env( initialize_stage_0_variables(
  cell.dat = cell.dat,
  scatter.plots = scatter.plots,
  transform.rawdata = transform.rawdata,
  qupath.separator = qupath.separator,
  cellular.segment.compartment = cellular.segment.compartment,
  cellular.segment.readout = cellular.segment.readout,
  asinh.cofactor = asinh.cofactor
), envir = .GlobalEnv)




###### Start spill-over correction routine ######         
if(correct.spill.over ==1){
          
  # Watch out! This routine will be executed exactly one time!
    cell.dat <- do.spill.over.corr(
    dat = cell.dat,
    spillover.matrix.path = spillover.matrix.path,
    donor.col.name = "acq.chnl" # the first column name with the donor isotopes
  )
}
          

# STAGE 0 needs to be run several times if you want batch normalization.
# in the first run we just aim to connect metadata as to find the permanent column numbers needed to define cellular.cols in the end
# the second run collects percentiles for batch normalization
# the third run applies the preset percentiles and normalizes the raw data
        
# normalization is not needed for the first pass
if(how.often.ran.STAGE0>0) {        
        
                
        
  ###### Start batch normalization routine ######      
  if(calculate.scaling.factors==1){

    batchnormalize.marker <- batch.normalization.matrix[seq(1,length(batch.normalization.matrix),2)]
    batchnormalize.percentile <- as.numeric(batch.normalization.matrix[seq(2,length(batch.normalization.matrix),2) ])

    set.seed(723451) # for reproducibility
    percentile_palette <- createPalette(length(batch.aligning.percentiles), c("#1BC912"), M=100000, prefix = "")
    names(percentile_palette) <- batch.aligning.percentiles   

    # Compared to Lev, we hop into exploration and normalization automatically.
    # how.often.ran.STAGE0 is now either 1 or 2

    if(  how.often.ran.STAGE0==1 ) { 
      # we need to ignore the percentiles in batch.normalization.matrix
      message("Batch normalization routine started: collecting percentiles of markers to align the signals...")
      # if we run batch normalization, we need more packages:
      library(flowCore)
      library(cowplot)
      library(ggridges)
      library(DescTools) # call the function "Closest" for scaling  

      quantiles_table_long <- data.frame()

      setwd(OutputDirectory)
      dir.create('Batch_Normalization')
      setwd('Batch_Normalization') 

      # in order to track the histogram lines from within ggplot, we need to assign unique colors pre-assigned to each batch:
      set.seed(723451) # for reproducibility
      Batches_palette <- createPalette(amount.of.TMAs, c("#0000ff"), M=100000, prefix = "")
      names(Batches_palette) <- all_images  

      # its easier to just subset a pos. ctr. dataset out now:  
      spleen.data.frame <- subset(cell.dat, Tissue %in% postive.control.Tissue.name )  

      # we silently do the transformation for the spleen as well,
      # but we brutally override whatever normalization was chosen:
      spleen.data.frame <- do.data.normalization(
        dat = spleen.data.frame,
        use.cols = cellular.cols,
        # Transform
        do.transform = transform.rawdata,
        cofactor = asinh.cofactor,
        # min-max normalize?
        do.minmax = FALSE, # minmax.norm,
        new.min = min.norm,
        new.max = max.norm,
        # z-score normalize?
        do.zscore = FALSE #zscore.norm
      )

      # from here on, if we plot, we need the asinh transformed columns:
      # which are pasted like this in the function: paste0(names(value), "_t",cofactor)

      batchnormalize.marker.transformed <- paste0(batchnormalize.marker, "_t",asinh.cofactor)

      # pushing the percentiles into a separate df or list is a cute idea, but the ggplot engine wont work like that:
      # ggplot wants the transformed data and the percentile values in the same df, and that plotting takes ages

      pb <- progress_bar$new(format = "[:bar] :percent [Collecting Percentiles | :eta]",
        total = amount.of.TMAs*length(batchnormalize.marker.transformed), #
        show_after=0, #allows to call it right way 
        current = "|",    # Current bar character
        clear = FALSE) # make it persist
      pb$tick(0) # call in the progress bar without any progress, just to show it

      # We push the percentiles into two datasets:
      for (batch in all_images ){

        #1: make an own percentile df that we use to pull out the values when calculating the scaling factors in a second run through:
        # this here is calculating the percentiles of all channels of the current TMA, while......
        temp.df <- as.data.frame( 
          cbind( 
            sapply(  subset(  spleen.data.frame, Batch==batch, select=batchnormalize.marker.transformed) , 
              function(x) quantile(x, probs=batch.aligning.percentiles )  ) , 
            Percentile=batch.aligning.percentiles,
            Batch=batch 
          )
        )

        # and then bind this baby to the 
        quantiles_table_long <- rbind(quantiles_table_long, temp.df )

        for (channel in batchnormalize.marker.transformed){

          # ......this here subsets into batch and channel. We need this second percentile calculation (not as list, but as stupid vector) for the ggplot routine:

          temp.signal <- as.vector(subset(spleen.data.frame, Batch==batch, select=channel))

          temp.percentiles <-  sapply( temp.signal , 
            function(x) quantile(x, probs=batch.aligning.percentiles )  ) 

          # 2: push the percentiles into the raw data for the upcoming plotting routine:
          spleen.data.frame$p60[spleen.data.frame$Batch==batch] <- temp.percentiles[1]
          spleen.data.frame$p80[spleen.data.frame$Batch==batch] <- temp.percentiles[2]
          spleen.data.frame$p85[spleen.data.frame$Batch==batch] <- temp.percentiles[3]
          spleen.data.frame$p90[spleen.data.frame$Batch==batch] <- temp.percentiles[4]
          spleen.data.frame$p95[spleen.data.frame$Batch==batch] <- temp.percentiles[5]
          spleen.data.frame$p99[spleen.data.frame$Batch==batch] <- temp.percentiles[6]

          pb$tick()

        }# end run along all channels and fill up the spleen.data.frame with the according percentiles
      }# end run along all batches to calculated percentiles of current channel

      # now that everything is stored in the raw data dataframe, now we plot, thats faster, I thinK:
      message("Plotting collected percentiles on all batches...")

      pb <- progress_bar$new(format = "[:bar] :percent [Plotting unscaled | :eta]",
        total = length(batchnormalize.marker.transformed), #
        show_after=0, #allows to call it right way 
        current = "|",    # Current bar character
        clear = FALSE) # make it persist
      pb$tick(0) # call in the progress bar without any progress, just to show it

      #https://stackoverflow.com/questions/71555598/can-geom-vline-be-connected-across-facet-grid

      for (channel in batchnormalize.marker.transformed){

        # we first create a tibble that is split by batch, the same we use for the facet. this allows to split the percentiles batch-wise later on: 
        percentiles.tibble.facetplottin <- spleen.data.frame %>%
          group_by(Batch) %>%
          summarize(p60 = unique(p60),
            p80 = unique(p80),
            p85 = unique(p85),
            p90 = unique(p90),
            p95 = unique(p95),
            p99 = unique(p99)
          )

        png(filename=paste0( sub("_.*", "", channel ),"_t",asinh.cofactor,"_0_unscaled.png"), width = 1600, height=length(unique(spleen.data.frame$Batch))*100)
        suppressMessages(
          print(  
            ggplot(spleen.data.frame  , aes( x= spleen.data.frame[[channel]] )  ) +
              facet_grid(Batch ~ . , scales = "free_y", switch = "y")+ # this brings the facet tag to the left
              geom_area( alpha = 0.15 ,   color = 'black', fill = 'grey',
                stat = "bin",
                binwidth = function(x) 2 * IQR(x) / (length(x)^(1/3)), 
                size = 0.5    )+
              #scale_color_manual(values = scico(amount.of.TMAs, begin = 0.12, palette = "roma"))+
              #scale_fill_manual(values = scico(amount.of.TMAs, begin = 0.12, palette = "roma"))+
              # geom_vline(data=percentiles.tibble.facetplottin, aes(xintercept = p60), colour = "#a361c7"   ) +
              #  geom_vline(data=percentiles.tibble.facetplottin, aes(xintercept = p80), colour = "#5ba962"   ) +
              #  geom_vline(data=percentiles.tibble.facetplottin, aes(xintercept = p85), colour = "#c75a87"   ) +
              #  geom_vline(data=percentiles.tibble.facetplottin, aes(xintercept = p90), colour = "#ab973d"   ) +
              #  geom_vline(data=percentiles.tibble.facetplottin, aes(xintercept = p95), colour = "#648ace"   ) +
              #  geom_vline(data=percentiles.tibble.facetplottin, aes(xintercept = p99), colour = "#cb6342"   )+
              geom_vline(data=percentiles.tibble.facetplottin, aes(xintercept = p60), colour = percentile_palette["0.6"]   ) +
              geom_vline(data=percentiles.tibble.facetplottin, aes(xintercept = p80), colour = percentile_palette["0.8"]   ) +
              geom_vline(data=percentiles.tibble.facetplottin, aes(xintercept = p85), colour = percentile_palette["0.85"]   ) +
              geom_vline(data=percentiles.tibble.facetplottin, aes(xintercept = p90), colour = percentile_palette["0.9"]   ) +
              geom_vline(data=percentiles.tibble.facetplottin, aes(xintercept = p95), colour = percentile_palette["0.95"]   ) +
              geom_vline(data=percentiles.tibble.facetplottin, aes(xintercept = p99), colour = percentile_palette["0.99"]   )+
              #geom_point(data=percentiles.tibble.facetplottin, mapping=aes(x=p60, y=0, group = 1), colour = "#a361c7"   ) +
              #geom_point(data=percentiles.tibble.facetplottin, mapping=aes(x=p80, y=0, group = 1), colour = "#5ba962"   ) +
              #geom_point(data=percentiles.tibble.facetplottin, mapping=aes(x=p85, y=0, group = 1), colour = "#c75a87"   ) +
              #geom_point(data=percentiles.tibble.facetplottin, mapping=aes(x=p90, y=0, group = 1), colour = "#ab973d"   ) +
              #geom_point(data=percentiles.tibble.facetplottin, mapping=aes(x=p95, y=0, group = 1), colour = "#648ace"   ) +
              #geom_point(data=percentiles.tibble.facetplottin, mapping=aes(x=p99, y=0, group = 1), color= "#cb6342" )+
              geom_point(data=percentiles.tibble.facetplottin, mapping=aes(x=p60, y=0, group = 1), colour = percentile_palette["0.6"]   ) +
              geom_point(data=percentiles.tibble.facetplottin, mapping=aes(x=p80, y=0, group = 1), colour = percentile_palette["0.8"]   ) +
              geom_point(data=percentiles.tibble.facetplottin, mapping=aes(x=p85, y=0, group = 1), colour = percentile_palette["0.85"]   ) +
              geom_point(data=percentiles.tibble.facetplottin, mapping=aes(x=p90, y=0, group = 1), colour = percentile_palette["0.9"]   ) +
              geom_point(data=percentiles.tibble.facetplottin, mapping=aes(x=p95, y=0, group = 1), colour = percentile_palette["0.95"]   ) +
              geom_point(data=percentiles.tibble.facetplottin, mapping=aes(x=p99, y=0, group = 1), color= percentile_palette["0.99"] )+
              theme_minimal() +
              theme(axis.title.y = element_blank(), #take "count" away
                #these are the facet settings:
                strip.text.x = element_text(size=25,  face="bold"),
                strip.background = element_rect(colour="black", fill="grey90", size=1, linetype="solid"),
                strip.text.y.left = element_text(angle = 0), # and once the tag is on the left, keep it horizontal
                panel.border=element_blank(), 
                axis.line=element_line(), # turn off the x-axis line on the bottom of the facet
                axis.title.x = element_text(size=20),
                axis.text.x = element_text(size=20),
                axis.text.y = element_text(size=10),
                panel.spacing = unit(1.2, "lines")   
              )+
              xlab(as.character(channel))
          )#end print
        )
        invisible(dev.off())
        pb$tick()
      }# end run along all channels and plot ridges

      # end second time pass STAGE0: Percentiles were ignored check which percentiles to align 
      # to switch to final STAGE0 setup, we do not use how.often.ran.STAGE0==2 but force the script into that setting without any chance of getting out
    } else if(  how.often.ran.STAGE0>1 ) {
      # begin second pass: use the provided percentiles given in batch.normalization.matrix to calculate scaling factors:
      setwd(OutputDirectory)
      setwd('Batch_Normalization')    

      # watch out, make sure the last column of batchnormalize.percentile is actually filled with numbers:
      if(any(is.na(batchnormalize.percentile )) == F ){

        message("Batch normalization routine continued: percentiles found, calculating scaling factors...")  

        scaled.cell.dat <- cell.dat

        scaled.spleen.data.frame <- subset(cell.dat, Tissue %in% postive.control.Tissue.name )  # this didnt do much yet. we gonna push scaled "_SF" cols in here and transform
        scaled.spleen.data.frame <- do.data.normalization(
          dat = scaled.spleen.data.frame,
          use.cols = cellular.cols,
          # Transform
          do.transform = transform.rawdata,
          cofactor = asinh.cofactor,
          # min-max normalize?
          do.minmax = FALSE, # minmax.norm,
          new.min = min.norm,
          new.max = max.norm,
          # z-score normalize?
          do.zscore = FALSE #zscore.norm
        )

        # the trick is that quantiles_table_long was rbound running through "batch in all_images". 
        # therefore, if we pull out temp.quantile.values in the upcoming loop, these values have the same order as all_images
        # so there is no need to cbind the batches everytime we build up the scaling.factor dataframe: the scaling factors follow the same order as in all_images:
        scaling.factors <- data.frame( Batch.ID = all_images   )

        # so we got percentiles set, lets see them. instead of running througth the names, lets run through the position of both vectors:

        # lets run through the markers again and calculate the SF, given the percentile we find in the matrix (or rather, the batchnormalize.percentile vector derived from that matrix )

        for(m in c(1:length(batchnormalize.percentile ))){ 

          #  batchnormalize.marker[m]  # watch out, thats raw signal, this is not part of  quantiles_table_long
          #  batchnormalize.marker.transformed[m]  # thats the one you wanna work on first for the SF
          #  batchnormalize.percentile[m] 

          # its important that we pull the TRANSFORMED Data here! We will calculate the SF on these data, but scale then the linear data!
          # temp.quantile.values <- deframe( subset(quantiles_tibble_long, Percentile==batchnormalize.percentile[m]  )[,batchnormalize.marker.transformed[m]]   )
          temp.quantile.values <- as.numeric( subset(quantiles_table_long, Percentile==batchnormalize.percentile[m]  )[,batchnormalize.marker.transformed[m] ] )  

          # find the closest value to the mean of all percentiles
          #https://stackoverflow.com/questions/43472234/fastest-way-to-find-nearest-value-in-vector
          anchor <- as.vector( Closest(x = temp.quantile.values, a = mean(temp.quantile.values)     ) )

          # Somesh warned me that if one has only two OTs to align, anchor will have two objects now and we gonna crash the script in the upcoming scaling.
          # so lets just make sure we pull the smaller of both values:
          anchor <- sort(anchor)[1]

          # match this anchor value to the OT nth memeber in temp.qunantiles (we do not )
          #match( anchor,  temp.quantile.values   )

          message( paste0("Scaling ", batchnormalize.marker[m], 
            " @ ", 100*batchnormalize.percentile[m], "th percentile: Closest to percentile-mean: ",  
            subset(quantiles_table_long, Percentile==batchnormalize.percentile[m]  )$Batch[ match( anchor,  temp.quantile.values   ) ])      )   

          # now we just use the quantiles of the transformed signals from above to calculate the scale factors:
          # rather than storing them all in separate df, we create a temporal one and then....
          temp.SF.values <- data.frame(
            #Batch = subset(quantiles_tibble_long, Percentile==batchnormalize.percentile[m]  )$Batch,
            Batch = subset(quantiles_table_long, Percentile %in% batchnormalize.percentile[m]  )$Batch, # this is equivalent with all_images, but I wanna be sure the order is the same
            SF = rep(anchor, length(temp.quantile.values)) /temp.quantile.values 
            #    ^ this rep() here is stupid, but without we would run into a "shorter vector is not a multiple of longer vector" warning which I hate to see in the console
          ) 

          # batch order and SF.values is always the same. I cross-checked that but see the comment when scaling.factors is initialized
          scaling.factors[, paste0( batchnormalize.marker[m] )] <- temp.SF.values[,2]

          # ...and then use classic dirty R to push a new column into the spleen.data.frame right away
          # we use the long spleen.data.frame$Batch vector to match every according SF in the temp.SF.values in one go, and use that scaling to create a new columN:

          temp.scaled.rawsignal <- spleen.data.frame[[ batchnormalize.marker[m] ]] * temp.SF.values[ match(spleen.data.frame$Batch, temp.SF.values$Batch   ) ,  ]$SF

          scaled.spleen.data.frame[, paste0( batchnormalize.marker[m], "_SF" )] <- temp.scaled.rawsignal

          #  WATCH OUT 
          # !!!! we are about to use the same column names like our raw data has!!!!  
          # as of v10, we do the scaling on a copy of cell.dat: scaled.cell.dat,
          # and we override the original raw data and push the scaled data into these columns!
          #  this has the advantage that we can feed that object to downstream engines (gating and plotting mostly) without the need to change their input column name structure

          # this time with use the million-cell long scaled.cell.dat$Batch vector to match the SF present in temp.SF.values:
          temp.scaled.rawsignal <- scaled.cell.dat[[ batchnormalize.marker[m] ]] * temp.SF.values[ match(scaled.cell.dat$Batch, temp.SF.values$Batch   ) ,  ]$SF
          scaled.cell.dat[, batchnormalize.marker[m] ] <- temp.scaled.rawsignal

          # quantiles_tibble_long <- as_tibble(quantiles_table_long) %>% type_convert()
          # filter(test, Percentile  == 0.95) 

          #library(tidyverse)

          #filter(quantiles_tibble_long, Percentile  == 0.95)

          #test %>% filter(test, Percentile  == 0.95)

          #quantiles_tibble_long %>% 
          #  filter(`Percentile`  = 0.95)

        } # end running with m along all set percentiles and markers

        # again, we block any min-max or zscale setting manually:
        scaled.spleen.data.frame <- do.data.normalization(
          dat = scaled.spleen.data.frame,
          use.cols = paste0(batchnormalize.marker, "_SF"),
          # Transform
          do.transform = transform.rawdata,
          cofactor = asinh.cofactor,
          # min-max normalize?
          do.minmax = FALSE, # minmax.norm,
          new.min = min.norm,
          new.max = max.norm,
          # z-score normalize?
          do.zscore = FALSE #zscore.norm
        )

        # And memorize which ones did not go into the scaling routine  
        # we set the difference based on the raw signal names, so we have to paste the processed tags in both cases: 
        # scaled columns get a _SF, the ones left out get their transformed tag back on:
        scaled.cellular.cols <- c( paste0(batchnormalize.marker, "_SF", "_t",asinh.cofactor) , paste0(setdiff(cellular.cols,batchnormalize.marker), "_t",asinh.cofactor) )       

        # and plot with the selected percentile:
        message("Plotting all scaled batches...")
        pb <- progress_bar$new(format = "[:bar] :percent [Plotting scaled | :eta]",
          total = length(scaled.cellular.cols), #
          show_after=0, #allows to call it right way 
          current = "|",    # Current bar character
          clear = FALSE) # make it persist
        pb$tick(0) # call in the progress bar without any progress, just to show it

        for (channel in scaled.cellular.cols){   

          png(filename=paste0( sub("_.*", "", channel ),"_t",asinh.cofactor,"_1_scaled.png"), width = 1600, height=length(unique(spleen.data.frame$Batch))*100)
          suppressMessages(
            print(  
              ggplot(scaled.spleen.data.frame  , aes( x= scaled.spleen.data.frame[[channel]] )  ) +
                facet_grid(Batch ~ . , scales = "free_y", switch = "y")+ # this brings the facet tag to the left
                geom_area( alpha = 0.15 ,   color = 'black', fill = 'grey',
                  stat = "bin",
                  binwidth = function(x) 2 * IQR(x) / (length(x)^(1/3)), 
                  size = 0.5    )+
                theme_minimal() +
                theme(axis.title.y = element_blank(), #take "count" away
                  #these are the facet settings:
                  strip.text.x = element_text(size=15,  face="bold"),
                  strip.background = element_rect(colour="black", fill="grey90", size=1, linetype="solid"),
                  strip.text.y.left = element_text(angle = 0), # and once the tag is on the left, keep it horizontal
                  panel.border=element_blank(), 
                  axis.line=element_line(), # turn off the x-axis line on the bottom of the facet
                  axis.title.x = element_text(size=20),
                  axis.text.x = element_text(size=20),
                  axis.text.y = element_text(size=10),
                  panel.spacing = unit(1.2, "lines")
                )+
                xlab(as.character(channel))
              #    .............
            )#end print
          )
          invisible(dev.off())
          pb$tick()  
        }# end run along all channels and plot ridges of scaled data

        # let us store away scaled.cell.dat
        setwd(OutputDataDirectory)
        qsave(scaled.cell.dat, "scaled.cell.dat.qs")
      }else{#protect from batch normalizing if at least one marker percentile is not defined NA stop("Transformation/Normalization engine stopped")  
        cat("\n\n\n----------------------------------- \n ERROR in STAGE 0 third run.   \n At least one channel has undefined percentiles to normalize on \n  1) Check your batch.normalization.matrix \n -> Re-run STAGE 0 to apply batch norm. \n-----------------------------------\n\n")# stop(), call. = FALSE)  
        stop()
      }
    }# end second pass bracket: percentiles are given in batch.normalization.matrix 

    #protect from running Levs code that is in adaptation
    if(run.experimental.code==1){

      # set up a tibble where we store all max values of every marker we will normalize
      spleen.normalization <- tibble(marker=character(), 
        batch=character(), 
        max.signal=numeric()
      )

      pb <- progress_bar$new(total = amount.of.TMAs,
        show_after=0, #allows to call it right way 
        clear = FALSE) # make it persist
      pb$tick(0) # call in the progress bar without any progress, just to show it

      # for (channel in batchnormalize.marker){
      for (batch in all_images ){

        # quantiles.list[[list.pos]] <-   t(   sapply(  subset(  spleen.data.frame, Batch==batch, select=batchnormalize.marker) , 
        #                                          function(x) quantile(x, probs=batch.correct.percentile )  )
        #      )

        temp.df <- as.data.frame( 
          cbind( 
            sapply(  subset(  spleen.data.frame, Batch==batch, select=batchnormalize.marker.transformed) , 
              function(x) quantile(x, probs=batch.correct.percentile )  ) , 
            Percentile=batch.correct.percentile,
            Batch=batch 
          )
        )

        # and then bind this baby to the 
        quantiles_table_long <- rbind(quantiles_table_long, temp.df )

        # list.pos <- list.pos+1
        pb$tick()
      }

      # names(quantiles.list) <- all_images 

      # now some tibble:
      quantiles.df.long <- dplyr::bind_rows(quantiles.list, .id = 'id')

      # now quick and dirty, cbind the list and add a batch name column:
      quantiles.df.long <-  do.call(rbind, quantiles.list) 
      quantiles.df.long <- cbind(quantiles.df.long, Batch= seq(1: (amount.of.TMAs* length(batchnormalize.marker)) )  ) 

      quantiles.df.long <- cbind(quantiles.df.long,  rep(all_images, each= length(batchnormalize.marker) ) )

      library(purrr)
      map_df(quantiles.list, .id="id")

      do.call(rbind, quantiles.list)
      unsplit(quantiles.list, quantiles.list$Batch)

      dplyr::bind_rows(quantiles.list, .id = 'id')

      data.table::rbindlist(quantiles.list, idcol = 'id')
      dplyr::bind_rows(unname(quantiles.list), .id = 'id')

      ggplot(spleen.data.frame, aes(x=spleen.data.frame[[m]], y=spleen.data.frame[,Batch]))+
        geom_density_ridges(fill='grey', scale = 1, alpha = 0.4, quantile_lines = T, quantiles = c(0.6, 0.8, 0.85, 0.90, 0.95, 0.99))+
        theme_cowplot()+
        theme(axis.title.x = element_text(size=35),
          axis.text.x = element_text(size=30),
          axis.text.y = element_text(size=20)
        )+
        xlab(as.character(m))

      # set up lyrics and progress bar for the GAM modelling:

      message( paste0("Starting batch normalization of ", length(batchnormalize.marker), " markers on ", amount.of.batches, " TMAs" )  )
      pb <- progress_bar$new(total = length(batchnormalize.marker)*amount.of.batches,
        show_after=0, #allows to call it right way 
        clear = FALSE) # make it persist
      pb$tick(0) # call in the progress bar without any progress, just to show it

      for(m in batchnormalize.marker ){

        # we make the pretty plot first and then do the work:
        #g <-
        ggplot(spleen.data.frame  , aes( x= spleen.data.frame[[m]] , fill = Batch , color = Batch )  ) +
          facet_grid(Batch ~ . , scales = "free_y")+ theme(panel.border=element_blank(), axis.line=element_line())+
          geom_area( alpha = 0.15 ,
            stat = "bin",
            binwidth = function(x) 2 * IQR(x) / (length(x)^(1/3)), 
            size = 0.5    )+
          scale_color_manual(values = scico(amount.of.TMAs, begin = 0.12, palette = "roma"))+
          scale_fill_manual(values = scico(amount.of.TMAs, begin = 0.12, palette = "roma"))+
          # scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),  labels = trans_format("log10", math_format(10^.x))) +
          # annotation_logticks(sides="b")   +
          theme_classic()+
          labs(title= paste0(sub("_.*", "", m ), " signal drift over all batches"  )    , 
            subtitle= paste0("Using ", nrow(spleen.data.frame), " cells of all spleen ROIs.") ,
            #caption="Created by M.Barone", 
            y="Amount of cells in bin", 
            x= sub("_.*", "", m )  
          )+
          theme(panel.background = element_rect(fill = "grey93", colour = "grey55", size = 0.5), # change 'colour' to black for informative axis
            axis.title.x=element_text(color="grey15", size=11),
            axis.title.y=element_text(color="grey15", size=11),
            axis.text=element_text(size=10),
            #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
            legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
            legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
            #legend.title=element_blank(),
            plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
            plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
          )

        filename <- paste0("SPLEEN_RawSignalDrift_byCondition_", sub("_.*", "", m ) , ".png")
        ggsave(filename,g, width = 30, height = 20, units = "cm",dpi = 1200)

        # now we gonna do the same, but this time we supply the unique colors to the histogram plotting engine: Batches_palette 
        # this looks terrible, but we can then track back the batch name via color code:
        g <-ggplot(spleen.data.frame  , aes( x= spleen.data.frame[[m]] , fill = Batch , color = Batch )  ) +
          geom_area( alpha = 0.15 ,
            stat = "bin",
            binwidth = function(x) 2 * IQR(x) / (length(x)^(1/3)), 
            size = 0.5    )+
          #    # scico::scale_fill_scico(palette = "roma", direction=-1)+ # bilbao
          #    scale_color_manual(values = scico(amount.of.TMAs, begin = 0.12, palette = "roma"))+
          #    scale_fill_manual(values = scico(amount.of.TMAs, begin = 0.12, palette = "roma"))+
          #    #facet_wrap(~Batch, scales = "free_y")+
          scale_color_manual(values = Batches_palette)+
          scale_fill_manual(values = Batches_palette)+  
          scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),  labels = trans_format("log10", math_format(10^.x))) +
          # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),         labels = trans_format("log10", math_format(10^.x))) +
          annotation_logticks(sides="b")   +
          theme_classic()+
          labs(title= paste0(sub("_.*", "", m ), " signal drift over all batches"  )    , 
            subtitle= paste0("Using ", nrow(spleen.data.frame), " cells of all spleen ROIs.") ,
            #caption="Created by M.Barone", 
            y="Amount of cells in bin", 
            x= sub("_.*", "", m )  
          )+
          theme(panel.background = element_rect(fill = "grey93", colour = "grey55", size = 0.5), # change 'colour' to black for informative axis
            axis.title.x=element_text(color="grey15", size=11),
            axis.title.y=element_text(color="grey15", size=11),
            axis.text=element_text(size=10),
            #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
            legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
            legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
            #legend.title=element_blank(),
            plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
            plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
          )

        # using the unique colors is not worth plotting, it looks ugly

        # extract the histogram lines from ggplot engine:
        histograms <-  ggplot_build(g)$data[[1]]

        #the extract of ggplot has only the assigned colors in, so we need to supply that object with TMA names
        histograms <- cbind(histograms, Batch = names(Batches_palette)[match(histograms$colour,Batches_palette)] )

        # for the current marker, we will do the fitting TMA by TMA.
        # there is a multiple GAM option, but I do it the slow way to store each variable away in the tibble:
        for(t in all_images){

          oneOThistograms <- subset(histograms, Batch %in% t)

          # ggplot(oneOThistograms  , aes( x= x , y = ncount)  ) +  geom_line( )

          #lets create a fit with gam
          k.start <- seq(65, 5, by = -10)

          for(k in k.start){

            model <- try(gam(
              ncount~
                # s(x),
                s(x, bs = "gp",  k=k), #  tp, gp gaussian 4 modes  bs = "gp",
              data=oneOThistograms,
              method = "ML" # REML
            ), silent = TRUE )

            #now the trick is to run GAM with many modes and decrease k until the fit works. 
            #then the following condition is met and we break out of the GAM model
            if (is.character(model) == FALSE) break 

          }

          plot_smooths(     model = model, series = x )

          pred <- predict.gam(model,se.fit=TRUE)
          uniquepred <-unique(pred$fit)
          uniquepred.SE <-unique(pred$se.fit)

          oneOThistograms <- cbind(oneOThistograms,pred=as.numeric(uniquepred) ) 
          oneOThistograms <- cbind(oneOThistograms,pred.SE=as.numeric(uniquepred.SE) )

          # and return the x value of the max, and drop all info into 

          #oneOThistograms[which.max(oneOThistograms$pred),]$x

          spleen.normalization <- spleen.normalization %>% add_row(
            marker = m,
            batch = t,
            max.signal = 10^oneOThistograms[which.max(oneOThistograms$pred),]$x
          )

          diff(sign(diff(oneOThistograms$pred)))==-2

          sign(diff(oneOThistograms$pred))

          g <- 
            ggplot(oneOThistograms  , aes( x= x , y=ncount, color = Batch )  ) +
            geom_line(colour="red", size=2 )+
            geom_line(data = oneOThistograms, aes(x=x, y=pred), colour="black", size=1.1, alpha=0.8)+
            #geom_line(data = oneOThistograms, aes(x=x, y=(pred+  (pred.SE)*sqrt(nrow(oneOThistograms))  )), colour="black", size=0.5, alpha=0.3)+
            #geom_line(data = oneOThistograms, aes(x=x, y=(pred- (pred.SE)*sqrt(nrow(oneOThistograms)) )), colour="black", size=0.5, alpha=0.3)+
            labs(title= paste("Smoothed Spleen signal of ", m), 
              subtitle=paste("TMA ID: ", t), 
              #caption="Created by M.Barone", 
              y="Normalized Spleen signal in TMA"#, 
              #x="Days into the experiment"
            )

          filename <- paste0(sub("_.*", "", m ),"_spleen_",t,"_GAMmodel.png")
          ggsave(filename,g, width = 30, height = 20, units = "cm",dpi = 600)

          # tick the progress bar before doing another TMA of said marker.
          pb$tick()
        } # end run along all TMAs with t
      }# end run along spleen.col markers and plot

      spleen.normalization <-   (spleen.normalization %>% group_by(marker) %>% mutate(SF =  max.signal[1] / max.signal  )) 
      ggplot(spleen.normalization.histograms  , aes( x= batch , y=SF, group= marker, color = marker )  ) +
        geom_line( size=2 )

      # initially I wanted to use mutate() to fish out OTs and apply the SF there, but thats not working yet. so lets do the slow way:

      temp.cell.dat <- data.frame()
      for(t in all_images){

        temp.subset.cell.dat <- subset(cell.dat, Batch %in% t) 
        temp.SF <- subset(spleen.normalization, batch %in% t) 

        for(b in batchnormalize.marker){

          # so this just overrides the marker raw signal without any further comment
          temp.subset.cell.dat[,c(b)] <- temp.subset.cell.dat[[b]] * temp.SF[grepl(b, temp.SF$marker ), ]$SF

        }# run along all batchnormalize.markers

        temp.cell.dat <- rbind(temp.cell.dat, temp.subset.cell.dat )

      }# run along all batches

      #this is nuts.....
      cell.dat <- temp.cell.dat

      rm(temp.cell.dat)
      rm(temp.subset.cell.dat)
      rm(temp.SF)
    } # end protect from experimental code exec

    if(run.experimental.code == 1){

      #spleen.normalization.histograms <- spleen.normalization
      #spleen.normalization %>% add_column( Normalized = NA )

      #spleen.normalization$Nor <- ave(spleen.normalization$max.signal, spleen.normalization$batch, FUN=function(x) x/max(x))   

      #spleen.normalization.histograms <-  (spleen.normalization %>% group_by(marker) %>% mutate(Nor = max.signal/max(max.signal))) 

      # lets test

      testcelldata <- data.frame(
        batch = c("OT11", "OT11", "OT11", "OT11",
          "OT2", "OT2", "OT2", "OT2",
          "OT13", "OT13", "OT13", "OT13"),
        ID = c("1", "2", "3", "4",
          "1", "2", "3", "4",
          "1", "2", "3", "4"),
        CD3_Nd143 = c(1, 1, 1, 1,
          2,2,2,2,
          3,3,3,3),
        CD4_Dy161 = c(1, 1, 1, 1,
          2,2,2,2,
          3,3,3,3)
      )

      spleen.normalization[1,]

      #testcelldata[spleen.normalization$marker[1]]*spleen.normalization$SF[1]

      testcelldata[[spleen.normalization$marker[1]]]

      #testcelldata %>% mutate_at( testcelldata[[spleen.normalization$marker[1]]] = testcelldata[[spleen.normalization$marker[1]]]  %>%   
      #                            replace( testcelldata$batch == noquote(spleen.normalization$batch[1]) ,  20   ) 
      #                      )

      testcelldata %>% rowwise() %>% mutate(
        ifelse(testcelldata$batch == noquote(spleen.normalization$batch[1]), 
          CD3_Nd143 = testcelldata[CD3_Nd143]*spleen.normalization$SF[1], 
          NULL)
      )

      #testcelldata %>% rowwise() %>% mutate(
      #  
      #  spleen.normalization$marker[1] = spleen.normalization$marker[1]  %>% 
      #  ifelse( testcelldata$batch == spleen.normalization$batch[1] ,  spleen.normalization$marker[1] *spleen.normalization$SF[1] ,        ))

    } # end protect running experimental code with dypl to mutate parts in cell.dat

  }#end batch normalization routine
     
} # protect from running in first pass      
      
    
} #end default stage 1: load data from first script and generate cell.dat anew 

if(rebuild.cell.data.anew == 0){

  #--------------------------------------
  # Set Input Directory
  # if we rebuild anew, the input directory is the QuPath_tsv_inputs folder
  #--------------------------------------

  setwd(PrimaryDirectory)
  setwd(paste0("output/",day.to.reload.data,"/Data/") )
  InputDirectory <- getwd()


  file_loading.Lyrics <-  paste("\n\n\n")
  file_loading.Lyrics <-  paste(file_loading.Lyrics, "_____________________________________________________________________________\n\n")

  # now there is a chance that you want to proceed with STAGE 2,
  # so lets load all files in Data directory:
  all.STAGE0.cell.dats <- list.files(pattern = '.qs')

  # as you can see, we load both cell.dats into the re-running code. 
  # We decide which one will enter into STAGE1 via the calculate.scaling.factors==1 condition.
  if("scaled.cell.dat.qs" %in% all.STAGE0.cell.dats){ scaled.cell.dat <- qread("scaled.cell.dat.qs")  
  file_loading.Lyrics <-  paste(file_loading.Lyrics, "   Loaded scaled.cell.dat -> proceed with STAGE 1 and calculate.scaling.factors=1\n",   sep="")}

  if("unscaled.cell.dat.qs" %in% all.STAGE0.cell.dats){ cell.dat <- qread("unscaled.cell.dat.qs")  
  file_loading.Lyrics <-  paste(file_loading.Lyrics, "   Loaded (un-scaled) cell.dat -> proceed with STAGE 1 and calculate.scaling.factors=0\n",   sep="")}

  if("clustered.gated.subsets.qs" %in% all.STAGE0.cell.dats){ clustered.gated.subsets <- qread("clustered.gated.subsets.qs")  
  file_loading.Lyrics <-  paste(file_loading.Lyrics, "   Loaded clustered.gated.subsets -> proceed with STAGE 3\n",   sep="")  }

  # gated.cell.dat from STAGE 3+ is ready-to-go: it is
  # gated, cluster collapsed and annotated, might contain masks
  all.STAGE3.files <- list.files(pattern = '.csv')

  if("all.cells.csv" %in% all.STAGE3.files){ all.cells  <- fread("all.cells.csv")  # this is not STAGE 3 produced, but we need it for late stages
  file_loading.Lyrics <-  paste(file_loading.Lyrics, "   Loaded all.cells -> proceed with STAGE 4+\n",   sep="")}

  if("gated.cell.dat.csv" %in% all.STAGE3.files){ gated.cell.dat <- fread("gated.cell.dat.csv")  
  file_loading.Lyrics <-  paste(file_loading.Lyrics, "   Loaded gated.cell.dat -> proceed with STAGE 4+\n",   sep="")}

  if("masked.gated.cell.dat.csv" %in% all.STAGE3.files){ masked.gated.cell.dat <- fread("masked.gated.cell.dat.csv")  
  file_loading.Lyrics <-  paste(file_loading.Lyrics, "   Loaded masked.gated.cell.dat -> proceed with STAGE 5+\n",   sep="")}       

  file_loading.Lyrics <-  paste(file_loading.Lyrics, "_____________________________________________________________________________\n\n")
  message("Loading backed up cellular data file(s) complete. Check the files and proceed accordingly:")
  cat( file_loading.Lyrics )

  # now that cell.dat is loaded, switch back to default input directory:      
  ### Set InputDirectory 
  setwd(PrimaryDirectory)
  setwd(paste0("QuPath_tsv_inputs/"))
  InputDirectory <- getwd()
          
  ###################### COPY-PASTE VARIABLES IF YOU DO NOT GENERATE CELL.DAT ANEW ######################
                  
  # we still need to know some variables that we would skip in STAGE2 and in the end of STAGE 3.
  # add them manually!! 
  if(exists("cell.dat")){  
  # now, assuming that gated.cell.dat is a child of cell.dat in that folder:
  # before doing anything, lets check if the cellular.cols vector 
  # is present in cell.dat:

  # we do not care for typos as we do if we rebuild cell.dat from scratch, but simply 
  # test if 

  if(!all(paste0(
    Markers2Process,
    qupath.separator,
    cellular.segment.compartment,
    qupath.separator,
    cellular.segment.readout ) %in% names(cell.dat)) ){

    cat("\n\n\n----------------------------------- \n ERROR in STAGE 0 initialization.   \n 
      The loaded cell.dat does not contain all markers to process \n  
      Please check the Markers2Process vector and the cell.dat columns \n 
      -> You are likely using an outdated cellular column list or cell.dat \n-----------------------------------\n\n") 
    stop()   

    
  } else {

    cellular.cols <- paste0(
    Markers2Process,
    qupath.separator,
    cellular.segment.compartment,
    qupath.separator,
    cellular.segment.readout )

    # all settings done. 
    # initialize following stage 0 varibles :
    # all_images, amount.of.TMAs, all_cores, amount.of.cores, markercomparision
    list2env( initialize_stage_0_variables(
      cell.dat = cell.dat,
      scatter.plots = scatter.plots,
      transform.rawdata = transform.rawdata,
      qupath.separator = qupath.separator,
      cellular.segment.compartment = cellular.segment.compartment,
      cellular.segment.readout = cellular.segment.readout,
      asinh.cofactor = asinh.cofactor
    ), envir = .GlobalEnv)

  }# end check if all markers are present in cell.dat
  }# END check cellular columns from cell.dat
          
  # STAGE 2 did the gating and dropped a list of lists, which tells us the amount of GATE subsets we had initially gated.
  # Watch out, dont gated.cell.dat$GATEsubset as a pointer, you might have deleted in there the trashbin gate already
  if(exists("clustered.gated.subsets")){          
  amountofGATEsubsets <- length(clustered.gated.subsets)-exclude.trash.bin.gate 
  #admittedly this is a bit stupid, since we could just count the subsets in gated.cell.dat, which is the more important dataframe here..
  }# END load clustered.gated.subsets-related variables that you miss since you proceed with STAGE3
          
  # STAGE 3 bound then only the first amountofGATEsubsets'th subsets together. You need to know a couple of palettes if you hop over STAGE3:

  if(exists("gated.cell.dat")){
  
  all.MC <- unique(gated.cell.dat$Collapsed_metacluster)
  amount.of.MC <- length(all.MC)  

  #     CREATE polychrome palette for the unannotated Phenograph MCs
  colRange <- unique(gated.cell.dat[["Collapsed_metacluster"]])
  colRange <- colRange[order(colRange)]
  #colRange <- as.character(colRange)
  set.seed(723451) # for reproducibility
  Phenograph_metacluster_palette <- createPalette(length(colRange), c("#ff0000"), M=100000, prefix = "")
  # so, the cool thing is that once we annotated these clusters, we can just do the same game again, overrride the palette and plot again.
  names(Phenograph_metacluster_palette) <- colRange



  # we re-create now annotated.Phenograph_metacluster_palette like we do after annotation:
  cluster.annots <- do.list.switch(cluster.annots)
  names(cluster.annots) <- c('Values', 'Annotated_metacluster')
  cluster.annots
  # we first need to adjust all.MC to incorporate also the annotation
  # the easy way is of course:
  all.annot.MC <- do.add.cols(as.data.frame(all.MC), 'all.MC', cluster.annots, 'Values')
  ### Whatever Metacluster was not annotated got now NA
  # in our case we dont want "other" for these, but we need the original numbers back in:
  all.annot.MC <- all.annot.MC %>% 
  mutate(Annotated_metacluster = coalesce(Annotated_metacluster,all.MC))


  # now we tweak the list to be accepted by do.add.cols, but this time we inject into gated.cell.dat:
  names(all.annot.MC) <- c('Values', 'Annotated_metacluster')

  annotated.Phenograph_metacluster_palette <- data.frame("Collapsed_metacluster"=names(Phenograph_metacluster_palette), "color"=Phenograph_metacluster_palette, row.names=NULL)
  annotated.Phenograph_metacluster_palette <- do.add.cols(annotated.Phenograph_metacluster_palette, 'Collapsed_metacluster', all.annot.MC, 'Values')
  # you just gotta love do.add.cols, no?

  annotated.Phenograph_metacluster_palette$Collapsed_metacluster <- as.numeric(as.character(annotated.Phenograph_metacluster_palette$Collapsed_metacluster))

  # if you annotated and then changed colors, you might want to manually change Phenograph_metacluster_palette after re-loading gated.cell.dat.







  #     CREATE polychrome palette for the used TMA batches
  colRange <- unique(gated.cell.dat[["Batch"]])
  colRange <- colRange[order(colRange)]
  #colRange <- as.character(colRange)
  set.seed(723451) # for reproducibility
  OTbatch_palette <- createPalette(length(colRange), c("#009E73"), M=100000, prefix = "")

  #swatch(OTbatch_palette) # wanna see it?
  names(OTbatch_palette) <- colRange   

  #     CREATE polychrome palette for the split conditions
  colRange <- unique(gated.cell.dat[[sample.col]])
  colRange <- colRange[order(colRange)]
  # we need a dirty patch, since these characters come in like this: "ND_1"   "ND_100" "ND_12"  and so on:
  colRange <- c(colRange[1],colRange[3:5],colRange[2],colRange[6:8]         )
  set.seed(723451) # for reproducibility
  Conditions_palette <- createPalette(length(colRange), c("#ff0000"), M=100000, prefix = "")

  #swatch(Conditions_palette) # wanna see it?
  # so, the cool thing is that once we annotated these clusters, we can just do the same game again, overrride the palette and plot again.
  names(Conditions_palette) <- colRange

  # update our plot.ROIs object since it still contains spleen and stuff:
  plot.ROIs  <- unique(gated.cell.dat$ROI) 
  amount.of.ROIs <- length(plot.ROIs)

  # this variable is initially created via cell.dat$Batch and then a second time end of STAGE 3 from gated.cell.dat like so:
  all_images <- sort(unique(gated.cell.dat$Batch)) 
  amount.of.TMAs <- length(all_images)


  # also create some vectors to run along:

  diet.vector <-  unique( gated.cell.dat$Diet[str_detect(gated.cell.dat$Diet, "ctr", negate=T)] )
  # with the weeks we sort them and remove the NA from the spleen sample:
  # the problem of the column is that its a character, so na.last does not work except when transforming them into numerals
  sampleweek.vector <-  sort( unique( suppressWarnings(as.numeric(gated.cell.dat$Sample_week)) ), na.last=NA ) 


  all.condtions <-   sort( unique(gated.cell.dat$Diet_week) )
  # and again, we need a dirty patch to get ND_100 back into position:
  all.condtions <- c(all.condtions[1],all.condtions[3:5],all.condtions[2],all.condtions[6:8]         )






  #     CREATE a palette for diet
  diet_palette <- c(
    "ND" = "#017DD5",
    "WD" = "#EB8B19"
  )

  # lets assign some cute colors and use it by calling scale_fill_manual(values = cond_palette)+
  cond_palette <- c(
    "ND_1"   = "#509CEF",
    "ND_12"  = "#2070C8",
    "ND_14"  = "#185EAB",
    "ND_24"  = "#0A4689",
    "ND_100" = "#042852",
    "WD_1"   = "#F48E63",
    "WD_12"  = "#CD5420",
    "WD_24"  = "#8F3008"
  )
          


  # watch out, in STAGE 1 you set the GATEsubset names together with the channels you used to cluster the data subset.
  # from this object, the code pulls the Names, from nowhere else. To get these GateSubet Names, load the object:
  # come back here, cluster.col




cluster.cols.GATEsubsets <- list(
  names(gated.cell.dat)[c(91,119,123,115,120,100,96,121,125,108,104,95,117,97,122,118)+4], 
  names(gated.cell.dat)[c(91,119,101,105,109,98,114,124,123,115,120,100,96,111,121,125,108,104,95,102,117,106,107,97,122,118)+4], # 88 is CD19
  names(gated.cell.dat)[c(91,119,101,105,109,98,114,124,123,115,120,100,96,111,121,125,108,104,95,102,117,106,107,97,122,118)+4],
  names(gated.cell.dat)[c(98,96,125,108,104,95,102,117,110,103,122,94,118,92)+4], 
  names(gated.cell.dat)[c(110)]
)

names(cluster.cols.GATEsubsets) <- GATEsubsetNames


cat("\n\n----------------------------------- \n Cluster channels for every subset loaded. \n -> Please check! \n\n")# stop(), call. = FALSE)  
print(cluster.cols.GATEsubsets)
cat("\n-----------------------------------\n")






}# END load gated.cell.dat-related variables that you miss since you proceed with STAGE4


















          
}#end alternate stage 0: read already produced cell.dat file from previous run      
      
      

      
      
# The lyrics in the end of STAGE 0 depend of what settings you chose:
# whether you rebuilt cell.dat from scratch or loaded a gated, clustered and collapsed gated.cell.dat
# and if you rebuild, STAGE 0 is run up to 3 times
      
      
if(rebuild.cell.data.anew == 1){      

  # for some reason I sometimes end this stage having a plotting dev open.
  # so if dev.list returns not a NULL; close it:
  if( !is.null(dev.list()) ){
  invisible(dev.off())
  }



  if(  how.often.ran.STAGE0==0 ) { 
    
    if( calculate.scaling.factors ==0 ) {  
      cat("\n\n\n----------------------------------- \n End of STAGE 0 frist run.   \n Spillover correction: ON/OFF \n Batch normalization: OFF \n  1) Drop markers not present in every experiment \n  2) Make sure channel names are consistent \n -> Proceed to STAGE 1 \n-----------------------------------\n\n")# stop(), call. = FALSE)  
    }   
    
    if(calculate.scaling.factors ==1 ) {  
      cat("\n\n\n----------------------------------- \n End of STAGE 0 frist run.   \n Spillover correction: ON/OFF \n Batch normalization: ON \n  1) Drop markers not present in every experiment \n  2) Make sure channel names are consistent \n -> Re-run to STAGE 0 to collect percentile signals \n-----------------------------------\n\n")# stop(), call. = FALSE)  
    }
    
    }
  if(  how.often.ran.STAGE0==1 ) { 
    
    if( calculate.scaling.factors ==0 ) {  
      cat("\n\n\n----------------------------------- \n End of STAGE 0 second run.   \n Spillover correction: ON/OFF \n Batch normalization: OFF \n  1) Drop markers not present in every experiment \n  2) Make sure channel names are consistent \n -> Proceed to STAGE 1 \n-----------------------------------\n\n")# stop(), call. = FALSE)  
    }   
    
    if(calculate.scaling.factors ==1 ) {  
      cat("\n\n\n----------------------------------- \n End of STAGE 0 second run.   \n Spillover correction: ON/OFF \n Batch normalization: ON \n 1) Check Batch_Normalization folder to align channel(s)  \n  2) Define channel-wise percentile in batch.normalization.matrix \n -> Re-run to STAGE 0 to apply batch norm. \n-----------------------------------\n\n")# stop(), call. = FALSE)  
    }
    
  }

  if(  how.often.ran.STAGE0>1 ) { 
    
    if( calculate.scaling.factors ==0 ) {  
      cat("\n\n\n----------------------------------- \n End of STAGE 0 third run.   \n Spillover correction: ON/OFF \n Batch normalization: OFF \n -> Proceed to STAGE 1 \n-----------------------------------\n\n")# stop(), call. = FALSE)  
    }   
    
    if(calculate.scaling.factors ==1 ) {  
      cat("\n\n\n----------------------------------- \n End of STAGE 0 third run.   \n Spillover correction: ON/OFF \n Batch normalization: ON \n  -> Proceed to STAGE 1 with scaled raw data \n-----------------------------------\n\n")# stop(), call. = FALSE)  
    }
    
  }



}else{
  
  if(exists("gated.cell.dat")){
cat("\n\n\n----------------------------------- \n End of STAGE 0. \n I found and loaded gated.cell.dat \n Proceed to STAGE 4 or farther! \n-----------------------------------\n\n")
  }else{
cat("\n\n\n----------------------------------- \n End of STAGE 0. \n I did not find gated.cell.dat \n Proceed with a listed file and the according STAGE \n-----------------------------------\n\n")
    
  }
}
      
  
} # end stage 0: read, drop bad ROIs and allocate consistent markers       







########################## STAGE 1 ##########################

# yes, we source that everytime anew, just in case we are fiddling in there I wanna be sure R loads the right functions
if(whichSTAGE > 0){
  setwd(PrimaryDirectory)
  source(paste0("utils_v",utils_version,".R"))
}

if(whichSTAGE == 1){
  

 
  
# in case we scaled, we will store away cell.dat and replace that object with scaled.cell.dat. 
if(calculate.scaling.factors==1){
  # override now the raw data:
  cell.dat <- scaled.cell.dat
  message("\n\n\n################################################\n# Scaled raw data are being used from here on! #\n################################################")
}
  
       
      
      

# just in case you need to remove a column now:      
#     if("outside.gate" %in% colnames(cell.dat)){
#     cell.dat <-  subset(cell.dat, select = -c(outside.gate) )
#    }
  
  

      

#------------------------------
# Flag cell segmentation penalty (square micrometers)
#------------------------------

# outliers in segmenation area and gated cells are flagged in the cell.dat object
# the flag is inverted, so 0 means active datapoint, 1 means outlier
    
# the code does not subset cell.dat, but will simply flag the cells
# for the area we need a range penalty:
  # it uses j <= i & i <= k  to flag "in range, 0"
cell.dat$outside.area <-mapply(function.range.penalty,i=cell.dat$Cell..Area.µm.2, j=cell.dat$area.highpass, k=cell.dat$area.lowpass  )



      
#------------------------------
# Gating cells by flagging:
# -> Transform marker signals
#------------------------------

# gating is done either on raw or transformed data,
# depending on whether transform.data is TRUE or FALSE
# if TRUE, we need to transform the data now using do.data.normalization:
if(transform.rawdata == TRUE){
  cell.dat <- do.data.normalization(
    dat = cell.dat,
    use.cols = cellular.cols,
    # Transform
    do.transform = transform.rawdata,
    cofactor = asinh.cofactor,
    # min-max normalize?
    do.minmax = FALSE,
    new.min = min.norm,
    new.max = max.norm,
    # z-score normalize?
    do.zscore = FALSE
  )
} # end transform.data == TRUE



#------------------------------
# Gating cells by flagging
# -> Threshold gates supplied in the metadata file
#------------------------------

# the plotting engine will will check if the marker on x/y axis is used as gate
# if so, it will add the gate line in the plot
# we will therefore collect all markers that are affected by a gate line and store them in gate_affected.markers
gate.col.expr <- c(".highpass", ".lowpass")   
threshold.gates <- names(cell.dat)[ grep( paste(gate.col.expr, collapse="|")   , colnames(cell.dat) ) ]
# this pulled any gate threshold, but also the area.highpass and area.lowpass columns, so drop these:
threshold.gates <- grep("area.", threshold.gates, invert=TRUE, value = TRUE)   
# if you re-run stage 1, you also pull out the already set flag columns
# remove them as well:
threshold.gates <- grep("outside.", threshold.gates, invert=TRUE, value = TRUE) 

# we prepare the gate_affected.markers vector either way:.after# for plotting, we need to know the channel AB affected by a gate:
gate_affected.markers <-    sub(".[^.]+$", "", threshold.gates ) # character 0 if no threshold.

# if you do not supply and thresholdgates (yet), gate_affected.markers stays an empty character vector for now and we skip the upcoming part: 
if (!identical(threshold.gates, character(0))) {

  #------------------------------
  # Flagging cells with threshold gates
  #------------------------------

  for (g in threshold.gates) {

    # HIGH-PASS GATE FLAGGING
    if (grepl(".highpass", g)) {
      # do not chop off the marker and then grep it from cell.dat!
      # you risk to fish more than one marker: CD3 grep will return CD3, CD31, CD321 and so on...
      # just remove the last dot and assemble the column name:

      temp.gated.marker <- lapply(
        sub(".[^.]+$", "", g),
        function(x) ifelse(
          transform.rawdata == TRUE,
          paste0(
            x,
            qupath.separator,
            cellular.segment.compartment,
            qupath.separator,
            cellular.segment.readout,
            "_t", asinh.cofactor
          ),
          paste0(
            x,
            qupath.separator,
            cellular.segment.compartment,
            qupath.separator,
            cellular.segment.readout
          )
        )
      )
      # we initialized it already, push the markes now in:
      gate_affected.markers <- c(gate_affected.markers, temp.gated.marker)

      # watch out, this column name structure is afterwards used to pull the flags, such as in "outside.",yMarker,".highpass". Caution if you wanna change that...
      cell.dat[, paste0("outside.", g)] <- mapply(
        function.lower.penalty,
        i = cell.dat[[temp.gated.marker]],
        j = as.numeric(cell.dat[[g]])
      )

      passed.cells <- table(cell.dat[[paste0("outside.", g)]])[1]
      message(
        g, " GATE: ", passed.cells, " cells (",
        round(passed.cells / nrow(cell.dat) * 100, 2),
        "%) passed high-pass gate of ", temp.gated.marker
      )

    } # end high-pass flagging

    # LOW-PASS GATE FLAGGING
    if (grepl(".lowpass", g)) {

      temp.gated.marker <- lapply(
        sub(".[^.]+$", "", g),
        function(x) ifelse(
          transform.rawdata == TRUE,
          paste0(
            x,
            qupath.separator,
            cellular.segment.compartment,
            qupath.separator,
            cellular.segment.readout,
            "_t", asinh.cofactor
          ),
          paste0(
            x,
            qupath.separator,
            cellular.segment.compartment,
            qupath.separator,
            cellular.segment.readout
          )
        )
      )

      # we initialized it already, push the markes now in:
      gate_affected.markers <- c(gate_affected.markers, temp.gated.marker)

      cell.dat[, paste0("outside.", g)] <- mapply(
        function.upper.penalty,
        i = cell.dat[[temp.gated.marker]],
        j = as.numeric(cell.dat[[g]])
      )

      passed.cells <- table(cell.dat[[paste0("outside.", g)]])[1]
      message(
        g, " GATE: ", passed.cells, " cells (",
        round(passed.cells / nrow(cell.dat) * 100, 2),
        "%) passed low-pass gate of ", temp.gated.marker
      )

    } # end low-pass flagging

  } # end run with g along all threshold gates and create the Hgate or Lgate flag columns in cell.dat

} # end flagging with threshold gates supplied in the metadata file



#--------------------------------
# MULTISEGMENT LINE GATE FLAGGING
#--------------------------------

# depending on whether on which side (right/left) the cell lays, it is flagged
# the x axis marker flags the cells:
# outside.xx.highpass flags the cells right of the line as 0
# outside.xx.lowpass flags the cells right of the line as 1
# the line can run diagonally TL to BR, or BL to TR, doesnt matter

# check the ranges. you must not have NA in there
# the polygon gate points must be outside the range of the data in both x and y direction.
#r$> range(cell.dat$CD45..Cell..Mean_t2)
#[1] 0.000000 8.787497
#r$> range(cell.dat$CD3e..Cell..Mean_t2)
#[1]  5.171347 10.179299
#r$> range(cell.dat$CD20..Cell..Mean_t2)
#[1] 0.00000 8.08536
#r$> range(cell.dat$E.cadherin..Cell..Mean_t2)
#[1] 0.000000 9.807788

polygon.gates <- list(
#do.not.use.polygon.gates <- list(

# Important: the format xx.highpass or xx.lowpass must contain the exact maker name xx
# xx is used to assemle the column name in cell.dat: use CD3e for example, not CD3
  
  
  "CD20.highpass"  = tribble(
    #x            #y   
    ~CD45, ~CD20,
    -1,  40 , # <- if you intend to fine-tune that graph, make sure that the scaling does not bring this corner point inside the data range!!
    1, 10,
    4,7,
    6,  6.25,
    8,6,
    10,5.9,
    1000,-1   # <- if you intend to fine-tune that graph, make sure that the scaling does not bring this corner point inside the data range!!
  ),
  
  
 
  "CD3e.highpass"  = tribble(
    #x            #y   
    ~CD45, ~CD3e,
   -1,  40 , 
    1, 10,
    5,9,
    6,  8.5,
    6.25, 7.8,
    7,7.6,
    8,7.5,
    10,7,
    1000,1
  ),




"CD45.highpass"  = tribble(
  #x            #y   
  ~CD45, ~E.cadherin,
   -1,  0 , 
    4,0.5,
    5,2,
    6,  3.5,
    6.5,6.5,
    7,6.9,
    7.5,10,
    8,12,
    1000,11
)
  
)# end polygon.gates list


# you might wanna fine-tune the position of the gate Image by Image
# if the gate line is off too far, its an indication 
# that you have titration batch effects. 
# then its worth batch-correction the data 
# you do it like this:
#finetune.CD45segmentline[  match ("OT11", scaling.factors$Batch.ID), "CD45_Er168"]  <- 0.6


# lets memorize which marker combinations are part of these gates:
if(exists("polygon.gates") ){

  # the marker inside each polygon gate markers are not in the proper column format yet.

  all.polygon.gate.markers <- c()
  for(g in names(polygon.gates)){
    g <- as.character(g) # just in case it is a factor
    # the polygon gate markers are not in the proper column format yet.
    # build up the correct column names:
    names(polygon.gates[[g]]) <-  lapply(names(polygon.gates[[g]]), function(x) ifelse( transform.rawdata == TRUE, 
    paste0(
          x,
          qupath.separator,
          cellular.segment.compartment,
          qupath.separator,
          cellular.segment.readout,
          "_t",asinh.cofactor
        ),
        paste0(
          x,
          qupath.separator,
          cellular.segment.compartment,
          qupath.separator,
          cellular.segment.readout
        )
    ))

    all.polygon.gate.markers <- c( all.polygon.gate.markers, names(polygon.gates[[g]]) )
    
  }# end run along the polygon.gates to find out which markers are used. then create a finetune matrix:

# any of the used makers in polygon.gates can be adjusted by batch to correct some titration artifacts
# for this we scale the polygon gate points by a scaling factor. Set that to 1 by default:
all.polygon.gate.markers <- unique(all.polygon.gate.markers)
# by default, a fine.tune list is created:
finetune.segmentline <- data.frame( Batch.ID = all_images   )
finetune.segmentline[, all.polygon.gate.markers ] <- 1

# you might wanna fine-tune the position of the gate Image by Image
# Warning: this is an indication for titration batch effects. 
# consider batch-correction the data instead since this here only affects gating, not the data itself 
# to fine-tune adapt the factor to move the line along each marker:
#finetune.segmentline[  match ("20250228_27X_bleaching trial (control)_Scan1.er.qptiff - resolution #1", scaling.factors$Batch.ID), "CD45..Cell..Mean_t2"]  <- 0.6


# flagging via polygons now:

# we will create multisegmented lines given the xy points in the polygon.gates list on the fly:
pb <- progress_bar$new(format = "[:bar] :percent [Setting gate flags | :eta]",
                       total =  length(names(polygon.gates)), #
                       show_after=0, #allows to call it right way 
                       current = "|",    # Current bar character
                       clear = TRUE) # make it persist
pb$tick(0) # call in the progress bar without any progress, just to show it  
  
for(g in names(polygon.gates)){
    
    
  f <- multisegmentGATEline( unlist(polygon.gates[[g]][2])  , unlist(polygon.gates[[g]][1])   )
  #f <- multisegmentGATEline( unlist(polygon.gates[[g]][1])  , unlist(polygon.gates[[g]][2])   )
    
    
  # before flagging, we need to make sure the defined line has the proper dependency along x or y:
  # if the multisegmented line acts on the marker on the x-axis:
  # in this case the line will check on the height y if the data point lays right hand of the line position calculated at height y:
  #   if( grepl( sub(".[^.]+$", "", g ), names(polygon.gates[[g]])[1] ) ){
  #     f <- multisegmentGATEline( unlist(polygon.gates[[g]][2])  , unlist(polygon.gates[[g]][1])   )
  #   }else if(  grepl( sub(".[^.]+$", "", g ), names(polygon.gates[[g]])[2] ) ){
  #     # in the other case the coordinates are simply swapped:
  #     f <- multisegmentGATEline( unlist(polygon.gates[[g]][1])  , unlist(polygon.gates[[g]][2])   )
  #   }
    
    
    #HIGH-PASS GATE FLAGGING
    if(  grepl(".highpass", g)   ){
      
      cell.dat[, paste0("outside.",g) ] <-  ifelse( cell.dat[[ names(polygon.gates[[g]])[1] ]] >   finetune.segmentline[  match ( cell.dat$Image, finetune.segmentline$Batch.ID), names(polygon.gates[[g]])[1]] * f(  cell.dat[[ names(polygon.gates[[g]])[2] ]] / finetune.segmentline[  match ( cell.dat$Image, finetune.segmentline$Batch.ID), names(polygon.gates[[g]])[2] ] )  , 0, 1 )
    }
    
    #LOW-PASS GATE FLAGGING
    if(  grepl(".lowpass", g)   ){
      cell.dat[, paste0("outside.",g) ] <-  ifelse( cell.dat[[ names(polygon.gates[[g]])[1] ]] >   finetune.segmentline[  match ( cell.dat$Image, finetune.segmentline$Batch.ID), names(polygon.gates[[g]])[1]] * f(  cell.dat[[ names(polygon.gates[[g]])[2] ]] / finetune.segmentline[  match ( cell.dat$Image, finetune.segmentline$Batch.ID), names(polygon.gates[[g]])[2] ] )  , 1, 0 )
    }
    
  
  
  
  
  # since this gate flag column did not come via metadata file, its affected marker does not make its way into the gate_affected.markers list and we need to add it now:

  # we need to pull the marker col name first.
  #strip off the gate tag and reconstruct the channel name:
  temp.gated.marker <-  lapply( sub(".[^.]+$", "", g )   , function(x) ifelse( transform.rawdata == TRUE, 
    paste0(
          x,
          qupath.separator,
          cellular.segment.compartment,
          qupath.separator,
          cellular.segment.readout,
          "_t",asinh.cofactor
        ),
        paste0(
          x,
          qupath.separator,
          cellular.segment.compartment,
          qupath.separator,
          cellular.segment.readout
        )
    ))

  # only after the channel is fully assigned, we bring it into gate_affected.markers
  gate_affected.markers <- c(gate_affected.markers, temp.gated.marker )
  
  pb$tick()

} # end run along the polygon.gates to set the flags
} # end flagging routine via polygons stored in polygon.gates




# here is the protected plot to check your polygon:
if(run.experimental.code==1){
  


#-----------------
# setup CD20 highpass
#-----------------


  # lets get closer# we want to see reasonable gates, so cut gated.cell.dat back to remove some hilarious outliers in xMarker and yMarker:
  xMarker<-"CD45..Cell..Mean_t2"
  yMarker<-"CD20..Cell..Mean_t2"
  colorMarker <- "CD3e..Cell..Mean_t2"

  
g <- ggplot( data=cell.dat ,
          aes( x= CD45..Cell..Mean_t2,   y= CD20..Cell..Mean_t2 , color = CD3e..Cell..Mean_t2 ) ) + # color =CD68..Mean as.factor(outside.CD3e.highpass)
    
    geom_point(#aes(color = Annotated_metacluster),
      alpha=1, 
      size=0.2 #, 
      # pch='.'
    )+

    geom_path(data = polygon.gates[[1]],  size = 1.1, colour = "#A100AB", alpha=1) +
    
    labs(title=paste0( "All cells") , 
         subtitle="Gating into B cells",
          x = xMarker,
          y = yMarker
    )+
    
 #   scale_color_manual(name = "Gate flag", 
#                       guide =  guide_legend(override.aes = list(size=4,alpha=1 )),
#                       values = c("#40B0A6", "#E1BE6A"), #colorfriendly. not this: c("#308441", "#F6CD7A")
#                       labels = c("in", "out") #0 blue in, 1 yellow out
#    )+  
  
  scico::scale_color_scico(palette = "roma", direction=-1, 
                          limits=force,
                          # midpoint = median,
                           #limits=c(lower.threshold,upper.threshold),
                          # na.value = "#8C0172",
                            name = colorMarker
  )+ # bilbao roma
  
   # zoom a bit in to get the outliers out of the pic 
  coord_cartesian(ylim=c(as.numeric(quantile(cell.dat[[yMarker]], 0.001,na.rm=TRUE)),as.numeric(quantile(cell.dat[[yMarker]], 0.9995,na.rm=TRUE))), 
                  xlim=c(as.numeric(quantile(cell.dat[[xMarker]], 0.001,na.rm=TRUE)),as.numeric(quantile(cell.dat[[xMarker]], 0.9995,na.rm=TRUE)))
                  )+

    scale_x_log10(breaks = minor_breaks, 
                  minor_breaks = minor_breaks#,
                  #labels = trans_format("log10", math_format(10^.x))
                  ) +
    scale_y_log10(breaks = minor_breaks, 
                  minor_breaks = minor_breaks#,
                  #labels = trans_format("log10", math_format(10^.x))
                  ) +
    
    stat_density_2d(#aes( fill = ..level..),        #fill = ..level..), 
      alpha = 0.2, 
      contour = TRUE,
      contour_var = "density",
      adjust=1.7,
      n = 500, # default 100, computationally intense, do not exceed that value too much..
      geom = "polygon", #geom: polygon or raster (equidistant lines or heatmap-ish picture)
      colour="white") +
    
    theme(panel.background = element_rect(fill = "grey95", colour = "grey55", linewidth = 0.5), # change 'colour' to black for informative axis
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=8),
          axis.text.x = element_text(angle = 90, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          legend.position="right",
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )
  
ggsave("CD20.jpg", g, width = 20, height = 16, units = "cm", dpi = 300)



g <- ggplot( data=cell.dat ,
          aes( x= CD45..Cell..Mean_t2,   y= CD20..Cell..Mean_t2 , color = as.factor(outside.CD20.highpass) ) ) + 
    
    geom_point(#aes(color = Annotated_metacluster),
      alpha=1, 
      size=0.2 #, 
      # pch='.'
    )+

    geom_path(data = polygon.gates[[1]],  size = 1.1, colour = "#A100AB", alpha=1) +
    
    labs(title=paste0( "All cells") , 
         subtitle="Gating into B cells",
          x = xMarker,
          y = yMarker
    )+
    
    scale_color_manual(name = "Gate flag", 
                       guide =  guide_legend(override.aes = list(size=4,alpha=1 )),
                      values = c("#40B0A6", "#E1BE6A"), #colorfriendly. not this: c("#308441", "#F6CD7A")
                       labels = c("in", "out") #0 blue in, 1 yellow out
    )+  
  

  
   # zoom a bit in to get the outliers out of the pic 
  coord_cartesian(ylim=c(as.numeric(quantile(cell.dat[[yMarker]], 0.001,na.rm=TRUE)),as.numeric(quantile(cell.dat[[yMarker]], 0.9995,na.rm=TRUE))), 
                  xlim=c(as.numeric(quantile(cell.dat[[xMarker]], 0.001,na.rm=TRUE)),as.numeric(quantile(cell.dat[[xMarker]], 0.9995,na.rm=TRUE)))
                  )+

    scale_x_log10(breaks = minor_breaks, 
                  minor_breaks = minor_breaks#,
                  #labels = trans_format("log10", math_format(10^.x))
                  ) +
    scale_y_log10(breaks = minor_breaks, 
                  minor_breaks = minor_breaks#,
                  #labels = trans_format("log10", math_format(10^.x))
                  ) +
    
    stat_density_2d(#aes( fill = ..level..),        #fill = ..level..), 
      alpha = 0.2, 
      contour = TRUE,
      contour_var = "density",
      adjust=1.7,
      n = 500, # default 100, computationally intense, do not exceed that value too much..
      geom = "polygon", #geom: polygon or raster (equidistant lines or heatmap-ish picture)
      colour="white") +
    
    theme(panel.background = element_rect(fill = "grey95", colour = "grey55", linewidth = 0.5), # change 'colour' to black for informative axis
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=8),
          axis.text.x = element_text(angle = 90, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          legend.position="right",
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )
  
ggsave("CD20_gated.jpg", g, width = 20, height = 16, units = "cm", dpi = 300)



#-----------------
# setup CD3e highpass
#-----------------


u <- list()
xMarker<-"CD45..Cell..Mean_t2"
yMarker<-"CD3e..Cell..Mean_t2"


colorMarker <- "CD8..Cell..Mean_t2"

lower.threshold <- as.numeric(quantile( cell.dat[,  ..colorMarker ], probs = c( 0.02),na.rm=TRUE ))
median <- as.numeric(quantile( cell.dat[,  ..colorMarker ], probs = c( 0.5),na.rm=TRUE ))
upper.threshold <- as.numeric(quantile( cell.dat[,  ..colorMarker ], probs = c( 0.98),na.rm=TRUE ))

u[[1]] <- ggplot( data=cell.dat ,
          aes( x= CD45..Cell..Mean_t2,   y= CD3e..Cell..Mean_t2 , color = CD8..Cell..Mean_t2 ) ) + # color =CD68..Mean as.factor(outside.CD3e.highpass)
    
    geom_point(#aes(color = Annotated_metacluster),
      alpha=1, 
      size=0.2 #, 
      # pch='.'
    )+

    geom_path(data = polygon.gates[[2]],  size = 1.1, colour = "#A100AB", alpha=1) +
    
    labs(title=paste0( "All cells") , 
         subtitle="Gating into T cells (colored by 2-98th precentile)",
          x = xMarker,
          y = yMarker
    )+
    
 #   scale_color_manual(name = "Gate flag", 
#                       guide =  guide_legend(override.aes = list(size=4,alpha=1 )),
#                       values = c("#40B0A6", "#E1BE6A"), #colorfriendly. not this: c("#308441", "#F6CD7A")
#                       labels = c("in", "out") #0 blue in, 1 yellow out
#    )+  
  
  scico::scale_color_scico(palette = "roma", direction=-1, 
                          oob=squish,
                           midpoint = median,
                           limits=c(lower.threshold,upper.threshold),
                          # na.value = "#8C0172",
                            name = colorMarker
  )+ # bilbao roma
  
   # zoom a bit in to get the outliers out of the pic 
  coord_cartesian(ylim=c(as.numeric(quantile(cell.dat[[yMarker]], 0.001,na.rm=TRUE)),as.numeric(quantile(cell.dat[[yMarker]], 0.9995,na.rm=TRUE))), 
                  xlim=c(as.numeric(quantile(cell.dat[[xMarker]], 0.001,na.rm=TRUE)),as.numeric(quantile(cell.dat[[xMarker]], 0.9995,na.rm=TRUE)))
                  )+

    scale_x_log10(breaks = minor_breaks, 
                  minor_breaks = minor_breaks#,
                  #labels = trans_format("log10", math_format(10^.x))
                  ) +
    scale_y_log10(breaks = minor_breaks, 
                  minor_breaks = minor_breaks#,
                  #labels = trans_format("log10", math_format(10^.x))
                  ) +
    
    stat_density_2d(#aes( fill = ..level..),        #fill = ..level..), 
      alpha = 0.2, 
      contour = TRUE,
      contour_var = "density",
      adjust=1.7,
      n = 500, # default 100, computationally intense, do not exceed that value too much..
      geom = "polygon", #geom: polygon or raster (equidistant lines or heatmap-ish picture)
      colour="white") +
    
    theme(panel.background = element_rect(fill = "grey95", colour = "grey55", linewidth = 0.5), # change 'colour' to black for informative axis
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=8),
          axis.text.x = element_text(angle = 90, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          legend.position="right",
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )
  


colorMarker <- "CD4..Cell..Mean_t2"
lower.threshold <- as.numeric(quantile( cell.dat[,  ..colorMarker ], probs = c( 0.02),na.rm=TRUE ))
median <- as.numeric(quantile( cell.dat[,  ..colorMarker ], probs = c( 0.5),na.rm=TRUE ))
upper.threshold <- as.numeric(quantile( cell.dat[,  ..colorMarker ], probs = c( 0.98),na.rm=TRUE ))


u[[2]] <- ggplot( data=cell.dat ,
          aes( x= CD45..Cell..Mean_t2,   y= CD3e..Cell..Mean_t2 , color = CD4..Cell..Mean_t2 ) ) + # color =CD68..Mean as.factor(outside.CD3e.highpass)
    
    geom_point(#aes(color = Annotated_metacluster),
      alpha=1, 
      size=0.2 #, 
      # pch='.'
    )+

    geom_path(data = polygon.gates[[2]],  size = 1.1, colour = "#A100AB", alpha=1) +
    
    labs(title=paste0( "All cells") , 
         subtitle="Gating into T cells (colored by 2-98th precentile)",
          x = xMarker,
          y = yMarker
    )+
    
 #   scale_color_manual(name = "Gate flag", 
#                       guide =  guide_legend(override.aes = list(size=4,alpha=1 )),
#                       values = c("#40B0A6", "#E1BE6A"), #colorfriendly. not this: c("#308441", "#F6CD7A")
#                       labels = c("in", "out") #0 blue in, 1 yellow out
#    )+  
  
  scico::scale_color_scico(palette = "roma", direction=-1, 
                          oob=squish,
                           midpoint = median,
                           limits=c(lower.threshold,upper.threshold),
                          # na.value = "#8C0172",
                            name = colorMarker
  )+ # bilbao roma
  
   # zoom a bit in to get the outliers out of the pic 
  coord_cartesian(ylim=c(as.numeric(quantile(cell.dat[[yMarker]], 0.001,na.rm=TRUE)),as.numeric(quantile(cell.dat[[yMarker]], 0.9995,na.rm=TRUE))), 
                  xlim=c(as.numeric(quantile(cell.dat[[xMarker]], 0.001,na.rm=TRUE)),as.numeric(quantile(cell.dat[[xMarker]], 0.9995,na.rm=TRUE)))
                  )+

    scale_x_log10(breaks = minor_breaks, 
                  minor_breaks = minor_breaks#,
                  #labels = trans_format("log10", math_format(10^.x))
                  ) +
    scale_y_log10(breaks = minor_breaks, 
                  minor_breaks = minor_breaks#,
                  #labels = trans_format("log10", math_format(10^.x))
                  ) +
    
    stat_density_2d(#aes( fill = ..level..),        #fill = ..level..), 
      alpha = 0.2, 
      contour = TRUE,
      contour_var = "density",
      adjust=1.7,
      n = 500, # default 100, computationally intense, do not exceed that value too much..
      geom = "polygon", #geom: polygon or raster (equidistant lines or heatmap-ish picture)
      colour="white") +
    
    theme(panel.background = element_rect(fill = "grey95", colour = "grey55", linewidth = 0.5), # change 'colour' to black for informative axis
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=8),
          axis.text.x = element_text(angle = 90, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          legend.position="right",
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )
  



  g <- do.call("arrangeGrob", c(u, ncol=2)) # use a grob to save the obj
ggsave("CD3.jpg", g, width = 40, height = 16, units = "cm", dpi = 300)



g <-  ggplot( data=cell.dat ,
          aes( x= CD45..Cell..Mean_t2,   y= CD3e..Cell..Mean_t2 , color = as.factor(outside.CD3e.highpass)  ) ) + # color =CD68..Mean as.factor(outside.CD3e.highpass)
    
    geom_point(#aes(color = Annotated_metacluster),
      alpha=1, 
      size=0.2 #, 
      # pch='.'
    )+

    geom_path(data = polygon.gates[[2]],  size = 1.1, colour = "#A100AB", alpha=1) +
    
    labs(title=paste0( "All cells") , 
         subtitle="Gating into T cells (colored by 2-98th precentile)",
          x = xMarker,
          y = yMarker
    )+
    
    scale_color_manual(name = "Gate flag", 
                       guide =  guide_legend(override.aes = list(size=4,alpha=1 )),
                       values = c("#40B0A6", "#E1BE6A"), #colorfriendly. not this: c("#308441", "#F6CD7A")
                       labels = c("in", "out") #0 blue in, 1 yellow out
    )+  
  

  
   # zoom a bit in to get the outliers out of the pic 
  coord_cartesian(ylim=c(as.numeric(quantile(cell.dat[[yMarker]], 0.001,na.rm=TRUE)),as.numeric(quantile(cell.dat[[yMarker]], 0.9995,na.rm=TRUE))), 
                  xlim=c(as.numeric(quantile(cell.dat[[xMarker]], 0.001,na.rm=TRUE)),as.numeric(quantile(cell.dat[[xMarker]], 0.9995,na.rm=TRUE)))
                  )+

    scale_x_log10(breaks = minor_breaks, 
                  minor_breaks = minor_breaks#,
                  #labels = trans_format("log10", math_format(10^.x))
                  ) +
    scale_y_log10(breaks = minor_breaks, 
                  minor_breaks = minor_breaks#,
                  #labels = trans_format("log10", math_format(10^.x))
                  ) +
    
    stat_density_2d(#aes( fill = ..level..),        #fill = ..level..), 
      alpha = 0.2, 
      contour = TRUE,
      contour_var = "density",
      adjust=1.7,
      n = 500, # default 100, computationally intense, do not exceed that value too much..
      geom = "polygon", #geom: polygon or raster (equidistant lines or heatmap-ish picture)
      colour="white") +
    
    theme(panel.background = element_rect(fill = "grey95", colour = "grey55", linewidth = 0.5), # change 'colour' to black for informative axis
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=8),
          axis.text.x = element_text(angle = 90, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          legend.position="right",
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )
  
ggsave("CD3_gated.jpg", g, width = 20, height = 16, units = "cm", dpi = 300)




#-----------------
# setup CD45 highpass
#-----------------

u <- list()
xMarker<-"CD45..Cell..Mean_t2"
yMarker<-"E.cadherin..Cell..Mean_t2"


colorMarker <- "CD68..Cell..Mean_t2"

lower.threshold <- as.numeric(quantile( cell.dat[,  ..colorMarker ], probs = c( 0.02),na.rm=TRUE ))
median <- as.numeric(quantile( cell.dat[,  ..colorMarker ], probs = c( 0.5),na.rm=TRUE ))
upper.threshold <- as.numeric(quantile( cell.dat[,  ..colorMarker ], probs = c( 0.98),na.rm=TRUE ))

u[[1]] <- ggplot( data=cell.dat ,
          aes( x= CD45..Cell..Mean_t2,   y= E.cadherin..Cell..Mean_t2 , color = CD68..Cell..Mean_t2 ) ) + 
    
    geom_point(#aes(color = Annotated_metacluster),
      alpha=1, 
      size=0.2 #, 
      # pch='.'
    )+

    geom_path(data = polygon.gates[[3]],  size = 1.1, colour = "#A100AB", alpha=1) +
    
    labs(title=paste0( "All cells") , 
         subtitle="Gating into remaining immune cells (colored by 2-98th precentile)",
          x = xMarker,
          y = yMarker
    )+
    
  scico::scale_color_scico(palette = "roma", direction=-1, 
                          oob=squish,
                           midpoint = median,
                           limits=c(lower.threshold,upper.threshold),
                          # na.value = "#8C0172",
                            name = colorMarker
  )+ # bilbao roma
  
   # zoom a bit in to get the outliers out of the pic 
  coord_cartesian(ylim=c(as.numeric(quantile(cell.dat[[yMarker]], 0.001,na.rm=TRUE)),as.numeric(quantile(cell.dat[[yMarker]], 0.9995,na.rm=TRUE))), 
                  xlim=c(as.numeric(quantile(cell.dat[[xMarker]], 0.001,na.rm=TRUE)),as.numeric(quantile(cell.dat[[xMarker]], 0.9995,na.rm=TRUE)))
                  )+

    scale_x_log10(breaks = minor_breaks, 
                  minor_breaks = minor_breaks#,
                  #labels = trans_format("log10", math_format(10^.x))
                  ) +
    scale_y_log10(breaks = minor_breaks, 
                  minor_breaks = minor_breaks#,
                  #labels = trans_format("log10", math_format(10^.x))
                  ) +
    
    stat_density_2d(#aes( fill = ..level..),        #fill = ..level..), 
      alpha = 0.2, 
      contour = TRUE,
      contour_var = "density",
      adjust=1.7,
      n = 500, # default 100, computationally intense, do not exceed that value too much..
      geom = "polygon", #geom: polygon or raster (equidistant lines or heatmap-ish picture)
      colour="white") +
    
    theme(panel.background = element_rect(fill = "grey95", colour = "grey55", linewidth = 0.5), # change 'colour' to black for informative axis
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=8),
          axis.text.x = element_text(angle = 90, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          legend.position="right",
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )
  


#yMarker<-"CD11b..Cell..Mean_t2"
colorMarker <- "CD11c..Cell..Mean_t2"

lower.threshold <- as.numeric(quantile( cell.dat[,  ..colorMarker ], probs = c( 0.02),na.rm=TRUE ))
median <- as.numeric(quantile( cell.dat[,  ..colorMarker ], probs = c( 0.5),na.rm=TRUE ))
upper.threshold <- as.numeric(quantile( cell.dat[,  ..colorMarker ], probs = c( 0.98),na.rm=TRUE ))


u[[2]] <- ggplot( data=cell.dat ,
          aes( x= CD45..Cell..Mean_t2,   y= E.cadherin..Cell..Mean_t2 , color = CD11c..Cell..Mean_t2 ) ) + # color =CD68..Mean as.factor(outside.CD3e.highpass)
    
    geom_point(#aes(color = Annotated_metacluster),
      alpha=1, 
      size=0.2 #, 
      # pch='.'
    )+

    geom_path(data = polygon.gates[[3]],  size = 1.1, colour = "#A100AB", alpha=1) +
    
    labs(title=paste0( "All cells") , 
         subtitle="Gating into remaining immune cells (colored by 2-98th precentile)",
          x = xMarker,
          y = yMarker
    )+
    
 #   scale_color_manual(name = "Gate flag", 
#                       guide =  guide_legend(override.aes = list(size=4,alpha=1 )),
#                       values = c("#40B0A6", "#E1BE6A"), #colorfriendly. not this: c("#308441", "#F6CD7A")
#                       labels = c("in", "out") #0 blue in, 1 yellow out
#    )+  
  
  scico::scale_color_scico(palette = "roma", direction=-1, 
                          oob=squish,
                           midpoint = median,
                           limits=c(lower.threshold,upper.threshold),
                          # na.value = "#8C0172",
                            name = colorMarker
  )+ # bilbao roma
  
   # zoom a bit in to get the outliers out of the pic 
  coord_cartesian(ylim=c(as.numeric(quantile(cell.dat[[yMarker]], 0.001,na.rm=TRUE)),as.numeric(quantile(cell.dat[[yMarker]], 0.9995,na.rm=TRUE))), 
                  xlim=c(as.numeric(quantile(cell.dat[[xMarker]], 0.001,na.rm=TRUE)),as.numeric(quantile(cell.dat[[xMarker]], 0.9995,na.rm=TRUE)))
                  )+

    scale_x_log10(breaks = minor_breaks, 
                  minor_breaks = minor_breaks#,
                  #labels = trans_format("log10", math_format(10^.x))
                  ) +
    scale_y_log10(breaks = minor_breaks, 
                  minor_breaks = minor_breaks#,
                  #labels = trans_format("log10", math_format(10^.x))
                  ) +
    
    stat_density_2d(#aes( fill = ..level..),        #fill = ..level..), 
      alpha = 0.2, 
      contour = TRUE,
      contour_var = "density",
      adjust=1.7,
      n = 500, # default 100, computationally intense, do not exceed that value too much..
      geom = "polygon", #geom: polygon or raster (equidistant lines or heatmap-ish picture)
      colour="white") +
    
    theme(panel.background = element_rect(fill = "grey95", colour = "grey55", linewidth = 0.5), # change 'colour' to black for informative axis
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=8),
          axis.text.x = element_text(angle = 90, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          legend.position="right",
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )
  



  g <- do.call("arrangeGrob", c(u, ncol=2)) # use a grob to save the obj
ggsave("CD45.jpg", g, width = 40, height = 16, units = "cm", dpi = 300)


g <-  ggplot( data=cell.dat ,
          aes( x= CD45..Cell..Mean_t2,   y= E.cadherin..Cell..Mean_t2 , color = as.factor(outside.CD45.highpass)  ) ) + # color =CD68..Mean as.factor(outside.CD3e.highpass)
    
    geom_point(#aes(color = Annotated_metacluster),
      alpha=1, 
      size=0.2 #, 
      # pch='.'
    )+

    geom_path(data = polygon.gates[[3]],  size = 1.1, colour = "#A100AB", alpha=1) +
    
    labs(title=paste0( "All cells") , 
         subtitle="Gating into T cells (colored by 2-98th precentile)",
          x = xMarker,
          y = yMarker
    )+
    
    scale_color_manual(name = "Gate flag", 
                       guide =  guide_legend(override.aes = list(size=4,alpha=1 )),
                       values = c("#40B0A6", "#E1BE6A"), #colorfriendly. not this: c("#308441", "#F6CD7A")
                       labels = c("in", "out") #0 blue in, 1 yellow out
    )+  
  

  
   # zoom a bit in to get the outliers out of the pic 
  coord_cartesian(ylim=c(as.numeric(quantile(cell.dat[[yMarker]], 0.001,na.rm=TRUE)),as.numeric(quantile(cell.dat[[yMarker]], 0.9995,na.rm=TRUE))), 
                  xlim=c(as.numeric(quantile(cell.dat[[xMarker]], 0.001,na.rm=TRUE)),as.numeric(quantile(cell.dat[[xMarker]], 0.9995,na.rm=TRUE)))
                  )+

    scale_x_log10(breaks = minor_breaks, 
                  minor_breaks = minor_breaks#,
                  #labels = trans_format("log10", math_format(10^.x))
                  ) +
    scale_y_log10(breaks = minor_breaks, 
                  minor_breaks = minor_breaks#,
                  #labels = trans_format("log10", math_format(10^.x))
                  ) +
    
    stat_density_2d(#aes( fill = ..level..),        #fill = ..level..), 
      alpha = 0.2, 
      contour = TRUE,
      contour_var = "density",
      adjust=1.7,
      n = 500, # default 100, computationally intense, do not exceed that value too much..
      geom = "polygon", #geom: polygon or raster (equidistant lines or heatmap-ish picture)
      colour="white") +
    
    theme(panel.background = element_rect(fill = "grey95", colour = "grey55", linewidth = 0.5), # change 'colour' to black for informative axis
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=8),
          axis.text.x = element_text(angle = 90, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          legend.position="right",
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )
  
ggsave("CD45_gated.jpg", g, width = 20, height = 16, units = "cm", dpi = 300)

















g <- ggplot( data=cell.dat ,
          aes( x= cell.dat[[xMarker]],   y= cell.dat[[yMarker]] , color = Object.cross.section ) ) + 
    
    geom_point(#aes(color = Annotated_metacluster),
      alpha=1, 
      size=0.2 #, 
      # pch='.'
    )+
    #geom_path(data = polygon.gates[[2]],  size = 1.1, colour = "#A100AB", alpha=1) +
    
    labs(title=paste0( "All cells") , 
         subtitle="Gating into T cells",
          x = xMarker,
          y = yMarker
    )+
    
    scale_color_manual(name = "Object cros-section", 
                       guide =  guide_legend(override.aes = list(size=4,alpha=1 )),
                       values = c( "#E1BE6A","#40B0A6")
                       #labels = c("in", "out") #0 blue in, 1 yellow out
    )+  
  

  
   # zoom a bit in to get the outliers out of the pic 
  coord_cartesian(ylim=c(as.numeric(quantile(cell.dat[[yMarker]], 0.001,na.rm=TRUE)),as.numeric(quantile(cell.dat[[yMarker]], 0.9995,na.rm=TRUE))), 
                  xlim=c(as.numeric(quantile(cell.dat[[xMarker]], 0.001,na.rm=TRUE)),as.numeric(quantile(cell.dat[[xMarker]], 0.9995,na.rm=TRUE)))
                  )+

    scale_x_log10(breaks = breaks, 
                  minor_breaks = minor_breaks,
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = breaks, 
                  minor_breaks = minor_breaks,
                  labels = trans_format("log10", math_format(10^.x))) +
    
    stat_density_2d(#aes( fill = ..level..),        #fill = ..level..), 
      alpha = 0.2, 
      contour = TRUE,
      contour_var = "density",
      adjust=1.7,
      n = 500, # default 100, computationally intense, do not exceed that value too much..
      geom = "polygon", #geom: polygon or raster (equidistant lines or heatmap-ish picture)
      colour="white") +
    
    theme(panel.background = element_rect(fill = "grey95", colour = "grey55", linewidth = 0.5), # change 'colour' to black for informative axis
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=8),
          #axis.text.x = element_text(angle = 45, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          legend.position="right",
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )
ggsave("Cell_profiles.jpg", g, width = 20, height = 16, units = "cm", dpi = 300)


g<- ggplot(subset(cell.dat) ,
             aes(x=Centroid.X.µm,y=Centroid.Y.µm,color=Object.cross.section,fill=Object.cross.section )) + #,color=as.factor(inside.lymphnode)
    geom_point(size=0.1, alpha=1) + 
    facet_wrap(~Image, ncol=4) + #, scales="free"
    theme_minimal() +
    scale_y_reverse()+
    coord_fixed() + 
    # extract Tissue_type from sample.meta for the current i:
    #labs(title= paste0("Core ", i, " - tissue: ", sample.meta[ sample.meta$TMA.core  == i,]$Tissue_type  )  )+
    
    scale_color_manual(name = "Object cross-section", 
                       guide =  guide_legend(override.aes = list(size=4,alpha=1 )),
                       values = c( "#E1BE6A","#40B0A6")
                       #labels = c("in", "out") #0 blue in, 1 yellow out
    )+  
      scale_fill_manual(name = "Object cross-section", 
                       values = c( "#E1BE6A","#40B0A6")
                       #labels = c("in", "out") #0 blue in, 1 yellow out
    )+  

    theme(panel.background = element_rect(fill = "grey10", colour = "grey15", linewidth = 0.5), # change 'colour' to black for informative axis
          panel.grid = element_line(color = "grey15"),
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=8),
          #axis.text.x = element_text(angle = 45, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          legend.position="left",
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )
ggsave("Cell_profiles.jpg", g, width = 50, height = 12, units = "cm", dpi = 1200)
  
ggplot( data=cell.dat ,
          aes( x= CD45..Cell..Mean_t2,   y= CD20..Cell..Mean_t2 , color = CD3e..Cell..Mean_t2 ) ) +  #   as.factor(outside.CD20.highpass)
  
    geom_point(#aes(color = Annotated_metacluster),
      alpha=1, 
      size=0.8 #, 
      # pch='.'
    )+
    #geom_path(data = polygon.gates[[1]],  size = 1.1, colour = "#A100AB", alpha=1) +
    
    labs(title=paste0( "All cells") , 
         subtitle="Gating into B cells"
    )+
    
    scale_color_manual(name = "Gate flag", 
                       guide =  guide_legend(override.aes = list(size=4,alpha=1 )),
                       values = c("#40B0A6", "#E1BE6A")#, #colorfriendly. not this: c("#308441", "#F6CD7A")
                       #labels = c("in", "out") #0 blue in, 1 yellow out
    )+ #, 
    
    scale_x_log10(breaks = breaks, 
                  minor_breaks = minor_breaks,
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = breaks, 
                  minor_breaks = minor_breaks,
                  labels = trans_format("log10", math_format(10^.x))) +
    
    
    stat_density_2d(#aes( fill = ..level..),        #fill = ..level..), 
      alpha = 0.2, 
      contour = TRUE,
      contour_var = "density",
      adjust=1.7,
      n = 500, # default 100, computationally intense, do not exceed that value too much..
      geom = "polygon", #geom: polygon or raster (equidistant lines or heatmap-ish picture)
      colour="white") +
    
    theme(panel.background = element_rect(fill = "grey95", colour = "grey55", linewidth = 0.5), # change 'colour' to black for informative axis
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=8),
          #axis.text.x = element_text(angle = 45, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          legend.position="right",
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )
  
ggsave("CD20_gate.jpg", g, width = 20, height = 16, units = "cm", dpi = 300)
  
  
  
 
  
  
  
  
  # since we already have a gating for B and T, we can take advantage now of that when we want to gate the remaining immune subset out.
# they are currently in G3 with the rest of the cells:


# we might need to clip the marker color:

xMarker <- "CD45..Mean"
yMarker <- "CD4..Mean"
colorMarker <- "Ki67..Mean"



percentile <- as.numeric(quantile( cell.dat$Ki67..Mean, # we compute the quantiles from non-ctr data
                                   probs = c( 0.02 , 0.5, 0.98),na.rm=TRUE )) # ,na.rm=TRUE can cause some weird behaviour of limits
lower.threshold <- percentile[1]
median <- percentile[2]
upper.threshold <- percentile[3]
  
#scico::scico(5, palette = 'roma')
#"#7E1700" "#B58B31" "#C0E9C2" "#389CC6" "#023198"
#scico::scico(7, palette = 'hawaii')
#"#8C0172" "#93384C" "#996330" "#9B951B" "#7FC55F" "#5EE2BB" "#B2F2FD"


g<- ggplot( data= cell.dat ,
        aes( x= CD45..Mean,   y= CD4..Mean , color =Ki67..Mean ) ) + 
  
  geom_point(#aes(color = Annotated_metacluster),
    alpha=1, 
    size=0.2 #, 
     #pch='.'
  )+
  #geom_path(data = polygon.gates[[2]],  size = 1.1, colour = "#A100AB", alpha=1) +
  
facet_grid(. ~ Image)+
  
  labs(title=paste0( "All cells") , 
       subtitle="Color gradient clipped at 2 and 98th percentile"
  )+
  
#  scale_color_manual(name = "Gate flag", 
#                     guide =  guide_legend(override.aes = list(size=4,alpha=1 )),
#                     values = c("#40B0A6", "#E1BE6A"), #colorfriendly. not this: c("#308441", "#F6CD7A")
#                     labels = c("in", "out") #0 blue in, 1 yellow out
#  )+ #, 
 scico::scale_color_scico(palette = "roma", direction=-1, 
                          midpoint = median,
                          limits=c(lower.threshold,upper.threshold),
                          na.value = "#8C0172"
                          )+ # bilbao roma
  
  #geom_path(data = polygon.gates[[3]],  size = 1.1, colour = "#A100AB", alpha=1) +
  
#  scale_color_gradient2(low="#5EE2BB", 
#                        mid="#9B951B",  
#                        high="#8C0172", 
#                        midpoint = median, #space='Lab' other values are depreciated
#                        limits=c(lower.threshold,upper.threshold), 
#                        #limits=c(quartiles[1],quartiles[3]), #oob=squish, # logical: squish outofbound vals into limits
#                        na.value = "#7E1700" )+
  
  coord_cartesian(ylim=c(as.numeric(quantile(cell.dat[[yMarker]], 0.005,na.rm=TRUE)),as.numeric(quantile(cell.dat[[yMarker]], 0.9998,na.rm=TRUE))), 
                  xlim=c(as.numeric(quantile(cell.dat[[xMarker]], 0.001,na.rm=TRUE)),as.numeric(quantile(cell.dat[[xMarker]], 0.9998,na.rm=TRUE)))
  )+
  
  scale_x_log10(breaks = breaks, 
                minor_breaks = minor_breaks,
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = breaks, 
                minor_breaks = minor_breaks,
                labels = trans_format("log10", math_format(10^.x))) +
  
  
  stat_density_2d(#aes( fill = ..level..),        #fill = ..level..), 
    alpha = 0.1, 
    contour = TRUE,
    contour_var = "density",
    adjust=1.7,
    n = 500, # default 100, computationally intense, do not exceed that value too much..
    geom = "polygon", #geom: polygon or raster (equidistant lines or heatmap-ish picture)
    colour="white") +
  
  theme(panel.background = element_rect(fill = "grey95", colour = "grey55", linewidth = 0.5), # change 'colour' to black for informative axis
        axis.title.x=element_text(color="grey15", size=11),
        axis.title.y=element_text(color="grey15", size=11),
        axis.text=element_text(size=8),
        #axis.text.x = element_text(angle = 45, hjust=1),
        #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
        legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
        legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
        #legend.title=element_blank(),
        legend.position="right",
        plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
        plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
  )
ggsave("CD45_CD4_Ki67_bySlide.jpg", g, width = 25, height = 16, units = "cm", dpi = 300)

g<- ggplot( data= cell.dat ,
            aes( x= CD45..Mean,   y= CD4..Mean , color =TMA.core ) ) + 
  
  geom_point(#aes(color = Annotated_metacluster),
    alpha=1, 
     size=0.2 #, 
    #pch='.'
  )+
  #geom_path(data = polygon.gates[[2]],  size = 1.1, colour = "#A100AB", alpha=1) +
  
  facet_grid(. ~ Image)+
  
  labs(title=paste0( "All cells") , 
       subtitle="split by both cores on each slide"
  )+
  
   scale_color_manual(name = "TMA core", 
                       guide =  guide_legend(override.aes = list(size=4,alpha=1 )),
                       values = c("#0ed2dc", "#d51ec7")#, #colorfriendly. not this: c("#308441", "#F6CD7A")
                       #labels = c("in", "out") #0 blue in, 1 yellow out
    )+ #, 
  #scico::scale_color_scico(palette = "roma", direction=-1, 
  #                         midpoint = median,
  #                         limits=c(lower.threshold,upper.threshold),
  #                         na.value = "#8C0172"
  #)+ # bilbao roma
  
  #geom_path(data = polygon.gates[[3]],  size = 1.1, colour = "#A100AB", alpha=1) +
  
  #  scale_color_gradient2(low="#5EE2BB", 
  #                        mid="#9B951B",  
  #                        high="#8C0172", 
  #                        midpoint = median, #space='Lab' other values are depreciated
  #                        limits=c(lower.threshold,upper.threshold), 
  #                        #limits=c(quartiles[1],quartiles[3]), #oob=squish, # logical: squish outofbound vals into limits
  #                        na.value = "#7E1700" )+
  
  coord_cartesian(ylim=c(as.numeric(quantile(cell.dat[[yMarker]], 0.005,na.rm=TRUE)),as.numeric(quantile(cell.dat[[yMarker]], 0.9998,na.rm=TRUE))), 
                  xlim=c(as.numeric(quantile(cell.dat[[xMarker]], 0.001,na.rm=TRUE)),as.numeric(quantile(cell.dat[[xMarker]], 0.9998,na.rm=TRUE)))
  )+
  
  scale_x_log10(breaks = breaks, 
                minor_breaks = minor_breaks,
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = breaks, 
                minor_breaks = minor_breaks,
                labels = trans_format("log10", math_format(10^.x))) +
  
  
  stat_density_2d(#aes( fill = ..level..),        #fill = ..level..), 
    alpha = 0.1, 
    contour = TRUE,
    contour_var = "density",
    adjust=1.7,
    n = 500, # default 100, computationally intense, do not exceed that value too much..
    geom = "polygon", #geom: polygon or raster (equidistant lines or heatmap-ish picture)
    colour="white") +
  
  theme(panel.background = element_rect(fill = "grey95", colour = "grey55", linewidth = 0.5), # change 'colour' to black for informative axis
        axis.title.x=element_text(color="grey15", size=11),
        axis.title.y=element_text(color="grey15", size=11),
        axis.text=element_text(size=8),
        #axis.text.x = element_text(angle = 45, hjust=1),
        #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
        legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
        legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
        #legend.title=element_blank(),
        legend.position="right",
        plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
        plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
  )
ggsave("CD45_CD4_core_bySlide.jpg", g, width = 25, height = 16, units = "cm", dpi = 300)


# Define bin width
binwidth <- 0.1


# Filter out invalid values for log10
cell.dat <- cell.dat %>%
  dplyr::filter(!is.na(CD45..Mean), CD45..Mean > 0) %>%
  dplyr::mutate(
    log_val = log10(CD45..Mean),
    bin = cut(
      log_val,
      breaks = seq(
        floor(min(log_val, na.rm = TRUE)),
        ceiling(max(log_val, na.rm = TRUE)),
        by = binwidth
      ),
      include.lowest = TRUE
    )
  )

bin_summary <- cell.dat %>%
  group_by(bin, Image) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(bin) %>%
  mutate(fraction = n / sum(n))

ggplot(bin_summary, aes(x = bin, y = fraction, fill = Image)) +
  geom_col(position = "fill", color = "black") +  # "fill" makes it proportional (0 to 1)
  labs(title = "Fraction of cell per condition CD45..Mean Bin",
       x = "CD45..Mean (log scale bins)",
       y = "Fraction",
       fill = "Condition") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(cell.dat, aes(x = CD45..Mean)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", alpha = 0.7) +
  scale_x_log10() +  # Logarithmic scale for x-axis
  labs(title = "Histogram of CD45..Mean Expression",
       x = "CD45..Mean (log scale)",
       y = "Count") +
  theme_minimal()

} # end run.experimental.code






###################### QC: gates and segmentation penalties ######################


# at this point it would make sense to have a look at gate_affected.markers. 
# until now, we just ignored plotting the gate lines if both markers were part of the scatter plot. which does not make much sense,
# it would be better to plot two plots, rather than plotting one scatter without gate lines...
all.gate.flags <- grep("outside.", names(cell.dat), value = TRUE ) 
all.gate.flags <- grep(".area", all.gate.flags, invert= TRUE, value = TRUE )

# we already got the gate-affected markers: gate_affected.markers
# and we already wrote out the threshold gate column names: threshold.gates

# what we miss is a way to run through all.gate.flags and find the threshold by which they were created.
# if the flag was created with a threshold, you would find a column missing "outside." :
all.marker.thresholds <-  sub("^.*?\\.", "", all.gate.flags) 
# to keep the start of the sentence you need to do it this way around: ^, 
#see https://stackoverflow.com/questions/72966851/r-how-to-extract-everything-after-first-occurance-of-a-dot

# outside.HLA.DR.highpass -> "outside.HLA.DR"
# this chops off outside. sub(".[^.]+$", "", all.gate.flags  ) 


if(plot.TMAwise.gates == 1){
  
  message( paste0("\nQC 1: Check Gates on TMA-wise scatter plots: ")  )  
  setwd(OutputDirectory)
  dir.create('1 - QC Gates on TMAs')
  setwd('1 - QC Gates on TMAs')

  # for the progess bar we need to differenciate whether we gonna loop through the gate flags at all or not:
  
  if( identical(all.gate.flags, character(0)) ){  
    message("There are no flags present in cell.dat. You have to set them either in the metadata or as multi-segmented lines")
    pb <- progress_bar$new(format = "[:bar] :percent [Plotting marker scatter | :eta]",
                           total = length(markercomparision)*length(all_images) , 
                           show_after=0, #allows to call it right way 
                           current = "|",    # Current bar character
                           clear = TRUE) # make it persist
    pb$tick(0) # call in the progress bar without any progress, just to show it
    
  }else{
    
    pb <- progress_bar$new(format = "[:bar] :percent [Plotting marker scatter | :eta]",
                           total =length(markercomparision)*length(all_images)*length(all.gate.flags)*2 , # *2 cause we run over the tickmark of the first condion in this if-else
                           show_after=0, #allows to call it right way 
                           current = "|",    # Current bar character
                           clear = T) # make it persist
    pb$tick(0) # call in the progress bar without any progress, just to show it
    
    
  }
  
  
  
for(m in c(1:length(markercomparision))){  
  
# scatter plotting routine is fully automated
# here, the scatter.plots object only contained the markers, but we created markercomparison on the fly with raw data, since at this point we dont have transformed signals:
  xMarker <-  markercomparision[[m]][1] 
  yMarker <-  markercomparision[[m]][2] 
  colorMarker <-  markercomparision[[m]][3]   
  
  
  
  # set up lyrics and progress bar for every marker anew:
  # set up lyrics and progress bar for every marker anew:
  if(!is.na(colorMarker)){
    # Compute the coloring of the dots with the aid of a third variable in the ROI 
    percentile <- as.numeric(quantile( subset(cell.dat, tissue.type.col %!in% postive.control.Tissue.name )[[markercomparision[[m]][3]]], # we compute the quantiles from non-ctr data
                                       probs = c( (1-percentile.colorlimits) , 0.5, percentile.colorlimits))) # ,na.rm=TRUE can cause some weird behaviour of limits
    lower.threshold <- percentile[1]
    median <- percentile[2]
    upper.threshold <- percentile[3]
    # message( paste0(m, " of ", length(markercomparision), ": ",    xMarker," against ",yMarker, ", colored by ",colorMarker )  )
  }else{# message( paste0(m, " of ", length(markercomparision), ": ",    xMarker," against ",yMarker, ", no third coloring marker chosen" )  )
      }
  
  
  l = list()
  list.pos <- 1
 
for(i in all_images ){
  
#you could also compute the color gradient OT-wise here: 
# Compute the coloring of the dots with the aid of a third variable in the ROI 
#  percentile <- as.numeric(quantile(temp.OTdata$Insulin_Pr141, 
#                                    probs = c( (1-percentile.colorlimits) , 0.5, percentile.colorlimits))) # ,na.rm=TRUE can cause some weird behaviour of limits
#  lower.threshold <- percentile[1]
#  median <- percentile[2]
#  upper.threshold <- percentile[3]


#watch out!! this only works if every ROI in TMA i has the same gate set!!

# and it only works if the Batch is subset by the variable "i" since plot.scatter.plot will use "i" to pull the correct line in scaling.factor if needed:


# for the first pass, no gates set, all.gate.flags is character(0) and our f loop does not start!! 
if( identical(all.gate.flags, character(0)) ){
  
  
  
  # in this case we just push the marker comparison right into the plotting engine and plot without flag colors neither gate lines:
  
  gate.printing.needed <- 0
  
#  l[[list.pos]] <- 
#    plot.scatter.gate.v2(indata=subset(cell.dat, Batch==i & Tissue %!in% postive.control.Tissue.name ) , 
#                       inx=markercomparision[[m]][1] , 
#                       iny=markercomparision[[m]][2], 
#                       incol= markercomparision[[m]][3],
#                       dot.size = 1.2, # default 1
#                       dot.alpha = 0.9, # default 0.2
#                       gate.flag = f,
#                       gate.line.alpha=1,
#                       passed.gate.color = "#308441",
#                       rejected.gate.color = "#F6CD7A"   # default "#EB8563"
#  )
#  list.pos <- list.pos+1

  
  l[[list.pos]] <-  
    
    plot.scatter.density.v2(indata=subset(cell.dat, Image==i & Tissue %!in% postive.control.Tissue.name ) , 
                            inx=markercomparision[[m]][1] , 
                            iny=markercomparision[[m]][2], 
                            incol= markercomparision[[m]][3],
                            dot.size = 0.08, # default 1
                            dot.alpha = 0.9, # default 0.2
                            gate.flag = f,
                            gate.line.alpha=1,
                            passed.gate.color = "#308441",
                            rejected.gate.color = "#F6CD7A",   # default "#EB8563"
                            titlestring=paste0(i," entire TMA incl. pos-ctr ",postive.control.Tissue.name)
                            #subtitlestring=paste0("Entire data set incl pos-ctr ", nrow(indata), " cells, incl pos-ctr")  # defaults to paste0("Plotted ", nrow(indata), " cells") 
    )
  
  list.pos <- list.pos+1
  
}else{
  
  # all.gate.flags has found flags, so lets plot accordingly:
  
  
  for(f in all.gate.flags ){
    
    # So: Rather than asking if the markercomparison are part of any gate, we do it the other way around: 
    # We run through all gates and check if they should be printed in the current marker comparison. 
    # This has the advantage that I can check the gate separately, and decide if it needs a simple hline/vline plotted, or a more complex multi-segmented line.
    # By just checking if the markercomparison is part of a flag, I dont know what method created that flag in the first place.
    
    # all.gate.flags are set and you can run through the flags and reconstruct the threshold origin (metadata or polygon)  
    temp.gated.marker <- sub(".[^.]+$", "",  sub("^.*?\\.", "", f)   ) 
    
    # check if the currently scanned gate is part of scatter  
    # the cool thing is that we can simply fall back to markercomparision[[m]][1] and markercomparision[[m]][2] which 
    # are the naked AB names:
    
    gate.printing.needed <- ifelse(  temp.gated.marker %in% c(scatter.plots[[m]][1],scatter.plots[[m]][2]) , 1, 0 )  
    
    # now we decide if the gate needs printing:
    if(gate.printing.needed==1 | gate.printing.needed==0){
      # check if the threshold comes from metadata or polygon
      is.metadata.threshold <- ifelse(  sub("^.*?\\.", "", f) %in% threshold.gates, 1, 0 ) # if 1, this needs a hline/vline, depending on where the marker is plotted
      
      is.highpass <-  ifelse(   grepl(".highpass", f)   , 1, 0 ) 
    

    
   l[[list.pos]] <-  
    
    plot.scatter.density.v2(indata=subset(cell.dat, Image==i & tissue.type.col %!in% postive.control.Tissue.name ) , 
      inx=markercomparision[[m]][1] , 
      iny=markercomparision[[m]][2], 
      incol= markercomparision[[m]][3],
      dot.size = 0.08, # default 1
      dot.alpha = 0.9, # default 0.2
      gate.flag = f,
      gate.line.alpha=1,
      passed.gate.color = "#308441",
      rejected.gate.color = "#F6CD7A"   # default "#EB8563" 
    )
   
   list.pos <- list.pos+1
  
    }# end gate.printing is needed
    
   
    pb$tick() # this is the ticker for "gate flags found, we need to check all gate drawings
    
  }#end run with f along all.gate.flags
  }# end plotting: all.gate.flags has found flags in cell.dat

  
  #this is the ticker if no gate flags were found and we never run through the gate thresholds to draw them
  # however, if we have flags set, the loop nevertheless runs over this tickmark. this causes the *2 when defining pb more up
  pb$tick()
  
    } # end run along all TMAs to assemble to plots of the current marker comparison, including color-coding the flagged cells
    
 # after all subset ROIs are in the grob, store it with the markercomposition in the file name: 
 # message("\nArranging the plots, be patient...\n")
  
  # we gonna try to arrange them on a 3x2 ratio
  # since we add rows rather than cols:
  
  g <- do.call("grid.arrange", c(l, nrow = ceiling( (list.pos-1)*2/6 ) ))
  filename <- paste0("QC_scatterplot_TMAwise_",xMarker,"vs",yMarker,"_color",colorMarker,"_noAreaGate_noMarkerGate.png")
  ggsave(filename,g, width = 3*18*ceiling( (list.pos-1)*2/6 )   , height = 18*ceiling( (list.pos-1)*2/6 ) , units = "cm",dpi = 350,limitsize = FALSE)
 
  invisible(dev.off())
  

  
   
}# end run along all marker combinations we need to plot per TMA i
  
} # END QC plotting: TMA-wise gate positions and flags

  
if(plot.segementation.areas == 1){
  
  setwd(OutputDirectory)
  dir.create('2 - QC Cell Segment Areas on TMAs')
  setwd('2 - QC Cell Segment Areas on TMAs')
  
  
#................QC: cellArea Outliers..................
# break the x axis string into two lines:

g <- 
  ggplot(data=cell.dat, aes(x=gsub("(.{15,}?)_", "\\1\n", cell.dat$Image), y= Cell..Area.µm.2,color=as.factor(outside.area) )) + #, colour=Batch 
  
geom_quasirandom( size = 2 ) + #alpha=0.8
 
  scale_color_manual(name = "Area outlier flag", values = c("#308441", "#F6CD7A"))+ #, 
  
  # if you did not do object classifier, 
  # this column will split by nucleated and non-nucleated cells:
  facet_wrap(~Classification)+ 
  theme_classic()+
  ylab( "Segmentation area [log10(1/px)]") + 
  labs(title=paste0("Segmentation areas of ",nrow(cell.dat)," cells" )  , 
       subtitle="Entire dataset, no gate" ,
       #caption="Created by M.Barone", 
       # y="Maximal Velocity [AU/sec]", 
         x="TMA"
  )+
  #facet_wrap(~Tissue)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  
  annotation_logticks(sides = "l",) + 
  
  theme(#axis.title.x = element_blank(), 
    #legend.title=element_text(size=13), # turn off with element_blank(),
    #legend.text=element_text(size=11),
    axis.text=element_text(size=13),#,axis.text.x=element_text(size=10),
    axis.title.y=element_text(size=12), axis.title.x=element_text(size=15),
    legend.position="right", 
    # legend.key.width=unit(1,"line")
    axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1),
    #these are the facet settings:
    strip.text.x = element_text(size=12,  face="bold"),
    strip.background = element_rect(colour="black", fill="grey90", size=1, linetype="solid"),
  ) 

message("\nQC 2: Plotting cell segmentation area outliers...")
filename <- "QC_SegementationAreas.pdf" 
ggsave(filename,g, width = 45, height = 20, units = "cm",dpi = 1200)






} # END QC plotting: segmentation areas




   

        

  

###################### Drop pos-ctr. tissue and segmentation area outlier for good ###################### 
# in the first advanced script, we kicked the pos-ctr. out for clustering. so lets do that as well:  
        
# This is an idea I brought back from shinyIMC:   
# lets create a second cell.dat that is gated and subsettable. like this, we can change things here and do not loose the big dataset:
               
              
cells.before.outlierdrop <- nrow(cell.dat)     
cells.area.outliers <- table(cell.dat$outside.area)[2]   
#in case none dropped, we just created an NA:
if( is.na(cells.area.outliers)) cells.area.outliers<-0




# drop segmentation outliers, pos control ROIs and rename object for upcoming gating:
#gated.cell.dat <- subset(cell.dat,  outside.area==0 & tissue.type.col %!in% postive.control.Tissue.name )     
gated.cell.dat <- subset(cell.dat,   tissue.type.col %!in% postive.control.Tissue.name )  
# store a list of all ROIs away so that it is accessible for the TMA-choice slider further down the page
plot.ROIs  <- unique(gated.cell.dat$TMA.core) 
        
cells.after.outlierdrop <- nrow(gated.cell.dat) 


#
        

message( paste0("\nCleared data from ",postive.control.Tissue.name," and segmentation area penalties:\n",

                cells.area.outliers, " cells (", round( 100-(cells.before.outlierdrop- cells.area.outliers)/cells.before.outlierdrop*100   ,3)   ,"%) dropped out due to segmentation area penalties\n",        
                
                
                cells.after.outlierdrop, " cells remain for gating (", round(  cells.after.outlierdrop / cells.before.outlierdrop *100   ,2), "% of the pre-filtered dataset)" 
                
               
                
)  )      
   
       
        

        


        

###################### SET UP GATING ROUTINE (second part) ######################

# Since we got the flags of the thresholds and the multi-segmented line all set up, we can push them all into our hierarchical.gating.routine function in one go

# Adapt this function to set the boolean gates of every subset:
# Set the conditions with care and cross-check on the plots created in QC

hierarchical.gating.routine <- function(outside.CD20.highpass , 
                                        outside.CD3.highpass , 
                                       outside.CD45.highpass 
                                         ){
  
  if( outside.CD20.highpass==0){
    #if gate condition is met, flag the cell for population 1:
    n <- 1
  }else if( outside.CD3.highpass==0 ){
    n <- 2
  }else if( outside.CD45.highpass==0 ){
    n <- 3
  }else{
    #non immune cells, CD19- CD45-
   n <- 4
  }#else{
    # trash bin subset: this one should be empty
  #  n <- 5
 # }
  
  return(n)
} # end function definition of the hierarchical gates


# mapply now passes these gate column flags to our boolean gating routine and deposits the GATEsubset numbers:        
gated.cell.dat$GATEsubset <-mapply(hierarchical.gating.routine,
                                   outside.CD20.highpass=gated.cell.dat$outside.CD20.highpass, 
                                   outside.CD3.highpass=gated.cell.dat$outside.CD3e.highpass,
                                   outside.CD45.highpass=gated.cell.dat$outside.CD45.highpass
                                  
)






# For the hierarchical gating and also the upcoming clustering, we simply gonna subset in a loop for the n GATEsubset(s) we stored in gated.cell.dat:
amountofGATEsubsets <- max(gated.cell.dat$GATEsubset)-exclude.trash.bin.gate  # pre-defining this here allows us to take out the "rest" subset





        
 
        
# this message is automated by now:
summedupcells<-0
message( "\n------ Checking if hierarical gates catch all cells ------")

        
for(s in 1:amountofGATEsubsets){
          message( paste0(  nrow(subset(gated.cell.dat,  GATEsubset==s )), " cells in subset ", s   ) )
          summedupcells <- summedupcells+nrow(subset(gated.cell.dat,  GATEsubset==s ))
        }
message( paste0("gated.cell.dat contains totally ", nrow(gated.cell.dat), " cells," ) )        
message( paste0("of which ", summedupcells,  " cells are gated into ",s," selected subsets"  ) )
      
        

        


if(plot.scatterplots.gatingsteps==1){

# We gonna run now two plotting routines: 
        # 1) pre-gated subsets
        # 2) gates subsets
    # the first is interesting to see if the gates fit at all, since this is how the dataset looks right before gating into a subset. 
    # So thats basically the bounced back cells from the prior gate drawn on the upcoming gate
      
    # the second is then a look at the actual composition of the gated subsets
              
        

###################### QC pre-gated subsets on TMAs ###################### 
  
  setwd(OutputDirectory)
  dir.create('3 - QC Gating steps')
  setwd('3 - QC Gating steps')
 
  
# so here we need first all three subsets, then subset number one out, then subset number 2 out
        
        # since we got the GATEsubset variable in gated.cell.dat, we can now simply run through it and check on the rawdata how the gates look TMA-wise:
        message("\nQC 3: Plotting the the cells affected by an upcoming gating step..." )  
        

        #this time we start with the full dataset...
        temp.OTdata <- gated.cell.dat
        
        #...and take twice a subset out, namely number one and number two.
        for(s in 0:(amountofGATEsubsets-1)  ){
          
          # the trick is to start at GATEsubset 0 which is nonexistent, so subsetting for NOT 0 does nothing the first time we run through
          # the next time we gonna take subset 1 out and store temp.OTdata away:
          temp.OTdata <- subset( temp.OTdata, GATEsubset!=s )
          message(paste0("\n------------ GATE step ",s, " -----------\n      Subsets present: "))
          cat( sort( unique( temp.OTdata$GATEsubset) ) )
          message(paste0("\n      ",nrow(temp.OTdata)," cells"   ))
          
          
          

          #if you put the list here, you will get every gating step separately. so you get all marker comparisons with every gate on one plot
          if (scatterplots.gatingsteps.packaging == "step") {
            l = list()
            list.pos <- 1
          }
          


          
          
          
          
          # the first plot will split the dataset into markers, and will store away a scatter and density plot into l. This are currently three marker combos, 6 plots per subset
          # it will do that for all three subsets and in the end will assemble the plot, so in total 6x3 plots stored in l
          
          for(m in c(1:length(markercomparision))){  
            
       
            
            #for the entire DS we need to re-calculate the color gradient:
            percentile <- as.numeric(quantile(temp.OTdata[[markercomparision[[m]][3]]], 
                                              probs = c( (1-percentile.colorlimits) , 0.5, percentile.colorlimits) ) # ,na.rm=TRUE can cause some weird behaviour of limits
                                     ) 
            lower.threshold <- percentile[1]
            median <- percentile[2]
            upper.threshold <- percentile[3]
            
            
            xMarker <-  markercomparision[[m]][1] 
            yMarker <-  markercomparision[[m]][2] 
            colorMarker <-  markercomparision[[m]][3]  
            
            # set up lyrics and progress bar for every marker anew:
            

            
              #if you put the list here, you will get every marker separately. so you get more than one plot per gating step
            if (scatterplots.gatingsteps.packaging == "marker") {
             l = list()
             list.pos <- 1
            }
            
            
            
            # now we also make the TMA wise plots with the subset:
            for(i in all_images ){
              
              #we further split the current GATEsubset into each TMA:
              temp.TMAwise.OTdata <- subset(temp.OTdata, Image==i )  
              
              # WATCH OUT. There is a good chance that for small GATE subsets, you will end up with some TMA that got no cells.
              # To protect from plotting 0 cells:
              
              if(nrow(temp.TMAwise.OTdata)>1  ){

                #you could also compute the color gradient OT-wise here: 
                # Compute the coloring of the dots with the aid of a third variable in the ROI 
                #  percentile <- as.numeric(quantile(temp.OTdata$Insulin_Pr141, 
                #                                    probs = c( (1-percentile.colorlimits) , 0.5, percentile.colorlimits))) # ,na.rm=TRUE can cause some weird behaviour of limits
                #  lower.threshold <- percentile[1]
                #  median <- percentile[2]
                #  upper.threshold <- percentile[3]
                
                
                #watch out!! this only works if every ROI in TMA i has the same gate set!!
                #plotting CD45 gates
               # t[[TMAlist.pos]] <-
              #    
              #    plot.scatter.gate(indata=temp.TMAwise.OTdata, 
              #                     inx=markercomparision[[m]][1] , 
               #                    iny=markercomparision[[m]][2], 
               #                    incol= markercomparision[[m]][3],
               #                    titlestring=paste0(i, ", about to gate into: G", (s+1) ) 
               #   )
               ## 
              #  TMAlist.pos <- TMAlist.pos+1
                
                
                for(f in all.gate.flags ){
                  
                  # So: Rather than asking if the markercomparison are part of any gate, we do it the other way around: 
                  # We run through all gates and check if they should be printed in the current marker comparison. 
                  # This has the advantage that I can check the gate separately, and decide if it needs a simple hline/vline plotted, or a more complex multi-segmented line.
                  # By just checking if the markercomparison is part of a flag, I dont know what method created that flag in the first place.
                  
                  # all.gate.flags are set and you can run through the flags and reconstruct the threshold origin (metadata or polygon)  
                  temp.gated.marker <- sub(".[^.]+$", "",  sub("^.*?\\.", "", f)   ) 
                  
                  # check if the currently scanned gate is part of scatter  
                  # the cool thing is that we can simply fall back to markercomparision[[m]][1] and markercomparision[[m]][2] which 
                  # are the naked AB names:
                  gate.printing.needed <- ifelse(  temp.gated.marker %in% c(scatter.plots[[m]][1],scatter.plots[[m]][2]) , 1, 0 ) 
                  
                  # now we decide if the gate needs printing:
                  if(gate.printing.needed==1){
                    # check if the threshold comes from metadata or polygon
                    is.metadata.threshold <- ifelse(  sub("^.*?\\.", "", f) %in% threshold.gates, 1, 0 ) # if 1, this needs a hline/vline, depending on where the marker is plotted
                    
                    is.highpass <-  ifelse(   grepl(".highpass", f)   , 1, 0 ) 
               
                
                
                l[[list.pos]] <-    plot.scatter.density.v2(indata=temp.TMAwise.OTdata, 
                                                                                   inx=markercomparision[[m]][1] , 
                                                                                   iny=markercomparision[[m]][2], 
                                                                                   incol= markercomparision[[m]][3],
                                                                                   dot.size = 0.08, # default NULL for faster plotting by drawing real dots. 
                                                                                   dot.alpha = 0.9, # default 0.2
                                                                                   gate.flag = f,
                                                                                   gate.line.alpha=1,
                                                                                   passed.gate.color = "#308441",
                                                                                   rejected.gate.color = "#F6CD7A",   # default "#EB8563"
                                                                                   titlestring = paste0("About to gate into: G", (s+1) )
                                                           )
                
                
                list.pos <- list.pos+1
                
                  }# end gate.printing is needed
                }#end run with f along all.gate.flags
                
                
                # END protect from running into plotting TMA-wise scatter and densities if not enough cells are present for yticks$start in temp.TMAwise.OTdata 
              }else{
                message(paste0(i, " had only ",  nrow(temp.TMAwise.OTdata), " cells -> skipped the plotting"))}
              
              
              
              
            } # end run along all OTs and plot the gates
            

            
            
            
            
       
              
              if (scatterplots.gatingsteps.packaging == "marker") {
                
                # after all subset ROIs are in the grob, store it with the markercomposition in the file name: 
                
                message(paste0("Assembling scatter plots of ",xMarker,"vs",yMarker, " in gating step ",s, " into one object") ) 
                
                g <- do.call("grid.arrange", c(l, nrow = ceiling( (list.pos-1)*2/6 ) ))
                filename <- paste0("QC_scatterplot_TMAwise_gatingSTEP_",s,"_",xMarker,"vs",yMarker,"_color",colorMarker,".png")
                ggsave(filename,g, width = 1.5*18*ceiling( (list.pos-1)*2/6 )   , height = 18*ceiling( (list.pos-1)*2/6 ) , units = "cm",dpi = 350,limitsize = FALSE)
                
                invisible(dev.off())

              }
            

            
            
            
            
      
            

            
          } # end run along m the markers
          
          if (scatterplots.gatingsteps.packaging == "step") {
            message(paste0("Assembling all scatter plots of gating step ",s, " into one object") ) 
            
            g <- do.call("grid.arrange", c(l, nrow = ceiling( (list.pos-1)*2/6 ) ))
            filename <- paste0("QC_scatterplot_gatingSTEP_",s,".png")
            ggsave(filename,g, width = 1.5*18*ceiling( (list.pos-1)*2/6 )   , height = 18*ceiling( (list.pos-1)*2/6 ) , units = "cm",dpi = 350,limitsize = FALSE)
            
            invisible(dev.off())
            
          }
          
        } # end subsetting the dataset for pre-gate 
        
        

        
        
        
if(run.experimental.code==1){
  # watch out, QC4 returns an error since the 
  
  #if( identical(all.gate.flags, character(0)) ){
    # in this case we just push the marker comparison right into the plotting engine and plot without flag colors neither gate lines:
   #    gate.printing.needed <- 0
  # ... is not set correctly

        
        
######################  QC gated subsets on TMAs ###################### 
 
        setwd(OutputDirectory)
        dir.create('4 - QC Gated subsets')
        setwd('4 - QC Gated subsets')   
     
            
# since we got the GATEsubset variable in gated.cell.dat, we can now simply run through it and check on the rawdata how the gates look TMA-wise:
message("\nQC 4: Plotting the gated subsets with their gates..." )  
        
list.pos <- 1
l = list()   
        
        
for(s in 1:amountofGATEsubsets){
          
temp.OTdata <- subset( gated.cell.dat, GATEsubset==s )
          
message(paste0("\n------------ Subset Number ",s, " -----------\n      ",GATEsubsetNames[s],"\n      ",nrow(temp.OTdata)," cells"   ))
         
 
     
     
# the first plot will split the dataset into markers, and will store away a scatter and density plot into l. This are currently three marker combos, 6 plots per subset
# it will do that for all three subsets and in the end will assemble the plot, so in total 6x3 plots stored in l

for(m in c(1:length(markercomparision))){  
  
  TMAlist.pos <- 1
  t<- list()
  
  #for the entire DS we need to re-calculate the color gradient:
  percentile <- as.numeric(quantile(temp.OTdata[[markercomparision[[m]][3]]], 
                                    probs = c( (1-percentile.colorlimits) , 0.5, percentile.colorlimits))) 
  lower.threshold <- percentile[1]
  median <- percentile[2]
  upper.threshold <- percentile[3]
  
  
  xMarker <-  markercomparision[[m]][1] 
  yMarker <-  markercomparision[[m]][2] 
  colorMarker <-  markercomparision[[m]][3]  
  

  for(f in all.gate.flags ){
    
    # So: Rather than asking if the markercomparison are part of any gate, we do it the other way around: 
    # We run through all gates and check if they should be printed in the current marker comparison. 
    # This has the advantage that I can check the gate separately, and decide if it needs a simple hline/vline plotted, or a more complex multi-segmented line.
    # By just checking if the markercomparison is part of a flag, I dont know what method created that flag in the first place.
    
    # all.gate.flags are set and you can run through the flags and reconstruct the threshold origin (metadata or polygon)  
    temp.gated.marker <- sub(".[^.]+$", "",  sub("^.*?\\.", "", f)   ) 
    
    # check if the currently scanned gate is part of scatter  
    gate.printing.needed <- ifelse(  temp.gated.marker %in% c(xMarker,yMarker) , 1, 0 )  
    
    # now we decide if the gate needs printing:
    if(gate.printing.needed==1){
      # check if the threshold comes from metadata or polygon
      is.metadata.threshold <- ifelse(  sub("^.*?\\.", "", f) %in% threshold.gates, 1, 0 ) # if 1, this needs a hline/vline, depending on where the marker is plotted
      
      is.highpass <-  ifelse(   grepl(".highpass", f)   , 1, 0 ) 
      
      
  
  l[[list.pos]] <-  plot.scatter.density.v2(indata=temp.OTdata, 
                          inx=markercomparision[[m]][1] , 
                          iny=markercomparision[[m]][2], 
                          incol= markercomparision[[m]][3],
                          dot.size = 0.08, # default 1
                          dot.alpha = 0.9, # default 0.2
                          gate.flag = f,
                          gate.line.alpha=1,
                          passed.gate.color = "#308441",
                          rejected.gate.color = "#F6CD7A",   # default "#EB8563"
                          titlestring=paste0("Subset in: ", GATEsubsetNames[s] ) 
  )
  
  list.pos <- list.pos+1
  
    }# end gate.printing is needed
  }#end run with f along all.gate.flags
  
  
  message( paste0(m, " of ", length(markercomparision), ": ",    xMarker," against ",yMarker, ", colored by ",colorMarker, " (plotting entire dataset)" )  )
  # g <- do.call("grid.arrange", c(l, ncol=2 ))
  # filename <- paste0("QC_1stGATE_scatterplot_non-spleenDS_",xMarker,"vs",yMarker,"_color",colorMarker,"_AreaGate_noMarkerGate.png")
  # ggsave(filename,g, width = 60, height = 30, units = "cm",dpi = 100,limitsize = FALSE)
  
  invisible(dev.off()) 
  
   # now we also make the TMA wise plots with the subset
  # however, if there is only one TMA, re-doing the plot we just did above does not make much sense, so protect:
  if(length(all_images)>1){
  
  for(i in all_images ){
    
    #we further split the current GATEsubset into each TMA:
    temp.TMAwise.OTdata <- subset(temp.OTdata, Image==i )  
    
    # WATCH OUT. There is a good chance that for small GATE subsets, you will end up with some TMA that got no cells.
    # To protect from plotting 0 cells:
    
    if(nrow(temp.TMAwise.OTdata)>1  ){
      
      
      
  
    #you could also compute the color gradient OT-wise here: 
    # Compute the coloring of the dots with the aid of a third variable in the ROI 
    #  percentile <- as.numeric(quantile(temp.OTdata$Insulin_Pr141, 
    #                                    probs = c( (1-percentile.colorlimits) , 0.5, percentile.colorlimits))) # ,na.rm=TRUE can cause some weird behaviour of limits
    #  lower.threshold <- percentile[1]
    #  median <- percentile[2]
    #  upper.threshold <- percentile[3]
    
    
    #watch out!! this only works if every ROI in TMA i has the same gate set!!
    #plotting CD45 gates
   t[[TMAlist.pos]] <- 
      plot.scatter.density.v2(indata=temp.TMAwise.OTdata, 
                                               inx=markercomparision[[m]][1] , 
                                               iny=markercomparision[[m]][2], 
                                               incol= markercomparision[[m]][3],
                                               dot.size = 0.08, # default 1
                                               dot.alpha = 0.9, # default 0.2
                                               gate.flag = f,
                                               gate.line.alpha=1,
                                               passed.gate.color = "#308441",
                                               rejected.gate.color = "#F6CD7A",   # default "#EB8563"
                                              #  titlestring=paste0("Subset in: ", GATEsubsetNames[s] ) # title string defaults to i as all_images object
    )
    
    TMAlist.pos <- TMAlist.pos+1
    
    
    
    
    # END protect from running into plotting TMA-wise scatter and densities if not enough cells are present for yticks$start in temp.TMAwise.OTdata 
    }else{
      message(paste0(i, " had only ",  nrow(temp.TMAwise.OTdata), " cells -> skipped the plotting"))}
    
    
    
    
  } # end run along all OTs and plot the gates
  
  # after all subset ROIs are in the grob, store it with the markercomposition in the file name: 
  
  message( paste0(m, " of ", length(markercomparision), ": ",    xMarker," against ",yMarker, ", colored by ",colorMarker, " (plotting TMA-wise)" )  )
  
  message(paste0("Assembling TMA-wise plots of subset ",GATEsubsetNames[s], "...") ) 
  
  g <- do.call("grid.arrange", c(t, nrow = ceiling( (TMAlist.pos-1)*2/6 ) ))
  filename <- paste0("QC_scatterplot_TMAwise_",GATEsubsetNames[s],"_",xMarker,"vs",yMarker,"_color",colorMarker,".png")
  ggsave(filename,g, width = 1.5*18*ceiling( (list.pos-1)*2/6 )   , height = 18*ceiling( (list.pos-1)*2/6 ) , units = "cm",dpi = 350,limitsize = FALSE)
  
  invisible(dev.off())
  
  

  
  }#end protecting from plotting TMA-wise if the dataset contains only one TMA in total
  
}# end run along all markers we need to plot per markercombo

    
} # end run with s through all gatesubsets in gated.cell.dat
        
        
#after the whole gating was simulated, lets plot them here:
message("Assembling hierarchical gating plot...")   

g <- do.call("grid.arrange", c(l, nrow = ceiling( (list.pos-1)*2/6 ) ))
filename <- paste0("QC_hierarchicalGATES_scatterplots.png")
ggsave(filename,g, width = 1.5*18*ceiling( (list.pos-1)*2/6 )   , height = 18*ceiling( (list.pos-1)*2/6 ) , units = "cm",dpi = 100,limitsize = FALSE)

invisible(dev.off())

}#protect from running QC4
  
        
        
        
rm(temp.OTdata)  
rm(temp.TMAwise.OTdata)
}#end plotting machine for the gating steps themselves ("About to gate into subset n")        

       

        
        
        
        
        
        
  


        
        
        
        
        
########## Only now that we gated and got all relevant cells left, we can do transformation and normalization:
        
        
        message("\nTranforming and normalizing gated.cell.dat..." )   
        
        
        # this function does all in one go. I skipped these long names as they annoyed me.
        # Function is fully defined in the STAGE 1 SETTINGS on top.
        gated.cell.dat <-  do.data.normalization(dat=gated.cell.dat, 
                                           use.cols=cellular.cols, 
                                           # Transform
                                           do.transform = transform.rawdata,
                                           cofactor = asinh.cofactor,
                                           # min-max normalize?
                                           do.minmax = minmax.norm,
                                           new.min = min.norm, 
                                           new.max = max.norm,
                                           # z-score normalize?
                                           do.zscore = zscore.norm
                                           ) 
        

# watch out, if you want to calculate cluster densities rel. to the total segementation area, 
#do NOT use gated.cell.dat of the later stages since this does not contain the trashbin gate anymore, possibly the majority of all cells on that ROI!!
# for that we need to create a separate dataframe:
# at this point in the code, gated.cell.dat and all.cells is the same dataframe (gated.cell.dat containing flags and transformed channels, but we are interested in Area)        
all.cells <- gated.cell.dat
# watch out, all.cells cannot be made if you re-load an old gated.cell.dat. you need to manually load cell.dat, kick ROIs out and make the outside.area flag from within STAGE1
setwd(OutputDataDirectory)
fwrite( all.cells , 'all.cells.csv')        
      
        


# time to prepare the marker list to cluster the subsets via markers.to.cluster.per.subset

# displaying the selected markers for each subset:       

message("Cellular cols available:")
print(as.matrix(names(gated.cell.dat)))
message("Cellular cols currently selected for clustering the subsets:")



# we first write out the GATEsubsetNames
GATEsubsetNames <- unlist(markers.to.cluster.per.subset[c(T, F)])


GATEsubset.Markers2Cluster <- markers.to.cluster.per.subset[c(F, T)]
names(GATEsubset.Markers2Cluster) <- GATEsubsetNames # why not...

# create the raw data columns:
# nope! we order the markers along global marker order first, then do this thing here:
#GATEsubset.Markers2Cluster <- lapply(GATEsubset.Markers2Cluster, function(x) ifelse(is.na(x), NA, paste0(x, qupath.separator,cellular.segment.readout)  ))
 


# we just use a simple for loop to build up a list that pulls the column names:
        cluster.cols.GATEsubsets <- list(  )
        for(g in 1:length(GATEsubset.Markers2Cluster) ){
          
          
          # if you have use.global.marker.order set, lets order them:
           if( use.global.marker.order == "yes"){
             
          # global.marker.order is just the markers, pre-order them....
           cluster.cols.GATEsubsets[[g]] <- global.marker.order[  global.marker.order  %in%  GATEsubset.Markers2Cluster[[g]]  ] 
           # ...and now create the raw data column name:
           cluster.cols.GATEsubsets[[g]] <- unlist( lapply(cluster.cols.GATEsubsets[[g]], function(x)  paste0(x, qupath.separator,cellular.segment.readout)  )) 
          
           # depending on the processing request, the column now gets suffixes:
           if(transform.rawdata==TRUE){cluster.cols.GATEsubsets[[g]] <- unlist( lapply(cluster.cols.GATEsubsets[[g]], function(x)  paste0(x, "_t",asinh.cofactor)  ))    }    
           if(minmax.norm==TRUE){  cluster.cols.GATEsubsets[[g]] <- unlist( lapply(cluster.cols.GATEsubsets[[g]], function(x)  paste0(x, "_m",min.norm,"m",max.norm)  ))    }      
           if(zscore.norm==TRUE){  cluster.cols.GATEsubsets[[g]] <- unlist( lapply(cluster.cols.GATEsubsets[[g]], function(x)  paste0(x,  "_z")  ))        }
           
         }else{
           cluster.cols.GATEsubsets[[g]] <- names(cell.dat %>% dplyr:: select(starts_with( paste0(GATEsubset.Markers2Cluster[[g]], "_") )))
           # depending on the processing request, the column now gets suffixes:
           if(transform.rawdata==TRUE){cluster.cols.GATEsubsets[[g]] <- paste0(cluster.cols.GATEsubsets[[g]], "_t",asinh.cofactor)  }
           if(minmax.norm==TRUE){  cluster.cols.GATEsubsets[[g]] <- paste0(cluster.cols.GATEsubsets[[g]], , "_m",min.norm,"m",max.norm)  }
           if(zscore.norm==TRUE){  cluster.cols.GATEsubsets[[g]] <- paste0(cluster.cols.GATEsubsets[[g]],  "_z")  }
           
        
         }
          
        }# end assembling cluster.cols.GATEsubsets list
        names(cluster.cols.GATEsubsets) <- GATEsubsetNames
        print(cluster.cols.GATEsubsets)
        
        

     
# in the end of STAGE1, we set a counter for the upcoming stage:        
how.often.ran.STAGE2 <- 0 


cat("\n\n\n----------------------------------- \n End of STAGE 1. \n Gate subsets defined: using gated.cell.dat \n Proceed to STAGE 2 for clustering \n-----------------------------------\n\n")


# clear from elements we dont need anymore:
rm(summedupcells)  
} # end stage 1: allocate the markers to cluster        
 


########################## STAGE 2 ##########################   

if(whichSTAGE == 2){ 
 

# From now on we wont do any raw data processing, so let us create a small string that is plotted in all upcoming graphs, to keep track of the processing data:
  data.source.lyrics <- c("Data source: ")
  
  if(correct.spill.over==1 ){
    data.source.lyrics <- paste0(data.source.lyrics, "spillover-corr., ")
  }
  
  if(calculate.scaling.factors==1 ){
    data.source.lyrics <- paste0(data.source.lyrics, "batch-norm., ")
  }
  
  if(transform.rawdata ){
    data.source.lyrics <- paste0(data.source.lyrics, "asinh_cf",asinh.cofactor)
  }
  if(minmax.norm){
    data.source.lyrics <- paste0(data.source.lyrics, ", minmax norm.")
  }
  if(zscore.norm){
    data.source.lyrics <- paste0(data.source.lyrics, ", z-score norm.")
  }
  
  
  

 if(run.old.code == 1){
  # ---- rerun assembly of cluster.cols.GATEsubsets
  # watch out! 
  # if you re-run STAGE2 after modifying the markerlist you use to cluster, the following block that ran in STAGE1 needs to be re-executed to update the cluster list:
  # from the STAGE 1 settings, we have already the col numbers we pull from gated.cell.dat
  # this list needs to be first sliced, and all odd entries thrown out
  GATEsubset.Markers2Cluster <- GATEsubset.Markers2Cluster[c(FALSE, TRUE)]
  names(GATEsubset.Markers2Cluster) <- GATEsubsetNames # why not...
  
  
  # before we load the markers to cluster, lets set up the global marker order:
  
  # if cellular.cols is already in the proper order:
  if(global.marker.order.is.ordered == "yes"){
    global.marker.order <- cellular.cols
  }else{
    # watch out, at this part in the code, you cannot just pull starts_with from cell.dat since the normalization and transformation already happend!!
    # so you need to order cellular.cols by global.marker.order
    
    
    
    cellular.cols %>% dplyr:: select(starts_with( paste0(global.marker.order, "_") ) )
    
    
    global.marker.order <- cellular.cols[ match(global.marker.order, cellular.cols) ]
    
    
    x[order(match(x,y))]
    
    
    
  }
  
  # now lets assemble to proper column name for the clustering:

  
  global.marker.order <- names(cell.dat %>% dplyr:: select(starts_with( paste0(global.marker.order, "_") )))
  
  if(transform.rawdata==TRUE){global.marker.order <- paste0(global.marker.order, "_t",asinh.cofactor)  }
  if(minmax.norm==TRUE){  global.marker.order <- paste0(global.marker.order, , "_m",min.norm,"m",max.norm)  }
  if(zscore.norm==TRUE){  global.marker.order <- paste0(global.marker.order,  "_z")  }
  
  
  
  
  
  # we just use a simple for loop to build up a list that pulls the column names:
  cluster.cols.GATEsubsets <- list(  )
  for(g in 1:length(GATEsubset.Markers2Cluster)){
    
    # if you have use.global.marker.order set, lets order them:
    if( use.global.marker.order == "yes"){
      cluster.cols.GATEsubsets[[g]] <- global.marker.order[  global.marker.order  %in%  names(gated.cell.dat)[GATEsubset.Markers2Cluster[[g]]] ] 
    }else{
      cluster.cols.GATEsubsets[[g]] <- names(gated.cell.dat)[GATEsubset.Markers2Cluster[[g]]]
    }
    
  }# end assembling cluster.cols.GATEsubsets list
  names(cluster.cols.GATEsubsets) <- GATEsubsetNames
  # ---- end rerun assembly of cluster.cols.GATEsubsets
  
 }#protect from running old code
  


  
  
    setwd(OutputDirectory)
  
  # here is the try out to return to the settings of STAGE1
  # rm("clustered.gated.subsets")
  # how.often.ran.STAGE2 <- 0

  
    
    
######----------  AFTER FIRST PASS CLUSTERING: Remove clusters according to the delete.clustersinsubsets  list defined in STAGE2 settings ----------######     

  if( exists("clustered.gated.subsets") ){
    
    # this means you ran STAGE2 already, and you might run it a second time to clear out clusters
    
    dir.create('STAGE2 Second Pass')
    setwd('STAGE2 Second Pass')
    Output.STAGE2.Directory <- getwd()
    
    
    # we need to make sure that we do not keep on deleting cluster numbers as we re-run and re-run this STAGE!
    # So after the first time, the counter will read 1, thats when we are allowed to delete numbers. 
    # Afterwards the counter will be higher and that will block continous deleting of said cluster numbers
    
    if("Phenograph_metacluster" %in% colnames(clustered.gated.subsets[[1]] ) & how.often.ran.STAGE2==1) {
      
    # watch out, this is risky stuff here, so I want the user to confirm the cluster delete here:
      if (confirm.cluster.delete.after.first.run.STAGE2 == "yes") {
      
      message("You are running STAGE 2 the second time\nI will delete the chosen clusters out of clustered.gated.subsets now:")
      
        for(x in 1:length(delete.clustersinsubsets) ) {
          temp.subset <- delete.clustersinsubsets[[x]][1]
          temp.clustertodelete <- delete.clustersinsubsets[[x]][2]
          clustered.gated.subsets[[ temp.subset  ]] <- subset( clustered.gated.subsets[[ temp.subset  ]]  , !(Phenograph_metacluster== temp.clustertodelete  ) ) 
          message(paste0("Removed cluster ", temp.clustertodelete, " from subset ",temp.subset ))
          rm(temp.subset ,temp.clustertodelete )
        }
        
      
      }# end protect deleting clusters if user did not confirm them
      
      # You might have chosen some markers to fish out bad cells using "untypical" markers. 
      # Now that the clusters are deleted, you might want these markers gone. In that case, re-define the cluster markers within here:
      
    #  cluster.cols.GATEsubsets <- list(
    #    names(gated.cell.dat)[c(78:83,85:93,96:98,102:104)], # 84 is CD19, insulin marker out! Now also 99 (GLP) and 105 (Nkx) goes
    #    names(gated.cell.dat)[c(78:83,85:93,96:98,102:104)], # 84 is CD19, insulin marker out! Now also 99 (GLP) and 105 (Nkx) goes
    #    names(gated.cell.dat)[c(78:83,85:93,96:98,102:104)] # 84 is CD19, insulin marker out! Now also 99 (GLP) and 105 (Nkx) goes
    #  )
      
      # and you also need to update the plotting order of course in that case!

      #marker.order 
      
      
      how.often.ran.STAGE2 <- how.often.ran.STAGE2+1 # this pushed the counter above 1 for sure
    }else{
      
      message("You are running STAGE 2 more than two times\n NO CLUSTERS will be deleted anymore")
      how.often.ran.STAGE2 <- how.often.ran.STAGE2+1
    }
    
    
    
    
  }else{
    
    # Welcome to STAGE2. clustered.gated.subsets does not exist, so lets get to work with UMAPs and HM:
    message("You are running STAGE 2 the first time\nTo delete clusters from subsets: use the delete.clustersinsubsets list set and confirm deletion in STAGE 2 SETTINGS.")
    dir.create('STAGE2 First Pass')
    setwd('STAGE2 First Pass')
    Output.STAGE2.Directory <- getwd()
    
    
     # split gated.cell.dat and rename object
  clustered.gated.subsets <- split(gated.cell.dat, f = gated.cell.dat$GATEsubset) 
  how.often.ran.STAGE2 <- how.often.ran.STAGE2+1
  
    
  }
  

  
  
  
# even thought we also split the "rest" or trashbin subset into that element, we not gonna call that since we run only until amountofGATEsubsets, which takes into account 
  #whether you want to also cluster this huge rest subset.
  
  
  ## Again, this still needs manual intervention.
  # if you ran STAGE 2 the first time, you might want to kick out clusters that are clearly misplaced cells


 
#-------  Clustering and dimensionality reduction engine start ------------ 
for(s in 1:amountofGATEsubsets  ){
  
  # so, while this split was a cool way forward, we cannot run here with s in names(clustered.gated.subsets), since the list names are characters.
  # better stay safe with a 1:amountofGATEsubsets variable that is numeric for sure and will be able to pull out the according things from other lists and vectors

 
  # load the cluster columns from our list for the according subset, and set the marker.order for that subset:
  cluster.cols <- cluster.cols.GATEsubsets[[s]] # use double brackets to subset into the vector and not just the list.
  marker.order <- cluster.cols.GATEsubsets[[s]] # depreciated since we took care of that: global.marker.order[global.marker.order %in% cluster.cols]
  

  #cluster.cols.GATEsubsets
  
  
  
  message(paste0("\n------------ Subset Number ",s, " -----------\n      ",GATEsubsetNames[s],"\n      ",nrow(clustered.gated.subsets[[s]])," cells"   ))
  message(paste0("Using the following markers to cluster:"   ))
  print(as.matrix(cluster.cols))
  
  

    if("Phenograph_metacluster" %in% colnames(clustered.gated.subsets[[s]] ) ) {
      
      message(  "\n-------   Clustering startup   --------\nWill remove ", max( clustered.gated.subsets[[s]]$Phenograph_metacluster)," clusters from clustered.gated.cell.dat before running"  )
      clustered.gated.subsets[[s]] <-  subset(clustered.gated.subsets[[s]], select = -c(Phenograph_metacluster) )
    }else{
      message(  "\n-------   Clustering startup   --------\nNo clusters present in clustered.gated.cell.dat before running"  )
    }
    
    
    
    # create the Phenograph clusters:
    
    
    if(cluster.engine=="Rphenoannoy"){
      set.seed(123)  
      Rphenograph_cluster <- Rphenoannoy(data=as.matrix(clustered.gated.subsets[[s]][, cluster.cols, with = FALSE]), k=phenograph.k)
      Rphenograph_cluster <- as.numeric(membership(Rphenograph_cluster[[2]]))
    }
    
    
    
    if(cluster.engine=="cytof"){
      set.seed(123)
      Rphenograph_cluster <- cytof_cluster(xdata = as.matrix(clustered.gated.subsets[[s]][, cluster.cols, with = FALSE]), 
                                           #Rphenograph_k = 30, # depreciated as of V1.4 
                                           method="Rphenograph")
    }
    
    
    
    
    clustered.gated.subsets[[s]][["Phenograph_metacluster"]] <- Rphenograph_cluster
    
    
    
    
 
    
  
# we gonna use the same code as further down when collapsing, so we will overwrite all of these objects again:
    #.................. Prepare colors for our clusters  ...............      
    
    #     CREATE polychrome palette for the unannotated Phenograph MCs
    colRange <- unique(clustered.gated.subsets[[s]][["Phenograph_metacluster"]])
    colRange <- colRange[order(colRange)]
    #colRange <- as.character(colRange)
    set.seed(723451) # for reproducibility
    Phenograph_metacluster_palette <- createPalette(length(colRange), c("#ff0000"), M=100000, prefix = "")
    # so, the cool thing is that once we annotated these clusters, we can just do the same game again, overrride the palette and plot again.
    names(Phenograph_metacluster_palette) <- colRange  
    
    
    
    
    #NEW:.................. MAKE COMPLEX HEATMAP  ...............    
    
    exp <- clustered.gated.subsets[[s]][, lapply(.SD, "mean", na.rm = TRUE), by='Phenograph_metacluster',  
                                        .SDcols = cluster.cols]
    
    ## new way to plot complex heatmap
    
    #z-normalize feature expression (iterate over columns with MARGIN=2)
    zscored.exp <- cbind(  exp[,1]  ,   apply(exp[,-1], scale, MARGIN = 2)  ) 
    
    # we gonna bring back only the signals into the matrix:
    zscored.exp.mat <- as.matrix(zscored.exp[,-1])
    
    # prepare column and row names:
    rownames(zscored.exp.mat) <- paste0(exp[[1]] )
    colnames(zscored.exp.mat) <- sub("_.*", "", colnames(zscored.exp.mat) )
    
    col_fun = colorRamp2(c(range(zscored.exp.mat)[1], 0, range(zscored.exp.mat)[2]), c("blue", "white", "red"))
    
    
    
    column_ha = HeatmapAnnotation(cluster = rownames(zscored.exp.mat),
                                  
                                  annotation_legend_param = list(
                                    cluster = list(
                                      ncol = 2, 
                                      title = "Cluster",
                                      title_position = "topcenter",
                                      at = names(Phenograph_metacluster_palette),
                                      grid_height = unit(0.02*length(marker.order), "cm"),
                                      grid_width = unit(0.02*length(marker.order), 'cm'),
                                      labels_gp = gpar(fontsize = 0.8*length(marker.order)),
                                      title_gp = gpar(fontsize = 0.8*length(marker.order))
                                    )
                                  ),
                                  
                                  col = list(cluster=Phenograph_metacluster_palette),
                                  
                                  na_col = "black", #white
                                  show_annotation_name = FALSE
    )
    
    hm <-   Heatmap( t( zscored.exp.mat ) ,
                     col = col_fun,
                     row_order = sub("_.*", "", marker.order ) , # this is the only that needs stripping of transformed_rescaled
                     cluster_columns = fh,
                     top_annotation = column_ha,
                     rect_gp = gpar(col = "white", lwd = 1),
                     heatmap_legend_param = list(
                       title = "z-score",
                       direction = 'horizontal',
                       title_position = "topcenter",
                       legend_width = unit(0.25*length(marker.order), "cm"),
                       grid_width = unit(0.02*length(marker.order), 'cm'),
                       labels_gp = gpar(fontsize = 0.8*length(marker.order)),
                       title_gp = gpar(fontsize = 0.8*length(marker.order))
                     )
                     
    )
    

    
 # gonna store the hm away as grob and unite it with the UMAP that is not created at this point:
    hm.grob = grid.grabExpr(
      draw(hm,
           heatmap_legend_side = 'top',
           row_sub_title_side = 'left',
           padding = unit(c(2, 2, 2, 10), "mm"))
    )    
    
    
    
    
  
  
  
  
  # ........ UMAP ........
  if(run.UMAP == 1){  

    umap.coords <- umap(as.matrix(clustered.gated.subsets[[s]][, ..cluster.cols]), 
                        #umap(as.matrix(clustered.gated.subsets[[s]][, cluster.cols, with = FALSE]), 
                        n_neighbors = 20, # 16. shinyIMC looked cool with 16 - 0.03. I had here 10 - 0.01 to really separate the clusters
                        min_dist = 0.08,  #0.08
                        metric = 'euclidean',
                        verbose= TRUE
                        ) #n_neighbors = 20, min_dist = 0.1, metric = 'euclidean'
    
    message("UMAP calculation done")
    
    
    # we push the umap coordinates right back into the list of dataframes:
    
    # watch out, if you source uwot before spectre, you need this layout:
    #colnames(umap.coords$layout)[1:2] <- paste0('UMAP', 1:2)
    #clustered.gated.subsets[[s]][["UMAP1"]] <- umap.coords$layout[, 1] #umap.coords[, 1]
    #clustered.gated.subsets[[s]][["UMAP2"]] <- umap.coords$layout[, 2] #umap.coords[, 2]
    
    #if you source uwot after, the umap.coords element looks different
    colnames(umap.coords)[1:2] <- paste0('UMAP', 1:2)
    clustered.gated.subsets[[s]][["UMAP1"]] <- umap.coords[, 1]
    clustered.gated.subsets[[s]][["UMAP2"]] <- umap.coords[, 2]
    
    
    

    # we have a temporary palette at hand:
    
    g <-  make.colour.plot( clustered.gated.subsets[[s]] , 'UMAP1', 'UMAP2', "Phenograph_metacluster", 'factor', 
                            titlestring = paste0("UMAP by unannotated Phenograph MC")  ,
                            subtitlestring = paste0("Subset: ", GATEsubsetNames[s], " (",  nrow(clustered.gated.subsets[[s]]), " cells plotted)" ),
                            point.alpha = 0.4, # lets not overdraw the UMAP and reduce alpha of the geom_points
                            dot.size = 0.9, # less alpha and smaller dots on the half-width-page graphs
                            colours = "polychrome",
                            polychromepalette = Phenograph_metacluster_palette,
                            add.label = TRUE, #col.axis are the Phenograph/FlowSOM metacluster numbers. Print these pls
                            legend.loc = "top",
                            legend.text.size = 8,
                            save.to.disk = F)
    # filename <- paste0(GATEsubsetNames[s],"_UMAP_unannotatedPhenoclusters.png")    
    #  ggsave(filename,g, width = 30, height = 20, units = "cm",dpi = 1200)
    
    
    
  }#end UMAP block  
  
  
 # now that the UMAP is stored as grob, lets created the graph:
  
  
  hm.umap.fig <-  ggarrange(hm.grob, g, widths = c(3,2) )  
  


  
  hm.umap.fig <-   annotate_figure(hm.umap.fig,
                                   top = text_grob(paste0("Subset: ", GATEsubsetNames[s]), color = "black", face = "bold", size = 17),
                                   bottom = text_grob(data.source.lyrics, color = "red",
                                                      hjust = 1, x = 1, face = "italic", size = 13),
                                   #left = text_grob("Figure arranged using ggpubr", color = "green", rot = 90),
                                   #right = text_grob(bquote("Superscript: ("*kg~NH[3]~ha^-1~yr^-1*")"), rot = 90),
                                   fig.lab = "Uncollapsed unannotated", fig.lab.face = "bold"
  ) 
  
  
  # make sure, the GATEsubset names dont contain funny characters:
dev.off()
  height = (length(marker.order))
  width = 1.6*height # its not 3/2 bang on, since there is a legend in between...
  
  if(min(width, height)<4){
    width = 10*width
    height = width
  }
  
  # now, dpi must not exceeed 50000px. adjust dpi so that we end up with a 15000px image on the long edge:
  dpi <- 6000/max(width, height)
  
  ggsave(  plot = hm.umap.fig, 
          filename = paste0(gsub("[:/\\\\<>|?*]", "_", GATEsubsetNames[s]),"_UMAP_Phenograph.jpg"), 
          width = width, 
          height = height ,
          limitsize = FALSE,
         dpi = dpi)
  
  

  
  
  

  
  
  
  
  
  
} # end run with s through all gatesubsets in gated.cell.dat


# after we did all the magic in the gated subsets, lets just store this list away:
  
  
  setwd(OutputDirectory)
  dir.create('Data')
  setwd('Data')
  OutputDataDirectory <- getwd()
  qsave(clustered.gated.subsets, "clustered.gated.subsets.qs")
  
  #before we mess around with collapsing clusters, lets backup the dataset
  bkup.clustered.gated.subsets <- clustered.gated.subsets
  #thanks.
  
  

  message("Clustering/DimRed done, gated subsets stored as backup as clustered.gated.subsets.qs in Data folder.")


  
  
  
  
  
  
  
######################  process the pos-ctr. if asked ######################   
  if(cell.analyis.spleen==1){


    message("Clustering/DimRed started for pos-ctr...")

    
    #spleen.data.frame   cluster.cols.SPLEENsubsets
    Rphenograph_cluster <- cytof_cluster(xdata = as.matrix( spleen.data.frame[, cluster.cols.SPLEENsubsets, with = FALSE] ), 
                                         Rphenograph_k = 30, # metak defaults to 30, which is too high IMHO
                                         method="Rphenograph")
    spleen.data.frame[["Phenograph_metacluster"]] <- Rphenograph_cluster
    
    #     CREATE polychrome palette for the unannotated Phenograph MCs
    colRange <- unique(spleen.data.frame[["Phenograph_metacluster"]])
    colRange <- colRange[order(colRange)]
    #colRange <- as.character(colRange)
    set.seed(723451) # for reproducibility
    spleen_Phenograph_metacluster_palette <- createPalette(length(colRange), c("#ff0000"), M=100000, prefix = "")
    # so, the cool thing is that once we annotated these clusters, we can just do the same game again, overrride the palette and plot again.
    names(spleen_Phenograph_metacluster_palette) <- colRange  
    
    #NEW:.................. MAKE COMPLEX HEATMAP  ...............    
    
    exp <- spleen.data.frame[, lapply(.SD, "mean", na.rm = TRUE), by='Phenograph_metacluster',  
                             .SDcols = cluster.cols]
    
    ## new way to plot complex heatmap
    
    #z-normalize feature expression (iterate over columns with MARGIN=2)
    zscored.exp <- cbind(  exp[,1]  ,   apply(exp[,-1], scale, MARGIN = 2)  ) 
    
    # we gonna bring back only the signals into the matrix:
    zscored.exp.mat <- as.matrix(zscored.exp[,-1])
    
    # prepare column and row names:
    rownames(zscored.exp.mat) <- paste0(exp[[1]] )
    colnames(zscored.exp.mat) <- sub("_.*", "", colnames(zscored.exp.mat) )
    
    col_fun = colorRamp2(c(range(zscored.exp.mat)[1], 0, range(zscored.exp.mat)[2]), c("blue", "white", "red"))
    
    column_ha = HeatmapAnnotation(cluster = rownames(zscored.exp.mat),
                                  
                                  annotation_legend_param = list(
                                    cluster = list(
                                      ncol = 2, 
                                      title = "Cluster",
                                      title_position = "topcenter",
                                      at = names(spleen_Phenograph_metacluster_palette),
                                      grid_height = unit(0.02*length(marker.order), "cm"),
                                      grid_width = unit(0.02*length(marker.order), 'cm'),
                                      labels_gp = gpar(fontsize = 0.8*length(marker.order)),
                                      title_gp = gpar(fontsize = 0.8*length(marker.order))
                                    )
                                  ),
                                  
                                  #col = list(cluster=hm_cols),
                                  col = list(cluster=spleen_Phenograph_metacluster_palette),
                                  
                                  na_col = "black", #white
                                  show_annotation_name = FALSE
    )
    
    hm <-   Heatmap( t( zscored.exp.mat ) ,
                     col = col_fun,
                     row_order = sub("_.*", "", marker.order ) , # this is the only that needs stripping of transformed_rescaled
                     cluster_columns = fh,
                     top_annotation = column_ha,
                     rect_gp = gpar(col = "white", lwd = 1),
                     heatmap_legend_param = list(
                       title = "z-score",
                       direction = 'horizontal',
                       title_position = "topcenter",
                       legend_width = unit(0.25*length(marker.order), "cm"),
                       grid_width = unit(0.02*length(marker.order), 'cm'),
                       labels_gp = gpar(fontsize = 0.8*length(marker.order)),
                       title_gp = gpar(fontsize = 0.8*length(marker.order))
                     )
                     
    )
    
    
    
    # gonna store the hm away as grob and unite it with the UMAP that is not created at this point:
    hm.grob = grid.grabExpr(
      draw(hm,
           heatmap_legend_side = 'top',
           row_sub_title_side = 'left',
           padding = unit(c(2, 2, 2, 10), "mm"))
    )    
    
    umap.coords <- umap(as.matrix( spleen.data.frame[, cluster.cols.SPLEENsubsets, with = FALSE] ), 
                        n_neighbors = 20, # 16. shinyIMC looked cool with 16 - 0.03. I had here 10 - 0.01 to really separate the clusters
                        min_dist = 0.08,  #0.08
                        metric = 'euclidean',
                        verbose= TRUE
    ) #n_neighbors = 20, min_dist = 0.1, metric = 'euclidean'
    
    colnames(umap.coords)[1:2] <- paste0('UMAP', 1:2)
    spleen.data.frame[["UMAP1"]] <- umap.coords[, 1]
    spleen.data.frame[["UMAP2"]] <- umap.coords[, 2]
    
    g <-  make.colour.plot( spleen.data.frame , 'UMAP1', 'UMAP2', "Phenograph_metacluster", 'factor', 
                            titlestring = paste0("UMAP by unannotated Phenograph MC")  ,
                            subtitlestring = paste0("Subset: spleen (",  nrow(spleen.data.frame), " cells plotted)" ),
                            point.alpha = 0.4, # lets not overdraw the UMAP and reduce alpha of the geom_points
                            dot.size = 0.9, # less alpha and smaller dots on the half-width-page graphs
                            colours = "polychrome",
                            polychromepalette = spleen_Phenograph_metacluster_palette,
                            add.label = TRUE, #col.axis are the Phenograph/FlowSOM metacluster numbers. Print these pls
                            legend.loc = "top",
                            save.to.disk = F)
    
    
    hm.umap.fig <-  ggarrange(hm.grob, g, widths = c(3,2) )  
    
    hm.umap.fig <-   annotate_figure(hm.umap.fig,
                                     top = text_grob(paste0("Spleen ROIs"), color = "black", face = "bold", size = 17),
                                     bottom = text_grob(data.source.lyrics, color = "red",
                                                        hjust = 1, x = 1, face = "italic", size = 13),
                                     #left = text_grob("Figure arranged using ggpubr", color = "green", rot = 90),
                                     #right = text_grob(bquote("Superscript: ("*kg~NH[3]~ha^-1~yr^-1*")"), rot = 90),
                                     fig.lab = "Uncollapsed unannotated", fig.lab.face = "bold"
    ) 
    
    
    
    
    png(paste0("Spleen_UMAP_Phenograph.png"), 
        width = 60*(ncol(zscored.exp.mat[,-1])), 
        height = 29*(length(marker.order)) )
    
    print(hm.umap.fig)
    
    invisible(dev.off())
  }#end analyse spleen as asked by cell.analyis.spleen   
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  



# this is a bracket that protects to run the old code without hierarchical gating. It runs gated.cell.dat as base
if(run.old.code == 1){



# ........ FlowSOM ........
  if(run.FlowSOM == 1){
        #this error: cannot take a sample larger than the population means you have not enough data for that dim
  
  if("FlowSOM_metacluster" %in% colnames(gated.cell.dat) ) {
    
    message(  "\n-------   FlowSOM   --------\nGonna remove ", max( gated.cell.dat$FlowSOM_metacluster)," metaclusters from isolated gated.cell.dat before running"  )
    gated.cell.dat <-  subset(gated.cell.dat, select = -c(FlowSOM_cluster) )
    gated.cell.dat <-  subset(gated.cell.dat, select = -c(FlowSOM_metacluster) )
  }else{
    message(  "\n-------   FlowSOM   --------\nNo metaclusters present in isolated gated.cell.dat before running"  )
  }
  

  gated.cell.dat <- run.flowsom(gated.cell.dat, cluster.cols, xdim = 20, ydim = 20, meta.k = 21) #meta.k = 15, you can use "auto", default 14x14
  gated.cell.dat
  
  
  ### Expression heatmap for FlowSOM 
  
  exp <- do.aggregate(gated.cell.dat, cluster.cols, 'FlowSOM_metacluster', 'mean') #cellular.cols
  make.pheatmap(exp, 'FlowSOM_metacluster', cluster.cols,
                transpose=TRUE, # F by default MARKERS are columns, CLUSTERS are rows
                file.name = paste0("Heatmap_FlowSOM_unannotated.png")
  ) #cellular.cols
  rm(exp)
  }#end FlowSOM block
  
  

# ........ Phenograph ........    
  
  if(run.Phenograph == 1){  
  if("Phenograph_metacluster" %in% colnames(gated.cell.dat) ) {
    
    message(  "\n-------   Phenograph   --------\nGonna remove ", max( gated.cell.dat$Phenograph_metacluster)," metaclusters from isolated gated.cell.dat before running"  )
    gated.cell.dat <-  subset(gated.cell.dat, select = -c(Phenograph_metacluster) )
  }else{
    message(  "\n-------   Phenograph   --------\nNo metaclusters present in isolated gated.cell.dat before running"  )
  }
  
  
  
        library(cytofkit) # rphenograph this maskes now Rphenograph function
        # create the Phenograph clusters:
        set.seed(84)
        Rphenograph_cluster <- cytof_cluster(xdata = as.matrix(gated.cell.dat[, cluster.cols, with = FALSE]), 
                                             Rphenograph_k = 20, # metak defaults to 30, which is too high IMHO
                                             method="Rphenograph")
        gated.cell.dat[["Phenograph_metacluster"]] <- Rphenograph_cluster
        
        
        
        
        ### Expression heatmap for Phenograph
        exp <- do.aggregate(gated.cell.dat, cluster.cols, 'Phenograph_metacluster', 'mean') #cellular.cols
        make.pheatmap(exp, 'Phenograph_metacluster', cluster.cols,
                      #dendrograms.sort = TRUE, #dendrograms.sort = FALSE  is default
                      transpose=TRUE, # F by default MARKERS are columns, CLUSTERS are rows
                      file.name = paste0("Phenograph_unannotated.png")
        ) #cellular.cols
        rm(exp)
        
  }#end Phenograph block      

        
        
# ................... dim reduction ...................           
        
# ........ tSNE ........
if(run.tSNE == 1){  
        gated.cell.dat <- run.fitsne(gated.cell.dat, cluster.cols, perplexity = 200) # org perplex 200, default 30
        gated.cell.dat
}#end tSNE block
        
        
# ........ UMAP ........
if(run.UMAP == 1){  
        library(uwot) # umap
        umap.coords <- umap(as.matrix(gated.cell.dat[, cluster.cols, with = FALSE]), 
                            n_neighbors = 16, # 20. shinyIMC looked cool with 16 - 0.03. I had here 10 - 0.01 to really separate the clusters
                            min_dist = 0.08,  #0.1
                            metric = 'euclidean') #n_neighbors = 20, min_dist = 0.1, metric = 'euclidean'
        
        message("UMAP calculation done")
        colnames(umap.coords)[1:2] <- paste0('UMAP', 1:2)
        
        gated.cell.dat[["UMAP1"]] <- umap.coords[, 1]
        gated.cell.dat[["UMAP2"]] <- umap.coords[, 2]
        
        
        
        

        # you might wanna see the two WD against ND on the map, right?
        
        g <-   make.colour.plot( gated.cell.dat , 'UMAP1', 'UMAP2', "Diet", 'factor', 
                                 titlestring = paste0("UMAP by sample and week")  ,
                                 subtitlestring = "Entire pancreatic dataset" ,
                                 point.alpha = 0.8, # lets not overdraw the UMAP and reduce alpha of the geom_points
                                 dot.size = 0.9, # less alpha and smaller dots on the half-width-page graphs
                                 add.label = TRUE, #col.axis are the Phenograph/FlowSOM metacluster numbers. Print these pls
                                 legend.loc = "top",
                                 save.to.disk = F)
        filename <- "UMAP_entirePancDS_WDvsND.png"
        ggsave(filename,g, width = 30, height = 20, units = "cm",dpi = 1200)
        
        
        
        
        
 
        
        
        
        
  
        
        #     CREATE polychrome palette for the unannotated Phenograph MCs
        colRange <- unique(gated.cell.dat[["Phenograph_metacluster"]])
        colRange <- colRange[order(colRange)]
        #colRange <- as.character(colRange)
        set.seed(723451) # for reproducibility
        Phenograph_metacluster_palette <- createPalette(length(colRange), c("#ff0000"), M=100000, prefix = "")
        # so, the cool thing is that once we annotated these clusters, we can just do the same game again, overrride the palette and plot again.
        names(Phenograph_metacluster_palette) <- colRange
        
        
        
        
        g <-  make.colour.plot( gated.cell.dat , 'UMAP1', 'UMAP2', "Phenograph_metacluster", 'factor', 
                                        titlestring = paste0("UMAP by unannotated Phenograph MC")  ,
                                        subtitlestring = "Entire pancreatic dataset" ,
                                        point.alpha = 0.4, # lets not overdraw the UMAP and reduce alpha of the geom_points
                                        dot.size = 0.9, # less alpha and smaller dots on the half-width-page graphs
                                        colours = "polychrome",
                                        polychromepalette = Phenograph_metacluster_palette,
                                        add.label = TRUE, #col.axis are the Phenograph/FlowSOM metacluster numbers. Print these pls
                                        legend.loc = "top",
                                        save.to.disk = F)
        filename <- "UMAP_unannotatedPhenoclusters.png"
        ggsave(filename,g, width = 30, height = 20, units = "cm",dpi = 1200)
        
        
        
}#end UMAP block       
        
  

  
        
        
        g <-          make.colour.plot( gated.cell.dat , 'UMAP1', 'UMAP2', sample.col, 'factor', 
                          titlestring = paste0("UMAP by sample and week")  ,
                          subtitlestring = "Entire pancreatic dataset" ,
                          point.alpha = 0.4, # lets not overdraw the UMAP and reduce alpha of the geom_points
                          dot.size = 0.9, # less alpha and smaller dots on the half-width-page graphs
                          colours = "polychrome",
                          polychromepalette = Conditions_palette,
                          add.label = TRUE, #col.axis are the Phenograph/FlowSOM metacluster numbers. Print these pls
                          legend.loc = "top",
                          save.to.disk = F)
        filename <- "UMAP_bySampleandWeek.png"
        ggsave(filename,g, width = 30, height = 20, units = "cm",dpi = 1200)
        
        
        

        
        
        
        g <- make.colour.plot( gated.cell.dat , 'UMAP1', 'UMAP2', "Diet", 'factor', 
                          titlestring = paste0("UMAP by conditions")  ,
                          subtitlestring = "Entire pancreatic dataset" ,
                          point.alpha = 0.2, # lets not overdraw the UMAP and reduce alpha of the geom_points
                          dot.size = 0.8, # less alpha and smaller dots on the half-width-page graphs
                          add.label = TRUE, #col.axis are the Phenograph/FlowSOM metacluster numbers. Print these pls
                          legend.loc = "top",
                          save.to.disk = F)
        filename <- "UMAP_byCondition.png"
        ggsave(filename,g, width = 30, height = 20, units = "cm",dpi = 1200)
        
        
        
        
        
        #     CREATE polychrome palette for the used OT batches
        colRange <- unique(gated.cell.dat[["Batch"]])
        colRange <- colRange[order(colRange)]
        #colRange <- as.character(colRange)
        set.seed(723451) # for reproducibility
        OTbatch_palette <- createPalette(length(colRange), c("#009E73"), M=100000, prefix = "")
        
        swatch(OTbatch_palette) # wanna see it?
        # so, the cool thing is that once we annotated these clusters, we can just do the same game again, overrride the palette and plot again.
        names(OTbatch_palette) <- colRange   
        
        
        
        
      g <-   make.colour.plot( gated.cell.dat[gated.cell.dat[['Tissue']] != 'spleen',] , 'UMAP1', 'UMAP2', "Batch", 'factor', 
                          title = paste0("UMAP by sample and week")  ,
                          point.alpha = 0.4, # lets not overdraw the UMAP and reduce alpha of the geom_points
                          dot.size = 0.9, # less alpha and smaller dots on the half-width-page graphs
                          colours = "polychrome",
                          polychromepalette = OTbatch_palette,
                          add.label = TRUE, #col.axis are the Phenograph/FlowSOM metacluster numbers. Print these pls
                          legend.loc = "top",
                          save.to.disk = F)
      filename <- "UMAP_Pancreas_batches.png"
      ggsave(filename,g, width = 30, height = 20, units = "cm",dpi = 1200)
      
      
      

        

        

       # watch out. even though we could use here the  Conditions_palette for consistency, I did not adapt make.multi.plot function yet. so these colors are different and taken from spectral I guess.
        make.multi.plot(gated.cell.dat, 'UMAP1', 'UMAP2', sample.col, 'ROI', col.type = 'factor', figure.title = 'Annotated metacluster by ROI')
        filename <- "UMAPmultiplot_ROIsise_annotatedMetacluster.png"
        ggsave(filename,g, width = 30, height = 20, units = "cm",dpi = 1200)
            
        
        
} # end run.old.code protection bracket      
    
    
    
    
    
    
    
    
    
    
  if( how.often.ran.STAGE2==1) {
    
    cat("\n\n\n----------------------------------- \n End of STAGE 2 FIRST PASS. \n If you need to clean out clusters, define them and re-run STAGE 2 \n If not proceed to STAGE 3 \n-----------------------------------\n\n") 
    
  }else if( how.often.ran.STAGE2==2) {
    
    cat("\n\n\n----------------------------------- \n End of STAGE 2 SECOND PASS. \n Defined clusters were deleted and re-calculated. \n Cluster-deleting is blocked from now on \n-----------------------------------\n\n") 
    
  }else{
    
    cat("\n\n\n----------------------------------- \n End of STAGE 2 THIRD+ PASS. \n You were re-running the same code \n Cluster-deleting is blocked \n-----------------------------------\n\n") 
  }
    
    
    
    
    
    
    
    
} # end stage 2: cluster and dim-reduction




########################## STAGE 3 ##########################
if(whichSTAGE == 3){

  
# there is a cryptic "out of bound" error that is hard to pin down if you dont know what the script does here, and its happening
# when collapse.levels is not defined for all GATEsubsets. in that case the for s loop along amountofGATEsubsets runs out of numbers.
# its a silly mistake when you adjust the gating strategy and end up with one more subset, then you need to adjust the collapse.levels vector with an additional number:
if(   ClusterCollapsing==1 &  !(length(collapse.levels)==amountofGATEsubsets)      ){
  cat("\n\n\n----------------------------------- \n ERROR in STAGE 3 initialization.   \n Collapse levels are not defined for all GATEsubsets \n  1) Check your collapse.levels vector and make sure it matches your subsets \n -> Re-run STAGE 3 to pool the subsets back together. \n-----------------------------------\n\n")# stop(), call. = FALSE)  
  stop()   
} # end STAGE3 if you did not define collapse levels for all subgates 

  
  setwd(OutputDirectory)
  

######################  Collapse clusters ###################### 


homogeneity.of.subsets <- data.frame( Collapse.level = integer(),
                                        Subset = integer(),
                                        RMSD = integer()    )

  for(s in 1:amountofGATEsubsets  ){

    # load the cluster columns from our list for the according subset, and set the marker.order for the subset:
    cluster.cols <- cluster.cols.GATEsubsets[[s]] # use double brackets to subset into the vector and not just the list.
    marker.order <- cluster.cols.GATEsubsets[[s]] # depreciated: global.marker.order[global.marker.order %in% cluster.cols]

if (collapse.levels[[s]] == "max"){
  
  maximal.amount.of.clusters <-  max(clustered.gated.subsets[[s]]$Phenograph_metacluster) 
  # you use 0 index, and hclust wants more or equal than 2 objects to cluster: 
  #thats why the max amount of collapses is 2 less than clusters we got:
  max.collapse.level <- maximal.amount.of.clusters-2
}else{
  max.collapse.level <- collapse.levels[[s]]
}

   

# for every of the n collapse levels, we define the collapsed clusters a new, so every next level will have more and more clusters to collapse

for(l in 0:max.collapse.level  ){
  
  
  

    
    
    # for sequential cluster number collapse, we need to do that only once in the beginning of the loop. After that, we recycle these numbers
    if(ClusterCollapsing==1 && l==0){
      
      # we move the original Metacluster numbers in a separate column
      clustered.gated.subsets[[s]]$Phenograph_metacluster -> clustered.gated.subsets[[s]]$Collapsed_metacluster
      #and use three-digit metacluster numbers so that we can merge the subsets later on:
      clustered.gated.subsets[[s]]$Collapsed_metacluster <-  100*s+clustered.gated.subsets[[s]][,Collapsed_metacluster]
    }
    


      
  
  
# we need to protect hclust from collapsing less than two clusters
# this holds especially for the subset(s) that have less clusters than the others: these need skipping while we collapse to "max"
    
if( length(unique(clustered.gated.subsets[[s]]$Collapsed_metacluster)) >= 2 ){
  
  
  
  # LEVEL 0 is uncollapsed, so we need another bracket here to skip l=0 from the collapse routine and just plot:  
  
  if( l > 0 ){
  
     message(paste0("\nCollapsing subset ", GATEsubsetNames[s], " to level ", l    ))
  
  
  #---- find the two closest clusters in clustered.gated.subsets[[s]] ---  
# the idea is from there https://slowkow.com/notes/pheatmap-tutorial/    
    
    exp <- do.aggregate(dat=clustered.gated.subsets[[s]], 
                        use.cols=cluster.cols, 
                        by='Collapsed_metacluster', 
                        func = 'mean') #
    

    # store away the metacluster numbers
    metacluster.numbers <- exp$Collapsed_metacluster
    # cut the metacluster column away and make matrix
    # WATCH OUT! For the plotting of these signals, we would make another z-score now
    # we do that for maximizing contrast on the heatmap. However, here we find the two closest clusters,
    # so we do not massage the data for that, but simply pass the original distribution on:
    exp <- as.matrix(exp[,-1])
    
    # give the rows their numbers back
    rownames(exp) <-metacluster.numbers
    # and since we wont plot, we dont chop off the _t_z underscores in the column names either...
    
    # calculate dendrogram
    # rather than hclust( dist( exp ) ) we now use the fast function fastcluster:
    mat_cluster_cols <- fh(exp)
    # use our fancy function to sort them along closeness
    mat_cluster_cols <- sort_hclust(mat_cluster_cols)
    
    # like this the two most close clusters are now left and are pulled as the first two numbers from this vector.
    # Watch out, you can pull the order, which is the nth elemnt, and you use that nth element to get the label out (which is our rownames)
    clusters.to.collapse <-     mat_cluster_cols$labels[     sort( mat_cluster_cols$order[1:2] ) ]
    
    # overwrite the bigger number with the smaller number into which we collapse
    #message(paste0("Collapsing the closest clusters ", clusters.to.collapse[2], " into ", clusters.to.collapse[1] ))
    clustered.gated.subsets[[s]]$Collapsed_metacluster[which( clustered.gated.subsets[[s]]$Collapsed_metacluster == clusters.to.collapse[2] )] <- clusters.to.collapse[1]
    
    
    # we now have to mutate the numbers so that we get a nice sequence again
    # lev does it differently: he moves the first number, and then every number above that gets moved one down. This is certainly less messy in the code, but you cannot
    # create these lists I did above to see which clusters are moved. I move the numbers all at once, he moves and re-numbers them sequenctially.
    
    
    # lets predefine which numbers are missing
    temp.cluster.numbers <- unique(clustered.gated.subsets[[s]]$Collapsed_metacluster)
    
    missing.numbers <- setdiff( min(temp.cluster.numbers):max(temp.cluster.numbers)  , temp.cluster.numbers) 
    
    while( length(missing.numbers) > 0 ){
      
      # now we need to make sure that the number one above missing.numbers[1] is actually part of the metaclusters:
      
      # we gonna find the cluster number that is adjacent to the missing cluster number. we use which to return the element and pull the cluster out.
      # this is the cluster we will re-assign to missing.numbers[1]:
      temp.cluster.number.to.move <- temp.cluster.numbers[    which(temp.cluster.numbers > missing.numbers[1] )[1] ]
      
      clustered.gated.subsets[[s]]$Collapsed_metacluster[which( clustered.gated.subsets[[s]]$Collapsed_metacluster == temp.cluster.number.to.move )] <- missing.numbers[1]
      
     # message(paste0("Re-assigned cluster ", temp.cluster.number.to.move, " to next smaller unused number ", missing.numbers[1]))
      
      # now we need to update the two vectors:
      temp.cluster.numbers <- unique(clustered.gated.subsets[[s]]$Collapsed_metacluster)
      missing.numbers <- setdiff( min(temp.cluster.numbers):max(temp.cluster.numbers)   , temp.cluster.numbers) 
      
      #message(paste0("Following cluster numbers are still missing:" ))
      #cat(missing.numbers)
      
      
    }# end while bracket
    
    
    
  }#only enter collapse routine if l passed level 0 
  # level 0 skipped the upper bracket but we need MAE/RMSD from the uncollapsed nevertheless, so here we go:
    
    
    
    
    # now we calculate homogeneity of the subset at this collapse level
    #it important that the cluster numbers come in last, since the calc.homogeneity function wants the cluster numbers in the last column. (and removes it for rmsd calculation)
    temp.cluster.cols.CollapsedClusters <- cbind( clustered.gated.subsets[[s]][, cluster.cols, with = FALSE] , 
                                                  clustered.gated.subsets[[s]][,Collapsed_metacluster]  )  
    
    #if you wanna see how the channel signals are present in that subset, here is a way to scan the range of every marker:
   # t(sapply(temp.cluster.cols.CollapsedClusters, range))
    
    
    
    
    # the idea to measure spread is RMSD, we measure accuracy of the cluster
    # with this function, we can return a single number that is calculated by the whole dataframe: it takes RMSD or MAE given what you put in the header
    temp.RMSD <- calc.homogeneity(m=temp.cluster.cols.CollapsedClusters  )
    homogeneity.of.subsets %>% add_row(
      Collapse.level = l,
      Subset = s,
      RMSD = as.numeric(temp.RMSD)
    ) -> homogeneity.of.subsets
    
    rm(temp.cluster.cols.CollapsedClusters)
    
    if (clusterhomogeneitymeasure == "MAE") {   message(paste0("-> Subset MAE: ", round( temp.RMSD,3)  ))  }
    if (clusterhomogeneitymeasure == "RMSD") {   message(paste0("-> Subset RMSD: ", round( temp.RMSD,3)  ))  }
    
    
    # I wanna see all these plots only if I go for max collapse, or if I am at the end of collapsing with l:
      if (collapse.levels[[s]] == "max" | (max.collapse.level==l & max.collapse.level>0 ) ){
      
        
    # watch out, if you let the code collapse the last two clusters, hclust that is used to plot the heatmap in the following function will cry.
    # so you gotta protect that here as well.
    
    if( length(unique(clustered.gated.subsets[[s]]$Collapsed_metacluster)) > 1 ){
      
      
#.................. Prepare colors for our clusters  ...............      
## we re-create a temporary palette for every subgate:
# watch out, this here CAN NOT be the final collapsed palette, since we did not merge the subsets into one!! 
# This palette will be created in a the end. so you need to re-run the heatmap again by hand to get the same colors like the final one!
      
      #     CREATE polychrome palette for the unannotated Phenograph MCs
      colRange <- unique(clustered.gated.subsets[[s]][["Collapsed_metacluster"]])
      colRange <- colRange[order(colRange)]
      #colRange <- as.character(colRange)
      set.seed(723451) # for reproducibility
      Phenograph_metacluster_palette <- createPalette(length(colRange), c("#ff0000"), M=100000, prefix = "")
      # so, the cool thing is that once we annotated these clusters, we can just do the same game again, overrride the palette and plot again.
      names(Phenograph_metacluster_palette) <- colRange
      
      
      
      
#.................. MAKE COMPLEX HEATMAP  ...............    

    exp <- clustered.gated.subsets[[s]][, lapply(.SD, "mean", na.rm = TRUE), by = "Collapsed_metacluster", 
                          .SDcols = cluster.cols]
    
    ## new way to plot complex heatmap
    
    #z-normalize feature expression (iterate over columns with MARGIN=2)
    zscored.exp <- cbind(  exp[,1]  ,   apply(exp[,-1], scale, MARGIN = 2)  ) 
    
    # we gonna bring back only the signals into the matrix:
    zscored.exp.mat <- as.matrix(zscored.exp[,-1])
   
    # prepare column and row names:
    rownames(zscored.exp.mat) <- paste0(exp[[1]] )
    colnames(zscored.exp.mat) <- sub("\\.\\..*", "", colnames(zscored.exp.mat))
    
    col_fun = colorRamp2(c(range(zscored.exp.mat)[1], 0, range(zscored.exp.mat)[2]), c("blue", "white", "red"))
    
    # so the trick is to initiate column_ha first with the alphabetical order:
    column_ha = HeatmapAnnotation(cluster = rownames(zscored.exp.mat),
                                  
                                  
                                  # this is the Cluster legend. the problem is the "at", which takes over the alphabetical order of the subset.Phenograph_mc_palette
                                  annotation_legend_param = list(
                                    cluster = list(
                                      ncol = 1, 
                                      title = "Cluster",
                                      #title_position = "topcenter",
                                      at = rownames(zscored.exp.mat) ,#rownames(zscored.exp.mat), # names(subset.Phenograph_metacluster_palette),
                                      grid_height = unit(0.02*length(marker.order), "cm"),
                                      grid_width = unit(0.04*length(marker.order), 'cm'),
                                      labels_gp = gpar(fontsize = 0.8*length(marker.order)),
                                      title_gp = gpar(fontsize = 0.8*length(marker.order))
                                    )
                                  ),
                                  
                                  #col = list(cluster=hm_cols),
                                  col = list(cluster=Phenograph_metacluster_palette),
                                  
                                  na_col = "black", #white
                                  show_annotation_name = T
    )
    
    # then initiate the map
    hm <-   Heatmap( t( zscored.exp.mat ) ,
                     col = col_fun,
                     row_order = sub("\\.\\..*", "", marker.order) , # this is the only that needs stripping of transformed_rescaled
                     cluster_columns = fh,
                     top_annotation = column_ha,
                     
                     show_column_names = F, # turn off the names of the clusters under the hm
                     
                     #column_names_rot = 45,
                     rect_gp = gpar(col = "white", lwd = 1),
                     
                     heatmap_legend_param = list(
                       title = "z-score",
                       # direction = 'horizontal',
                       #title_position = "topcenter", # this is the title z-score
                       at=c(-2,-1,0,1,2,3,4),
                       title_position = "lefttop-rot",
                       legend_width = unit(0.25*length(marker.order), "cm"),
                       legend_height = unit(0.3*length(marker.order), "cm"),
                       grid_width = unit(0.02*length(marker.order), 'cm'),
                       labels_gp = gpar(fontsize = 0.6*length(marker.order)),
                       title_gp = gpar(fontsize = 0.8*length(marker.order))
                     )
                     
    )
    
    
    #then bring back column_ha
    
    column_ha = HeatmapAnnotation(cluster = rownames(zscored.exp.mat),
                                  
                                  
                                  # this is the Cluster legend. the problem is the "at", which takes over the alphabetical order of the subset.Phenograph_mc_palette
                                  annotation_legend_param = list(
                                    cluster = list(
                                      ncol = 1, 
                                      title = "Cluster",
                                      #title_position = "topcenter",
                                      at = rownames(zscored.exp.mat)[column_order(hm)] ,#rownames(zscored.exp.mat), # names(subset.Phenograph_metacluster_palette),
                                      grid_height = unit(0.02*length(marker.order), "cm"),
                                      grid_width = unit(0.04*length(marker.order), 'cm'),
                                      labels_gp = gpar(fontsize = 0.8*length(marker.order)),
                                      title_gp = gpar(fontsize = 0.8*length(marker.order))
                                    )
                                  ),
                                  
                                  #col = list(cluster=hm_cols),
                                  col = list(cluster=Phenograph_metacluster_palette),
                                  
                                  na_col = "black", #white
                                  show_annotation_name = T
    )
    
    
    
    
    # then initiate the map again with the new column_ha order:
    hm <-   Heatmap( t( zscored.exp.mat ) ,
                     col = col_fun,
                     row_order = sub("\\.\\..*", "", marker.order), # this is the only that needs stripping of the qupath string
                     cluster_columns = fh,
                     top_annotation = column_ha,
                     
                     show_column_names = F, # turn off the names of the clusters under the hm
                     
                     #column_names_rot = 45,
                     rect_gp = gpar(col = "white", lwd = 1),
                     
                     heatmap_legend_param = list(
                       title = "z-score",
                       # direction = 'horizontal',
                       #title_position = "topcenter", # this is the title z-score
                       at=c(-2,-1,0,1,2,3,4),
                       title_position = "lefttop-rot",
                       legend_width = unit(0.25*length(marker.order), "cm"),
                       legend_height = unit(0.3*length(marker.order), "cm"),
                       grid_width = unit(0.02*length(marker.order), 'cm'),
                       labels_gp = gpar(fontsize = 0.6*length(marker.order)),
                       title_gp = gpar(fontsize = 0.8*length(marker.order))
                     )
                     
    ) 
    
    
    # I forgot why we remove the zscored.exp here:
    rm(zscored.exp)   
    
    
# save every heatmap away:
#    png(filename = paste0(GATEsubsetNames[s],"_Heatmap_collapsed_",l,".png"), 
#        width = 35*(ncol(zscored.exp.mat[,-1])), 
#        height = 29*(length(marker.order)) ) 
#    
#    draw(hm,
#         heatmap_legend_side = 'top',
#         row_sub_title_side = 'left',
#         padding = unit(c(2, 2, 2, 10), "mm"))
#    
#    invisible(dev.off())
#    
#
#    rm(exp)
#    rm(zscored.exp)
#    rm(zscored.exp.mat)

    
    hm.grob = grid.grabExpr(
          draw(hm,
               heatmap_legend_side = 'top',
               row_sub_title_side = 'left',
               padding = unit(c(2, 2, 2, 10), "mm"))
    )    
    
        
    

#.................. MAKE UMAP  ...............    
    
    #we need different title strings given RMSD or MAE measure:
    paste0("UMAP by Phenograph clusters. Collapsed to level ", l,", RMSD: ", round(temp.RMSD,3)) 
    if (clusterhomogeneitymeasure == "MAE") {  subtitlestring <- paste0("MAE: ", round(temp.RMSD,3), " (",  nrow(clustered.gated.subsets[[s]]), " cells plotted)")   }
    if (clusterhomogeneitymeasure == "RMSD") {  subtitlestring <- paste0("RMSD: ", round(temp.RMSD,3), " (",  nrow(clustered.gated.subsets[[s]]), " cells plotted)")   }
    

    g <-  make.colour.plot.adapted( dat=clustered.gated.subsets[[s]] , x.axis = 'UMAP1', y.axis ='UMAP2', 
                                    col.axis ="Collapsed_metacluster", 
                                    col.type = 'factor', 
                                    add.label = TRUE, #coloring comes via Annotated_metaclusters. Print these pls
                                    titlestring =  paste0("UMAP collapsed to level ", l)  ,
                                    subtitlestring = subtitlestring,
                                    point.alpha = 0.8, # lets not overdraw the UMAP and reduce alpha of the geom_points
                                    dot.size = 1.6, # less alpha and smaller dots on the half-width-page graphs
                                    colours = "polychrome",
                                    polychromepalette = Phenograph_metacluster_palette,
                                    repel.label.size = 3,
                                    label.label.repelforce = 30, # repulsion between overlapping text labels. default 30 to push them apart
                                    label.datapoint.pullforce = 0.8, # attraction between label and datapoint. default 0.8 to help release the label
                                    min.line.length.tolabel = 0, # minimal distance under which the segment line is not drawn anymore: 0 draw all, Inf turn off even if far
                                    legend.loc = "none",
                                    legend.text.size = 15, # default 18
                                    save.to.disk = F)
    
  #  filename <- paste0(GATEsubsetNames[s],"_UMAP_unannotatedPhenoclusters_collapse_level",l,".png")    
  #  ggsave(filename,g, width = 30, height = 20, units = "cm",dpi = 1200)
    
    
  
    hm.umap.fig <-  ggarrange(hm.grob, g, widths = c(6,5) )  
    
  
  if(l==0){
    
    collapse.lyrics = "Uncollapsed dataset"
    
  }else{
    
    collapse.lyrics = paste0("Collapsed ", clusters.to.collapse[2], " into ", clusters.to.collapse[1], " (everything above ",clusters.to.collapse[2]," is newly assigned)" )
  }
  
  
  hm.umap.fig <-   annotate_figure(hm.umap.fig,
                  top = text_grob(paste0("Subset: ", GATEsubsetNames[s]), color = "black", face = "bold", size = 17),
                  bottom = text_grob(data.source.lyrics, color = "red",
                                     hjust = 1, x = 1, face = "italic", size = 13),
                  left = text_grob(collapse.lyrics, color = "black", rot = 90),
                  #right = text_grob(bquote("Superscript: ("*kg~NH[3]~ha^-1~yr^-1*")"), rot = 90),
                  fig.lab = "Automatic collapse", fig.lab.face = "bold"
  ) 
   
  
 
  
       png(paste0(GATEsubsetNames[s],"_collapsed-",l,".png"), 
        width = 60*(ncol(zscored.exp.mat[,-1])), 
        height = 29*(length(marker.order)) )
  

       ggsave(paste0(gsub("[:/\\\\<>|?*]", "_", GATEsubsetNames[s]),"_collapsed-",l,".pdf"), hm.umap.fig, width = 15, height = 10, dpi=400) # use the ratio to make the HM squared
       
       
       height = (length(marker.order))
       width = 3*height # its not 3/2 bang on, since there is a legend in between...
       
       if(min(width, height)<4){
         width = 10*width
         height = width
       }
       
       # now, dpi must not exceeed 50000px. adjust dpi so that we end up with a 15000px image on the long edge:
       dpi <- 15000/max(width, height)
       
       ggsave(  plot = hm.umap.fig, 
                filename = paste0(gsub("[:/\\\\<>|?*]", "_", GATEsubsetNames[s]),"_collapsed-",l,".pdf"), 
                width = width, 
                height = height ,
                limitsize = FALSE,
                dpi = dpi)



    
    }# protect make.pheatmap hclust to plot one cluster
    
    
    
    
    
    
    
    
      }#only plot if you either reach the end of collapsing or its set to max
    
   
    
}# end protect hclust from collapsing one single cluster.
    
    
    
}# end processing all GATEsubsets at different collapsing levels running along l    
}# end run through the GATEsubsets with s and allocate
  
  
  

  




#invisible(dev.off())
g <- ggplot(homogeneity.of.subsets, aes(x=Collapse.level , y=RMSD, group=factor(Subset) )  ) +
  geom_line(aes(color=factor(Subset)))+
  geom_point(aes(color=factor(Subset)), size=2.5)+
  theme_classic()+
  labs(title=paste0("Homogeneity of subsets CD19+ and CD19-CD45+" )  , 
       subtitle=paste0(clusterhomogeneitymeasure, " summed up over all markers in all clusters"),
       #caption="Created by M.Barone", 
       y=paste0("Total normalized ",clusterhomogeneitymeasure," [AU]"), 
       x="Cluster Collapse level"
  )

filename <- paste0( "Homogeneity_" ,clusterhomogeneitymeasure,"_maxCollapseLevel",l,".png")    
ggsave(filename,g, width = 20, height = 20, units = "cm",dpi = 1200)




  

        
        
        
# now for the moment we dont need more in that script, so its important to store away our clustered.gated.subsets as normal csv file.
# For this, we gonna have to drop the trash bin subset if this was not used:



######################  RE-UNITE gated subsets to one dataframe: gated.cell.dat. ###################### 
gated.cell.dat <- do.call(rbind, clustered.gated.subsets[c(1:amountofGATEsubsets)])
        




###################### ---- Create the polychrome color palettes now for good --- ###################### 
# wanna see a palette: swatch(xxx_palette)  

#     CREATE polychrome palette for the unannotated Phenograph MCs
colRange <- unique(gated.cell.dat[["Collapsed_metacluster"]])
colRange <- colRange[order(colRange)]
#colRange <- as.character(colRange)
set.seed(723451) # for reproducibility
Phenograph_metacluster_palette <- createPalette(length(colRange), c("#ff0000"), M=100000, prefix = "")
names(Phenograph_metacluster_palette) <- colRange



#     CREATE polychrome palette for the used TMA batches
colRange <- unique(gated.cell.dat[["Image"]])
colRange <- colRange[order(colRange)]
#colRange <- as.character(colRange)
set.seed(723451) # for reproducibility
TMA_palette <- createPalette(length(colRange), c("#009E73"), M=100000, prefix = "")
names(TMA_palette) <- colRange   

all_images <- sort(unique(gated.cell.dat$Image)) 
amount.of.TMAs <- length(all_images)



#     CREATE polychrome palette for the Patient.ID
colRange <- unique(gated.cell.dat[["Patient.ID"]])  
colRange <- colRange[order(colRange)]
set.seed(723451) # for reproducibility
Patient.ID_palette <- createPalette(length(colRange), c("#32c5a1"), M=100000, prefix = "")
names(Patient.ID_palette) <- colRange
#swatch(Patient.ID_palette)


#     CREATE polychrome palette for the cell classes
colRange <- unique(gated.cell.dat[["Classification"]]) 
colRange <- colRange[order(colRange)]
set.seed(723451) # for reproducibility
class_palette <- createPalette(length(colRange), c("#e77b21"), M=100000, prefix = "")
names(class_palette) <- colRange
#swatch(class_palette)






# the next two are not needed for the IM project:
#     CREATE a palette for diet
tissue.type_palette <- c("Healthy" = "#017DD5", 
                  "Tumor" = "#EB8B19")




# update our plot.ROIs object since it still contains spleen and stuff:
plot.ROIs  <- sort(unique(gated.cell.dat$ROI) )
amount.of.ROIs <- length(plot.ROIs)






all.MC <- unique(gated.cell.dat$Collapsed_metacluster)
amount.of.MC <- length(all.MC)





if(plot.global.heatmap.all.subsets==1){

### Expression heatmap for Phenograph
# watch out, this cluster.cols is the selection for the very last subset that ran before. Maybe it makes sense to use here marker order itself to use all 

exp <- gated.cell.dat[, lapply(.SD, "mean", na.rm = TRUE), by = "Collapsed_metacluster", 
    .SDcols = global.marker.order]


#z-normalize feature expression (iterate over columns with MARGIN=2)
zscored.exp <- cbind(  exp[,1]  ,   apply(exp[,-1], scale, MARGIN = 2)  ) 

# we gonna bring back only the signals into the matrix:
zscored.exp.mat <- as.matrix(zscored.exp[,-1])

# prepare column and row names:
rownames(zscored.exp.mat) <- paste0(exp[[1]] )
colnames(zscored.exp.mat) <- sub("\\.\\..*", "", colnames(zscored.exp.mat))

col_fun = colorRamp2(c(range(zscored.exp.mat)[1], 0, range(zscored.exp.mat)[2]), c("blue", "white", "red"))


column_ha = HeatmapAnnotation(cluster = rownames(zscored.exp.mat),
                              
                              annotation_legend_param = list(
                                cluster = list(
                                  ncol = 2, 
                                  title = "Cluster",
                                  title_position = "topcenter",
                                  at = names(Phenograph_metacluster_palette),
                                  grid_height = unit(0.02*length(global.marker.order), "cm"),
                                  grid_width = unit(0.02*length(global.marker.order), 'cm'),
                                  labels_gp = gpar(fontsize = 0.8*length(global.marker.order)),
                                  title_gp = gpar(fontsize = 0.8*length(global.marker.order))
                                )
                              ),
                              
                              #col = list(cluster=hm_cols),
                              col = list(cluster=Phenograph_metacluster_palette),
                              
                              na_col = "black", #white
                              show_annotation_name = FALSE
)


hm <-   Heatmap( t( zscored.exp.mat ) ,
                 col = col_fun,
                 row_order = sub("_.*", "", global.marker.order ) , # this is the only that needs stripping of transformed_rescaled
                 cluster_columns = fh,
                 top_annotation = column_ha,
                 rect_gp = gpar(col = "white", lwd = 1),
                 heatmap_legend_param = list(
                   title = "z-score",
                   direction = 'horizontal',
                   title_position = "topcenter",
                   legend_width = unit(0.25*length(global.marker.order), "cm"),
                   grid_width = unit(0.02*length(global.marker.order), 'cm'),
                   labels_gp = gpar(fontsize = 0.8*length(global.marker.order)),
                   title_gp = gpar(fontsize = 0.8*length(global.marker.order))
                 )
                 
              )


setwd(OutputDirectory)
invisible(dev.off())
png(filename = paste0('Heatmap_CollapsedPooled.png'), width = 35*(ncol(zscored.exp.mat[,-1])), height = 29*(length(global.marker.order)) ) #, width = 65*(ncol(zscored.exp.mat[,-1])), height = 58*(length(marker.order)) )

draw(hm,
     heatmap_legend_side = 'top',
     row_sub_title_side = 'left',
     padding = unit(c(2, 2, 2, 10), "mm"))

invisible(dev.off())


rm(exp)
rm(zscored.exp)
}# end only plot the heatmap over the entire dataset.







######################  Annotate clusters ###################### 

#maybe its a good help to know how many MC we are dealing with. We use all.MC to move the un-annotated cluster names back into the Annotation



#.................. do it by hand  ..............       

if(use.my.annotation.list==TRUE){ 

  
  
# now, before we do anything here, lets just check if we already got an annotation in gated.cell.dat.
# This most probably happens because you re-loaded data from a previous run, or just re-annotate the clusters in the same session.
# we gonna play safe and just remove it for good:
  
  if("Annotated_metacluster" %in% colnames(gated.cell.dat) ) {
    
    message(  "\n-------   Manual Annotation   --------\nGonna remove Annotated_metacluster from isolated gated.cell.dat before running"  )
    gated.cell.dat <-  subset(gated.cell.dat, select = -c(Annotated_metacluster) )
  }else{
    message(  "\n-------   Manual Annotation   --------\nNo manual annotation present in isolated gated.cell.dat. Running..."  )
  }
  
  
  # work with cluster.annots that is provided in the STAGE 3 SETTINGS:     
  cluster.annots <- do.list.switch(cluster.annots)
  names(cluster.annots) <- c('Values', 'Annotated_metacluster')
  cluster.annots
  
  
  
  # Rather than injecting the annotation of cluster.annots right into the gated.cell.dat dataframe, lets do it differently:
  
  # we first need to adjust all.MC to incorporate also the annotation
  # the easy way is of course:
  all.annot.MC <- do.add.cols(as.data.frame(all.MC), 'all.MC', cluster.annots, 'Values')
  ### Whatever Metacluster was not annotated got now NA
  # in our case we dont want "other" for these, but we need the original numbers back in:
  all.annot.MC <- all.annot.MC %>% 
    mutate(Annotated_metacluster = coalesce(Annotated_metacluster,all.MC))
  
  
  # now we tweak the list to be accepted by do.add.cols, but this time we inject into gated.cell.dat:
  names(all.annot.MC) <- c('Values', 'Annotated_metacluster')
  gated.cell.dat <- do.add.cols(gated.cell.dat, 'Collapsed_metacluster', all.annot.MC, 'Values')
  
  
  # thats it.
  # override the last status of gated.cell.dat (STAGE 3 at least) in your data folder:
  setwd(OutputDataDirectory)
  fwrite( gated.cell.dat , 'gated.cell.dat.csv')
  
  
  # last thing to do is to create a color palette for the annotated clusters:
  # watch out, the make.color.plot engine does not accept a full-range Phenograph_metacluster_palette if you just plot a subset of it.
  # you most definetly want to also plot the subsets, and for that we need the three digit number of the clusters, which we would overwrite
  # if we just matched the names of Phenograph_metacluster_palette with the all.annot.MC$Values entries.
  # so we need a more complex thing than just a character string:
  
  
  annotated.Phenograph_metacluster_palette <- data.frame("Collapsed_metacluster"=names(Phenograph_metacluster_palette), "color"=Phenograph_metacluster_palette, row.names=NULL)
  annotated.Phenograph_metacluster_palette <- do.add.cols(annotated.Phenograph_metacluster_palette, 'Collapsed_metacluster', all.annot.MC, 'Values')
  # you just gotta love do.add.cols, no?
  
  annotated.Phenograph_metacluster_palette$Collapsed_metacluster <- as.numeric(as.character(annotated.Phenograph_metacluster_palette$Collapsed_metacluster))
  
# only now I would tweak the colors:
  # Birigt doesnt like two closely colored clusters:
  # unknown and CD127+ pDCs
  # Ki67hi ILC-like and neutrophils

    annotated.Phenograph_metacluster_palette <-  annotated.Phenograph_metacluster_palette %>%
    # other is not interesting and gets a boring af gray:
    mutate(color=replace(color, str_detect(Annotated_metacluster, "other"), "#AFAFAF")) %>% 
    # the others Birgit didnt like:
    mutate(color=replace(color, str_detect(Annotated_metacluster, "neutrophils"), "#B71C00")) %>%
    mutate(color=replace(color, str_detect(Annotated_metacluster, "Tregs"), "#0072B2")) %>%  # color-blind friendly
    # somesh has the macs in green:
    mutate(color=replace(color, str_detect(Annotated_metacluster, "F4.80lo macrophages"), "#E69F00")) %>% # color-blind friendly
    mutate(color=replace(color, str_detect(Annotated_metacluster, "F4.80- macrophages"), "#09A3CC")) %>% # color-blind friendly
    mutate(color=replace(color, str_detect(Annotated_metacluster, "F4.80hi macrophages"), "#D55E00")) %>% # color-blind friendly
    # and use a strong color for our activated effectors:
    mutate(color=replace(color, str_detect(Annotated_metacluster, "activ. effector-like CD8+"), "#FF2D00")) %>% # color-blind friendly adapt from CC79A7
    
    # since we moved them all around, we also need to adapt some similarly colored clusters in the UMAP
    mutate(color=replace(color, str_detect(Annotated_metacluster, "MHCII lo macrophages 1"), "#009E73")) %>% # they get a similar blue 1CCFFF like macs 2
    mutate(color=replace(color, str_detect(Annotated_metacluster, "MHCII lo macrophages 2"), "#2CC79D")) %>% 
    
    mutate(color=replace(color, str_detect(Annotated_metacluster, "CD4- pDCs"), "#503AEE")) %>%  # they are too close to our F480 neg macs. they get a similar A8B2FD like the CD4pos pDCS
    
    mutate(color=replace(color, str_detect(Annotated_metacluster, "NK cells"), "#BE1F96")) %>% # this also changed the string of the CD16lo!!!
    mutate(color=replace(color, str_detect(Annotated_metacluster, "CD16lo NK cells"), "#B672A5"))  # move the NKs closer to each other for good
  
    
    
    # end protect from running manual annotation  
}else if(use.my.annotation.list==FALSE){
  # in that case, lets just copy-paste Collapsed into the annotated column  
  message(  "\n-------   No Annotation   --------\nCopying the collapsed metaclusters to annotation column..."  ) 
  
  gated.cell.dat$Annotated_metacluster <- gated.cell.dat$Collapsed_metacluster
  
 # we still need the palette, but this one is just the copy-paste of Collapsed_metacluster:
  annotated.Phenograph_metacluster_palette <- data.frame("Collapsed_metacluster"=names(Phenograph_metacluster_palette), "color"=Phenograph_metacluster_palette, row.names=NULL)
  annotated.Phenograph_metacluster_palette$Collapsed_metacluster <- as.numeric(as.character(annotated.Phenograph_metacluster_palette$Collapsed_metacluster))
  annotated.Phenograph_metacluster_palette$Annotated_metacluster <- annotated.Phenograph_metacluster_palette$Collapsed_metacluster
}  

 
  # you might wanna plot the UMAP now anew using the annotated clusters
  # for the heatmap that shall use Annotated_metacluster now, we need to supply the channels used for clustering
  # this object was already loaded: cluster.cols.GATEsubsets
  
  setwd(OutputDirectory)
  
  
  for(s in 1:amountofGATEsubsets  ){
  #for(s in 2:2  ){
    
message("\nRe-plotting Heatmap/UMAP of ", GATEsubsetNames[s], ", this time annotated" )  
    
    
    # load the cluster columns from our list for the according subset:
    cluster.cols <- cluster.cols.GATEsubsets[[s]] # use double brackets to subset into the vector and not just the list.  
    marker.order <- cluster.cols.GATEsubsets[[s]] # depreciated: global.marker.order[global.marker.order %in% cluster.cols]
    
    # watch out, the make.color.plot engine does not accept a full-range Phenograph_metacluster_palette if you just plot a subset of it.
    # And if you plot the heatmap with the full range palette, they get all listed in the legend. So lets create a temporary Palette subset:
    #subset.Phenograph_metacluster_palette <- Phenograph_metacluster_palette[ which( 100*(s) < names(annotated.Phenograph_metacluster_palette) & names(annotated.Phenograph_metacluster_palette) < 100*(s+1)  ) ]
    #names(subset.Phenograph_metacluster_palette) <- cluster.annots$Annotated_metacluster[match( names(subset.Phenograph_metacluster_palette), cluster.annots$Values )]
    #subset.Phenograph_metacluster_palette <- subset.Phenograph_metacluster_palette[!duplicated(names(subset.Phenograph_metacluster_palette))]
    
    

 # we gotta create a character vector again:
    #new way:
    
    subset.Phenograph_metacluster_palette <-   subset(annotated.Phenograph_metacluster_palette, Annotated_metacluster %in% unique(subset( gated.cell.dat, GATEsubset==s)$Annotated_metacluster) ) %>%
      pull( color,Annotated_metacluster  ) # pulling color out and naming it Annotated_metacluster
    
    
    #old way depreciated.
  #  subset.Phenograph_metacluster_palette <- subset(annotated.Phenograph_metacluster_palette, Collapsed_metacluster > 100*(s) &  Collapsed_metacluster < 100*(s+1)   ) 
    
   #  #temp.pal <- subset.Phenograph_metacluster_palette[,color]# for some reason, this ends up empty when not annotating with a list. its not the missing tibble, I checked that
  #  temp.pal <- subset.Phenograph_metacluster_palette$color
  #  names(temp.pal ) <-subset.Phenograph_metacluster_palette$Annotated_metacluster
  #  subset.Phenograph_metacluster_palette <- temp.pal
  #  rm(temp.pal)
    
  #  subset.Phenograph_metacluster_palette <- subset.Phenograph_metacluster_palette[ unique(names(subset.Phenograph_metacluster_palette)) ] 
    # watch out, this here: subset.Phenograph_metacluster_palette[!duplicated(subset.Phenograph_metacluster_palette)] kicks off "NK cells" since the string is present in "CD16lo NK cells"
    

    
    # ........ Phenograph ........    
    
    if(run.Phenograph == 1){  
      
      #NEW:.................. MAKE COMPLEX HEATMAP  ...............    
      
      exp <- subset( gated.cell.dat, GATEsubset==s)[, lapply(.SD, "mean", na.rm = TRUE), by='Annotated_metacluster',  
                                                          .SDcols = cluster.cols]
      
      ## new way to plot complex heatmap
      
      #z-normalize feature expression (iterate over columns with MARGIN=2)
      zscored.exp <- cbind(  exp[,1]  ,   apply(exp[,-1], scale, MARGIN = 2)  ) 
      
      # we gonna bring back only the signals into the matrix:
      zscored.exp.mat <- as.matrix(zscored.exp[,-1])
      
      # prepare column and row names:
      rownames(zscored.exp.mat) <- paste0(exp[[1]] )
      colnames(zscored.exp.mat) <-  sub("\\.\\..*", "", colnames(zscored.exp.mat))
      
      
      col_fun = colorRamp2(c(range(zscored.exp.mat)[1], 0, range(zscored.exp.mat)[2]), c("blue", "white", "red"))
      
      
      
      # so the trick is to initiate column_ha first with the alphabetical order:
      column_ha = HeatmapAnnotation(cluster = rownames(zscored.exp.mat),
                                    
                                    
                                    # this is the Cluster legend. the problem is the "at", which takes over the alphabetical order of the subset.Phenograph_mc_palette
                                    annotation_legend_param = list(
                                      cluster = list(
                                        ncol = 1, 
                                        title = "Cluster",
                                        #title_position = "topcenter",
                                        at = rownames(zscored.exp.mat) ,#rownames(zscored.exp.mat), # names(subset.Phenograph_metacluster_palette),
                                        grid_height = unit(0.02*length(marker.order), "cm"),
                                        grid_width = unit(0.04*length(marker.order), 'cm'),
                                        labels_gp = gpar(fontsize = 0.8*length(marker.order)),
                                        title_gp = gpar(fontsize = 0.8*length(marker.order))
                                      )
                                    ),
                                    
                                    #col = list(cluster=hm_cols),
                                    col = list(cluster=subset.Phenograph_metacluster_palette),
                                    
                                    na_col = "black", #white
                                    show_annotation_name = T
      )
      
      # then initiate the map
      hm <-   Heatmap( t( zscored.exp.mat ) ,
                       col = col_fun,
                       row_order = sub("\\.\\..*", "", marker.order) , # this is the only that needs stripping of transformed_rescaled
                       cluster_columns = fh,
                       top_annotation = column_ha,
                       
                       show_column_names = F, # turn off the names of the clusters under the hm
                       
                       #column_names_rot = 45,
                       rect_gp = gpar(col = "white", lwd = 1),
                       
                       heatmap_legend_param = list(
                         title = "z-score",
                         # direction = 'horizontal',
                         #title_position = "topcenter", # this is the title z-score
                         at=c(-2,-1,0,1,2,3,4),
                         title_position = "lefttop-rot",
                         legend_width = unit(0.25*length(marker.order), "cm"),
                         legend_height = unit(0.3*length(marker.order), "cm"),
                         grid_width = unit(0.02*length(marker.order), 'cm'),
                         labels_gp = gpar(fontsize = 0.6*length(marker.order)),
                         title_gp = gpar(fontsize = 0.8*length(marker.order))
                       )
                       
      )
      
      
      #then bring back column_ha
      
      column_ha = HeatmapAnnotation(cluster = rownames(zscored.exp.mat),
                                    
                                    
                                    # this is the Cluster legend. the problem is the "at", which takes over the alphabetical order of the subset.Phenograph_mc_palette
                                    annotation_legend_param = list(
                                      cluster = list(
                                        ncol = 1, 
                                        title = "Cluster",
                                        #title_position = "topcenter",
                                        at = rownames(zscored.exp.mat)[column_order(hm)] ,#rownames(zscored.exp.mat), # names(subset.Phenograph_metacluster_palette),
                                        grid_height = unit(0.02*length(marker.order), "cm"),
                                        grid_width = unit(0.04*length(marker.order), 'cm'),
                                        labels_gp = gpar(fontsize = 0.8*length(marker.order)),
                                        title_gp = gpar(fontsize = 0.8*length(marker.order))
                                      )
                                    ),
                                    
                                    #col = list(cluster=hm_cols),
                                    col = list(cluster=subset.Phenograph_metacluster_palette),
                                    
                                    na_col = "black", #white
                                    show_annotation_name = T
      )
      
      

      
      # then initiate the map again with the new column_ha order:
      hm <-   Heatmap( t( zscored.exp.mat ) ,
                       col = col_fun,
                       row_order = sub("\\.\\..*", "", marker.order), # this is the only that needs stripping of the qupath string
                       cluster_columns = fh,
                       top_annotation = column_ha,
                       
                       show_column_names = F, # turn off the names of the clusters under the hm
                      
                       #column_names_rot = 45,
                       rect_gp = gpar(col = "white", lwd = 1),
                       
                       heatmap_legend_param = list(
                         title = "z-score",
                         # direction = 'horizontal',
                         #title_position = "topcenter", # this is the title z-score
                         at=c(-2,-1,0,1,2,3,4),
                         title_position = "lefttop-rot",
                         legend_width = unit(0.25*length(marker.order), "cm"),
                         legend_height = unit(0.3*length(marker.order), "cm"),
                         grid_width = unit(0.02*length(marker.order), 'cm'),
                         labels_gp = gpar(fontsize = 0.6*length(marker.order)),
                         title_gp = gpar(fontsize = 0.8*length(marker.order))
                       )
                       
      ) 
      
      # gonna store the hm away as grob and unite it with the UMAP that is not created at this point:
      hm.grob = grid.grabExpr(
        draw(hm,
             heatmap_legend_side = 'left',
             row_sub_title_side = 'left',
             padding = unit(c(20, 2, 2, 2), "mm") # use padding to make the HM squared
        )
      )    
      
      
      
      #Save hm alone LAIA
      #png(paste0(GATEsubsetNames[s],"_Phenograph_annotated.png"), 
      # width = 40*(ncol(zscored.exp.mat[,-1])), 
      # height = 40*(length(marker.order)) )
      
      # grid.draw(hm.grob)
      
      #  invisible(dev.off())  
      
      
    }#end Phenograph block   
    
    
    
   
    # ........ UMAP ........
    # prepare UMAP subtitle string:
    subtitlestring <- paste0("Subset: ", GATEsubsetNames[s], " (",  nrow( subset( gated.cell.dat, GATEsubset==s) ), " annotated cells plotted)" ) 
    if(run.UMAP == 1){  
      g <-  make.colour.plot.adapted( dat=subset( gated.cell.dat, GATEsubset==s) , x.axis = 'UMAP1', y.axis ='UMAP2', 
                                      col.axis ="Annotated_metacluster", 
                                      col.type = 'factor', 
                                      add.label = TRUE, #coloring comes via Annotated_metaclusters. Print these pls
                                      titlestring = paste0("UMAP by annotated clusters")  ,
                                      subtitlestring = subtitlestring,
                                      point.alpha = 0.8, # lets not overdraw the UMAP and reduce alpha of the geom_points
                                      dot.size = 1.6, # less alpha and smaller dots on the half-width-page graphs
                                      colours = "polychrome",
                                      polychromepalette = subset.Phenograph_metacluster_palette,
                                      repel.label.size = 3,
                                      label.label.repelforce = 30, # repulsion between overlapping text labels. default 30 to push them apart
                                      label.datapoint.pullforce = 0.8, # attraction between label and datapoint. default 0.8 to help release the label
                                      min.line.length.tolabel = 0, # minimal distance under which the segment line is not drawn anymore: 0 draw all, Inf turn off even if far
                                      legend.loc = "none",
                                      legend.text.size = 15, # default 18
                                      save.to.disk = F)
      
      
      g2 <-  make.colour.plot.adapted( dat=subset( gated.cell.dat, GATEsubset==s) , x.axis = 'UMAP1', y.axis ='UMAP2', 
                                      col.axis ="Annotated_metacluster", 
                                      col.type = 'factor', 
                                      add.label = TRUE, #coloring comes via Annotated_metaclusters. Print these pls
                                      titlestring = paste0("UMAP by annotated clusters")  ,
                                      subtitlestring = subtitlestring,
                                      point.alpha = 1, # lets not overdraw the UMAP and reduce alpha of the geom_points
                                      dot.size = 0.05, # less alpha and smaller dots on the half-width-page graphs
                                      colours = "polychrome",
                                      polychromepalette = subset.Phenograph_metacluster_palette,
                                      repel.label.size = 4,
                                      label.label.repelforce = 30, # repulsion between overlapping text labels. default 30 to push them apart
                                      label.datapoint.pullforce = 0.4, # attraction between label and datapoint. default 0.8 to help release the label
                                      min.line.length.tolabel = 0, # minimal distance under which the segment line is not drawn anymore: 0 draw all, Inf turn off even if far
                                      legend.loc = "right",
                                      legend.text.size = 15, # default 18
                                      save.to.disk = F)
      
      ggsave(paste0(gsub("[:/\\\\<>|?*]", "_", GATEsubsetNames[s]),"_UMAP.pdf"), g2, width = 7, height = 6, dpi=400) # use the ratio to make the HM squared
      
      
      
      
      
      
      
      
      
      
      
      
      rm(g2)
      
    }#end UMAP block  
    
    
    hm.umap.fig <-  ggarrange(hm.grob, g, widths = c(6,5) )  
    
    #Save UMAP alone LAIA
    #filename <- paste0(GATEsubsetNames[s],"_UMAP_annotated.png") # watch out, make sure the pdf does not contain any .csv since you gonna pull these files next time you run the script!
    #ggsave(filename,g, width = 30, height = 15, units = "cm",dpi = 600)
    
    
    hm.umap.fig <-   annotate_figure(hm.umap.fig,
                                     top = text_grob(paste0("Subset: ", GATEsubsetNames[s]), color = "black", face = "bold", size = 17),
                                     bottom = text_grob(data.source.lyrics, color = "red",
                                                        hjust = 1, x = 1, face = "italic", size = 13),
                                     #left = text_grob("Figure arranged using ggpubr", color = "green", rot = 90),
                                     #right = text_grob(bquote("Superscript: ("*kg~NH[3]~ha^-1~yr^-1*")"), rot = 90),
                                     fig.lab = "Collapsed, annotated", fig.lab.face = "bold"
    ) 
    
    

    ggsave(paste0(gsub("[:/\\\\<>|?*]", "_", GATEsubsetNames[s]),"_Phenograph_annotated.pdf"), hm.grob, width = 15, height = 10, dpi=400) # use the ratio to make the HM squared
    ggsave(paste0(gsub("[:/\\\\<>|?*]", "_", GATEsubsetNames[s]),"_UMAP_Phenograph_annotated.pdf"), hm.umap.fig, width = 15, height = 10, dpi=400) # use the ratio to make the HM squared

    
    height = (length(marker.order))
    width = 3*height # its not 3/2 bang on, since there is a legend in between...
    
    if(min(width, height)<4){
      width = 10*width
      height = width
    }
    
    # now, dpi must not exceeed 50000px. adjust dpi so that we end up with a 15000px image on the long edge:
    dpi <- 15000/max(width, height)
    
    ggsave(  plot = hm.umap.fig, 
             filename = paste0(gsub("[:/\\\\<>|?*]", "_", GATEsubsetNames[s]),"_UMAP_Phenograph_annotated_design2.pdf"), 
             width = width, 
             height = height ,
             limitsize = FALSE,
             dpi = dpi)
    

  } # end run with s through all gatesubsets in gated.cell.dat
  
  
  





############# Plot Phenograph&UMAP now for good with the global colors set #############
# now that Phenograph_metacluster_palette contains all subsets, we need to re-plot the Gate-wise HM and UMAPs.
# but this time we use gated.cell.dat, the collapsed set as source and the numbers from within:

# gated.cell.dat, Collapsed_metacluster colored by Phenograph_metacluster_palette

setwd(OutputDirectory)


for(s in 1:amountofGATEsubsets  ){
  
  # load the cluster columns from our list for the according subset:
  cluster.cols <- cluster.cols.GATEsubsets[[s]] # use double brackets to subset into the vector and not just the list.  
  marker.order <- cluster.cols.GATEsubsets[[s]] # depreciated: global.marker.order[global.marker.order %in% cluster.cols]
  
  # watch out, the make.color.plot engine does not accept a full-range Phenograph_metacluster_palette if you just plot a subset of it.
  # And if you plot the heatmap with the full range palette, they get all listed in the legend. So lets create a temporary Palette subset:
  subset.Phenograph_metacluster_palette <- Phenograph_metacluster_palette[ which( 100*(s) < names(Phenograph_metacluster_palette) & names(Phenograph_metacluster_palette) < 100*(s+1)  ) ]
  
  
  
  # ........ Phenograph ........    
  
  if(run.Phenograph == 1){  
    
    #NEW:.................. MAKE COMPLEX HEATMAP  ...............    
    
    exp <- subset( gated.cell.dat, GATEsubset==s)[, lapply(.SD, "mean", na.rm = TRUE), by='Collapsed_metacluster',  
                                        .SDcols = cluster.cols]
    
    ## new way to plot complex heatmap
    
    #z-normalize feature expression (iterate over columns with MARGIN=2)
    zscored.exp <- cbind(  exp[,1]  ,   apply(exp[,-1], scale, MARGIN = 2)  ) 
    
    # we gonna bring back only the signals into the matrix:
    zscored.exp.mat <- as.matrix(zscored.exp[,-1])
    
    # prepare column and row names:
    rownames(zscored.exp.mat) <- paste0(exp[[1]] )
    colnames(zscored.exp.mat) <-  sub("\\.\\..*", "", colnames(zscored.exp.mat))
    
    col_fun = colorRamp2(c(range(zscored.exp.mat)[1], 0, range(zscored.exp.mat)[2]), c("blue", "white", "red"))
    
    
    # this needs adaptation later on. take markers off the heat map when plotting the subgates. use with a good argumentation, this is a tad dangerous
    # when plotting the final colors on HM, we can take markers out that would disturb the interpretation:
    
   # zscored.exp.mat <- subset(zscored.exp.mat, select=-c(get("NKx6.1"),get("GLP.1R")))
    
    # also take them out of the current marker.order:
    #marker.order.noIsletMarkers <-  grep(paste0( c("NKx6.1", "GLP.1R" )   , collapse = "|") , marker.order, invert = TRUE, value = TRUE)
    marker.order.noIsletMarkers <-marker.order
    
    
    
    
    column_ha = HeatmapAnnotation(cluster = rownames(zscored.exp.mat),
                                  
                                  annotation_legend_param = list(
                                    cluster = list(
                                      ncol = 2, 
                                      title = "Cluster",
                                      title_position = "topcenter",
                                      at = names(subset.Phenograph_metacluster_palette),
                                      grid_height = unit(0.02*length(marker.order.noIsletMarkers), "cm"),
                                      grid_width = unit(0.02*length(marker.order.noIsletMarkers), 'cm'),
                                      labels_gp = gpar(fontsize = 0.8*length(marker.order.noIsletMarkers)),
                                      title_gp = gpar(fontsize = 0.8*length(marker.order.noIsletMarkers))
                                    )
                                  ),
                                  
                                  #col = list(cluster=hm_cols),
                                  col = list(cluster=subset.Phenograph_metacluster_palette),
                                  
                                  na_col = "black", #white
                                  show_annotation_name = FALSE
    )
    
    
    
    hm <-   Heatmap( t( zscored.exp.mat ) ,
                     col = col_fun,
                     row_order = sub("\\.\\..*", "", marker.order) , # this is the only that needs stripping of transformed_rescaled
                     cluster_columns = fh,
                     top_annotation = column_ha,
                     rect_gp = gpar(col = "white", lwd = 1),
                     heatmap_legend_param = list(
                       title = "z-score",
                       direction = 'horizontal',
                       title_position = "topcenter",
                       legend_width = unit(0.25*length(marker.order.noIsletMarkers), "cm"),
                       grid_width = unit(0.02*length(marker.order.noIsletMarkers), 'cm'),
                       labels_gp = gpar(fontsize = 0.8*length(marker.order.noIsletMarkers)),
                       title_gp = gpar(fontsize = 0.8*length(marker.order.noIsletMarkers))
                     )
                     
    )
    
   
    
    # gonna store the hm away as grob and unite it with the UMAP that is not created at this point:
    hm.grob = grid.grabExpr(
      draw(hm,
           heatmap_legend_side = 'top',
           row_sub_title_side = 'left',
           padding = unit(c(2, 2, 2, 10), "mm"))
    )    
    
    
    
  }#end Phenograph block   
  
  
  # ........ UMAP ........
  if(run.UMAP == 1){  
  
    subtitlestring <- paste0("Subset: ", GATEsubsetNames[s], " (",  nrow(clustered.gated.subsets[[s]]), " cells plotted)" )   
    
    
  g <-  make.colour.plot( clustered.gated.subsets[[s]] , 'UMAP1', 'UMAP2', "Collapsed_metacluster", 'factor', 
                          titlestring = paste0("UMAP by collapsed Phenograph MC")  ,
                          subtitlestring = subtitlestring,
                          point.alpha = 0.4, # lets not overdraw the UMAP and reduce alpha of the geom_points
                          dot.size = 0.9, # less alpha and smaller dots on the half-width-page graphs
                          colours = "polychrome",
                          polychromepalette = subset.Phenograph_metacluster_palette,
                          add.label = TRUE, #col.axis are the Phenograph/FlowSOM metacluster numbers. Print these pls
                          legend.loc = "top",
                          legend.text.size = 15, # default 18
                          save.to.disk = F)
}#end UMAP block  


hm.umap.fig <-  ggarrange(hm.grob, g, widths = c(3,2) )  




hm.umap.fig <-   annotate_figure(hm.umap.fig,
                                 top = text_grob(paste0("Subset: ", GATEsubsetNames[s]), color = "black", face = "bold", size = 17),
                                 bottom = text_grob(data.source.lyrics, color = "red",
                                                    hjust = 1, x = 1, face = "italic", size = 13),
                                 #left = text_grob("Figure arranged using ggpubr", color = "green", rot = 90),
                                 #right = text_grob(bquote("Superscript: ("*kg~NH[3]~ha^-1~yr^-1*")"), rot = 90),
                                 fig.lab = "Collapsed", fig.lab.face = "bold"
) 




png(paste0(GATEsubsetNames[s],"_UMAP_Phenograph.png"), 
    width = 60*(ncol(zscored.exp.mat[,-1])), 
    height = 29*(length(marker.order.noIsletMarkers)) )

print(hm.umap.fig)

invisible(dev.off())


} # end run with s through all gatesubsets in gated.cell.dat

# end reprint and overwrite HM UMAP per subset



# once the clusters are annotated and the individual subsets are united, plot each cluster on the scatter plots of marker.comparison markers,
if(plot.cluster_wise.scatterplots==1){


setwd(OutputDirectory)
dir.create('5 - QC Cluster distribution at gates')
setwd('5 - QC Cluster distribution at gates')   

# since we got the GATEsubset variable in gated.cell.dat, we can now simply run through it and check where the clusters are relative to the gates:
message("\nQC 5: Plotting every cluster in its subsets with their gates..." ) 



for(s in 1:amountofGATEsubsets){
  
  
setwd(OutputDirectory)  
setwd('5 - QC Cluster distribution at gates')   
dir.create(GATEsubsetNames[s])  
setwd(GATEsubsetNames[s])    
  
  
# we create an assembled plot for every subset:  
list.pos <- 1
l = list()    
  
  temp.OTdata <- subset( gated.cell.dat, GATEsubset==s )
  
  message(paste0("\n------------ Subset Number ",s, " -----------\n      ",GATEsubsetNames[s],"\n      ",nrow(temp.OTdata)," cells"   ))
  
  
  pb <- progress_bar$new(format = "[:bar] :percent [Plotting clusters by marker expression | :eta]",
                         total = length(unique(temp.OTdata$Annotated_metacluster))*length(markercomparision)   , #
                         show_after=0, #allows to call it right way 
                         current = "|",    # Current bar character
                         clear = FALSE) # make it persist
  pb$tick(0) # call in the progress bar without any progress, just to show it

for(c in unique(temp.OTdata$Annotated_metacluster) ){
  

    subset.Phenograph_metacluster_palette <-   subset(annotated.Phenograph_metacluster_palette, Collapsed_metacluster > 100*(s) &  Collapsed_metacluster < 100*(s+1)   )  %>%
      pull( color,Annotated_metacluster  ) # pulling color out and naming it Annotated_metacluster
  # and then just replace the remaining clusters with gray:
    subset.Phenograph_metacluster_palette[ names(subset.Phenograph_metacluster_palette) %!in% c ] <- "#B8B8B8"
    
    
  
  for(m in c(1:length(markercomparision))){    
  
    xMarker <- sub("_[^_]+$", "", markercomparision[[m]][1] )
    yMarker <- sub("_[^_]+$", "", markercomparision[[m]][2] )

    
g <-  ggplot(temp.OTdata, 
         aes_string(x = markercomparision[[m]][1],       y = markercomparision[[m]][2] ), 
        
  )+
  geom_point( aes(colour = as.factor(temp.OTdata$Annotated_metacluster)),  alpha=0.8,size=1)   +
      
      
  geom_point( data=subset(temp.OTdata, Annotated_metacluster == c), color= subset.Phenograph_metacluster_palette[c],  alpha=1,size=3 )   + 
  
      scale_color_manual( values = subset.Phenograph_metacluster_palette,
                          name="Cluster")+  
      scale_x_log10(breaks = breaks, 
                    minor_breaks = minor_breaks,
                    labels = trans_format("log10", math_format(10^.x))) +
      
      scale_y_log10(breaks = breaks, 
                    minor_breaks = minor_breaks,
                    labels = trans_format("log10", math_format(10^.x))) +

      annotation_logticks()   +
      ylab( yMarker ) + 
      xlab( xMarker ) +  
      coord_fixed() + 
      guides(colour = guide_legend(override.aes = list(size=2, alpha=1),  #we also gonna override the dot size how its depicted in the legend
                                   byrow=TRUE) )+ # and I want the elements listed row-by-row
      
      theme_minimal()  +
      theme(#axis.title.x = element_blank(), 
        legend.title=element_text(size=13), # turn off with element_blank(),
        legend.text=element_text(size=11),
        axis.text=element_text(size=13),#,axis.text.x=element_text(size=10),
        axis.title.y=element_text(size=12), axis.title.x=element_text(size=12),
         legend.position="top", 
        legend.key.width=unit(1,"line"),
        #aspect.ratio=1 #force ggplot on squared plotting # turn off if you force fixed coords of course
      )



filename <- paste0("QC_scatterplot_Cluster_",c,"_",xMarker,"vs",yMarker,".jpg")
ggsave(filename,g, width = 30, height = 20, units = "cm",dpi = 600)
pb$tick()  
  
  
  }# end run along all markers we need to plot per markercombo 
} # end run with c along all clusters of the current gatesubset to plot 

} # end run with s through all gatesubsets in gated.cell.dat
}# end QC plot each cluster highlighted on the marker comparison scatter plots  
  
  
  
  
  
  
  






setwd(OutputDataDirectory)
 fwrite( gated.cell.dat , 'gated.cell.dat.csv')

        
 
 
cat("\n\n\n----------------------------------- \n End of STAGE 3. \n Clusters collapsed and subsets re-united into gated.cell.dat \n Proceed to STAGE 4 or farther \n-----------------------------------\n\n") 
 
 
 
           
        
} # end stage 3: annotate and collapse clusters

        
########################## STAGE 4 ##########################        
if(whichSTAGE == 4){  
  

# we gonna use the collapsed and bound back togeter gated.cell.dat and not the clustered.gated.subsets[[s]]
# before we keep on plotting cellular data, we need to mask the lymph nodes if there were any.
  
# since we dont know which of both imageJ masks are provided, lets load the additional columns into gated.cell.dat  
  
##########################  Use EDM matrices to calculate cluster-wise distances to cellular signals ##########################   
  
  
  
  setwd(OutputDirectory)
  dir.create('FijiMask_Plots')
  setwd('FijiMask_Plots')
  FijiMasks <- getwd()
  

#-------------------------- LOAD lymph nodes --------------------------  
  
if(mask.lymphnodes==1){

  # we try to mask the lymph nodes with a mask created with a second imagej script:
  #now we gonna switch over to the imagej output folder:
  setwd(lymphnode.binarymask.path)
  # not sure if this is needed at all:
  binarymask.list <- list.files(pattern = "\\.txt$") # this ensures that only .txt is listed, and not .ttxt or so
  
  # we not gonna run along the EDM.list, since this contains ALL ROIs, also the positive controls spleen ROIs. 
  # Instead we use the ROIs present in gated.cell.dat: plot.ROIs  
  maxEDMwaringlyrics <- c("--------- Quality Control WARNING ---------")
  
  

  
  
  # initiate a dataframe to store our distances as we run through the ROIs
  temp.df <- data.frame()
  temp.df2 <- data.frame()
  lymphnode.area.data <- data.frame()
  

    
    if (lymphnode.mask.type == "BM") {
      
      
pb <- progress_bar$new(format = "[:bar] :percent [Locating lymph nodes | :eta]",
                             total = length(plot.ROIs), #
                             show_after=0, #allows to call it right way 
                             current = "|",    # Current bar character
                             clear = FALSE) # make it persist
pb$tick(0) # call in the progress bar without any progress, just to show it
      
      
      
    for(r in plot.ROIs){  
    
    # we subset the cells of the current ROI out: 
    temp.ROI <- subset(gated.cell.dat, ROI %in% r) 
    # we also need to do that with all.cells since we need the info which cells were all inside.lymphnode when all cells were still there:
    temp.ROI2 <- subset(all.cells,  ROI  %in% r) 
    
    # in order to match the cell into its lymph node mask, we need to round their position to an integer
    # lets just create two more columns in the subset
    temp.ROI$x.integer <- round(temp.ROI$x,0)
    temp.ROI$y.integer <- round(temp.ROI$y,0)
    
    temp.ROI2$x.integer <- round(temp.ROI2$x,0)
    temp.ROI2$y.integer <- round(temp.ROI2$y,0)
    
    # lets pre-load the according lymph node binary mask for that specific ROI:
    setwd(lymphnode.binarymask.path)
    lymphnode.mask <- read.table(file=paste0(r,"_BM.txt"), 
                                 head=FALSE,
                                 sep="\t" 
    ) 
    
    # I wanna be sure to track the ones with a lot of mask area:
    mask.pixels <- sum(rowSums(lymphnode.mask==255))
    total.pixels <- nrow(lymphnode.mask)*ncol(lymphnode.mask)
    percent.mask <- round( 100*mask.pixels/total.pixels , 1)
    
    # only consider the well-behaving masks. for the moment we take also here the value that is given for the islets.
    if(percent.mask < max.mask.coverage.percent){
      
      # for density calculations on these ROIs, we need to substract the area of the lymph node
      # by the time we calculate this, we will be missing the lymphnode cells in the dataset,
      # so better to create the area data set here..
      lymphnode.area.data <- rbind(lymphnode.area.data, data.frame(ROI = r,
                                   lymphnode.area = mask.pixels
      ))
      
      # we gonna put some warnings in, since we are running in autopilot here without graphical feedback:
      
      #give me a warning if a ROI has a lot of coverage nevertheless:
      if(percent.mask > 5){  
        maxEDMwaringlyrics <- c( maxEDMwaringlyrics, paste0(r, " has a lyph node coverage of ", percent.mask , "%"    ))
      }
      
      
      #now we use the binary.mask.scan to find lymphnode coordinates on the ROIs:
      temp.ROI$inside.lymphnode <- mapply( function(x,y) lymphnode.mask[y,x] , x=temp.ROI$x.integer, y=temp.ROI$y.integer)/255
      
      temp.ROI2$inside.lymphnode <- mapply( function(x,y) lymphnode.mask[y,x] , x=temp.ROI2$x.integer, y=temp.ROI2$y.integer)/255
      
      
      if(plot.lymphnodeflags.on.ROIs==1){
      setwd(FijiMasks)
      plot.ratio <-  ncol(lymphnode.mask) /nrow(lymphnode.mask) 
      
      g <-  ggplot( temp.ROI ,aes(x=x.integer,y=y.integer,color=as.factor(inside.lymphnode) )) +
        geom_point(size=2.5, alpha=0.8) + 
        theme_minimal() +
        scale_y_reverse()+
        #scale_color_gradient(low='gray',high='blue3')
        scale_color_manual(values = c( "0" = "gray",
                                       "1" ="blue3"))+
        labs(color = "Lymph node flag")+
        labs(title= paste("ROI ID ", r) , 
             subtitle=paste0(nrow(temp.ROI), " cells plotted"),
             y="y [px]", 
             x= "x [px]" 
        )+
        theme(#axis.title.x = element_blank(), 
          legend.title=element_text(size=13), # turn off with element_blank(),
          legend.text=element_text(size=11),
          axis.text=element_text(size=13),#,axis.text.x=element_text(size=10),
          axis.title.y=element_text(size=12), axis.title.x=element_text(size=12),
          legend.position="top", 
          legend.key.width=unit(1,"line"),
          aspect.ratio=1 #force ggplot on squared plotting # turn off if you force fixed coords of course
        ) 
      
      plot(g)
      
      filename <- paste0(r,"_lymphnode_flag.jpg")
      ggsave(filename,g, width = 20*plot.ratio, height = 20, units = "cm",dpi = 600)
      }#plot lymphnode flags on ROIs if asked to
      
      # end protect from running into masks with an unrealistically high islet coverage 
    }else{
      maxEDMwaringlyrics <- c( maxEDMwaringlyrics, paste0(" ! ",r, " was ignored: lymph node coverage  (", percent.mask , "%) above given threshold"    ))
      # either way we need to create the lymph node column for the upcoming rbind
      temp.ROI$inside.lymphnode <- NA
      temp.ROI2$inside.lymphnode <- NA
    }
  
  
  
    #this is slow, but for the moment Im lazy to find a fast mapply way to process the flags and plot...
    temp.df <- rbind(temp.df,temp.ROI  )  # gated.cell.dat
    temp.df2 <- rbind(temp.df2,temp.ROI2  ) # all.cells
    pb$tick()
    
  } # end run along all ROIs in gated.cell.dat with r
  
  # in the end move temp.df back to gated.cell.dat and remove the other
  gated.cell.dat <- temp.df
  all.cells <- temp.df2
  rm(temp.df,temp.df2)
  
  #only print the warning if we actually ran into weird lymph node binary masks
  if(length(maxEDMwaringlyrics)>1){ print(maxEDMwaringlyrics) }
  }#end process lymphnode.mask.type being Binary Mask
  
  
  if (lymphnode.mask.type == "EDM") {
    
    
    pb <- progress_bar$new(format = "[:bar] :percent [Calc. euclid. dist. to lymph nodes | :eta]",
                           total = length(plot.ROIs), #
                           show_after=0, #allows to call it right way 
                           current = "|",    # Current bar character
                           clear = FALSE) # make it persist
    pb$tick(0) # call in the progress bar without any progress, just to show it
    
    
    for(r in plot.ROIs){
    # the beginning of EDM-type is similar to Binary Mask processing:
    
    # we subset the cells of the current ROI out: 
    temp.ROI <- subset(gated.cell.dat, ROI %in% r) 
    # we also need to do that with all.cells since we need the info which cells were all inside.lymphnode when all cells were still there:
    temp.ROI2 <- subset(all.cells,  ROI  %in% r) 
    
    # in order to match the cell into its lymph node mask, we need to round their position to an integer
    # lets just create two more columns in the subset
    temp.ROI$x.integer <- round(temp.ROI$x,0)
    temp.ROI$y.integer <- round(temp.ROI$y,0)
    
    temp.ROI2$x.integer <- round(temp.ROI2$x,0)
    temp.ROI2$y.integer <- round(temp.ROI2$y,0)
    
    # lets pre-load the according lymph node binary mask for that specific ROI:
    setwd(lymphnode.binarymask.path)
    lymphnode.mask  <- read.table(file=paste0(r,"_EDM.txt"), 
                                 head=FALSE,
                                 sep="\t" 
    ) 
    

    
    
    #protect from analyzing ROIs that have no islets at all. For reasons I dont understand, an empty mask is filled with the number 46341. 
    # Change this number in the STAGE settings in imageJ.emptymask if needed.
    if(max(lymphnode.mask)!=0 && max(lymphnode.mask)!=imageJ.emptymask  ){
      
      
      # I wanna be sure to track the ones with a lot of mask area:
      mask.pixels <- sum(rowSums(lymphnode.mask==0))
      total.pixels <- nrow(lymphnode.mask)*ncol(lymphnode.mask)
      percent.mask <- round( 100*mask.pixels/total.pixels , 1)
      short.edge <- min(nrow(lymphnode.mask),ncol(lymphnode.mask) )
    
    
    # only consider the well-behaving masks. for the moment we take also here the value that is given for the islets.
    if(percent.mask < max.mask.coverage.percent){
      
      # for density calculations on these ROIs, we need to substract the area of the lymph node
      # by the time we calculate this, we will be missing the lymphnode cells in the dataset,
      # so better to create the area data set here..
      lymphnode.area.data <- rbind(lymphnode.area.data, data.frame(ROI = r,
                                                                   lymphnode.area = mask.pixels
      ))
      
      # we gonna put some warnings in, since we are running in autopilot here without graphical feedback:
      
      #give me a warning if a ROI has a lot of coverage nevertheless:
      if(percent.mask > 5){  
        maxEDMwaringlyrics <- c( maxEDMwaringlyrics, paste0(r, " has a lymph node coverage of ", percent.mask , "%"    ))
      }
      
      
      # give me a warning if there is a ton of small islets across the ROI so that the max distance on the whole ROI is half the short edge of the image
      if(max(lymphnode.mask) < 0.5*short.edge ){
        maxEDMwaringlyrics <- c( maxEDMwaringlyrics, paste0(r, " has a max distance between lymph nodes of only ", max(lymphnode.mask) , "px (",round( 100*max(lymphnode.mask)/short.edge   ,1),"% of short edge)"    ))
      }
      

      
      #now we can use the defined function and a beautiful R function to calculate distance from both xy columns in one go:
      temp.ROI$lymph.node.distance <- extract.distance(mask=lymphnode.mask, 
                                                       x=temp.ROI$x.integer, 
                                                       y=temp.ROI$y.integer)
      
      temp.ROI2$lymph.node.distance <- extract.distance(mask=lymphnode.mask, 
                                                       x=temp.ROI2$x.integer, 
                                                       y=temp.ROI2$y.integer)

      
      # the binary mask just produces in-out-flag. so lets do that as well here: 
      temp.ROI$inside.lymphnode <-  ifelse( temp.ROI$lymph.node.distance <= (0+lymphnode.extension)   , 1, 0)
      temp.ROI2$inside.lymphnode <-  ifelse( temp.ROI2$lymph.node.distance <= (0+lymphnode.extension)   , 1, 0)
      

            
      
      if(plot.lymphnodeflags.on.ROIs==1){
        setwd(FijiMasks)
        plot.ratio <-  ncol(lymphnode.mask) /nrow(lymphnode.mask) 
        
        g <-  ggplot( temp.ROI ,aes(x=x.integer,y=y.integer,color=as.factor(inside.lymphnode) )) +
          geom_point(size=1.8, alpha=1) + 
          geom_point(data=subset(temp.ROI, lymph.node.distance>0 & lymph.node.distance<= (0+lymphnode.extension) ), color="blue3",  size=1.8, alpha=1) +
          theme_minimal() +
          scale_y_reverse()+
          scale_color_manual(breaks = c("0", "1"), 
                             values=c("gray50", "red"),
                             labels = c("out", "in"),
                             name="Lymph node flag")+
          
          
          
          labs(title= paste("ROI ID ", r) , 
               subtitle=paste0(nrow(temp.ROI), " cells plotted"),
               y="y [px]", 
               x= "x [px]" 
          )+
          theme(#axis.title.x = element_blank(), 
            legend.title=element_text(size=13), # turn off with element_blank(),
            legend.text=element_text(size=11),
            axis.text=element_text(size=13),#,axis.text.x=element_text(size=10),
            axis.title.y=element_text(size=12), axis.title.x=element_text(size=12),
            legend.position="top", 
            legend.key.width=unit(1,"line"),
            aspect.ratio=1 #force ggplot on squared plotting # turn off if you force fixed coords of course
          ) 
        
        plot(g)
        
        filename <- paste0(r,"_lymphnode_flag.jpg")
        ggsave(filename,g, width = 20*plot.ratio, height = 20, units = "cm",dpi = 600)
      }#plot lymphnode flags on ROIs if asked to
      
      # end protect from running into masks with an unrealistically high islet coverage 
    }else{
      maxEDMwaringlyrics <- c( maxEDMwaringlyrics, paste0(" ! ",r, " was ignored: lymph node coverage  (", percent.mask , "%) above given threshold"    ))
      # either way we need to create the lymph node column for the upcoming rbind
      temp.ROI$inside.lymphnode <- NA
      temp.ROI2$inside.lymphnode <- NA
    }
    
    }else{# end protecting analyzing ROIs that have no lymph node at all. For the distance, we use a distance high enough to not be caught:
      temp.ROI$lymph.node.distance <- imageJ.emptymask
      temp.ROI2$lymph.node.distance <- imageJ.emptymask
      temp.ROI$inside.lymphnode <- 0
      temp.ROI2$inside.lymphnode <- 0
    }
    
    #this is slow, but for the moment Im lazy to find a fast mapply way to process the flags and plot...
    temp.df <- rbind(temp.df,temp.ROI  )  # gated.cell.dat
    temp.df2 <- rbind(temp.df2,temp.ROI2  ) # all.cells
    pb$tick()
    
  } # end run along all ROIs in gated.cell.dat with r
  
  # in the end move temp.df back to gated.cell.dat and remove the other
  gated.cell.dat <- temp.df
  all.cells <- temp.df2
  rm(temp.df,temp.df2)
  
  #only print the warning if we actually ran into weird lymph node binary masks
  if(length(maxEDMwaringlyrics)>1){ print(maxEDMwaringlyrics) }
    
    
    
    
    
  }#end process lymphnode.mask.type being EDM
  
  
  
  
  
}# end processing lymphnode masks
  
  
  
  
  
if(run.EDM.calculation==1){
  
#-------------------------- LOAD Distance map  --------------------------
    
    
    #now we gonna switch over to the imagej output folder:
    setwd(EDMpath)
    # not sure if this is needed at all:
    EDM.list <- list.files(pattern = "\\.txt$") # this ensures that only .txt is listed, and not .ttxt or so
    
    
   if(exists('Measurements.csv')) fiji.measurement <- fread('Measurements.csv')
    
    
    maxEDMwaringlyrics <- c("--------- Quality Control WARNING ---------")
    
    
    pb <- progress_bar$new(format = "[:bar] :percent [Calc. euclid. distances to mask | :eta]",
                           total = length(plot.ROIs), #
                           show_after=0, #allows to call it right way 
                           current = "|",    # Current bar character
                           clear = FALSE) # make it persist
    pb$tick(0) # call in the progress bar without any progress, just to show it   
    
    # initiate a dataframe to store our distances as we run through the ROIs
    temp.df <- data.frame()
    temp.df2 <- data.frame()
    
    
    for(r in plot.ROIs){
      
      # we subset the cells of the current ROI out: 
      temp.ROI <- subset(gated.cell.dat, ROI %in% r)  
      
      # we also need to do that with all.cells since we need the info which cells were all inside.lymphnode when all cells were still there:
      temp.ROI2 <- subset(all.cells,  ROI  %in% r) 
      
      # in order to match the cell to its distance on the EDM file, we need to round their position to an integer
      # lets just create two more columns in the subset
      temp.ROI$x.integer <- round(temp.ROI$x,0)
      temp.ROI$y.integer <- round(temp.ROI$y,0)
    
      temp.ROI2$x.integer <- round(temp.ROI2$x,0)
      temp.ROI2$y.integer <- round(temp.ROI2$y,0)
    
      # lets pre-load the according EDM for that specific ROI:
      setwd(EDMpath)
      
      
      #names(gated.cell.dat)[names(gated.cell.dat) == 'Ablation ID'] <- 'Ablation.ID'
      

      
      # load the EDM file for the current ROI if r is matched, else load it using "Ablation ID" column in gated.cell.dat
      if( paste0(r,"_EDM.txt") %in% EDM.list ){
        
        EDM <- read.table(file=paste0(r,"_EDM.txt"), 
                          head=FALSE,
                          sep="\t"     )
      
  }else if( paste0( unique(temp.ROI$Ablation.ID),"_EDM.txt")  %in% EDM.list ) {
        EDM <- read.table(file= paste0( unique(temp.ROI$Ablation.ID),"_EDM.txt"), 
                          head=FALSE,
                          sep="\t" 
        )
      }else{
        print(paste0("EDM file for ",r," not found. Skipping this ROI."))
        next
      }
      
      
      #protect from analyzing ROIs that have no islets at all. For reasons I dont understand, an empty mask is filled with the number 46341. 
      # Change this number in the STAGE settings in imageJ.emptymask if needed.
      if(max(EDM)!=0 && max(EDM)!=imageJ.emptymask  ){
        
        
        # I wanna be sure to track the ones with a lot of mask area:
        mask.pixels <- sum(rowSums(EDM==0))
        total.pixels <- nrow(EDM)*ncol(EDM)
        percent.mask <- round( 100*mask.pixels/total.pixels , 1)
        short.edge <- min(nrow(EDM),ncol(EDM) )
        
        # only consider the well-behaving masks, these have typically 0.4-3% islet coverage, not more.
        if(percent.mask < max.mask.coverage.percent){
          
          # we gonna put some warnings in, since we are running in autopilot here without graphical feedback:
          
          #give me a warning if a ROI has a lot of coverage nevertheless:
          if(percent.mask > 5){  
            maxEDMwaringlyrics <- c( maxEDMwaringlyrics, paste0(r, " has a mask coverage of ", percent.mask , "%"    ))
          }
          
          # give me a warning if there is a ton of small islets across the ROI so that the max distance on the whole ROI is half the short edge of the image
          if(max(EDM) < 0.5*short.edge ){
            maxEDMwaringlyrics <- c( maxEDMwaringlyrics, paste0(r, " has a max distance between islets of only ", max(EDM) , "px (",round( 100*max(EDM)/short.edge   ,1),"% of short edge)"    ))
          }
          
          
          #now we can use the defined function and a beautiful R function to calculate distance from both xy columns in one go:
          #temp.ROI$distance <-mapply(extract.distance, x=temp.ROI$x.integer, y=temp.ROI$y.integer)
          
          
          #scaling.factors[, paste0( batchnormalize.marker[m] )] 
          
          
          temp.ROI[, paste0(EDM.name, ".distance")] <- extract.distance(mask=EDM, 
                                x=temp.ROI$x.integer, 
                                y=temp.ROI$y.integer)
          
          
          temp.ROI2[, paste0(EDM.name, ".distance")] <- extract.distance(mask=EDM, 
                                                                        x=temp.ROI2$x.integer, 
                                                                        y=temp.ROI2$y.integer)
          
          # now, be aware that we dont include an option to define mask extension here. 
          # your maks has to sit precisely anyway, if you want to draw conclusions, so there should not be any need to fiddle with the mask here.
         #temp.ROI[, paste0("inside.",EDM.name )] <-  ifelse( temp.ROI[, paste0(EDM.name, ".distance")] == 0 , 1, 0) this fails since ifelse asks for a vector
          temp.ROI[, paste0("inside.",EDM.name )] <-  ifelse( temp.ROI[[ paste0(EDM.name, ".distance")]] == 0 , 1, 0)
          temp.ROI2[, paste0("inside.",EDM.name )] <-  ifelse( temp.ROI2[[ paste0(EDM.name, ".distance")]] == 0 , 1, 0)
          
          # we also define a peri-mask area if asked:
          if(islet.rim.width>0){
            
            temp.ROI[, paste0("inside.rim.",EDM.name )] <-  ifelse( temp.ROI[[ paste0(EDM.name, ".distance")]] > 0 & temp.ROI[[ paste0(EDM.name, ".distance")]] <= islet.rim.width , 1, 0)  
            temp.ROI2[, paste0("inside.rim.",EDM.name )] <-  ifelse( temp.ROI2[[ paste0(EDM.name, ".distance")]] > 0 & temp.ROI2[[ paste0(EDM.name, ".distance")]] <= islet.rim.width , 1, 0) 
          }
          
          
          
          
          # end protect from running into masks with an unrealistically high islet coverage 
        }else{
          maxEDMwaringlyrics <- c( maxEDMwaringlyrics, paste0(" ! ",r, " was ignored: islet mask coverage  (", percent.mask , "%) above given threshold"    ))
          # either way we need to create the distance column for the upcoming rbind
          temp.ROI[, paste0(EDM.name, ".distance")] <- NA
          temp.ROI[, paste0("inside.",EDM.name )] <- NA
          temp.ROI2[, paste0(EDM.name, ".distance")] <- NA
          temp.ROI2[, paste0("inside.",EDM.name )] <- NA
          
          if(islet.rim.width>0){  
            temp.ROI[, paste0("inside.rim.",EDM.name )] <-  NA 
            temp.ROI2[, paste0("inside.rim.",EDM.name )] <-  NA
          }
          
        }
      }else{# end protecting analyzing ROIs that have no islet at all
        temp.ROI[, paste0(EDM.name, ".distance")] <- NA
        temp.ROI[, paste0("inside.",EDM.name )] <- NA
        temp.ROI2[, paste0(EDM.name, ".distance")] <- NA
        temp.ROI2[, paste0("inside.",EDM.name )] <- NA
        if(islet.rim.width>0){  
          temp.ROI[, paste0("inside.rim.",EDM.name )] <-  NA
          temp.ROI2[, paste0("inside.rim.",EDM.name )] <-  NA
          }
      }
      
      
      
      
      
      
      #this is slow, but for the moment Im lazy to find a fast mapply way to process the flags and plot...
      temp.df <- rbind(temp.df,temp.ROI  )
      temp.df2 <- rbind(temp.df2,temp.ROI2  )
      pb$tick()
      
    } # end run along all ROIs in gated.cell.dat with r
    
    # in the end move temp.df back to gated.cell.dat and remove the other
    gated.cell.dat <- temp.df
    all.cells <- temp.df2
    rm(temp.df,temp.df2)
    
    #only print the warning if we actually ran into weird lymph node binary masks
    if(length(maxEDMwaringlyrics)>1){ print(maxEDMwaringlyrics) }
      
      
   # IM_20221222_OT8_001 IM_20221222_OT8_s0_p4_r2_a2_ac
    setwd(OutputDirectory)
    
    dir.create('5 - QC Masks on ROIs')
    setwd('5 - QC Masks on ROIs')
    
    
    for(r in plot.ROIs){
      
    # this is just a temporary plot that needs to be removed later:
    g <- ggplot( subset(all.cells, ROI %in% r) ,aes(x=x.integer,y=y.integer,color=as.factor(inside.rim.TE) )) +
      geom_point(size=1.2, alpha=1) + 
      theme_minimal() +
      scale_y_reverse()+
      #scale_color_gradient(low='gray',high='blue3')
      # scale_color_manual(values = c( "0" = "gray",
      #                               "1" ="blue3"))+

      
       #  values=c("#8C8C8C", "#EB9718", "#14D81E"), # 814165 is mask outline
      scale_color_manual(breaks = c("0", "1"), 
                         values=c("#8C8C8C", "#EB9718"),
                         labels = c("out", "in"),
                         name="TE rim flag")+
      
            geom_point(data=subset(all.cells, ROI %in% r & inside.TE==1), color="red",  size=1.2, alpha=1) +
      
      
      labs(title= paste("ROI ID ", r)  , 
           subtitle= paste0(nrow(subset(all.cells, ROI %in% r)), " cells plotted. Rim width: ",islet.rim.width, "px"),
           y="y [px]", 
           x= "x [px]" 
      )+
      theme(#axis.title.x = element_blank(), 
        legend.title=element_text(size=13), # turn off with element_blank(),
        legend.text=element_text(size=11),
        axis.text=element_text(size=13),#,axis.text.x=element_text(size=10),
        axis.title.y=element_text(size=12), axis.title.x=element_text(size=12),
        legend.position="right", 
        legend.key.width=unit(1,"line"),
        aspect.ratio=1 #force ggplot on squared plotting # turn off if you force fixed coords of course
      ) 
    
    plot(g)
    
    filename <- paste0(r,"_",islet.rim.width, "pxRIM.jpg")
    ggsave(filename,g, width = 30, height = 22, units = "cm",dpi = 600)
    
    }
    
    
    
    
    
    ggplot( subset(all.cells, ROI %in% c("IM_20221222_OT6_006")) ,aes(x=x.integer,y=y.integer,color=CK18_Yb172_t2 )) +
      geom_point(size=1.1, alpha=1) + 
      theme_minimal() +
      scale_y_reverse()
    
    
    
    
    
    
  }# end calulating distances to EDMask
  

  
  
  # ok, we are getting there. We store away the gated.cell.dat anew  
  setwd(OutputDataDirectory)
  fwrite( gated.cell.dat , 'gated.cell.dat.csv') 
# and we do that for all.cells as well:
  fwrite( all.cells , 'all.cells.csv')  
  message("Gated.cell.dat and all.cells updated in ./Data folder")
  


#  [1] "--------- Quality Control WARNING ---------"                          
#  [2] "BeCeFa_20201109_OT11_s0_p7_r4_a4_ac has a lyph node coverage of 10.8%"
#  [3] "BeCeFa_20210611_OT18_s0_p5_r1_a1_ac has a lyph node coverage of 6.5%" 
#  [4] "BeCeFa_20211022_OT23_s0_p6_r5_a5_ac has a lyph node coverage of 6.8%" 
#  [5] "BeCeFa_20211022_OT23_s0_p6_r6_a6_ac has a lyph node coverage of 9.8%" 
#  [6] "BeCeFa_20220128_OT30_s0_p4_r4_a4_ac has a lyph node coverage of 5.3%" 
  
  
  #  [1] "--------- Quality Control WARNING ---------"                                               
  #[2] "BeCeFa_20210312_OT13_s0_p3_r1_a1_ac has a max distance of only 959px (47.5% of short edge)"
  #[3] "BeCeFa_20210423_OT15_s0_p5_r3_a3_ac has a mask coverage of 43.8%"                          
  #[4] "BeCeFa_20210423_OT15_s0_p5_r3_a3_ac has a max distance of only 426px (38.6% of short edge)"
  #[5] "BeCeFa_20210423_OT15_s0_p5_r4_a4_ac has a mask coverage of 5.2%"                           
  ##[6] "BeCeFa_20210609_OT17_s0_p5_r2_a2_ac has a mask coverage of 25.8%"                          
  #[7] "BeCeFa_20210609_OT17_s0_p5_r2_a2_ac has a max distance of only 320px (30% of short edge)"  
  #[8] "BeCeFa_20210609_OT17_s0_p5_r5_a5_ac has a max distance of only 987px (46.3% of short edge)"
  ##[9] "BeCeFa_20210611_OT18_s0_p5_r1_a1_ac has a mask coverage of 17.1%"                          
  ##[10] "BeCeFa_20210611_OT18_s0_p5_r1_a1_ac has a max distance of only 470px (26.3% of short edge)"
  #[11] "BeCeFa_20210611_OT18_s0_p5_r5_a5_ac has a mask coverage of 7.9%"                           
  #[12] "BeCeFa_20210614_OT19_s0_p5_r3_a3_ac has a max distance of only 735px (47.2% of short edge)"
  #[13] "BeCeFa_20210614_OT19_s0_p5_r4_a4_ac has a max distance of only 752px (37.6% of short edge)"
  #[14] "BeCeFa_20211015_OT22_s0_p7_r4_a4_ac has a mask coverage of 5.3%"                           
  #[15] "BeCeFa_20211015_OT22_s0_p7_r4_a4_ac has a max distance of only 421px (22.1% of short edge)"
  #[16] "BeCeFa_20211022_OT23_s0_p6_r3_a3_ac has a mask coverage of 15.5%"                          
  #[17] "BeCeFa_20211022_OT23_s0_p6_r3_a3_ac has a max distance of only 348px (24.1% of short edge)"
  #[18] "BeCeFa_20211022_OT23_s0_p6_r5_a5_ac has a mask coverage of 11.2%"                          
  #[19] "BeCeFa_20211022_OT23_s0_p6_r5_a5_ac has a max distance of only 352px (20% of short edge)"  
  #[20] "BeCeFa_20220128_OT30_s0_p4_r1_a1_ac has a max distance of only 433px (38.7% of short edge)"
  #[21] "BeCeFa_20210318_OT14_s0_p5_r3_a3_ac has a max distance of only 882px (43.2% of short edge)"
  #[22] "BeCeFa_20211029_OT24_s0_p6_r4_a4_ac has a mask coverage of 27.1%"                          
  #[23] "BeCeFa_20211029_OT24_s0_p6_r4_a4_ac has a max distance of only 438px (28.9% of short edge)"  
  
  
  

  
if(mask.lymphnodes==1){
  
  if(lymphnode.extension>0){
    subtitlestring <- paste0("Lymph node mask extended by ",lymphnode.extension, "px")
  }else{
    subtitlestring <- paste0("Lymph node mask not extended")
  }
  
  
  # in most cases, you made the binary mask to get rid of cells on the ROIs. 
  # So instead of constant subsetting into inside.lymphnode==0, lets just create a gated.cell.dat with these cells masked:
    masked.gated.cell.dat <- subset(gated.cell.dat, inside.lymphnode==0)
    
    # you gotta put some order into the ROI ID, if not the ROI-wise plots get all mixed up along the x-axis:
    
    masked.gated.cell.dat$ROI <- factor(masked.gated.cell.dat$ROI,      # Reordering group factor levels
                                           levels = plot.ROIs )
    
    
    # and store it away:
    setwd(OutputDataDirectory)
    fwrite( masked.gated.cell.dat , 'masked.gated.cell.dat.csv') 
    
    
    lymphnode.frequency.dat<-data.frame()
    for(r in plot.ROIs){
      
     temp.df <- data.frame(as.data.frame(table( subset(gated.cell.dat, ROI==r)[,inside.lymphnode])))
     temp.df <- cbind(temp.df, ROI=r)
     temp.df <- cbind(temp.df, Batch=unique(subset(gated.cell.dat, ROI==r)$Batch))
     
     lymphnode.frequency.dat <- rbind(lymphnode.frequency.dat, temp.df)
     
    }
    
    names(lymphnode.frequency.dat) <- c("lymphnode.flag", "Count","ROI", "Batch")
    
    
    lymphnode.frequency.dat$ROI <-  gsub(".*OT","",lymphnode.frequency.dat$ROI )
    lymphnode.frequency.dat$ROI <- paste0("OT",lymphnode.frequency.dat$ROI )
 
    
    setwd(OutputDirectory)
    # and Im particularly interested in how many cells are masked on each TMA!
    g <- ggplot(lymphnode.frequency.dat, aes(x=factor(ROI), y = Count, fill=lymphnode.flag )) + 
      #facet_wrap(~Cluster )+
      geom_bar(stat="identity", position = position_stack(reverse = T)) +
      scale_fill_manual(breaks = c("0", "1"), 
                        values=c("gray50", "blue3"),
                        labels = c("out", "in"),
                        name="Lymphnode flag")+
    
  
      labs(title="Contribution of lymphnodes to total segemented cells on ROIs" , 
           subtitle=subtitlestring, 
           #caption="Created by M.Barone", 
           y="Cell count on ROI", 
           x=""
      )+
      theme_light() +
      theme(#panel.background = element_rect(fill = "grey90", colour = "grey55", size = 0.5), # change 'colour' to black for informative axis
            axis.title.x=element_text(color="grey15", size=11),
            axis.title.y=element_text(color="grey15", size=11),
            axis.text=element_text(size=10),
            #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
            legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
            legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
            #legend.title=element_blank(),
            plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
            plot.subtitle = element_text(color = "Black", size = 12, hjust = 0),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
      )
    
    plot(g)
    
    filename <- paste0("LymphnodeContribution_perROI_",lymphnode.extension,"pxLymphmaskextend.jpg")
    ggsave(filename,g, width = 50, height = 30, units = "cm",dpi = 1200)
    
    
    
    
    
    
    # if we got the mask in, lets plot the relative contribution of every cluster on the ROIs:
    

   
    
    for(s in 1:amountofGATEsubsets  ){
    #  for(s in 2:2  ){ 
      
      if(lymphnode.extension>0){
        subtitlestring <- paste0( nrow(subset(masked.gated.cell.dat, GATEsubset == s   )), " non-lymph node cells in subset. Lymph node mask extended by ",lymphnode.extension, "px")
      }else{
        subtitlestring <- subtitlestring <- paste0( nrow(subset(masked.gated.cell.dat, GATEsubset == s   )), " non-lymph node cells in subset. Lymph node mask not extended")
      }
    
    # we gotta create a character vector again:
    subset.Phenograph_metacluster_palette <- subset(annotated.Phenograph_metacluster_palette, Collapsed_metacluster > 100*(s) &  Collapsed_metacluster < 100*(s+1)   ) 
    
    temp.pal <- subset.Phenograph_metacluster_palette[,color] 
    names(temp.pal ) <-subset.Phenograph_metacluster_palette$Annotated_metacluster
    subset.Phenograph_metacluster_palette <- temp.pal
    rm(temp.pal)
    subset.Phenograph_metacluster_palette <-  subset.Phenograph_metacluster_palette[!duplicated(subset.Phenograph_metacluster_palette)]
    
     # the ROI ID tag is too long, so lets just cut it back:
    short.ROI.tags <- sort(unique(subset(masked.gated.cell.dat, GATEsubset ==s   )$ROI))
   # short.ROI.tags <- factor(short.ROI.tags,      # Reordering group factor levels
  #                           levels = plot.ROIs )
    
    short.ROI.tags <-     paste0("OT", gsub(".*OT","", short.ROI.tags ) )
    

    
    
  g<-  ggplot(subset(masked.gated.cell.dat, GATEsubset ==s   ), aes(x=ROI, y = Annotated_metacluster, fill=Annotated_metacluster )) + 
      #facet_wrap(~Cluster )+
      geom_bar(stat="identity", position = position_stack(reverse = T))+
      scale_fill_manual( values = subset.Phenograph_metacluster_palette )+ # push the polychrome palette in

     scale_x_discrete(labels= short.ROI.tags) +
      labs(title=paste0("Subset ", GATEsubsetNames[s], ", non-lymphnode cells: Contribution of clusters on ROIs ")   , 
            subtitle=subtitlestring, 
           #caption="Created by M.Barone", 
           y="Cell count on ROI", 
           x=""
      )+
      theme_light() +
      theme(#panel.background = element_rect(fill = "grey90", colour = "grey55", size = 0.5), # change 'colour' to black for informative axis
        axis.title.x=element_text(color="grey15", size=11),
        axis.title.y=element_text(color="grey15", size=11),
        axis.text=element_text(size=8),
        axis.text.y=element_blank(),
        #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
        legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
        legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
        #legend.title=element_blank(),
        plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
        plot.subtitle = element_text(color = "Black", size = 12, hjust = 0),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
      )
    
        plot(g)
    
    filename <- paste0(GATEsubsetNames[s],"_non-lymphnode_ClusterComposition_perROI_",lymphnode.extension,"pxLymphmaskextend.png")  
    ggsave(filename,g, width = 50, height = 40, units = "cm",dpi = 1200)
    
      }#end run with s along the gate subsets
 
    
  } # end some plotting if lymphnode masking was turned on
    
    
}#end stage 4: adding binary and distance measures into gated.cell.dat
        
########################## STAGE 5 ##########################         
if(whichSTAGE == 5){  
  
# STAGE 5: further QC and descriptive plots 
  
  

  setwd(OutputDirectory) 
  

  
  
###### Split UMAP by sample and ROI and TMA to check if there is one subset in data that provides disproportionally much to the UMAP ###### 
  
if(plot.DimRed.by.subset == 1){
    
    
    setwd(OutputDirectory) 
    dir.create('6 - QC UMAPs by sample')
    setwd('6 - QC UMAPs by sample')

    pb <- progress_bar$new(format = "[:bar] :percent [Plot. sample-wise UMAP per subset | :eta]",
                           total = amountofGATEsubsets * length(tissue.type.vector), #
                           show_after=0, #allows to call it right way 
                           current = "|",    # Current bar character
                           clear = FALSE) # make it persist
    pb$tick(0) # call in the progress bar without any progress, just to show it  
  
    i <- 1
    u <- list()
    for(d in tissue.type.vector){
      
      # I wanna have the UMAPs separated, and not an overlay over all subsets. 
      for(s in 1:amountofGATEsubsets  ){
        titlestring = paste0("UMAP by collapsed Phenograph MC") 
        subtitlestring <- paste0("Subset: ", GATEsubsetNames[s], " (",  nrow( subset(gated.cell.dat, TumorType == d & GATEsubset == s )  ), " cells plotted)" ) 
        
        u[[i]] <- make.colour.plot.adapted( subset(gated.cell.dat, TumorType == d & GATEsubset == s ) , x.axis = 'UMAP1', y.axis ='UMAP2', 
                                    col.axis = sample.col, # sample.col is Patient_ID
                                    col.type = 'factor', 
                                    titlestring = paste0("UMAP subset: Tissue type ",d, " by patient")  ,
                                    subtitlestring = subtitlestring,
                                    point.alpha = 0.5, # lets not overdraw the UMAP and reduce alpha of the geom_points
                                    dot.size = 1.5, # less alpha and smaller dots on the half-width-page graphs
                                    colours = "polychrome",
                                    polychromepalette = Conditions_palette,
                                    repel.label.size = 5,
                                    label.label.repelforce = 50, # repulsion between overlapping text labels. default 30 to push them apart
                                    label.datapoint.pullforce = 0.1, # attraction between label and datapoint. default 0.8 to help release the label
                                    min.line.length.tolabel = 0, # minimal distance under which the segment line is not drawn anymore: 0 draw all, Inf turn off even if far
                                    add.label = TRUE, #col.axis are the Phenograph/FlowSOM metacluster numbers. Print these pls
                                    legend.loc = "right", # in this plotting routine, turn legend off!!!!
                                    legend.text.size = 15, # default 18
                                    save.to.disk = F)
        
        i<- i+1 
        pb$tick()
      } # split along subsets with s
    } # end run along d for both diets
    
   # leg <- get_legend(u[[8]])
    
   # # g <- do.call("arrangeGrob", common.legend = TRUE, legend="bottom", c(u, ncol=amountofGATEsubsets)) # use a grob to save the obj
    
   # #g <-  ggarrange(plotlist=u, widths = c(amountofGATEsubsets,1), common.legend = TRUE, legend="right")
   # #g <-  ggarrange(plotlist=u, widths = c(amountofGATEsubsets,1) )
    #g <-  ggarrange(plotlist=u, legend="right", legend.grob = leg, nrow = 2 )
    
    
    g <- do.call("arrangeGrob", c(u, nrow=2)) # use a grob to save the obj
    filename <- "UMAP_byTissuetype_coloredbyPatientID.jpg" # watch out, make sure the pdf does not contain any .csv since you gonna pull these files next time you run the script!
    ggsave(filename,g, width = 55, height = 30, units = "cm",dpi = 800)
  }#end plot UMAP by subset
    
 
if(plot.DimRed.by.ROI == 1){
      
      setwd(OutputDirectory)
      dir.create('6 - QC UMAPs by ROI')
      setwd('6 - QC UMAPs by ROI')
      
      pb <- progress_bar$new(format = "[:bar] :percent [Plot. ROI-wise UMAP per subset | :eta]",
                             total = amountofGATEsubsets * amount.of.ROIs, #
                             show_after=0, #allows to call it right way 
                             current = "|",    # Current bar character
                             clear = FALSE) # make it persist
      pb$tick(0) # call in the progress bar without any progress, just to show it  
      
      u <- list()
      for(r in plot.ROIs){
        
        
        # I wanna have the UMAPs separated, and not an overlay over all subsets. 
        for(s in 1:amountofGATEsubsets  ){
          
          # same game as before: make.color.plot engine does not accept a full-range Phenograph_metacluster_palette if you just plot a subset of it.
          # we gotta create a character vector again:    
          subset.Phenograph_metacluster_palette <- subset(annotated.Phenograph_metacluster_palette, Collapsed_metacluster > 100*(s) &  Collapsed_metacluster < 100*(s+1)   ) 
          temp.pal <- subset.Phenograph_metacluster_palette[,color] 
          names(temp.pal ) <-subset.Phenograph_metacluster_palette$Annotated_metacluster
          subset.Phenograph_metacluster_palette <- temp.pal
          rm(temp.pal)
          
          subset.Phenograph_metacluster_palette <-  subset.Phenograph_metacluster_palette[!duplicated(subset.Phenograph_metacluster_palette)] 
          # ^ this can now be pushed into the plotting engine:
          
          
          u[[s]] <- make.colour.plot.adapted( subset(gated.cell.dat, ROI == r & GATEsubset == s ) , x.axis = 'UMAP1', y.axis ='UMAP2', 
                                      col.axis ="Annotated_metacluster", 
                                      col.type = 'factor',
                                      add.label = TRUE, #coloring comes via Annotated_metaclusters. Print these pls
                                      titlestring = paste0("Subset ", GATEsubsetNames[s]) ,
                                      # subtitlestring = NA,
                                      point.alpha = 0.4, # lets not overdraw the UMAP and reduce alpha of the geom_points
                                      dot.size = 5, # less alpha and smaller dots on the half-width-page graphs
                                      colours = "polychrome",
                                      polychromepalette = subset.Phenograph_metacluster_palette,
                                      repel.label.size = 3,
                                      label.label.repelforce = 30, # repulsion between overlapping text labels. default 30 to push them apart
                                      label.datapoint.pullforce = 0.8, # attraction between label and datapoint. default 0.8 to help release the label
                                      min.line.length.tolabel = 0, # minimal distance under which the segment line is not drawn anymore: 0 draw all, Inf turn off even if far
                                      legend.loc = "right",
                                      legend.text.size = 15, # default 18
                                      save.to.disk = F)
          
          pb$tick()
        }# end run along s for the both subsets
        
        leg <- get_legend(u[[1]])
        
        #g <-  ggarrange(plotlist=u, legend="right", legend.grob = leg )
        g <-  ggarrange(plotlist=u )
        
        g <-   annotate_figure(g,
                               top = text_grob(paste0("ROI ID: ",r), color = "black", face = "bold", size = 17),
                               bottom = text_grob(data.source.lyrics, color = "red",
                                                  hjust = 1, x = 1, face = "italic", size = 13),
                               #left = text_grob("Figure arranged using ggpubr", color = "green", rot = 90),
                               #right = text_grob(bquote("Superscript: ("*kg~NH[3]~ha^-1~yr^-1*")"), rot = 90),
                               fig.lab = paste0("Condition on ROI: ",gated.cell.dat$Diet_week[ grep( r , gated.cell.dat$ROI)[1]]), fig.lab.face = "bold"
        ) 
        
        
        
        png(paste0(r,"_subsetUMAPs.png"), 
            width = 1200, 
            height = 600,
        )
        print(g)
        invisible(dev.off()) 
        
        
        
      }# run along all ROIs with r and store the grobs
      
      
      
    }#end plot UMAP per ROI
    
    
if(plot.DimRed.by.TMA == 1){
      
      
      setwd(OutputDirectory) 
      dir.create('6 - QC UMAPs by TMA')
      setwd('6 - QC UMAPs by TMA')
      
      pb <- progress_bar$new(format = "[:bar] :percent [Plot. TMA-wise UMAP per subset | :eta]",
                             total = amountofGATEsubsets * amount.of.TMAs, #
                             show_after=0, #allows to call it right way 
                             current = "|",    # Current bar character
                             clear = FALSE) # make it persist
      pb$tick(0) # call in the progress bar without any progress, just to show it  
      
      
      u <- list()
      for(r in all_images){
        
        #temp.OTdata  <- subset(gated.cell.dat, ROI == r )  
        
        # I wanna have the UMAPs separated, and not an overlay over all subsets. 
        for(s in 1:amountofGATEsubsets  ){
          
          # same game as before: make.color.plot engine does not accept a full-range Phenograph_metacluster_palette if you just plot a subset of it.
          # we gotta create a character vector again:    
          subset.Phenograph_metacluster_palette <- subset(annotated.Phenograph_metacluster_palette, Collapsed_metacluster > 100*(s) &  Collapsed_metacluster < 100*(s+1)   ) 
          temp.pal <- subset.Phenograph_metacluster_palette[,color] 
          names(temp.pal ) <-subset.Phenograph_metacluster_palette$Annotated_metacluster
          subset.Phenograph_metacluster_palette <- temp.pal
          rm(temp.pal)
          
          subset.Phenograph_metacluster_palette <-  subset.Phenograph_metacluster_palette[!duplicated(subset.Phenograph_metacluster_palette)] 
          # ^ this can now be pushed into the plotting engine:
          
          u[[s]] <- make.colour.plot.adapted( subset(gated.cell.dat, Batch == r & GATEsubset == s ) , x.axis = 'UMAP1', y.axis ='UMAP2', 
                                      col.axis ="Annotated_metacluster", 
                                      col.type = 'factor',
                                      add.label = TRUE, #coloring comes via Annotated_metaclusters. Print these pls
                                      titlestring = paste0("Subset ", GATEsubsetNames[s]) ,
                                      # subtitlestring = NA,
                                      point.alpha = 0.4, # lets not overdraw the UMAP and reduce alpha of the geom_points
                                      dot.size = 3.5, # less alpha and smaller dots on the half-width-page graphs
                                      colours = "polychrome",
                                      polychromepalette = subset.Phenograph_metacluster_palette,
                                      repel.label.size = 3,
                                      label.label.repelforce = 30, # repulsion between overlapping text labels. default 30 to push them apart
                                      label.datapoint.pullforce = 0.8, # attraction between label and datapoint. default 0.8 to help release the label
                                      min.line.length.tolabel = 0, # minimal distance under which the segment line is not drawn anymore: 0 draw all, Inf turn off even if far
                                      legend.loc = "right",
                                      legend.text.size = 15, # default 18
                                      save.to.disk = F)
          
          pb$tick() 
        }# end run along s for the both subsets
        
        
        
        #leg <- get_legend(u[[1]])
        #g <-  ggarrange(plotlist=u, legend="right", legend.grob = leg )
        g <-  ggarrange(plotlist=u )
        
        # there is a good chance you got multiple condtions on that ROI:
        
        if( length( unique(gated.cell.dat$Diet_week[ grep( r , gated.cell.dat$ROI)])  )==1   ){
          
          g <-   annotate_figure(g,
                                 top = text_grob(paste0("TMA ID: ",r), color = "black", face = "bold", size = 17),
                                 bottom = text_grob(data.source.lyrics, color = "red",
                                                    hjust = 1, x = 1, face = "italic", size = 13),
                                 #left = text_grob("Figure arranged using ggpubr", color = "green", rot = 90),
                                 #right = text_grob(bquote("Superscript: ("*kg~NH[3]~ha^-1~yr^-1*")"), rot = 90),
                                 fig.lab = paste0("Condition on ROI: ",gated.cell.dat$Diet_week[ grep( r , gated.cell.dat$ROI)[1]]), fig.lab.face = "bold"
          )
          
        }else{
          
          g <-   annotate_figure(g,
                                 top = text_grob(paste0("TMA ID: ",r), color = "black", face = "bold", size = 17),
                                 bottom = text_grob(data.source.lyrics, color = "red",
                                                    hjust = 1, x = 1, face = "italic", size = 13),
                                 #left = text_grob("Figure arranged using ggpubr", color = "green", rot = 90),
                                 #right = text_grob(bquote("Superscript: ("*kg~NH[3]~ha^-1~yr^-1*")"), rot = 90),
                                 fig.lab = paste0("Multiple conditions on TMA!"), fig.lab.face = "bold"
          ) 
        }
        
        
        png(paste0(r,"_subsetUMAPs.png"), 
            width = 1200, 
            height = 600,
        )
        
        print(g)
        invisible(dev.off()) 
        
      }# run along all ROIs with r and store the grobs
      
      
      
    }#end plot UMAP per TMA    
    
  
  
 
###### Spatial plotting of cluster data ###### 

  
# Read in spatial data object if one of the two plots are asked for:
  
if( plot.each.subset.clusters.on.ROIs==1 | plot.all.subset.clusters.on.ROIs==1    ){
  
  
  setwd(InputDirectory)

    #if you wanna run this part of the code just shortly, you can switch to the WD like so:
   #setwd("D:/IMI/IMC/Western-diet/Analysis/spectre/Output 1 - add masks/20220601/Data")


  # you might need to re-run STAGE 0 to check the batch normalization graphs, in that case:
  
  # protect from re-loading a huge file:
  if( !exists("spatial.dat") ){
    #  list.files(getwd(), '.qs')
    message("\nLoading spatial data file (spatial.dat.qs)...")
    timeStart<-Sys.time()
    spatial.dat <- qread('spatial.dat.qs') 
    timeEnd<-Sys.time()
    
    message(paste0("Loading spatial.dat.qs completed in ", round(difftime(timeEnd, timeStart, units='mins'),2), " min") )
}# END just load spatial.dat if it doesnt exist
  
  
  ### Prepare for plotting party in a separate folder:
    setwd(OutputDirectory)
  dir.create('7 - QC Clusters on ROIs')
  setwd('7 - QC Clusters on ROIs')  
  
  
}# end only load spatial.dat if the plotting engine is asked for  
    
#plot.ROIs.rest <- plot.ROIs[71:88]

if(plot.each.subset.clusters.on.ROIs == 1){
  
  pb <- progress_bar$new(format = "[:bar] :percent [Plotting clusters of each subset on ROIs | :eta]",
                         total = amountofGATEsubsets * amount.of.ROIs, #
                         show_after=0, #allows to call it right way 
                         current = "|",    # Current bar character
                         clear = FALSE) # make it persist
  pb$tick(0) # call in the progress bar without any progress, just to show it   
  
  
  for(s in 1:amountofGATEsubsets){
  #  for(s in 2:amountofGATEsubsets){
    
    setwd(OutputDirectory)
    setwd('7 - QC Clusters on ROIs')
    dir.create(paste0('Subset ', s))
    setwd(paste0('Subset ', s))
    
  #  plot.ROIsrest <- plot.ROIs[27:length(plot.ROIs)]
    # for the paper we are only showing these 7 ROIs:
   # plot.ROIsrest <- plot.ROIs[c(7,5,18,76,48,20,62)]
    
    
   # for(i in plot.ROIsrest ){ #plot.ROIsrest
    for(i in plot.ROIs ){ #plot.ROIsrest
      
      temp <- gated.cell.dat[gated.cell.dat[['ROI']] == i & GATEsubset == s,]
      
      
         subset.Phenograph_metacluster_palette <-   subset(annotated.Phenograph_metacluster_palette, Annotated_metacluster %in% unique(temp$Annotated_metacluster) ) %>%
       pull( color,Annotated_metacluster  ) # pulling color out and naming it Annotated_metacluster
      
      
      
      
      
      #temp <- masked.gated.cell.dat[masked.gated.cell.dat[['ROI']] == i & GATEsubset == s,]
      
      # we are only interested in spatial plots of 4 clusters:
     # temp <- subset(temp, Annotated_metacluster %in% c("activ. effector-like CD8+",
    #                                                    "F4.80lo macrophages",
    #                                                    "F4.80- macrophages",
     #                                                   "F4.80hi macrophages"))
      
      
      # and now we only bring out the 4 colors:
     #   subset.Phenograph_metacluster_palette <-   subset(annotated.Phenograph_metacluster_palette, Annotated_metacluster %in% c("activ. effector-like CD8+",
    #                                                                                                                             "F4.80lo macrophages",
    #                                                                                                                             "F4.80- macrophages",
     #                                                                                                                            "F4.80hi macrophages") ) %>%
     # pull( color,Annotated_metacluster  ) # pulling color out and naming it Annotated_metacluster
      
         subtitlelyrics <- paste0(GATEsubsetNames[s], ". Patient ID: ",unique(temp$Patient_ID), " - ", unique(temp$TumorType), " tissue"  )

      
      #v1.0, and already adapted: 
      make.spatial.plot.adapted(dat=spatial.dat, 
                                image.roi= i, 
                                image.channel='aSMA_Yb174', #  Insulin_Pr141 in adv spat 1 add masks, this is used via variable base DNA2_Ir193. Here it makes more sense to see the islets..
                                mask.outlines = 'cell.mask', # and this via variable mask. 
                                cell.dat = temp,
                                cell.col = 'Annotated_metacluster', 
                                cell.col.type = 'factor',
                                dot.size = 3, #4 be aware that size will be overridden in the legend to 4 
                                dot.alpha= 1, # be aware that alpha will be overridden in the legend to 1 
                                cell.colours = "polychrome", # this is how the segmentation cell coordinate dots are colored
                                polychromepalette = subset.Phenograph_metacluster_palette,
                                amountofdiscretcols= 4, # wild guess amount.of.MC ,
                                image.min.threshold = 0.03, #0.04 lower col onset percentile #defaut 0-0.5-1, alternat. 0.025-0.5-0.975
                                image.med.threshold = 0.50, # mid point col percentile with:
                                markersignal.mid.col = "#12010A", #2F031B watch out, most of the pixels are background or dont express DNA, so the median is aaall the waaay down."black",##2F031B # this is default dark blue "#2B739B". here we change to mix into the mask col: so its a slighly darker reddish purple
                                image.max.threshold = 0.94, # high point col percentile with:
                                markersignal.high.col = "#EEE4E9", # this is default blue "#56B4E9". here we change to mix into the mask col: so its a slighly darker reddish purple B71E73
                                image.mask.colour = "#814165", # default CC79A7
                                titlestring= paste0("Ablation ID: ",i  ),
                                subtitlestring= subtitlelyrics,
                                legend.loc = "right",
                                plot.width = 18,
                                plot.height = 12,
                                file.extension = ".pdf"
      )
      
      pb$tick()
      
    }# end run along all ROIs
  }# end run along all subsets
  
  
  
}# end plot clusters of each subset on ROIs

if(plot.all.subset.clusters.on.ROIs == 1){
  
  
  setwd(OutputDirectory)
  setwd('7 - QC Clusters on ROIs')
  dir.create(paste0('All Subsets'))
  setwd(paste0('All Subsets'))
  
  pb <- progress_bar$new(format = "[:bar] :percent [Plotting clusters of all subsets on ROIs | :eta]",
                         total = amount.of.ROIs, #
                         show_after=0, #allows to call it right way 
                         current = "|",    # Current bar character
                         clear = FALSE) # make it persist
  pb$tick(0) # call in the progress bar without any progress, just to show it 
  

    for(i in plot.ROIs ){
      
        temp <- gated.cell.dat[gated.cell.dat[['ROI']] == i,]
        
        
        subset.Phenograph_metacluster_palette <-   subset(annotated.Phenograph_metacluster_palette, Annotated_metacluster %in% unique(temp$Annotated_metacluster) ) %>%
          pull( color,Annotated_metacluster  ) # pulling color out and naming it Annotated_metacluster
        
        
        
        
        subtitlelyrics <- paste0("Patient ID: ",unique(temp$Patient_ID), " - ", unique(temp$TumorType), " tissue"  )
 
        #v1.0, and already adapted: 
        make.spatial.plot.adapted(dat=spatial.dat, 
                          image.roi= i, 
                          image.channel='aSMA_Yb174', #  Insulin_Pr141 in adv spat 1 add masks, this is used via variable base DNA2_Ir193. Here it makes more sense to see the islets..
                          mask.outlines = 'cell.mask', # and this via variable mask. 
                          cell.dat = temp,
                          cell.col = 'Annotated_metacluster', 
                          cell.col.type = 'factor',
                          dot.size = 1.2, # be aware that size will be overridden in the legend to 4 
                          dot.alpha= 1, # be aware that alpha will be overridden in the legend to 1 
                          cell.colours = "polychrome", # this is how the segmentation cell coordinate dots are colored
                          polychromepalette = subset.Phenograph_metacluster_palette,
                          amountofdiscretcols= 4, # wild guess amount.of.MC 
                          image.min.threshold = 0.03, #0.04 lower col onset percentile #defaut 0-0.5-1, alternat. 0.025-0.5-0.975
                          image.med.threshold = 0.50, # mid point col percentile with:
                          markersignal.mid.col = "#12010A", #2F031B watch out, most of the pixels are background or dont express DNA, so the median is aaall the waaay down."black",##2F031B # this is default dark blue "#2B739B". here we change to mix into the mask col: so its a slighly darker reddish purple
                          image.max.threshold = 0.92, # high point col percentile with:
                          markersignal.high.col = "#EEE4E9", # this is default blue "#56B4E9". here we change to mix into the mask col: so its a slighly darker reddish purple B71E73
                          image.mask.colour = "#814165", # default CC79A7
                          titlestring= paste0("Ablation ID: ",i ),
                          subtitlestring= subtitlelyrics,
                          legend.loc = "right",
                          plot.width = 18,
                          plot.height = 12,
                          file.extension = ".pdf"
                           )
  
        pb$tick()
        
    }


}# end plot clusters of all subsets on ROIs



    
}#end stage 5: do additional QC plots 
  



########################## STAGE 6 ########################## 
# STAGE 6: Plot stuff very project-specific stuff. this needs adaptation for sure. pack anything into the run.old.code bracket and bring stuff out to adapt
if(whichSTAGE == 6){  
  
  # supply stuff we skipped along the stages:
  
  # From now on we wont do any raw data processing, so let us create a small string that is plotted in all upcoming graphs, to keep track of the processing data:
  data.source.lyrics <- c("Data source: ")
  
  if(correct.spill.over==1 ){
    data.source.lyrics <- paste0(data.source.lyrics, "spillover-corr., ")
  }
  
  if(calculate.scaling.factors==1 ){
    data.source.lyrics <- paste0(data.source.lyrics, "batch-norm., ")
  }
  
  if(transform.rawdata ){
    data.source.lyrics <- paste0(data.source.lyrics, "asinh_cf",asinh.cofactor)
  }
  if(minmax.norm){
    data.source.lyrics <- paste0(data.source.lyrics, ", minmax norm.")
  }
  if(zscore.norm){
    data.source.lyrics <- paste0(data.source.lyrics, ", z-score norm.")
  }
  
  

  
  # lets assign some cute colors and use it by calling scale_fill_manual(values = cond_palette)+
  if(run.old.code==T){
  tumorStage_grouped_palette <- c("I" = "#EEA800",
                                  "II & III" = "#C5582B", 
                                  "IVA & IVB" = "#6F2405"
                            
  )
  
  
  longtimesurvival_palette <- c("0" = "#712D06",
                                  "1" = "#29A50B"
                                  
  )
  
  longtimesurvival_palette_string <- c("<36months" = "#712D06",
                                ">36months" = "#29A50B"
                                
  )
  

  
  tissue_palette <- c("ROI" = "#482703",
                      "tumor" = "#BA400D",
                      "stroma" = "#AF7E09"
                                
  )
  }# some Felsenstein Project color palettes
  
  
  
  ##----------- Phenium --------------
  #     CREATE polychrome palette for the cell classes
  colRange <- unique(PCFdata[["Classification"]]) 
  colRange <- colRange[order(colRange)]
  set.seed(723451) # for reproducibility
  class_palette <- createPalette(length(colRange), c("#e77b21"), M=100000, prefix = "")
  names(class_palette) <- colRange
  #swatch(class_palette)
  
  
  
  percentile <- as.numeric(quantile( all.cells$CD4..Cell..Mean, # we compute the quantiles from non-ctr data
                                     probs = c( 0.05 , 0.5, 0.95),na.rm=TRUE )) # ,na.rm=TRUE can cause some weird behaviour of limits
  lower.threshold <- percentile[1]
  median <- percentile[2]
  upper.threshold <- percentile[3]
  
  
  u<-list()
  u[[1]]<-  ggplot( all.cells ,
          aes(x=Centroid.X.µm,y=Centroid.Y.µm,color=CD4..Cell..Mean,fill=CD4..Cell..Mean )) + #,color=as.factor(inside.lymphnode)
    geom_point(size=0.2, alpha=1) + 
    
    theme_minimal() +
    scale_y_reverse()+
    coord_fixed() + 
    # extract Tissue.type from sample.meta for the current i:
    labs(title= paste0("PCF: CD4 Protein Expression"  ) ,
         subtitle = paste0("signal clipped at 5th and 95th") 
    )+
    
    
    scico::scale_color_scico(palette = "roma", direction=-1, 
                             midpoint = median,
                             limits=c(lower.threshold,upper.threshold),
                             na.value = "#8C0172"
    )+ # bilbao roma
    
    scico::scale_fill_scico(palette = "roma", direction=-1, 
                             midpoint = median,
                             limits=c(lower.threshold,upper.threshold),
                             na.value = "#8C0172"
    )+ # bilbao roma
    
    theme(panel.background = element_rect(fill = "grey10", colour = "grey15", linewidth = 0.5), # change 'colour' to black for informative axis
          panel.grid = element_line(color = "grey15"),
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=8),
          #axis.text.x = element_text(angle = 45, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          legend.position="left",
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )
  
  

percentile <- as.numeric(quantile( all.cells$CD3D, # we compute the quantiles from non-ctr data
                                   probs = c( 0.01 , 0.5, 0.85),na.rm=TRUE )) # ,na.rm=TRUE can cause some weird behaviour of limits
lower.threshold <- percentile[1]
median <- percentile[2]
upper.threshold <- percentile[3]


u[[2]]<-  ggplot( all.cells ,
          aes(x=Centroid.X.µm,y=Centroid.Y.µm,color=CD3D,fill=CD3D )) + #,color=as.factor(inside.lymphnode)
    geom_point(size=0.2, alpha=1) + 
    
    theme_minimal() +
    scale_y_reverse()+
    coord_fixed() + 
    # extract Tissue.type from sample.meta for the current i:
    labs(title= paste0("Xenium: CD3D gene count"  )   ,
         subtitle = paste0("signal clipped at 1st and 85th") 
    )+
    
    scico::scale_color_scico(palette = "roma", direction=-1, 
                             midpoint = median,
                             limits=c(lower.threshold,upper.threshold),
                             na.value = "#8C0172"
    )+ # bilbao roma
    
    scico::scale_fill_scico(palette = "roma", direction=-1, 
                            midpoint = median,
                            limits=c(lower.threshold,upper.threshold),
                            na.value = "#8C0172"
    )+ # bilbao roma
    
  theme(panel.background = element_rect(fill = "grey10", colour = "grey15", linewidth = 0.5), # change 'colour' to black for informative axis
        panel.grid = element_line(color = "grey15"),
        axis.title.x=element_text(color="grey15", size=11),
        axis.title.y=element_text(color="grey15", size=11),
        axis.text=element_text(size=8),
        #axis.text.x = element_text(angle = 45, hjust=1),
        #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
        legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
        legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
        #legend.title=element_blank(),
        legend.position="left",
        plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
        plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
  )
  
  
  
  g <- do.call("arrangeGrob", c(u, ncol=2)) # use a grob to save the obj
  filename <- "Phenium_protCD4_geneCD3D_sidebyside.jpg" # watch out, make sure the pdf does not contain any .csv since you gonna pull these files next time you run the script!
  ggsave(filename,g, width = 55, height = 30, units = "cm",dpi = 800)
  
  

  
  percentile <- as.numeric(quantile( all.cells$CD4..Cell..Mean, # we compute the quantiles from non-ctr data
                                     probs = c( 0.01 , 0.5, 0.85),na.rm=TRUE )) # ,na.rm=TRUE can cause some weird behaviour of limits
  lower.threshold <- percentile[1]
  median <- percentile[2]
  upper.threshold <- percentile[3]
  
  
  
  g <- ggplot(all.cells, aes(x=Ki67..Cell..Mean, y=CD3E, color=CD4..Cell..Mean)) + # shape=Correlating_MC,
    geom_point(size=1)+

    scico::scale_color_scico(name = "CD4 (PCF)",
                             palette = "roma", direction=-1, 
                             midpoint = median,
                             limits=c(lower.threshold,upper.threshold),
                             na.value = "#8C0172"
    )+ # bilbao roma
    
    #geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
    #stat_cor(aes(color = Classification), label.x = 3)+
    
    scale_x_log10(  breaks = trans_breaks("log10", function(x) 10^x),  labels = trans_format("log10", math_format(10^.x))) +
   # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides="b")   +
    
    
    #coord_equal() +
    labs(#title=paste0("Size distribution of gated cells in ",d, " and ",w," weeks" )  , 
      title=paste0("CD3E gene expression vs Ki67 Expression " )  , 
      y="CD3E gene (Xenium)", 
      x="Ki67 (PCF)"
    )+
    theme(panel.background = element_rect(fill = "grey90", colour = "grey55", size = 0.5), # change 'colour' to black for informative axis
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=8),
          #axis.text.x = element_text(angle = 45, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          legend.position="right",
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )
  
  ggsave("Correlationplot_CD3EvsKi67.png",g, width = 25, height = 18, units = "cm",dpi = 800)
  
  
  
  
  
  
  
  
  
  ##---------- Adrian's Poster Plots -----------------  
  
  # for the cluster-wise analyis we need to pre-define groups of clusters.
  # we dont have an annotation yet, so these are the collapsed numbers:
  CD4.clusters <- c(211,208,201,212,203,205,202,210)
  CD8.clusters <- c(215,213,203,206,217,216,214,209,207,204)
  

  

  
  setwd(OutputDirectory)
  dir.create(paste0('Adrians_poster'))
  setwd(paste0('Adrians_poster'))
  
  i <- "C-3"  
temp.df <- subset( gated.cell.dat, GATEsubset==2 & TMA.core==i)


#lets order the temp.df:

temp.df$Collapsed_metacluster <-  factor(temp.df$Collapsed_metacluster,      # Reordering group factor levels
                                         levels = unique(c( CD4.clusters,CD8.clusters)) )


 make.colour.plot.adapted( dat=subset( gated.cell.dat, GATEsubset==2 & TMA.core==i) , x.axis = 'UMAP1', y.axis ='UMAP2', 
                                   col.axis ="Classification", 
                                   col.type = 'factor', 
                                   add.label = TRUE, #coloring comes via Annotated_metaclusters. Print these pls
                                   titlestring = paste0("UMAP by annotated clusters")  ,
                                   subtitlestring = subtitlestring,
                                   point.alpha = 1, # lets not overdraw the UMAP and reduce alpha of the geom_points
                                   dot.size = 1, # less alpha and smaller dots on the half-width-page graphs
                                   colours = "polychrome",
                                   polychromepalette = class_palette,
                                   repel.label.size = 4,
                                   label.label.repelforce = 30, # repulsion between overlapping text labels. default 30 to push them apart
                                   label.datapoint.pullforce = 0.4, # attraction between label and datapoint. default 0.8 to help release the label
                                   min.line.length.tolabel = 0, # minimal distance under which the segment line is not drawn anymore: 0 draw all, Inf turn off even if far
                                   legend.loc = "right",
                                   legend.text.size = 15, # default 18
                                   save.to.disk = F)
  
  #ggsave(paste0(gsub("[:/\\\\<>|?*]", "_", GATEsubsetNames[s]),"_UMAP.pdf"), g, width = 7, height = 6, dpi=400) # use the ratio to make the HM squared
  
  

  
 
g<- temp.df %>%
   # Get the counts
   group_by(Collapsed_metacluster,Classification) %>%
   summarise(count = n()) %>%
   
   ggplot(. , aes(fill=Classification, y=count, x=as.factor(Collapsed_metacluster))) + 
   geom_bar(position="fill", stat="identity")+
   scale_fill_manual( values = class_palette  ,
                      limits = force 
   )+
  labs(title= paste0("Core ", i, " - tissue: ", sample.meta[ sample.meta$TMA.core  == i,]$Tissue.type  )   ,
       subtitle= paste0("Cluster Frequencies in un-altered vs. metaplasia" ) ,
       x="Cluster", y="Frequency")+
 
  
  theme(panel.background = element_rect(fill = "grey99", colour = "grey99", linewidth = 0.5), # change 'colour' to black for informative axis
        panel.grid = element_line(color = "grey99"),
        axis.title.x=element_text(color="grey15", size=11),
        axis.title.y=element_text(color="grey15", size=11),
        axis.text=element_text(size=12),
        #axis.text.x = element_text(angle = 45, hjust=1),
        #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
        legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
        legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
        #legend.title=element_blank(),
        legend.position="left",
        plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
        plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
  )
  ggsave("Tcells_core_C3_by_tissue.pdf",g, width = 25, height = 10, units = "cm",dpi = 1200)
   
   
  
  
  temp.df.CD8 <- subset( temp.df, GATEsubset==2 & TMA.core==i & Collapsed_metacluster %in% CD8.clusters)
  
  
  temp.df.CD8$Collapsed_metacluster <-  factor(temp.df.CD8$Collapsed_metacluster,      # Reordering group factor levels
                                           levels = CD8.clusters )
  
g<-temp.df.CD8 %>%
    # Get the counts
    group_by(Classification,Collapsed_metacluster) %>%
  
    summarise(count = n()) %>%
   group_by(Classification)%>%
   mutate(percent= count/sum(count))%>%
    
    ggplot(. , aes(fill=Collapsed_metacluster, y=percent,  x=as.factor(Classification))) + 
   geom_bar(stat = "identity", width = 0.5)+

    scale_fill_manual( values = Phenograph_metacluster_palette  ,
                       limits = force 
    )+
    labs(title= paste0("Core ", i, " - tissue: ", sample.meta[ sample.meta$TMA.core  == i,]$Tissue.type  )   ,
         subtitle= paste0("CD8 cluster Frequencies in un-altered vs. metaplasia" ) ,
         x="Tissue", y="Cluster Frequency [per tissue]")+
    
    
    theme(panel.background = element_rect(fill = "grey99", colour = "grey99", linewidth = 0.5), # change 'colour' to black for informative axis
          panel.grid = element_line(color = "grey99"),
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=12),
          #axis.text.x = element_text(angle = 45, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          legend.position="left",
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )
 g
  ggsave("CD8Tcell_frequencies_core_C3_by_tissue.pdf",g, width = 25, height = 15, units = "cm",dpi = 1200)
  
  
  
  
  
  
  
  

   
g<- ggplot( subset( all.cells ,TMA.core==i) ,
        aes(x=Centroid.X.µm,y=Centroid.Y.µm,color=Classification,fill=Classification )) + #,color=as.factor(inside.lymphnode)
   geom_point(size=0.2, alpha=1) + 
   
   theme_minimal() +
   scale_y_reverse()+
   coord_fixed() + 
   # extract Tissue.type from sample.meta for the current i:
   labs(title= paste0("Core ", i, " - tissue: ", sample.meta[ sample.meta$TMA.core  == i,]$Tissue.type  )   
   )+
   
   
   scale_colour_manual(values = class_palette,
                       guide =  guide_legend(override.aes = list(size=6,alpha=1 )))+
   scale_fill_manual(values = class_palette)+
   
   theme(panel.background = element_rect(fill = "grey10", colour = "grey15", linewidth = 0.5), # change 'colour' to black for informative axis
         panel.grid = element_line(color = "grey15"),
         axis.title.x=element_text(color="grey15", size=11),
         axis.title.y=element_text(color="grey15", size=11),
         axis.text=element_text(size=8),
         #axis.text.x = element_text(angle = 45, hjust=1),
         #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
         legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
         legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
         #legend.title=element_blank(),
         legend.position="left",
         plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
         plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
   )
  
g<- ggsave(paste0("Core_", i, "_tissue_", sample.meta[ sample.meta$TMA.core  == i,]$Tissue.type ,".pdf" ),g, width = 30, height = 20, units = "cm",dpi = 1200)
  
  
  
  
 g<- ggplot( subset( all.cells ,TMA.core==i) ,
         aes(x=Centroid.X.µm,y=Centroid.Y.µm )) + #,color=as.factor(inside.lymphnode)
   geom_point(aes(color="#888888",fill="#888888"),
              size=0.2, alpha=1) + 
   
   
   
   geom_point(data= subset( gated.cell.dat, GATEsubset==2 & TMA.core=="C-3"),
              aes(color=as.factor(Collapsed_metacluster),fill=as.factor(Collapsed_metacluster)),
                size=2, alpha=1) + 
   
       scale_colour_manual(values = Phenograph_metacluster_palette,
                        guide =  guide_legend(override.aes = list(size=6,alpha=1 )))+
    scale_fill_manual(values = Phenograph_metacluster_palette)+
   
   theme_minimal() +
   scale_y_reverse()+
   coord_fixed() + 
   # extract Tissue.type from sample.meta for the current i:
   labs(title= paste0("Core ", i, " - tissue: ", sample.meta[ sample.meta$TMA.core  == i,]$Tissue.type  )   
   )+
   

   
   theme(panel.background = element_rect(fill = "grey10", colour = "grey15", linewidth = 0.5), # change 'colour' to black for informative axis
         panel.grid = element_line(color = "grey15"),
         axis.title.x=element_text(color="grey15", size=11),
         axis.title.y=element_text(color="grey15", size=11),
         axis.text=element_text(size=8),
         #axis.text.x = element_text(angle = 45, hjust=1),
         #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
         legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
         legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
         #legend.title=element_blank(),
         legend.position="left",
         plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
         plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
   )
  
 g<- ggsave(paste0("Core_", i, "_tissue_", sample.meta[ sample.meta$TMA.core  == i,]$Tissue.type ,"_G2_clusters.pdf" ),g, width = 30, height = 20, units = "cm",dpi = 1200)
  
  
  
  
  
  
  
  
  
  
  
  if(run.old.code==1){
  
  # Isis really just wants the CD45+ subset to be analyzed with the heatmap in Fig3 of the manuscript. 
  
  # subset gated.cell.dat for CD45 and re-create that heatmap there (its code that comes from v16 code and old runs)
  
  #Welche Immunzellen exprimieren am meisten PD-L1 ?
  # > Als Output möchten wir die reduzierte Heatmap (wie besprochen) > vermutlich ILC Zellen (evntl. NK Zellen) die PD-L1 + Zellen
  #> Bitte zusätzlich auch ein Heatmap, von nicht Immunzellen (CD45-negativen und CK18-negativen Zellen) > das bezieht sich auf eine spezifische Frage der Reviewer ob Stromazellen nicht evnt auch PD-L1 exprimieren
  
  
  ## ------ generate temp.df the immune subset by gating CD45+ -------
  
  #birgits idee: CD45+ vorgaten, clustern, dann schauen welche sind PDL1+  und dann alle immunzell marker zum clustern nehmen (ohne aSMA und whatever)
  temp.df <- subset(gated.cell.dat, outside.CD45.highpass == 0 )
  

  
  # how about we bring in the transfomed data and then z-score on the subset of the CD45+ gated cells?
  # rather than using cluster.cols with _t2_z
  # starting with CA19.9 finishning with Tim3_t2
 # cluster.cols <-names(gated.cell.dat)[c(101:108,111:112,114:116,118,121,123,125:128,130:134)+0]
  
  # actually, Fig3 does not have all markers in, just these here:
  
  cluster.cols <-names(gated.cell.dat)[c(108,112,116,106,105,121,133,123,128,107,118,130,115,102,125,114,104,131)+0]
  
  
  # remove everything after the last underscore in global.marker since it contains the _z
  #marker.order <- sub("_[^_]+$", "", global.marker.order)[sub("_[^_]+$", "", global.marker.order) %in% cluster.cols]
  
  cluster.cols <- paste0(cluster.cols, "_z")
  marker.order <- global.marker.order[global.marker.order %in% cluster.cols]
  # if you still want to bring in the z-scored data, you can do it like this:
  
  
  # ........ Phenograph ........    
  
  
  if("Phenograph_metacluster" %in% colnames( temp.df ) ) {
    
    message(  "\n-------   Phenograph   --------\nGonna remove ", max( temp.df$Phenograph_metacluster)," metaclusters from isolated gated.cell.dat before running"  )
    temp.df <-  subset( temp.df, select = -c(Phenograph_metacluster) )
  }else{
    message(  "\n-------   Phenograph   --------\nNo metaclusters present in isolated gated.cell.dat before running"  )
  }
  
  
  
  
  # create the Phenograph clusters:
  
  phenograph.k <- 26 # thats how much we have in fig3b
    set.seed(123)  
    Rphenograph_cluster <- Rphenoannoy(data=as.matrix(temp.df[, cluster.cols, with = FALSE]), k=phenograph.k)
    Rphenograph_cluster <- as.numeric(membership(Rphenograph_cluster[[2]]))
  
  

  temp.df[["Phenograph_metacluster"]] <- Rphenograph_cluster
  
  



  
  #NEW:.................. MAKE COMPLEX HEATMAP  ...............    
  
  exp <- temp.df[, lapply(.SD, "mean", na.rm = TRUE), by='Phenograph_metacluster',  
                 .SDcols = cluster.cols]
  
  ## new way to plot complex heatmap
  
  #z-normalize feature expression (iterate over columns with MARGIN=2)
  zscored.exp <- cbind(  exp[,1]  ,   apply(exp[,-1], scale, MARGIN = 2)  ) 
  
  # we gonna bring back only the signals into the matrix:
  zscored.exp.mat <- as.matrix(zscored.exp[,-1])
  
  # prepare column and row names:
  rownames(zscored.exp.mat) <- paste0(exp[[1]] )
  colnames(zscored.exp.mat) <- sub("_.*", "", colnames(zscored.exp.mat) )
  
  col_fun = colorRamp2(c(range(zscored.exp.mat)[1], 0, range(zscored.exp.mat)[2]), c("blue", "white", "red"))
  
  
  # this needs adaptation later on. take markers off the heat map when plotting the subgates. use with a good argumentation, this is a tad dangerous
  # when plotting the final colors on HM, lets take the two islet markers out, they disturb the interpretation:
  
  # zscored.exp.mat <- subset(zscored.exp.mat, select=-c(get("NKx6.1"),get("GLP.1R")))
  
  # also take them out of the current marker.order:
  #marker.order.noIsletMarkers <-  grep(paste0( c("NKx6.1", "GLP.1R" )   , collapse = "|") , marker.order, invert = TRUE, value = TRUE)
  marker.order.noIsletMarkers <-marker.order
  
  
  #     CREATE polychrome palette for the unannotated Phenograph MCs
  colRange <- unique(temp.df[["Phenograph_metacluster"]])
  colRange <- colRange[order(colRange)]
  #colRange <- as.character(colRange)
  set.seed(723451) # for reproducibility
  subset.Phenograph_metacluster_palette <- createPalette(length(colRange), c("#ff0000"), M=100000, prefix = "")
  # so, the cool thing is that once we annotated these clusters, we can just do the same game again, overrride the palette and plot again.
  names(subset.Phenograph_metacluster_palette) <- colRange  
  
  
  # before proceeeding we need to assemble this into a dataframe like annotated.Phenograph_metacluster_palette
  CD45_palette <- data.frame("Collapsed_metacluster"=names(subset.Phenograph_metacluster_palette), "color"=subset.Phenograph_metacluster_palette, row.names=NULL)
  
  # so the trick is to initiate column_ha first with the alphabetical order:
  column_ha = HeatmapAnnotation(cluster = rownames(zscored.exp.mat),
                                
                                annotation_legend_param = list(
                                  cluster = list(
                                    ncol = 2, 
                                    title = "Cluster",
                                    title_position = "topcenter",
                                    at = names(subset.Phenograph_metacluster_palette),
                                    grid_height = unit(0.02*length(marker.order.noIsletMarkers), "cm"),
                                    grid_width = unit(0.02*length(marker.order.noIsletMarkers), 'cm'),
                                    labels_gp = gpar(fontsize = 0.8*length(marker.order.noIsletMarkers)),
                                    title_gp = gpar(fontsize = 0.8*length(marker.order.noIsletMarkers))
                                  )
                                ),
                                
                                #col = list(cluster=hm_cols),
                                col = list(cluster=subset.Phenograph_metacluster_palette),
                                
                                na_col = "black", #white
                                show_annotation_name = FALSE
  )
  
  
  # then initiate the map
  hm <-   Heatmap( t( zscored.exp.mat ) ,
                   col = col_fun,
                   row_order = sub("_.*", "", marker.order.noIsletMarkers ) , # this is the only that needs stripping of transformed_rescaled
                   cluster_columns = fh,
                   top_annotation = column_ha,
                   rect_gp = gpar(col = "white", lwd = 1),
                   heatmap_legend_param = list(
                     title = "z-score",
                     direction = 'horizontal',
                     title_position = "topcenter",
                     legend_width = unit(0.25*length(marker.order.noIsletMarkers), "cm"),
                     grid_width = unit(0.02*length(marker.order.noIsletMarkers), 'cm'),
                     labels_gp = gpar(fontsize = 0.8*length(marker.order.noIsletMarkers)),
                     title_gp = gpar(fontsize = 0.8*length(marker.order.noIsletMarkers))
                   )
                   
  )
  
  #then bring back column_ha
  column_ha = HeatmapAnnotation(cluster = rownames(zscored.exp.mat),
                                
                                
                                # this is the Cluster legend. the problem is the "at", which takes over the alphabetical order of the subset.Phenograph_mc_palette
                                annotation_legend_param = list(
                                  cluster = list(
                                    ncol = 1, 
                                    title = "Cluster",
                                    #title_position = "topcenter",
                                    at = rownames(zscored.exp.mat)[column_order(hm)] ,#rownames(zscored.exp.mat), # names(subset.Phenograph_metacluster_palette),
                                    grid_height = unit(0.02*length(marker.order), "cm"),
                                    grid_width = unit(0.04*length(marker.order), 'cm'),
                                    labels_gp = gpar(fontsize = 0.8*length(marker.order)),
                                    title_gp = gpar(fontsize = 0.8*length(marker.order))
                                  )
                                ),
                                
                                #col = list(cluster=hm_cols),
                                col = list(cluster=subset.Phenograph_metacluster_palette),
                                
                                na_col = "black", #white
                                show_annotation_name = T
  )
  
  
  # then initiate the map again with the new column_ha order:
  hm <-   Heatmap( t( zscored.exp.mat ) ,
                   col = col_fun,
                   row_order = sub("_.*", "", marker.order ) , # this is the only that needs stripping of transformed_rescaled
                   cluster_columns = fh,
                   top_annotation = column_ha,
                   rect_gp = gpar(col = "white", lwd = 1),
                   heatmap_legend_param = list(
                     title = "z-score",
                     #direction = 'horizontal',
                     title_position = "lefttop-rot", #topcenter
                     legend_width = unit(0.25*length(marker.order), "cm"),
                     legend_height = unit(0.3*length(marker.order), "cm"),
                     grid_width = unit(0.02*length(marker.order), 'cm'),
                     labels_gp = gpar(fontsize = 0.6*length(marker.order)),
                     title_gp = gpar(fontsize = 0.8*length(marker.order))
                   )
  )
  
  
  
  
  
  
  
  
  
  
  
  
  
  # gonna store the hm away as grob and unite it with the UMAP that is not created at this point:
  hm.grob = grid.grabExpr(
    draw(hm,
         heatmap_legend_side = 'left',
         row_sub_title_side = 'left',
         padding = unit(c(2, 2, 2, 10), "mm"))
  )    
  
  
  
  
  # we gonna have to make a new UMAP:
  
  umap.coords <- umap(as.matrix(temp.df[, cluster.cols, with = FALSE]), 
                      n_neighbors = 20, # 16. shinyIMC looked cool with 16 - 0.03. I had here 10 - 0.01 to really separate the clusters
                      min_dist = 0.08,  #0.08
                      metric = 'euclidean',
                      verbose= TRUE
  ) #n_neighbors = 20, min_dist = 0.1, metric = 'euclidean'
  
  
  
  
  # we push the umap coordinates right back into the list of dataframes:
  
  #if you source uwot after, the umap.coords element looks different
  colnames(umap.coords)[1:2] <- paste0('UMAP', 1:2)
  temp.df[["UMAP1"]] <- umap.coords[, 1]
  temp.df[["UMAP2"]] <- umap.coords[, 2]
  
  
  subtitlestring <- paste0("Gated into all CD45+ cells (",  nrow(temp.df), " cells plotted)" )   
  
  
  g <-  make.colour.plot.adapted( temp.df , 'UMAP1', 'UMAP2', "Phenograph_metacluster", 'factor', 
                                  titlestring = paste0("UMAP by Phenograph_metaclusters")  ,
                                  subtitlestring = subtitlestring,
                                  point.alpha = 0.4, # lets not overdraw the UMAP and reduce alpha of the geom_points
                                  dot.size = 0.9, # less alpha and smaller dots on the half-width-page graphs
                                  colours = "polychrome",
                                  polychromepalette = subset.Phenograph_metacluster_palette,
                                  add.label = TRUE, #col.axis are the Phenograph/FlowSOM metacluster numbers. Print these pls
                                  legend.loc = "none",
                                  legend.text.size = 15, # default 18
                                  save.to.disk = F)
  
  
  
  
  hm.umap.fig <-  ggarrange(hm.grob, g, widths = c(6,4) )  # 3,2
  
  hm.umap.fig <-   annotate_figure(hm.umap.fig,
                                   top = text_grob(paste0("CD45+ subset of entire dataset"), color = "black", face = "bold", size = 17),
                                   bottom = text_grob(data.source.lyrics, color = "red",
                                                      hjust = 1, x = 1, face = "italic", size = 13),
                                   #left = text_grob("Figure arranged using ggpubr", color = "green", rot = 90),
                                   #right = text_grob(bquote("Superscript: ("*kg~NH[3]~ha^-1~yr^-1*")"), rot = 90),
                                   fig.lab = "Collapsed", fig.lab.face = "bold"
  ) 
  
  
  #"CD45+ subset of entire dataset_UMAP_Phenograph.png"
  
  ggsave( "CD45+ subset of entire dataset_UMAP_Phenograph.pdf", hm.umap.fig, width = 15, height = 8, dpi=600) # 4 Marker w11 h4
  
  
  
  
  ## ----- QC CD45 gate posision across entire data set gated.cell.dat: ----------------- 
  g <- 
    ggplot( data=gated.cell.dat ,
            aes( x= CD45_Sm152_t2,   y= CD3_Er170_t2 , color = as.factor(outside.CD45.highpass) ) ) + 
    
    geom_point(#aes(color = Annotated_metacluster),
      alpha=1, 
      size=0.8 #, 
      # pch='.'
    )+
    
    labs(title=paste0( "Entire dataset") , 
         subtitle="Gate postion for CD45+ cells", 
         # caption=test.stat.lyrics, 
         y="CD3", 
         x="CD45"
    )+
    
    scale_color_manual(name = "Gate flag", 
                       guide =  guide_legend(override.aes = list(size=4,alpha=1 )),
                       values = c("#40B0A6", "#E1BE6A"), #colorfriendly. not this: c("#308441", "#F6CD7A")
                       labels = c("in", "out") #0 blue in, 1 yellow out
    )+ #, 
    
    scale_x_log10(breaks = breaks, minor_breaks = minor_breaks, labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = breaks, minor_breaks = minor_breaks, labels = trans_format("log10", math_format(10^.x))) +
    
    
    stat_density_2d(#aes( fill = ..level..),        #fill = ..level..), 
      alpha = 0.2, 
      contour = TRUE,
      contour_var = "density",
      adjust=1.5,
      n = 220, # default 100, computationally intense, do not exceed that value too much..
      geom = "polygon", #geom: polygon or raster (equidistant lines or heatmap-ish picture)
      colour="white") +
    
    theme(panel.background = element_rect(fill = "grey95", colour = "grey55", linewidth = 0.5), # change 'colour' to black for informative axis
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=8),
          #axis.text.x = element_text(angle = 45, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          legend.position="right",
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )
  
  
  ggsave("CD45_gate_check.pdf",g, width = 30, height = 20, units = "cm",dpi = 600) 
  
  
  
  
  
  
## ----- re-group UICC metadata:  pool II and III and IVa and IVb into two groups -----------------
  
  temp.df <- temp.df %>%
    mutate(UICC_grouped = case_when(
      UICC_Staging %in% c("I")  ~ "I",
      UICC_Staging %in% c("II","III")  ~ "II & III",
      UICC_Staging %in% c("IVA","IVB")  ~ "IVA & IVB"
    ))
  
  
  temp.df <- temp.df %>%
    mutate(tissue.region = case_when(
      inside.TE == 1  ~ "tumor",
      inside.TE == 0  ~ "stroma"
    ))
  
  
  
##----------------------- setup frequency.data.table: DENSITIES and COMPOSITION calculation ROI-wise --------------------------------------  
  
  
  pb <- progress_bar$new(format = "[:bar] :percent [Cluster-wise calculations CD45+ subset | :eta]",
                         total = length(plot.ROIs)*3, # we run through the ROIs with three regions
                         show_after=0, #allows to call it right way 
                         current = "|",    # Current bar character
                         clear = FALSE) # make it persist
  pb$tick(0) # call in the progress bar without any progress, just to show it  
  
  
  # currently we do that only for subset 2....
  

  
  frequency.data.table <- data.frame()
  
  for(i in plot.ROIs ){
    
    temp.ROI <- subset(temp.df,  ROI== i) 
    
    # watch out, the triple gating pushes some ROIs outside this gate, such as "IM_20221222_OT6_s0_p4_r6_a6_ac":
    if(nrow(temp.ROI)>0){
      
      ROI.table <-  as.data.frame(table(temp.ROI$Phenograph_metacluster))
      names(ROI.table) <- c("MC", "Cells")
      
      
      # the absent ones we usually feed via a all.annot.MC object, but we dont have that in this project, so lets use annotated.Phenograph_metacluster_palette
      
      temp.absent <- data.frame( "MC" = setdiff( unique(temp.df$Phenograph_metacluster), ROI.table$MC),
                                 "Cells"= rep(0, length(   setdiff( unique(temp.df$Phenograph_metacluster), ROI.table$MC) )      ) 
      )
      
      
      # at that point, MC is factor, and rbinding the numeric MC of the absent cluster would kill the absent cluster names
      ROI.table$MC <- as.character(ROI.table$MC)
      temp.absent$MC <- as.character(temp.absent$MC)
      
      ROI.table <- rbind(ROI.table,temp.absent)
      rm(temp.absent)
      
      
      for(a in c("ROI","tumor","stroma") ){
        
        
        if(a == "ROI"){
          
          temp.totalSegmentationArea <- sum(subset(all.cells, ROI== i  )$Area)
          total.CD45 <- nrow(subset(temp.ROI, outside.CD45.highpass == 0  )) 
          
          # you need to subset for the current area:
          temp.table <- ROI.table
          
          
          
        }else if(a == "tumor"){
          
          temp.totalSegmentationArea <- sum(subset(all.cells, ROI== i & inside.TE == 1  )$Area)
          total.CD45 <- nrow(subset(temp.ROI, outside.CD45.highpass == 0 & inside.TE == 1  ))
          
          
          # you need to subset for the current area:
          temp.table <- data.frame(
            "MC" = ROI.table$MC,
            "Cells"= sapply(ROI.table$MC, function(x) sum(temp.ROI$Phenograph_metacluster == x & temp.ROI$inside.TE == 1  , na.rm=T ))
          )
          
          
          
        }else if(a == "stroma"){
          
          temp.totalSegmentationArea <- sum(subset(all.cells, ROI== i & inside.TE == 0  )$Area)
          total.CD45 <- nrow(subset(temp.ROI, outside.CD45.highpass == 0 & inside.TE == 0  ))
          
          # you need to subset for the current area:
          temp.table <- data.frame(
            "MC" = ROI.table$MC,
            "Cells"= sapply(ROI.table$MC, function(x) sum(temp.ROI$Phenograph_metacluster == x & temp.ROI$inside.TE == 0  , na.rm=T ))
          )
          
        }
        
        
        # I know this is seems like a stupid way to build up that frequency.data.table. But it allows to easily add more measurments while keeping the three areas separated:
        temp.table <- cbind(temp.table,  data.frame(
          
          "ROI"= i, 
          "Ref_Area"= a, 
          # patient-related meta
          "Patient_ID"= unique(temp.ROI$Patient_ID),  
          "age"= unique(temp.ROI$Age_at_diagnose_surgery_date),  
          # tumor-related meta
          "UICC_Stage"= unique(temp.ROI$UICC_Staging),
          "UICC_Stage_grouped"= unique(temp.ROI$UICC_grouped),
          # survival-related meta
          "monitored_survival_months" = unique(temp.ROI$Overall_Survival_months), # watch out, even the un-dead have here different months...
          "longtime_survivor"= unique(temp.ROI$Survival_more36months),
          "died" = unique(temp.ROI$Dead),
          
          
          "density_mm2" = 1e+06 * temp.table$Cells/temp.totalSegmentationArea,
          
          "percentofCD45" = 100 * temp.table$Cells/total.CD45 # be aware that we run through that cbind three times, and "total.CD4" here is the sum of the local area summed up!
          
        )
        )
        
        
        
        frequency.data.table <- rbind(frequency.data.table, temp.table)
        pb$tick()
      }# loop with a trhough all three areas and assemble the list
    }#protect from processing ROIs that do not have any cell in the current gate
  }# run through all ROIs
  
  

  # mutate longtime_survivor in longtime_survivor to <36months if 0, and >=36months if 1
  frequency.data.table <- 
    frequency.data.table %>%
    mutate(longtime_survivor_string = ifelse(longtime_survivor == 0, "<36months", ">36months"  ))
  
  
  
  # prepare the strip for the upcoming plotting:
  facet.labs <- unique(temp.df$Phenograph_metacluster)
  names(facet.labs) <- unique(temp.df$Phenograph_metacluster) 
  
  
  strip <- strip_themed(background_x = elem_list_rect(fill = CD45_palette$color))
  
  
  
  
  ## ------------------------ CD45_densities_byTumorSTAGE ------------------------
  
  # we bring here the new code from analyzing Laias IMC data in a statistaclly more clean way:
  # define a global statistic to keep in all plots:
  test.method <- "wilcox.test"# "wilcox.test", "t.test"
  adj.p <- "BH"# "BH" #"none"
  test.stat.lyrics <- paste0(test.method,", p.adj: ",adj.p )
  
  
  
  ## ----------- prepare normalizations across different conditions ---------- 
  
 
  
  
  # running variable i, y coordinate, split by conditon:
  normalizations.to.do <- list( c("density_mm2", "UICC_Stage_grouped"), 
                                c("density_mm2", "longtime_survivor_string"),
                                c("density_mm2", "Ref_Area"),
                                
                                
                                c("percentofCD45", "UICC_Stage_grouped"), 
                                c("percentofCD45", "longtime_survivor_string"),
                                c("percentofCD45", "Ref_Area")
                                
                                
                                )
  
  color.palette.to.use <- list(  tumorStage_grouped_palette,
                             longtimesurvival_palette_string,
                             tissue_palette,
                             
                             
                             tumorStage_grouped_palette,
                             longtimesurvival_palette_string,
                             tissue_palette
                             
                             )
  
  
  titlestrings <- c(  "CD45+ cluster densities by UICC stage",
                      "CD45+ cluster densities by long-time survival",
                      "CD45+ cluster densities by tumor area",
                      
                      "CD45+ cluster frequencies by UICC stage",
                      "CD45+ cluster frequencies by long-time survival",
                      "CD45+ cluster frequencies by tumor area"
  )
  
  yaxisstrings <- c( "Cluster density [cells * mm^-2]",
                     "Cluster density [cells * mm^-2]",
                     "Cluster density [cells * mm^-2 local area]",
                     
                     "Cluster frequency [% of CD45+ cells]",
                     "Cluster frequency [% of CD45+ cells]",
                     "Cluster frequency [% of CD45+ cells]"
                     
                     
  )
  
  # will be pasted to .pdf
  filestrings <- c("CD45_densities_byTumorSTAGE",
                   "CD45_densities_bySurvival",
                   "CD45_densities_byTumorArea",
                   
                   "CD45_frequencies_byTumorSTAGE",
                   "CD45_frequencies_bySurvival",
                   "CD45_frequencies_byTumorArea"
  )
  
  
  
## ------ for plotting, bring some order into the x coordinates -------
  
  
  
  for(x in c(1:6)){
    
    i <- normalizations.to.do[[x]][1] # x coordinates
    j <- normalizations.to.do[[x]][2] # y coordinates, not absolutely needed down in the code
    
    temp.palette <- color.palette.to.use[[x]]
    
    titlestring <- titlestrings[x]
    filestring <- filestrings[x]
    yaxisstring <- yaxisstrings[x]
    
    
    
    
    ## STATISTICAL TESTING start, for every cluster separately, and for every comparison of interest:
    stats_table <- data.frame()
    
    for(m in unique(temp.df$Phenograph_metacluster) ){
      
      # in this loop we need to create the formula string beforehand, since compare_means(formula= ... ) is stupid and doesnt recgnize the running variable i, 
      # instead it searches for a column called i, or whatever you try. so lets it like this:
      
      formula <- paste(i, "~", j)
      
      
      
      
      compp <- compare_means(as.formula(formula),
                             data = subset(frequency.data.table, MC %in% m),
                             method = test.method,
                             paired = F,
                             p.adjust.method = adj.p
      )   
      
      
      # compp <- compp[compp$p.adj<0.05,]
      compp <- compp[compp$p<0.05,]
      
      if(nrow(compp)>0){
        compp <- cbind(compp,"MC"=m)
        compp <- cbind(compp,"y"= seq(2000, length=nrow(compp), by=70)  )
        stats_table <- rbind(stats_table,compp)
      }
    }
    
    # for geom_significant to work, we need to move some stuff around in stats_table
    if(nrow(stats_table) >0){
      stats_table <- stats_table %>% 
        dplyr::rename(
          measure = .y.,
          xmin = group1,
          xmax = group2
        )
      
      
      
      #---------- here is the snipplet to filter only for interesting comparisons. not needed in this case here ------------ 
      #  filter.combinations <- filter.combinations.invert <- as.data.frame(outer(unique(frequenciesofinterestingClusters$Annotated_metacluster)
      #                                                                           ,c( "Healthy","Tumor"), paste, sep="_"))
      #  
      #  names(filter.combinations) <- c("xmin","xmax")
      #  # now, there is a chance that in the stats_table, H and T are inverted, so lets put that into filter.combinations as well:
      #  names(filter.combinations.invert) <- c("xmax","xmin")
      #  # then bring the inverted back into our filter combos:
      #  filter.combinations <-rbind(filter.combinations,filter.combinations.invert)
      #  rm(filter.combinations.invert)
      #.-......................................................................................................................  
      
      
      # we gonna re-set the y positions after filtering:
      stats_table_filtered_yAdjust <- data.frame()
      
      for(m in unique( stats_table$MC ) ){
        
        stats_temp <- subset(stats_table, MC %in% m)
        
        #for(m in unique(stats_table$Annotated_metacluster)){
        #  stats_temp <- subset(stats_table, Annotated_metacluster %in% m)
        
        
        # max.data.value <- max( subset(interstingClusterDensity.perROI.data.table, MC %in% m)[,i] , na.rm=T)
        max.data.value <- max( as.numeric(quantile( subset(frequency.data.table,MC %in% m)[,i], probs = 0.99, na.rm=T )))
        
        
        stats_temp$y <- seq(max.data.value , 
                            length=nrow(stats_temp), by= max( subset(frequency.data.table, MC %in% m)[,i] , na.rm=T) * 0.2 )  
        
        
        # watch out, Inf is passing through here:
        #bracketoffset <- 0.4*max( interstingClusterDensity.perROI.data.table[, i], na.rm=T )
        # better use is.finite and pass the vector like that:
        #bracketoffset <- 0.4*max(  interstingClusterDensity.perROI.data.table[, i][is.finite(interstingClusterDensity.perROI.data.table[, i])]  )
        #print(bracketoffset)                                          
        
        stats_table_filtered_yAdjust <- rbind(stats_table_filtered_yAdjust, stats_temp)
        
        
      }# end run with m along the annotated mc in the stats_table
      
      
      # not needed at this point:
      # prepare the columns we use on X axis to have the proper order:
      #interstingClusterDensity.permouse.data.table$Cond <- factor(interstingClusterDensity.permouse.data.table$Cond,  levels = all.condtions )
      #interstingClusterDensity.permouse.data.table$MC <- factor(interstingClusterDensity.permouse.data.table$MC,  levels = cluster.order )
      stats_table_filtered_yAdjust$MC <- factor(stats_table_filtered_yAdjust$MC,  levels = unique(temp.df$Phenograph_metacluster) )
      
      
      #only process the statistics if you got some in the data in stats_table 
    }else{
      # if the stats_table remained empty since everything was NS, we still need to reset the y-adjust since we absolutely need to protect from drawing stars from the previous run:
      stats_table_filtered_yAdjust <- stats_table # both are empty now
    }
    
    
    # do we need to level that guy beforehand?
    frequency.data.table$MC <- factor(frequency.data.table$MC,  levels = unique(temp.df$Phenograph_metacluster) )
    # yes we do have to
    
    
    ## plotting the facet start:
    
    g<- 
      ggplot(subset(frequency.data.table, MC %in% unique(temp.df$Phenograph_metacluster))   , 
             aes(x = as.factor( subset(frequency.data.table, MC %in% unique(temp.df$Phenograph_metacluster))[,j] ), 
                 y = subset(frequency.data.table, MC %in% unique(temp.df$Phenograph_metacluster))[,i], 
                 color=as.factor( subset(frequency.data.table, MC %in% unique(temp.df$Phenograph_metacluster))[,j] )
             )
      )+
      facet_wrap2(.~MC,
                  strip = strip,
                  labeller = labeller(source = facet.labs),
                  nrow = 5, scales="free_y"
      )+
      geom_violin(trim=T)+
  geom_quasirandom(alpha=1,size=2) +
      scale_colour_manual( values = temp.palette , limits = force  )+
      
      # significance brackets are cut off, expand axis:
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
      
      
      labs(title=titlestring , 
           #subtitle=subtitlestring, 
           caption=test.stat.lyrics, 
           y=yaxisstring, 
           x=""
      )+
      
      theme(panel.background = element_rect(fill = "grey90", colour = "grey55", linewidth = 0.5), # change 'colour' to black for informative axis
            axis.title.x=element_text(color="grey15", size=11),
            axis.title.y=element_text(color="grey15", size=11),
            axis.text=element_text(size=8),
            #axis.text.x = element_text(angle = 45, hjust=1),
            #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
            legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
            legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
            #legend.title=element_blank(),
            legend.position="none",
            plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
            plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
      )
    # only add brackets if there are significances above threshold
    if(nrow(stats_table_filtered_yAdjust) >0){
      g <- g+
        geom_signif(
          data = stats_table_filtered_yAdjust,
          aes(xmin = as.factor(xmin), 
              xmax = as.factor(xmax), 
              annotations = p.signif,#p.format,#p.signif, 
              y_position = y),
          textsize = 6, 
          tip_length=0.03,
          color= "black",
          vjust = 0.5,
          manual = TRUE
        )
    }#only process the statistics if you got some in the data
    
    plot(g)
    ggsave(paste0(filestring,".pdf"), g, width = 20, height = 10, dpi=600)
    
    
    #table(temp.df$Phenograph_metacluster)/nrow(temp.df)*100 
    
    
    
    
  }#end run with x through all normalizations
  
  
  
  
  
  
  
  
  ## ------------------------ CD45_composion_byTumorSTAGE ------------------------
  
  i <- "percentofCD45"
  
  stats_table <- data.frame()
  
  for(m in unique(temp.df$Phenograph_metacluster)){
    
    # in this loop we need to create the formula string beforehand, since compare_means(formula= ... ) is stupid and doesnt recgnize the running variable i, 
    # instead it searches for a column called i, or whatever you try. so lets it like this:
    
    formula <- paste(i, "~ UICC_Stage_grouped")
    
    
    
    
    compp <- compare_means(as.formula(formula),
                           data = subset(frequency.data.table, MC %in% m),
                           method = test.method,
                           paired = F,
                           p.adjust.method = adj.p
    )   
    
    
    # compp <- compp[compp$p.adj<0.05,]
    compp <- compp[compp$p<0.05,]
    
    if(nrow(compp)>0){
      compp <- cbind(compp,"MC"=m)
      compp <- cbind(compp,"y"= seq(2000, length=nrow(compp), by=70)  )
      stats_table <- rbind(stats_table,compp)
    }
  }
  
  # for geom_significant to work, we need to move some stuff around in stats_table
  if(nrow(stats_table) >0){
    stats_table <- stats_table %>% 
      dplyr::rename(
        measure = .y.,
        xmin = group1,
        xmax = group2
      )
    
    
    
    #---------- here is the snipplet to filter only for interesting comparisons. not needed in this case here ------------ 
    #  filter.combinations <- filter.combinations.invert <- as.data.frame(outer(unique(frequenciesofinterestingClusters$Annotated_metacluster)
    #                                                                           ,c( "Healthy","Tumor"), paste, sep="_"))
    #  
    #  names(filter.combinations) <- c("xmin","xmax")
    #  # now, there is a chance that in the stats_table, H and T are inverted, so lets put that into filter.combinations as well:
    #  names(filter.combinations.invert) <- c("xmax","xmin")
    #  # then bring the inverted back into our filter combos:
    #  filter.combinations <-rbind(filter.combinations,filter.combinations.invert)
    #  rm(filter.combinations.invert)
    #.-......................................................................................................................  
    
    
    # we gonna re-set the y positions after filtering:
    stats_table_filtered_yAdjust <- data.frame()
    
    for(m in unique( stats_table$MC ) ){
      
      stats_temp <- subset(stats_table, MC %in% m)
      
      #for(m in unique(stats_table$Annotated_metacluster)){
      #  stats_temp <- subset(stats_table, Annotated_metacluster %in% m)
      
      
      # max.data.value <- max( subset(interstingClusterDensity.perROI.data.table, MC %in% m)[,i] , na.rm=T)
      max.data.value <- max( as.numeric(quantile( subset(frequency.data.table,MC %in% m)[,i], probs = 0.99, na.rm=T )))
      
      
      stats_temp$y <- seq(max.data.value , 
                          length=nrow(stats_temp), by= max( subset(frequency.data.table, MC %in% m)[,i] , na.rm=T) * 0.2 )  
      
      
      # watch out, Inf is passing through here:
      #bracketoffset <- 0.4*max( interstingClusterDensity.perROI.data.table[, i], na.rm=T )
      # better use is.finite and pass the vector like that:
      #bracketoffset <- 0.4*max(  interstingClusterDensity.perROI.data.table[, i][is.finite(interstingClusterDensity.perROI.data.table[, i])]  )
      #print(bracketoffset)                                          
      
      stats_table_filtered_yAdjust <- rbind(stats_table_filtered_yAdjust, stats_temp)
      
      
    }# end run with m along the annotated mc in the stats_table
    
    
    # not needed at this point:
    # prepare the columns we use on X axis to have the proper order:
    #interstingClusterDensity.permouse.data.table$Cond <- factor(interstingClusterDensity.permouse.data.table$Cond,  levels = all.condtions )
    #interstingClusterDensity.permouse.data.table$MC <- factor(interstingClusterDensity.permouse.data.table$MC,  levels = cluster.order )
    stats_table_filtered_yAdjust$MC <- factor(stats_table_filtered_yAdjust$MC,  levels = unique(temp.df$Phenograph_metacluster) )
    
    
    #only process the statistics if you got some in the data in stats_table 
  }else{
    # if the stats_table remained empty since everything was NS, we still need to reset the y-adjust since we absolutely need to protect from drawing stars from the previous run:
    stats_table_filtered_yAdjust <- stats_table # both are empty now
  }
  
  
  
  g<- ggplot(subset(frequency.data.table, MC %in% unique(temp.df$Phenograph_metacluster) )   , 
             aes(x = as.factor(UICC_Stage_grouped), 
                 y = percentofCD45,
                 color=as.factor(UICC_Stage_grouped)
             )
  )+
    facet_wrap2(.~MC,
                strip = strip,
                labeller = labeller(source = facet.labs),
                nrow = 5, scales="free_y"
    )+
    geom_violin(trim=T)+
geom_quasirandom(alpha=1,size=2) +
    scale_colour_manual( values = tumorStage_grouped_palette, limits = force  )+
    
    # significance brackets are cut off, expand axis:
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
    
    
    labs(title="CD45+ cluster composition on ROIs" , 
         subtitle="Split by Tumor Stage Groups", 
         caption=test.stat.lyrics, 
         y="Cluster frequency [% of all CD45+ (threshold)]", 
         x=""
    )+
    
    theme(panel.background = element_rect(fill = "grey90", colour = "grey55", linewidth = 0.5), # change 'colour' to black for informative axis
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=8),
          #axis.text.x = element_text(angle = 45, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          legend.position="none",
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )
  # only add brackets if there are significances above threshold
  if(nrow(stats_table_filtered_yAdjust) >0){
    g <- g+
      geom_signif(
        data = stats_table_filtered_yAdjust,
        aes(xmin = as.factor(xmin), 
            xmax = as.factor(xmax), 
            annotations = p.signif,#p.format,#p.signif, 
            y_position = y),
        textsize = 6, 
        tip_length=0.03,
        color= "black",
        vjust = 0.5,
        manual = TRUE
      )
  }#only process the statistics if you got some in the data
  
  plot(g)
  ggsave("CD45_frequencies_byTumorSTAGE.pdf", g, width = 10, height = 6, dpi=600)
  
  
  
  
  
  

  ##---------- Cluster-wise analysis subset 2 ----------------- 
  s<-2
  
  setwd(OutputDirectory)
  dir.create('Plots')
  setwd('Plots')
  dir.create(paste0('Cluster-wise analysis ',GATEsubsetNames[s]))
  setwd(paste0('Cluster-wise analysis ',GATEsubsetNames[s]))


  
  
  
  # the density calc is ready, but the composition rel. to CD4+ needs a CD4+ flag since Birgit was not convinced to pool the CD4+ clusters. 
  # the reason for her being hesitant is that we got some rather big CD4CD8 doublepos clusters (two clusters, actually) so:
  # we used to set it like so:   
  #gated.cell.dat[, paste0("outside.CD4.highpass") ] <-  ifelse(  (gated.cell.dat[[ "CD4_Gd156_t2" ]] >  0.21 & gated.cell.dat[[ "CD3_Er170_t2" ]] >  0.1   )   , 0, 1 )
  # but we got a fancy CD3 gate already, so let use it:
  #### SET CD3-CD4 GATE ######  
  #depreciated, we use a multisegmented line:
  #gated.cell.dat[, paste0("outside.CD4.highpass") ] <-  ifelse(  (gated.cell.dat[[ "CD4_Gd156_t2" ]] >  0.21 & gated.cell.dat$outside.CD3.highpass == 0   )   , 0, 1 )
  

  
  
  # check:
  s=2
 g <- 
    ggplot( data=subset(gated.cell.dat, GATEsubset == s) ,
               aes( x= CD8_Dy162,   y= CD4_Gd156 , color = as.factor(outside.CD8.lowpass) ) ) + 
    
    geom_point(#aes(color = Annotated_metacluster),
      alpha=1, 
      size=0.8 #, 
      # pch='.'
    )+
      
      labs(title=paste0( "Gated subset ",GATEsubsetNames[s]) , 
           subtitle="Manually set CD4 gate to normalize to total CD4+ cells", 
           # caption=test.stat.lyrics, 
           y="CD4", 
           x="CD8"
      )+
    
    scale_color_manual(name = "Gate flag", 
                       guide =  guide_legend(override.aes = list(size=4,alpha=1 )),
                       values = c("#40B0A6", "#E1BE6A"), #colorfriendly. not this: c("#308441", "#F6CD7A")
                       labels = c("in", "out") #0 blue in, 1 yellow out
    )+ #, 
    
    scale_x_log10(breaks = breaks, minor_breaks = minor_breaks, labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = breaks, minor_breaks = minor_breaks, labels = trans_format("log10", math_format(10^.x))) +
    
    
    stat_density_2d(#aes( fill = ..level..),        #fill = ..level..), 
      alpha = 0.2, 
      contour = TRUE,
      contour_var = "density",
      adjust=1.0,
      n = 220, # default 100, computationally intense, do not exceed that value too much..
      geom = "polygon", #geom: polygon or raster (equidistant lines or heatmap-ish picture)
      colour="white") +
    
    theme(panel.background = element_rect(fill = "grey95", colour = "grey55", linewidth = 0.5), # change 'colour' to black for informative axis
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=8),
          #axis.text.x = element_text(angle = 45, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          legend.position="right",
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )
  
  
  ggsave("CD4_gate_check.jpg",g, width = 30, height = 20, units = "cm",dpi = 600) 
  
  
  
  
  # for grouping the stages, we need another column to pool II and III and IVa and IVb:
  
  gated.cell.dat <- gated.cell.dat %>%
    mutate(UICC_grouped = case_when(
      UICC_Staging %in% c("I")  ~ "I",
      UICC_Staging %in% c("II","III")  ~ "II & III",
      UICC_Staging %in% c("IVA","IVB")  ~ "IVA & IVB"
    ))
  

  gated.cell.dat <- gated.cell.dat %>%
    mutate(tissue.region = case_when(
      inside.TE == 1  ~ "tumor",
      inside.TE == 0  ~ "stroma"
    ))
  
  # > unique(subset(temp.df, UICC_Staging == "I")$UICC_grouped)
  #  [1] 1
  #  > unique(subset(temp.df, UICC_Staging == "II")$UICC_grouped)
  #  [1] 2
  #  > unique(subset(temp.df, UICC_Staging == "III")$UICC_grouped)
  #  [1] 2
  #  > unique(subset(temp.df, UICC_Staging == "IVA")$UICC_grouped)
  #  [1] 3
  #  > unique(subset(temp.df, UICC_Staging == "IVB")$UICC_grouped)
  #  [1] 3
  
  
  
  
  
  #----------------------- DENSITIES and COMPOSITION calculation ROI-wise --------------------------------------  
  
  
  pb <- progress_bar$new(format = "[:bar] :percent [Cluster-wise calculations | :eta]",
                         total = length(plot.ROIs)*3, # we run through the ROIs with three regions
                         show_after=0, #allows to call it right way 
                         current = "|",    # Current bar character
                         clear = FALSE) # make it persist
  pb$tick(0) # call in the progress bar without any progress, just to show it  
  
  
  # currently we do that only for subset 2....
  
  temp.df <- subset(gated.cell.dat, GATEsubset == s)
  
  frequency.data.table <- data.frame()
  
  for(i in plot.ROIs ){
    
    temp.ROI <- subset(gated.cell.dat, GATEsubset == s & ROI== i) 
    
    # watch out, the triple gating pushes some ROIs outside this gate, such as "IM_20221222_OT6_s0_p4_r6_a6_ac":
    if(nrow(temp.ROI)>0){
    
    ROI.table <-  as.data.frame(table(temp.ROI$Annotated_metacluster))
    names(ROI.table) <- c("MC", "Cells")

    
    # the absent ones we usually feed via a all.annot.MC object, but we dont have that in this project, so lets use annotated.Phenograph_metacluster_palette

    temp.absent <- data.frame( "MC" = setdiff(subset(annotated.Phenograph_metacluster_palette, Collapsed_metacluster > 100*(s) & Collapsed_metacluster < 100*(s+1))$Annotated_metacluster, ROI.table$MC),
                               "Cells"= rep(0, length(   setdiff(subset(annotated.Phenograph_metacluster_palette, Collapsed_metacluster > 100*(s) & Collapsed_metacluster < 100*(s+1))$Annotated_metacluster, ROI.table$MC) )      ) 
    )
    
    
    # at that point, MC is factor, and rbinding the numeric MC of the absent cluster would kill the absent cluster names
    ROI.table$MC <- as.character(ROI.table$MC)
    temp.absent$MC <- as.character(temp.absent$MC)
    
    ROI.table <- rbind(ROI.table,temp.absent)
    rm(temp.absent)
    
    
    for(a in c("ROI","tumor","stroma") ){
      
      
      if(a == "ROI"){
        
        temp.totalSegmentationArea <- sum(subset(all.cells, ROI== i  )$Area)
        total.CD4 <- nrow(subset(temp.ROI, outside.CD8.lowpass == 0  )) # CD4 highpass 0 equals to CD8 lowpass 0
        
       # you need to subset for the current area:
        temp.table <- ROI.table
        

        
      }else if(a == "tumor"){
        
        temp.totalSegmentationArea <- sum(subset(all.cells, ROI== i & inside.TE == 1  )$Area)
        total.CD4 <- nrow(subset(temp.ROI, outside.CD8.lowpass == 0  & inside.TE == 1  ))
        
        
        # you need to subset for the current area:
        temp.table <- data.frame(
                                 "MC" = ROI.table$MC,
                                 "Cells"= sapply(ROI.table$MC, function(x) sum(temp.ROI$Annotated_metacluster == x & temp.ROI$inside.TE == 1  , na.rm=T ))
                                 )
        
       
        
      }else if(a == "stroma"){
        
        temp.totalSegmentationArea <- sum(subset(all.cells, ROI== i & inside.TE == 0  )$Area)
        total.CD4 <- nrow(subset(temp.ROI, outside.CD8.lowpass == 0  & inside.TE == 0  ))
        
        # you need to subset for the current area:
        temp.table <- data.frame(
          "MC" = ROI.table$MC,
          "Cells"= sapply(ROI.table$MC, function(x) sum(temp.ROI$Annotated_metacluster == x & temp.ROI$inside.TE == 0  , na.rm=T ))
        )
        
      }
      
      
      # I know this is seems like a stupid way to build up that frequency.data.table. But it allows to easily add more measurments while keeping the three areas separated:
      temp.table <- cbind(temp.table,  data.frame(
        
        "ROI"= i, 
        "Ref_Area"= a, 
        # patient-related meta
        "Patient_ID"= unique(temp.ROI$Patient_ID),  
        "age"= unique(temp.ROI$Age_at_diagnose_surgery_date),  
        # tumor-related meta
        "UICC_Stage"= unique(temp.ROI$UICC_Staging),
        "UICC_Stage_grouped"= unique(temp.ROI$UICC_grouped),
        # survival-related meta
        "monitored_survival_months" = unique(temp.ROI$Overall_Survival_months), # watch out, even the un-dead have here different months...
        "longtime_survivor"= unique(temp.ROI$Survival_more36months),
        "died" = unique(temp.ROI$Dead),
        
        
         "density_mm2" = 1e+06 * temp.table$Cells/temp.totalSegmentationArea,

        "percentofCD4" = 100 * temp.table$Cells/total.CD4 # be aware that we run through that cbind three times, and "total.CD4" here is the sum of the local area summed up!

      )
      )
      
      
      
      frequency.data.table <- rbind(frequency.data.table, temp.table)
      pb$tick()
    }# loop with a trhough all three areas and assemble the list
    }#protect from processing ROIs that do not have any cell in the current gate
  }# run through all ROIs
  
  
  
  

  
  
  
  
  
  
  
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  

  
  
  
  # predictor variables to split the cluster properties:
  # tumor stage (gruppe 1: stage I, gruppe 2: stages II&III, gruppe 3: stage IVa / IVb
  #longterm survivors >36monaten
  #stroma vs tumor lokalisation
  
  
  
  
  normalizations.to.do <- c("density_mm2"#,  
                           # "percentofCD4.MOUSE.rimEXTislet",
                           # "percentofCD4.MOUSE.outsiderimEXTislet"
                            )
  
  
  
  titlestrings <- c(  paste0( "Gated subset ",GATEsubsetNames[s])#,
                      #paste0(islet.rim.width,"px-extended islet areas. Mouse-wise summed up"),
                     # paste0("Outside ",islet.rim.width,"px-extended islet areas. Mouse-wise summed up")
  )
  
  
  
  if(run.experimental.code==1){
  # the plotting will be homogenous by pre-defining the plotting engine as such:
  # since we wanna call ggplot multiple times over multiple areas on the ROIs, we need no for-loop, that fails, but some lapply and purr
  library(purrr)
  ## it gets even sexier than a for loop and call the plotting engine over column element as y:
  make.violin.multiplot.v2 = function (column, legendtitle, yaxisstring, data) {
    

ggplot( data, 
        aes_string(x = data$Cond, 
                    y = data[[ column ]], 
                    colour=data$MC 
        )
)+
  facet_wrap2(.~MC, 
              strip = strip,
              labeller = labeller(source = facet.labs),
              nrow = 1,
              scales="free_y")+
  geom_violin(trim=T)+
  geom_quasirandom(alpha=1,size=2) +
  scale_colour_manual( values = subset.Phenograph_metacluster_palette, limits = force  )+
  
  labs(title=column , 
        #subtitle=subtitlestring, 
        # caption=test.stat.lyrics, 
        y=yaxisstring, 
        x="",
        color=legendtitle
  )+
  theme(panel.background = element_rect(fill = "grey90", colour = "grey55", linewidth = 0.5), # change 'colour' to black for informative axis
        axis.title.x=element_text(color="grey15", size=11),
        axis.title.y=element_text(color="grey15", size=11),
        axis.text=element_text(size=8),
        axis.text.x = element_text(angle = 45, hjust=1),
        #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
        legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
        legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
        #legend.title=element_blank(),
        legend.position="right",
        plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
        plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
  )

}


  
  
  u <- lapply(normalizations.to.do ,make.violin.multiplot.v2, 
              yaxisstring="Frequency [% of CD4+]", 
              legendtitle="CD4+ cluster",  
              data=frequency.data.table) 
  
  # this does not allow you to push the titles from the list - yet. 
  # lets purr them into the list, by element:
  u <- map2(u, titlestrings, ~ .x + ggtitle(.y))
  }#experimental code of purring the plots in a loop
  
  

  # only right before plotting we kick the CD8s out for the composition CD4 of CD4
  
  # prepare the strip for the upcoming plotting:
  facet.labs <- CD4.clusters
  names(facet.labs) <- CD4.clusters 
  
  
  strip <- strip_themed(background_x = elem_list_rect(fill = subset(annotated.Phenograph_metacluster_palette, Annotated_metacluster %in% CD4.clusters )$color))
  
  
## ------------------------ CD4_densities_byTumorSTAGE ------------------------
  
  # we bring here the new code from analyzing Laias IMC data in a statistaclly more clean way:
  # define a global statistic to keep in all plots:
  test.method <- "wilcox.test"# "wilcox.test", "t.test"
  adj.p <- "BH"# "BH" #"none"
  test.stat.lyrics <- paste0(test.method,", p.adj: ",adj.p )
  
  
  # for every cluster we need to do that separately:
  
  i <- "density_mm2"
  
  stats_table <- data.frame()
  
  for(m in CD4.clusters){
    
    # in this loop we need to create the formula string beforehand, since compare_means(formula= ... ) is stupid and doesnt recgnize the running variable i, 
    # instead it searches for a column called i, or whatever you try. so lets it like this:
    
    formula <- paste(i, "~ UICC_Stage_grouped")
    

    
    
    compp <- compare_means(as.formula(formula),
                           data = subset(frequency.data.table, MC %in% m),
                           method = test.method,
                           paired = F,
                           p.adjust.method = adj.p
    )   
    
    
    # compp <- compp[compp$p.adj<0.05,]
    compp <- compp[compp$p<0.05,]
    
    if(nrow(compp)>0){
      compp <- cbind(compp,"MC"=m)
      compp <- cbind(compp,"y"= seq(2000, length=nrow(compp), by=70)  )
      stats_table <- rbind(stats_table,compp)
    }
  }
  
  # for geom_significant to work, we need to move some stuff around in stats_table
  if(nrow(stats_table) >0){
    stats_table <- stats_table %>% 
      dplyr::rename(
        measure = .y.,
        xmin = group1,
        xmax = group2
      )
  
    

    #---------- here is the snipplet to filter only for interesting comparisons. not needed in this case here ------------ 
  #  filter.combinations <- filter.combinations.invert <- as.data.frame(outer(unique(frequenciesofinterestingClusters$Annotated_metacluster)
  #                                                                           ,c( "Healthy","Tumor"), paste, sep="_"))
  #  
  #  names(filter.combinations) <- c("xmin","xmax")
  #  # now, there is a chance that in the stats_table, H and T are inverted, so lets put that into filter.combinations as well:
  #  names(filter.combinations.invert) <- c("xmax","xmin")
  #  # then bring the inverted back into our filter combos:
  #  filter.combinations <-rbind(filter.combinations,filter.combinations.invert)
  #  rm(filter.combinations.invert)
  #.-......................................................................................................................  
   
    
     # we gonna re-set the y positions after filtering:
    stats_table_filtered_yAdjust <- data.frame()
    
    for(m in unique( stats_table$MC ) ){
      
      stats_temp <- subset(stats_table, MC %in% m)
      
      #for(m in unique(stats_table$Annotated_metacluster)){
      #  stats_temp <- subset(stats_table, Annotated_metacluster %in% m)
      
      
      # max.data.value <- max( subset(interstingClusterDensity.perROI.data.table, MC %in% m)[,i] , na.rm=T)
      max.data.value <- max( as.numeric(quantile( subset(frequency.data.table,MC %in% m)[,i], probs = 0.98, na.rm=T )))
      
      
      stats_temp$y <- seq(max.data.value , 
                          length=nrow(stats_temp), by= max( subset(frequency.data.table, MC %in% m)[,i] , na.rm=T) * 0.3 )  
      
      
      # watch out, Inf is passing through here:
      #bracketoffset <- 0.4*max( interstingClusterDensity.perROI.data.table[, i], na.rm=T )
      # better use is.finite and pass the vector like that:
      #bracketoffset <- 0.4*max(  interstingClusterDensity.perROI.data.table[, i][is.finite(interstingClusterDensity.perROI.data.table[, i])]  )
      #print(bracketoffset)                                          
      
      stats_table_filtered_yAdjust <- rbind(stats_table_filtered_yAdjust, stats_temp)
      
      
    }# end run with m along the annotated mc in the stats_table
    
    
    # not needed at this point:
    # prepare the columns we use on X axis to have the proper order:
    #interstingClusterDensity.permouse.data.table$Cond <- factor(interstingClusterDensity.permouse.data.table$Cond,  levels = all.condtions )
    #interstingClusterDensity.permouse.data.table$MC <- factor(interstingClusterDensity.permouse.data.table$MC,  levels = cluster.order )
     stats_table_filtered_yAdjust$MC <- factor(stats_table_filtered_yAdjust$MC,  levels = CD4.clusters )
    
    
     #only process the statistics if you got some in the data in stats_table 
  }else{
    # if the stats_table remained empty since everything was NS, we still need to reset the y-adjust since we absolutely need to protect from drawing stars from the previous run:
    stats_table_filtered_yAdjust <- stats_table # both are empty now
  }
  
  
  
 g<- ggplot(subset(frequency.data.table, MC %in% CD4.clusters)   , 
          aes(x = as.factor(UICC_Stage_grouped), 
                     y = density_mm2, # density_mm2 percentofCD4
              color=as.factor(UICC_Stage_grouped)
          )
  )+
    facet_wrap2(.~MC,
                strip = strip,
                labeller = labeller(source = facet.labs),
                nrow = 5, scales="free_y"
                )+
    geom_violin(trim=T)+
    geom_quasirandom(alpha=1,size=2) +
    scale_colour_manual( values = tumorStage_grouped_palette, limits = force  )+
   
   # significance brackets are cut off, expand axis:
   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
   
   
   labs(title="Cluster Densities on ROIs" , 
        subtitle="Split by Tumor Stage Groups", 
        caption=test.stat.lyrics, 
        y="Cluster density [cells/mm^2]", 
        x=""
   )+

    theme(panel.background = element_rect(fill = "grey90", colour = "grey55", linewidth = 0.5), # change 'colour' to black for informative axis
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=8),
          #axis.text.x = element_text(angle = 45, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          legend.position="none",
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )
  # only add brackets if there are significances above threshold
  if(nrow(stats_table_filtered_yAdjust) >0){
    g <- g+
      geom_signif(
        data = stats_table_filtered_yAdjust,
        aes(xmin = as.factor(xmin), 
            xmax = as.factor(xmax), 
            annotations = p.signif,#p.format,#p.signif, 
            y_position = y),
        textsize = 6, 
        tip_length=0.03,
        color= "black",
        vjust = 0.5,
        manual = TRUE
      )
  }#only process the statistics if you got some in the data
  
  plot(g)
  ggsave("CD4_densities_byTumorSTAGE.pdf", g, width = 10, height = 6, dpi=600)
  
  
  
  
  
  
  
## ------------------------ CD4_frequencies_byTumorSTAGE ------------------------
  
  # we bring here the new code from analyzing Laias IMC data in a statistaclly more clean way:
  # define a global statistic to keep in all plots:
  test.method <- "wilcox.test"# "wilcox.test", "t.test"
  adj.p <- "BH"# "BH" #"none"
  test.stat.lyrics <- paste0(test.method,", p.adj: ",adj.p )
  
  
  
  
  # for every cluster we need to do that separately:
  
  i <- "percentofCD4"
  
  stats_table <- data.frame()
  
  for(m in CD4.clusters){
    
    # in this loop we need to create the formula string beforehand, since compare_means(formula= ... ) is stupid and doesnt recgnize the running variable i, 
    # instead it searches for a column called i, or whatever you try. so lets it like this:
    
    formula <- paste(i, "~ UICC_Stage_grouped")
    
    
    
    
    compp <- compare_means(as.formula(formula),
                           data = subset(frequency.data.table, MC %in% m),
                           method = test.method,
                           paired = F,
                           p.adjust.method = adj.p
    )   
    
    
    # compp <- compp[compp$p.adj<0.05,]
    compp <- compp[compp$p<0.05,]
    
    if(nrow(compp)>0){
      compp <- cbind(compp,"MC"=m)
      compp <- cbind(compp,"y"= seq(2000, length=nrow(compp), by=70)  )
      stats_table <- rbind(stats_table,compp)
    }
  }
  
  # for geom_significant to work, we need to move some stuff around in stats_table
  if(nrow(stats_table) >0){
    stats_table <- stats_table %>% 
      dplyr::rename(
        measure = .y.,
        xmin = group1,
        xmax = group2
      )
    
  
    
    #---------- here is the snipplet to filter only for interesting comparisons. not needed in this case here ------------ 
    #  filter.combinations <- filter.combinations.invert <- as.data.frame(outer(unique(frequenciesofinterestingClusters$Annotated_metacluster)
    #                                                                           ,c( "Healthy","Tumor"), paste, sep="_"))
    #  
    #  names(filter.combinations) <- c("xmin","xmax")
    #  # now, there is a chance that in the stats_table, H and T are inverted, so lets put that into filter.combinations as well:
    #  names(filter.combinations.invert) <- c("xmax","xmin")
    #  # then bring the inverted back into our filter combos:
    #  filter.combinations <-rbind(filter.combinations,filter.combinations.invert)
    #  rm(filter.combinations.invert)
    #.-......................................................................................................................  
    
    
    # we gonna re-set the y positions after filtering:
    stats_table_filtered_yAdjust <- data.frame()
    
    for(m in unique( stats_table$MC ) ){
      
      stats_temp <- subset(stats_table, MC %in% m)
      
      #for(m in unique(stats_table$Annotated_metacluster)){
      #  stats_temp <- subset(stats_table, Annotated_metacluster %in% m)
      
      
      # max.data.value <- max( subset(interstingClusterDensity.perROI.data.table, MC %in% m)[,i] , na.rm=T)
      max.data.value <- max( as.numeric(quantile( subset(frequency.data.table,MC %in% m)[,i], probs = 0.98, na.rm=T )))
      
      
      stats_temp$y <- seq(max.data.value , 
                          length=nrow(stats_temp), by= max( subset(frequency.data.table, MC %in% m)[,i] , na.rm=T) * 0.3 )  
      
      
      # watch out, Inf is passing through here:
      #bracketoffset <- 0.4*max( interstingClusterDensity.perROI.data.table[, i], na.rm=T )
      # better use is.finite and pass the vector like that:
      #bracketoffset <- 0.4*max(  interstingClusterDensity.perROI.data.table[, i][is.finite(interstingClusterDensity.perROI.data.table[, i])]  )
      #print(bracketoffset)                                          
      
      stats_table_filtered_yAdjust <- rbind(stats_table_filtered_yAdjust, stats_temp)
      
      
    }# end run with m along the annotated mc in the stats_table
    
    
    # not needed at this point:
    # prepare the columns we use on X axis to have the proper order:
    #interstingClusterDensity.permouse.data.table$Cond <- factor(interstingClusterDensity.permouse.data.table$Cond,  levels = all.condtions )
    #interstingClusterDensity.permouse.data.table$MC <- factor(interstingClusterDensity.permouse.data.table$MC,  levels = cluster.order )
    stats_table_filtered_yAdjust$MC <- factor(stats_table_filtered_yAdjust$MC,  levels = CD4.clusters )
    
    
    #only process the statistics if you got some in the data in stats_table 
  }else{
    # if the stats_table remained empty since everything was NS, we still need to reset the y-adjust since we absolutely need to protect from drawing stars from the previous run:
    stats_table_filtered_yAdjust <- stats_table # both are empty now
  }
  
  
  
  g<- ggplot(subset(frequency.data.table, MC %in% CD4.clusters)   , 
             aes(x = as.factor(UICC_Stage_grouped), 
                 y = percentofCD4,
                 color=as.factor(UICC_Stage_grouped)
             )
  )+
    facet_wrap2(.~MC,
                strip = strip,
                labeller = labeller(source = facet.labs),
                nrow = 5, scales="free_y"
    )+
    geom_violin(trim=T)+
    geom_quasirandom(alpha=1,size=2) +
    scale_colour_manual( values = tumorStage_grouped_palette, limits = force  )+
    
    # significance brackets are cut off, expand axis:
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
    
    
    labs(title="CD4+ cluster composition on ROIs" , 
         subtitle="Split by Tumor Stage Groups", 
         caption=test.stat.lyrics, 
         y="Cluster frequency [% of all CD4+ (threshold)]", 
         x=""
    )+
    
    theme(panel.background = element_rect(fill = "grey90", colour = "grey55", linewidth = 0.5), # change 'colour' to black for informative axis
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=8),
          #axis.text.x = element_text(angle = 45, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          legend.position="none",
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )
  # only add brackets if there are significances above threshold
  if(nrow(stats_table_filtered_yAdjust) >0){
    g <- g+
      geom_signif(
        data = stats_table_filtered_yAdjust,
        aes(xmin = as.factor(xmin), 
            xmax = as.factor(xmax), 
            annotations = p.signif,#p.format,#p.signif, 
            y_position = y),
        textsize = 6, 
        tip_length=0.03,
        color= "black",
        vjust = 0.5,
        manual = TRUE
      )
  }#only process the statistics if you got some in the data
  
  plot(g)
  ggsave("CD4_frequencies_byTumorSTAGE.pdf", g, width = 10, height = 6, dpi=600)
  
  
  
## ------------------------ CD4_frequencies_byTumorSTAGE_byLongTimeSurvival ------------------------
# I wanna see that cluster-wise reaction to tumor stage now also by survival, long time survival e.g.:  
  
  # we bring here the new code from analyzing Laias IMC data in a statistaclly more clean way:
  # define a global statistic to keep in all plots:
  test.method <- "wilcox.test"# "wilcox.test", "t.test"
  adj.p <- "BH"# "BH" #"none"
  test.stat.lyrics <- paste0(test.method,", p.adj: ",adj.p )
  
  
  
  
  # for every cluster we need to do that separately:
  
  i <- "percentofCD4"
  
  stats_table <- data.frame()
  
  for(m in CD4.clusters){
    
    # in this loop we need to create the formula string beforehand, since compare_means(formula= ... ) is stupid and doesnt recgnize the running variable i, 
    # instead it searches for a column called i, or whatever you try. so lets it like this:
    
    formula <- paste(i, "~ UICC_Stage_grouped")
    
    
    # we further need to subset the clusters for longtime_survivor and run the compp like so:
    for(l in c(0,1)){
     
      
    
    
    compp <- compare_means(as.formula(formula),
                           data = subset(frequency.data.table, MC %in% m & longtime_survivor == l ),
                           method = test.method,
                           paired = F,
                           p.adjust.method = adj.p
    )   
    
    
    # compp <- compp[compp$p.adj<0.05,]
    compp <- compp[compp$p<0.05,]
    
    if(nrow(compp)>0){
      compp <- cbind(compp,"MC"=m)
      compp <- cbind(compp,"longtime_survivor"=l)
      compp <- cbind(compp,"y"= seq(0, length=nrow(compp), by=0)  )
      stats_table <- rbind(stats_table,compp)
    }
    
    
    
    }# further subset each cluster by long-time survival 
  }#run through all clusters sep and calc statistics
  
  # for geom_significant to work, we need to move some stuff around in stats_table
  if(nrow(stats_table) >0){
    stats_table <- stats_table %>% 
      dplyr::rename(
        measure = .y.,
        xmin = group1,
        xmax = group2
      )
    
    
    
    #---------- here is the snipplet to filter only for interesting comparisons. not needed in this case here ------------ 
    #  filter.combinations <- filter.combinations.invert <- as.data.frame(outer(unique(frequenciesofinterestingClusters$Annotated_metacluster)
    #                                                                           ,c( "Healthy","Tumor"), paste, sep="_"))
    #  
    #  names(filter.combinations) <- c("xmin","xmax")
    #  # now, there is a chance that in the stats_table, H and T are inverted, so lets put that into filter.combinations as well:
    #  names(filter.combinations.invert) <- c("xmax","xmin")
    #  # then bring the inverted back into our filter combos:
    #  filter.combinations <-rbind(filter.combinations,filter.combinations.invert)
    #  rm(filter.combinations.invert)
    #.-......................................................................................................................  
    
    
    # we gonna re-set the y positions after filtering:
    stats_table_filtered_yAdjust <- data.frame()
    
    for(m in unique( stats_table$MC ) ){
      
      stats_temp <- subset(stats_table, MC %in% m)
      
      #for(m in unique(stats_table$Annotated_metacluster)){
      #  stats_temp <- subset(stats_table, Annotated_metacluster %in% m)
      
      
      # max.data.value <- max( subset(interstingClusterDensity.perROI.data.table, MC %in% m)[,i] , na.rm=T)
      max.data.value <- max( as.numeric(quantile( subset(frequency.data.table,MC %in% m)[,i], probs = 0.98, na.rm=T )))
      
      
      stats_temp$y <- seq(max.data.value , 
                          length=nrow(stats_temp), by= max( subset(frequency.data.table, MC %in% m)[,i] , na.rm=T) * 0.3 )  
      
      
      # watch out, Inf is passing through here:
      #bracketoffset <- 0.4*max( interstingClusterDensity.perROI.data.table[, i], na.rm=T )
      # better use is.finite and pass the vector like that:
      #bracketoffset <- 0.4*max(  interstingClusterDensity.perROI.data.table[, i][is.finite(interstingClusterDensity.perROI.data.table[, i])]  )
      #print(bracketoffset)                                          
      
      stats_table_filtered_yAdjust <- rbind(stats_table_filtered_yAdjust, stats_temp)
      
      
    }# end run with m along the annotated mc in the stats_table
    
    
    # not needed at this point:
    # prepare the columns we use on X axis to have the proper order:
    #interstingClusterDensity.permouse.data.table$Cond <- factor(interstingClusterDensity.permouse.data.table$Cond,  levels = all.condtions )
    #interstingClusterDensity.permouse.data.table$MC <- factor(interstingClusterDensity.permouse.data.table$MC,  levels = cluster.order )
    stats_table_filtered_yAdjust$MC <- factor(stats_table_filtered_yAdjust$MC,  levels = CD4.clusters )
    
    
    #only process the statistics if you got some in the data in stats_table 
  }else{
    # if the stats_table remained empty since everything was NS, we still need to reset the y-adjust since we absolutely need to protect from drawing stars from the previous run:
    stats_table_filtered_yAdjust <- stats_table # both are empty now
  }
  
 
# I dont want to change the underlying data, so lets try another approach to properly tag long/short time survival in the facet:
  
  temp.labeller = labeller(MC = facet.labs, # keep that one 
                      longtime_survivor = c("0"="<36m","1"=">36m" ))
  
  
  
  
  
  g<- ggplot(subset(frequency.data.table, MC %in% CD4.clusters)   , 
             aes(x = as.factor(UICC_Stage_grouped), 
                 y = percentofCD4,
                 color=as.factor(UICC_Stage_grouped)
             )
  )+
    facet_wrap2(longtime_survivor~MC,
                strip = strip,
                labeller=temp.labeller,
                #labeller = labeller(source = facet.labs),
                nrow = 2, scales="free_y"
    )+
    geom_violin(trim=T)+
    geom_quasirandom(alpha=1,size=2) +
    scale_colour_manual( values = tumorStage_grouped_palette, limits = force  )+
    
    # significance brackets are cut off, expand axis:
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
    
    
    labs(title="CD4+ cluster composition on ROIs" , 
         subtitle="Split by Tumor Stage Groups", 
         caption=test.stat.lyrics, 
         y="Cluster frequency [% of all CD4+ (threshold)]", 
         x=""
    )+
    
    theme(panel.background = element_rect(fill = "grey90", colour = "grey55", linewidth = 0.5), # change 'colour' to black for informative axis
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=8),
          #axis.text.x = element_text(angle = 45, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          legend.position="none",
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )
  # only add brackets if there are significances above threshold
  if(nrow(stats_table_filtered_yAdjust) >0){
    g <- g+
      geom_signif(
        data = stats_table_filtered_yAdjust,
        aes(xmin = as.factor(xmin), 
            xmax = as.factor(xmax), 
            annotations = p.signif,#p.format,#p.signif, 
            y_position = y),
        textsize = 6, 
        tip_length=0.03,
        color= "black",
        vjust = 0.5,
        manual = TRUE
      )
  }#only process the statistics if you got some in the data
  
  plot(g)
  ggsave("CD4_frequencies_byTumorSTAGE_byLongTimeSurvival.pdf", g, width = 16, height = 6, dpi=600)
  
  
  
  
  
  ## ------------------------ CD4_densities_byTumorSTAGE_byLongTimeSurvival ------------------------
  # I wanna see that cluster-wise reaction to tumor stage now also by survival, long time survival e.g.:  
  
  # we bring here the new code from analyzing Laias IMC data in a statistaclly more clean way:
  # define a global statistic to keep in all plots:
  test.method <- "wilcox.test"# "wilcox.test", "t.test"
  adj.p <- "BH"# "BH" #"none"
  test.stat.lyrics <- paste0(test.method,", p.adj: ",adj.p )
  
  
  
  
  # for every cluster we need to do that separately:
  
  i <- "density_mm2"
  
  stats_table <- data.frame()
  
  for(m in CD4.clusters){
    
    # in this loop we need to create the formula string beforehand, since compare_means(formula= ... ) is stupid and doesnt recgnize the running variable i, 
    # instead it searches for a column called i, or whatever you try. so lets it like this:
    
    formula <- paste(i, "~ UICC_Stage_grouped")
    
    
    # we further need to subset the clusters for longtime_survivor and run the compp like so:
    for(l in c(0,1)){
      
      
      
      
      compp <- compare_means(as.formula(formula),
                             data = subset(frequency.data.table, MC %in% m & longtime_survivor == l ),
                             method = test.method,
                             paired = F,
                             p.adjust.method = adj.p
      )   
      
      
      # compp <- compp[compp$p.adj<0.05,]
      compp <- compp[compp$p<0.05,]
      
      if(nrow(compp)>0){
        compp <- cbind(compp,"MC"=m)
        compp <- cbind(compp,"longtime_survivor"=l)
        compp <- cbind(compp,"y"= seq(0, length=nrow(compp), by=0)  )
        stats_table <- rbind(stats_table,compp)
      }
      
      
      
    }# further subset each cluster by long-time survival 
  }#run through all clusters sep and calc statistics
  
  # for geom_significant to work, we need to move some stuff around in stats_table
  if(nrow(stats_table) >0){
    stats_table <- stats_table %>% 
      dplyr::rename(
        measure = .y.,
        xmin = group1,
        xmax = group2
      )
    
    
    
    #---------- here is the snipplet to filter only for interesting comparisons. not needed in this case here ------------ 
    #  filter.combinations <- filter.combinations.invert <- as.data.frame(outer(unique(frequenciesofinterestingClusters$Annotated_metacluster)
    #                                                                           ,c( "Healthy","Tumor"), paste, sep="_"))
    #  
    #  names(filter.combinations) <- c("xmin","xmax")
    #  # now, there is a chance that in the stats_table, H and T are inverted, so lets put that into filter.combinations as well:
    #  names(filter.combinations.invert) <- c("xmax","xmin")
    #  # then bring the inverted back into our filter combos:
    #  filter.combinations <-rbind(filter.combinations,filter.combinations.invert)
    #  rm(filter.combinations.invert)
    #.-......................................................................................................................  
    
    
    # we gonna re-set the y positions after filtering:
    stats_table_filtered_yAdjust <- data.frame()
    
    for(m in unique( stats_table$MC ) ){
      
      stats_temp <- subset(stats_table, MC %in% m)
      
      #for(m in unique(stats_table$Annotated_metacluster)){
      #  stats_temp <- subset(stats_table, Annotated_metacluster %in% m)
      
      
      # max.data.value <- max( subset(interstingClusterDensity.perROI.data.table, MC %in% m)[,i] , na.rm=T)
      max.data.value <- max( as.numeric(quantile( subset(frequency.data.table,MC %in% m)[,i], probs = 0.98, na.rm=T )))
      
      
      stats_temp$y <- seq(max.data.value , 
                          length=nrow(stats_temp), by= max( subset(frequency.data.table, MC %in% m)[,i] , na.rm=T) * 0.3 )  
      
      
      # watch out, Inf is passing through here:
      #bracketoffset <- 0.4*max( interstingClusterDensity.perROI.data.table[, i], na.rm=T )
      # better use is.finite and pass the vector like that:
      #bracketoffset <- 0.4*max(  interstingClusterDensity.perROI.data.table[, i][is.finite(interstingClusterDensity.perROI.data.table[, i])]  )
      #print(bracketoffset)                                          
      
      stats_table_filtered_yAdjust <- rbind(stats_table_filtered_yAdjust, stats_temp)
      
      
    }# end run with m along the annotated mc in the stats_table
    
    
    # not needed at this point:
    # prepare the columns we use on X axis to have the proper order:
    #interstingClusterDensity.permouse.data.table$Cond <- factor(interstingClusterDensity.permouse.data.table$Cond,  levels = all.condtions )
    #interstingClusterDensity.permouse.data.table$MC <- factor(interstingClusterDensity.permouse.data.table$MC,  levels = cluster.order )
    stats_table_filtered_yAdjust$MC <- factor(stats_table_filtered_yAdjust$MC,  levels = CD4.clusters )
    
    
    #only process the statistics if you got some in the data in stats_table 
  }else{
    # if the stats_table remained empty since everything was NS, we still need to reset the y-adjust since we absolutely need to protect from drawing stars from the previous run:
    stats_table_filtered_yAdjust <- stats_table # both are empty now
  }
  
  
  # I dont want to change the underlying data, so lets try another approach to properly tag long/short time survival in the facet:
  
  temp.labeller = labeller(MC = facet.labs, # keep that one 
                           longtime_survivor = c("0"="<36m","1"=">36m" ))
  
  
  
  
  
  g<- ggplot(subset(frequency.data.table, MC %in% CD4.clusters)   , 
             aes(x = as.factor(UICC_Stage_grouped), 
                 y = density_mm2,
                 color=as.factor(UICC_Stage_grouped)
             )
  )+
    facet_wrap2(longtime_survivor~MC,
                strip = strip,
                labeller=temp.labeller,
                #labeller = labeller(source = facet.labs),
                nrow = 2, scales="free_y"
    )+
    geom_violin(trim=T)+
    geom_quasirandom(alpha=1,size=2) +
    scale_colour_manual( values = tumorStage_grouped_palette, limits = force  )+
    
    # significance brackets are cut off, expand axis:
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
    
    
    labs(title="CD4+ cluster densities on ROIs" , 
         subtitle="Split by Tumor Stage Groups", 
         caption=test.stat.lyrics, 
         y="Cluster density [cells/mm^2]", 
         x=""
    )+
    
    theme(panel.background = element_rect(fill = "grey90", colour = "grey55", linewidth = 0.5), # change 'colour' to black for informative axis
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=8),
          #axis.text.x = element_text(angle = 45, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          legend.position="none",
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )
  # only add brackets if there are significances above threshold
  if(nrow(stats_table_filtered_yAdjust) >0){
    g <- g+
      geom_signif(
        data = stats_table_filtered_yAdjust,
        aes(xmin = as.factor(xmin), 
            xmax = as.factor(xmax), 
            annotations = p.signif,#p.format,#p.signif, 
            y_position = y),
        textsize = 6, 
        tip_length=0.03,
        color= "black",
        vjust = 0.5,
        manual = TRUE
      )
  }#only process the statistics if you got some in the data
  
  plot(g)
  ggsave("CD4_densities_byTumorSTAGE_byLongTimeSurvival.pdf", g, width = 16, height = 6, dpi=600)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #
  
  ## ------------------------ CD4_densities_bylongtimesurvival ------------------------
  
  # we bring here the new code from analyzing Laias IMC data in a statistaclly more clean way:
  # define a global statistic to keep in all plots:
  test.method <- "wilcox.test"# "wilcox.test", "t.test"
  adj.p <- "BH"# "BH" #"none"
  test.stat.lyrics <- paste0(test.method,", p.adj: ",adj.p )
  
  
  # for every cluster we need to do that separately:
  
  i <- "density_mm2"
  
  stats_table <- data.frame()
  
  for(m in CD4.clusters){
    
    # in this loop we need to create the formula string beforehand, since compare_means(formula= ... ) is stupid and doesnt recgnize the running variable i, 
    # instead it searches for a column called i, or whatever you try. so lets it like this:
    
    formula <- paste(i, "~ longtime_survivor")
    
    
    
    
    compp <- compare_means(as.formula(formula),
                           data = subset(frequency.data.table, MC %in% m),
                           method = test.method,
                           paired = F,
                           p.adjust.method = adj.p
    )   
    
    
    # compp <- compp[compp$p.adj<0.05,]
    compp <- compp[compp$p<0.05,]
    
    if(nrow(compp)>0){
      compp <- cbind(compp,"MC"=m)
      compp <- cbind(compp,"y"= seq(2000, length=nrow(compp), by=70)  )
      stats_table <- rbind(stats_table,compp)
    }
  }
  
  # for geom_significant to work, we need to move some stuff around in stats_table
  if(nrow(stats_table) >0){
    stats_table <- stats_table %>% 
      dplyr::rename(
        measure = .y.,
        xmin = group1,
        xmax = group2
      )
    
    
    
    #---------- here is the snipplet to filter only for interesting comparisons. not needed in this case here ------------ 
    #  filter.combinations <- filter.combinations.invert <- as.data.frame(outer(unique(frequenciesofinterestingClusters$Annotated_metacluster)
    #                                                                           ,c( "Healthy","Tumor"), paste, sep="_"))
    #  
    #  names(filter.combinations) <- c("xmin","xmax")
    #  # now, there is a chance that in the stats_table, H and T are inverted, so lets put that into filter.combinations as well:
    #  names(filter.combinations.invert) <- c("xmax","xmin")
    #  # then bring the inverted back into our filter combos:
    #  filter.combinations <-rbind(filter.combinations,filter.combinations.invert)
    #  rm(filter.combinations.invert)
    #.-......................................................................................................................  
    
    
    # we gonna re-set the y positions after filtering:
    stats_table_filtered_yAdjust <- data.frame()
    
    for(m in unique( stats_table$MC ) ){
      
      stats_temp <- subset(stats_table, MC %in% m)
      
      #for(m in unique(stats_table$Annotated_metacluster)){
      #  stats_temp <- subset(stats_table, Annotated_metacluster %in% m)
      
      
       max.data.value <- max( subset(frequency.data.table, MC %in% m)[,i] , na.rm=T)
     # max.data.value <- max( as.numeric(quantile( subset(frequency.data.table,MC %in% m)[,i], probs = 0.98, na.rm=T )))
      
      
      stats_temp$y <- seq(max.data.value , 
                          length=nrow(stats_temp), by= max( subset(frequency.data.table, MC %in% m)[,i] , na.rm=T) * 0.3 )  
      
      
      # watch out, Inf is passing through here:
      #bracketoffset <- 0.4*max( interstingClusterDensity.perROI.data.table[, i], na.rm=T )
      # better use is.finite and pass the vector like that:
      #bracketoffset <- 0.4*max(  interstingClusterDensity.perROI.data.table[, i][is.finite(interstingClusterDensity.perROI.data.table[, i])]  )
      #print(bracketoffset)                                          
      
      stats_table_filtered_yAdjust <- rbind(stats_table_filtered_yAdjust, stats_temp)
      
      
    }# end run with m along the annotated mc in the stats_table
    
    
    # not needed at this point:
    # prepare the columns we use on X axis to have the proper order:
    #interstingClusterDensity.permouse.data.table$Cond <- factor(interstingClusterDensity.permouse.data.table$Cond,  levels = all.condtions )
    #interstingClusterDensity.permouse.data.table$MC <- factor(interstingClusterDensity.permouse.data.table$MC,  levels = cluster.order )
    stats_table_filtered_yAdjust$MC <- factor(stats_table_filtered_yAdjust$MC,  levels = CD4.clusters )
    
    
    #only process the statistics if you got some in the data in stats_table 
  }else{
    # if the stats_table remained empty since everything was NS, we still need to reset the y-adjust since we absolutely need to protect from drawing stars from the previous run:
    stats_table_filtered_yAdjust <- stats_table # both are empty now
  }
  
  
  
  g<- ggplot(subset(frequency.data.table, MC %in% CD4.clusters)   , 
             aes(x = as.factor(longtime_survivor), 
                 y = density_mm2, # density_mm2 percentofCD4
                 color=as.factor(longtime_survivor)
             )
  )+
    facet_wrap2(.~MC,
                strip = strip,
                labeller = labeller(source = facet.labs),
                nrow = 5, scales="free_y"
    )+
    geom_violin(trim=T)+
    geom_quasirandom(alpha=1,size=2) +
    scale_colour_manual( values = longtimesurvival_palette, limits = force  )+
    
    # significance brackets are cut off, expand axis:
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
    
    
    labs(title="CD4+ cluster densities on ROIs" , 
         subtitle="Split by long-time suvival (>36months)", 
         caption=test.stat.lyrics, 
         y="Cluster density [cells/mm^2]", 
         x=""
    )+
    
    theme(panel.background = element_rect(fill = "grey90", colour = "grey55", linewidth = 0.5), # change 'colour' to black for informative axis
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=8),
          #axis.text.x = element_text(angle = 45, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          legend.position="none",
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )
  # only add brackets if there are significances above threshold
  if(nrow(stats_table_filtered_yAdjust) >0){
    g <- g+
      geom_signif(
        data = stats_table_filtered_yAdjust,
        aes(xmin = as.factor(xmin), 
            xmax = as.factor(xmax), 
            annotations = p.signif,#p.format,#p.signif, 
            y_position = y),
        textsize = 6, 
        tip_length=0.03,
        color= "black",
        vjust = 0.5,
        manual = TRUE
      )
  }#only process the statistics if you got some in the data
  
  plot(g)
  ggsave("CD4_densities_bylongtimesurvival.pdf", g, width = 10, height = 6, dpi=600)
  
  
  
  
  
  
  
  
  
  
  ## ------------------------ CD4_frequencies_bylongtimesurvival ------------------------
  
  # we bring here the new code from analyzing Laias IMC data in a statistaclly more clean way:
  # define a global statistic to keep in all plots:
  test.method <- "wilcox.test"# "wilcox.test", "t.test"
  adj.p <- "BH"# "BH" #"none"
  test.stat.lyrics <- paste0(test.method,", p.adj: ",adj.p )
  
  
  # for every cluster we need to do that separately:
  
  i <- "percentofCD4"
  
  stats_table <- data.frame()
  
  for(m in CD4.clusters){
    
    # in this loop we need to create the formula string beforehand, since compare_means(formula= ... ) is stupid and doesnt recgnize the running variable i, 
    # instead it searches for a column called i, or whatever you try. so lets it like this:
    
    formula <- paste(i, "~ longtime_survivor")
    
    
    
    
    compp <- compare_means(as.formula(formula),
                           data = subset(frequency.data.table, MC %in% m),
                           method = test.method,
                           paired = F,
                           p.adjust.method = adj.p
    )   
    
    
    # compp <- compp[compp$p.adj<0.05,]
    compp <- compp[compp$p<0.05,]
    
    if(nrow(compp)>0){
      compp <- cbind(compp,"MC"=m)
      compp <- cbind(compp,"y"= seq(2000, length=nrow(compp), by=70)  )
      stats_table <- rbind(stats_table,compp)
    }
  }
  
  # for geom_significant to work, we need to move some stuff around in stats_table
  if(nrow(stats_table) >0){
    stats_table <- stats_table %>% 
      dplyr::rename(
        measure = .y.,
        xmin = group1,
        xmax = group2
      )
    
    
    
    #---------- here is the snipplet to filter only for interesting comparisons. not needed in this case here ------------ 
    #  filter.combinations <- filter.combinations.invert <- as.data.frame(outer(unique(frequenciesofinterestingClusters$Annotated_metacluster)
    #                                                                           ,c( "Healthy","Tumor"), paste, sep="_"))
    #  
    #  names(filter.combinations) <- c("xmin","xmax")
    #  # now, there is a chance that in the stats_table, H and T are inverted, so lets put that into filter.combinations as well:
    #  names(filter.combinations.invert) <- c("xmax","xmin")
    #  # then bring the inverted back into our filter combos:
    #  filter.combinations <-rbind(filter.combinations,filter.combinations.invert)
    #  rm(filter.combinations.invert)
    #.-......................................................................................................................  
    
    
    # we gonna re-set the y positions after filtering:
    stats_table_filtered_yAdjust <- data.frame()
    
    for(m in unique( stats_table$MC ) ){
      
      stats_temp <- subset(stats_table, MC %in% m)
      
      #for(m in unique(stats_table$Annotated_metacluster)){
      #  stats_temp <- subset(stats_table, Annotated_metacluster %in% m)
      
      
     # max.data.value <- max( subset(frequency.data.table, MC %in% m)[,i] , na.rm=T)
       max.data.value <- max( as.numeric(quantile( subset(frequency.data.table,MC %in% m)[,i], probs = 0.98, na.rm=T )))
      
      
      stats_temp$y <- seq(max.data.value , 
                          length=nrow(stats_temp), by= max( subset(frequency.data.table, MC %in% m)[,i] , na.rm=T) * 0.3 )  
      
      
      # watch out, Inf is passing through here:
      #bracketoffset <- 0.4*max( interstingClusterDensity.perROI.data.table[, i], na.rm=T )
      # better use is.finite and pass the vector like that:
      #bracketoffset <- 0.4*max(  interstingClusterDensity.perROI.data.table[, i][is.finite(interstingClusterDensity.perROI.data.table[, i])]  )
      #print(bracketoffset)                                          
      
      stats_table_filtered_yAdjust <- rbind(stats_table_filtered_yAdjust, stats_temp)
      
      
    }# end run with m along the annotated mc in the stats_table
    
    
    # not needed at this point:
    # prepare the columns we use on X axis to have the proper order:
    #interstingClusterDensity.permouse.data.table$Cond <- factor(interstingClusterDensity.permouse.data.table$Cond,  levels = all.condtions )
    #interstingClusterDensity.permouse.data.table$MC <- factor(interstingClusterDensity.permouse.data.table$MC,  levels = cluster.order )
    stats_table_filtered_yAdjust$MC <- factor(stats_table_filtered_yAdjust$MC,  levels = CD4.clusters )
    
    
    #only process the statistics if you got some in the data in stats_table 
  }else{
    # if the stats_table remained empty since everything was NS, we still need to reset the y-adjust since we absolutely need to protect from drawing stars from the previous run:
    stats_table_filtered_yAdjust <- stats_table # both are empty now
  }
  
  
  
  g<- ggplot(subset(frequency.data.table, MC %in% CD4.clusters)   , 
             aes(x = as.factor(longtime_survivor), 
                  y = percentofCD4,
                 color=as.factor(longtime_survivor)
             )
  )+
    facet_wrap2(.~MC,
                strip = strip,
                labeller = labeller(source = facet.labs),
                nrow = 5, scales="free_y"
    )+
    geom_violin(trim=T)+
    geom_quasirandom(alpha=1,size=2) +
    scale_colour_manual( values = longtimesurvival_palette, limits = force  )+
    
    # significance brackets are cut off, expand axis:
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
    

    labs(title="CD4+ cluster composition on ROIs" , 
         subtitle="Split by long-time suvival (>36months)",
         caption=test.stat.lyrics, 
         y="Cluster frequency [% of all CD4+ (threshold)]", 
         x=""
    )+
    
    
    
    theme(panel.background = element_rect(fill = "grey90", colour = "grey55", linewidth = 0.5), # change 'colour' to black for informative axis
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=8),
          #axis.text.x = element_text(angle = 45, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          legend.position="none",
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )
  # only add brackets if there are significances above threshold
  if(nrow(stats_table_filtered_yAdjust) >0){
    g <- g+
      geom_signif(
        data = stats_table_filtered_yAdjust,
        aes(xmin = as.factor(xmin), 
            xmax = as.factor(xmax), 
            annotations = p.signif,#p.format,#p.signif, 
            y_position = y),
        textsize = 6, 
        tip_length=0.03,
        color= "black",
        vjust = 0.5,
        manual = TRUE
      )
  }#only process the statistics if you got some in the data
  
  plot(g)
  ggsave("CD4_frequencies_bylongtimesurvival.pdf", g, width = 10, height = 6, dpi=600)
  
  
  
  
  
  
  ## ------------------------ CD4_frequencies_byTissue ------------------------
  
  # we bring here the new code from analyzing Laias IMC data in a statistaclly more clean way:
  # define a global statistic to keep in all plots:
  test.method <- "wilcox.test"# "wilcox.test", "t.test"
  adj.p <- "BH"# "BH" #"none"
  test.stat.lyrics <- paste0(test.method,", p.adj: ",adj.p )
  
  
  # for every cluster we need to do that separately:
  
  i <- "percentofCD4"
  
  stats_table <- data.frame()
  
  for(m in CD4.clusters){
    
    # in this loop we need to create the formula string beforehand, since compare_means(formula= ... ) is stupid and doesnt recgnize the running variable i, 
    # instead it searches for a column called i, or whatever you try. so lets it like this:
    
    formula <- paste(i, "~ Ref_Area")
    
    
    
    
    compp <- compare_means(as.formula(formula),
                           data = subset(frequency.data.table, MC %in% m),
                           method = test.method,
                           paired = F,
                           p.adjust.method = adj.p
    )   
    
    
    # compp <- compp[compp$p.adj<0.05,]
    compp <- compp[compp$p<0.05,]
    
    if(nrow(compp)>0){
      compp <- cbind(compp,"MC"=m)
      compp <- cbind(compp,"y"= seq(2000, length=nrow(compp), by=70)  )
      stats_table <- rbind(stats_table,compp)
    }
  }
  
  # for geom_significant to work, we need to move some stuff around in stats_table
  if(nrow(stats_table) >0){
    stats_table <- stats_table %>% 
      dplyr::rename(
        measure = .y.,
        xmin = group1,
        xmax = group2
      )
    
    
    
    #---------- here is the snipplet to filter only for interesting comparisons. not needed in this case here ------------ 
    #  filter.combinations <- filter.combinations.invert <- as.data.frame(outer(unique(frequenciesofinterestingClusters$Annotated_metacluster)
    #                                                                           ,c( "Healthy","Tumor"), paste, sep="_"))
    #  
    #  names(filter.combinations) <- c("xmin","xmax")
    #  # now, there is a chance that in the stats_table, H and T are inverted, so lets put that into filter.combinations as well:
    #  names(filter.combinations.invert) <- c("xmax","xmin")
    #  # then bring the inverted back into our filter combos:
    #  filter.combinations <-rbind(filter.combinations,filter.combinations.invert)
    #  rm(filter.combinations.invert)
    #.-......................................................................................................................  
    
    
    # we gonna re-set the y positions after filtering:
    stats_table_filtered_yAdjust <- data.frame()
    
    for(m in unique( stats_table$MC ) ){
      
      stats_temp <- subset(stats_table, MC %in% m)
      
      #for(m in unique(stats_table$Annotated_metacluster)){
      #  stats_temp <- subset(stats_table, Annotated_metacluster %in% m)
      
      
      # max.data.value <- max( subset(frequency.data.table, MC %in% m)[,i] , na.rm=T)
      max.data.value <- max( as.numeric(quantile( subset(frequency.data.table,MC %in% m)[,i], probs = 0.98, na.rm=T )))
      
      
      stats_temp$y <- seq(max.data.value , 
                          length=nrow(stats_temp), by= max( subset(frequency.data.table, MC %in% m)[,i] , na.rm=T) * 0.3 )  
      
      
      # watch out, Inf is passing through here:
      #bracketoffset <- 0.4*max( interstingClusterDensity.perROI.data.table[, i], na.rm=T )
      # better use is.finite and pass the vector like that:
      #bracketoffset <- 0.4*max(  interstingClusterDensity.perROI.data.table[, i][is.finite(interstingClusterDensity.perROI.data.table[, i])]  )
      #print(bracketoffset)                                          
      
      stats_table_filtered_yAdjust <- rbind(stats_table_filtered_yAdjust, stats_temp)
      
      
    }# end run with m along the annotated mc in the stats_table
    
    
    # not needed at this point:
    # prepare the columns we use on X axis to have the proper order:
    #interstingClusterDensity.permouse.data.table$Cond <- factor(interstingClusterDensity.permouse.data.table$Cond,  levels = all.condtions )
    #interstingClusterDensity.permouse.data.table$MC <- factor(interstingClusterDensity.permouse.data.table$MC,  levels = cluster.order )
    stats_table_filtered_yAdjust$MC <- factor(stats_table_filtered_yAdjust$MC,  levels = CD4.clusters )
    
    
    #only process the statistics if you got some in the data in stats_table 
  }else{
    # if the stats_table remained empty since everything was NS, we still need to reset the y-adjust since we absolutely need to protect from drawing stars from the previous run:
    stats_table_filtered_yAdjust <- stats_table # both are empty now
  }
  
  
  frequency.data.table$Ref_Area <- factor(frequency.data.table$Ref_Area,  levels = c("ROI" ,   "stroma",  "tumor") )
  
  
  g<- ggplot(subset(frequency.data.table, MC %in% CD4.clusters)   , 
             aes(x = Ref_Area, 
                 y = percentofCD4,
                 color=Ref_Area
             )
  )+
    facet_wrap2(longtime_survivor~MC,
                strip = strip,
                labeller = labeller(source = facet.labs),
                nrow = 2, scales="free_y"
    )+
    geom_violin(trim=T)+
    geom_quasirandom(alpha=1,size=2) +
    scale_colour_manual( values = tissue_palette, limits = force  )+
    
    # significance brackets are cut off, expand axis:
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
    
    
    labs(title="CD4+ cluster composition on ROIs" , 
         subtitle="Split by tumor mask",
         caption=test.stat.lyrics, 
         y="Cluster frequency [% of all CD4+ (threshold)]", 
         x=""
    )+
    
    
    
    theme(panel.background = element_rect(fill = "grey90", colour = "grey55", linewidth = 0.5), # change 'colour' to black for informative axis
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=8),
          #axis.text.x = element_text(angle = 45, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          legend.position="none",
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )
  # only add brackets if there are significances above threshold
  if(nrow(stats_table_filtered_yAdjust) >0){
    g <- g+
      geom_signif(
        data = stats_table_filtered_yAdjust,
        aes(xmin = as.factor(xmin), 
            xmax = as.factor(xmax), 
            annotations = p.signif,#p.format,#p.signif, 
            y_position = y),
        textsize = 6, 
        tip_length=0.03,
        color= "black",
        vjust = 0.5,
        manual = TRUE
      )
  }#only process the statistics if you got some in the data
  
  plot(g)
  #ggsave("CD4_frequencies_byTissue.pdf", g, width = 10, height = 6, dpi=600)
  ggsave("CD4_frequencies_byTissue_bylongtimesurvival.pdf", g, width = 15, height = 6, dpi=600)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # can we do some regressions with the survival in months?
  
  

  
  
  
  
g<-  ggplot(subset(frequency.data.table, MC %in% CD4.clusters)   , 
         aes(x = monitored_survival_months, 
             y = density_mm2, #percentofCD4,
             color=UICC_Stage_grouped
         )
  )+
    facet_wrap2(died~MC,
                strip = strip,
                labeller = labeller(source = facet.labs),
                nrow = 2, scales="free_y"
    )+

    geom_point(alpha=1,size=2) +
    geom_smooth(method = "lm", formula = 'y ~ x',se = F) +
    scale_colour_manual( values = tumorStage_grouped_palette, limits = force  )
  
  
  ggsave("CD4_densities_as_funct_of_survivaltime_byCluster_bydeath.pdf", g, width = 15, height = 6, dpi=600)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
 
 
  
  
  
  ## ------------------------ CD4_densities_byTissue ------------------------
  
  # we bring here the new code from analyzing Laias IMC data in a statistaclly more clean way:
  # define a global statistic to keep in all plots:
  test.method <- "wilcox.test"# "wilcox.test", "t.test"
  adj.p <- "BH"# "BH" #"none"
  test.stat.lyrics <- paste0(test.method,", p.adj: ",adj.p )
  
  
  # for every cluster we need to do that separately:
  
  i <- "density_mm2"
  
  stats_table <- data.frame()
  
  for(m in CD4.clusters){
    
    # in this loop we need to create the formula string beforehand, since compare_means(formula= ... ) is stupid and doesnt recgnize the running variable i, 
    # instead it searches for a column called i, or whatever you try. so lets it like this:
    
    formula <- paste(i, "~ Ref_Area")
    
    
    
    
    compp <- compare_means(as.formula(formula),
                           data = subset(frequency.data.table, MC %in% m),
                           method = test.method,
                           paired = F,
                           p.adjust.method = adj.p
    )   
    
    
    # compp <- compp[compp$p.adj<0.05,]
    compp <- compp[compp$p<0.05,]
    
    if(nrow(compp)>0){
      compp <- cbind(compp,"MC"=m)
      compp <- cbind(compp,"y"= seq(2000, length=nrow(compp), by=70)  )
      stats_table <- rbind(stats_table,compp)
    }
  }
  
  # for geom_significant to work, we need to move some stuff around in stats_table
  if(nrow(stats_table) >0){
    stats_table <- stats_table %>% 
      dplyr::rename(
        measure = .y.,
        xmin = group1,
        xmax = group2
      )
    
    
    
    #---------- here is the snipplet to filter only for interesting comparisons. not needed in this case here ------------ 
    #  filter.combinations <- filter.combinations.invert <- as.data.frame(outer(unique(frequenciesofinterestingClusters$Annotated_metacluster)
    #                                                                           ,c( "Healthy","Tumor"), paste, sep="_"))
    #  
    #  names(filter.combinations) <- c("xmin","xmax")
    #  # now, there is a chance that in the stats_table, H and T are inverted, so lets put that into filter.combinations as well:
    #  names(filter.combinations.invert) <- c("xmax","xmin")
    #  # then bring the inverted back into our filter combos:
    #  filter.combinations <-rbind(filter.combinations,filter.combinations.invert)
    #  rm(filter.combinations.invert)
    #.-......................................................................................................................  
    
    
    # we gonna re-set the y positions after filtering:
    stats_table_filtered_yAdjust <- data.frame()
    
    for(m in unique( stats_table$MC ) ){
      
      stats_temp <- subset(stats_table, MC %in% m)
      
      #for(m in unique(stats_table$Annotated_metacluster)){
      #  stats_temp <- subset(stats_table, Annotated_metacluster %in% m)
      
      
      # max.data.value <- max( subset(frequency.data.table, MC %in% m)[,i] , na.rm=T)
      max.data.value <- max( as.numeric(quantile( subset(frequency.data.table,MC %in% m)[,i], probs = 0.98, na.rm=T )))
      
      
      stats_temp$y <- seq(max.data.value , 
                          length=nrow(stats_temp), by= max( subset(frequency.data.table, MC %in% m)[,i] , na.rm=T) * 0.3 )  
      
      
      # watch out, Inf is passing through here:
      #bracketoffset <- 0.4*max( interstingClusterDensity.perROI.data.table[, i], na.rm=T )
      # better use is.finite and pass the vector like that:
      #bracketoffset <- 0.4*max(  interstingClusterDensity.perROI.data.table[, i][is.finite(interstingClusterDensity.perROI.data.table[, i])]  )
      #print(bracketoffset)                                          
      
      stats_table_filtered_yAdjust <- rbind(stats_table_filtered_yAdjust, stats_temp)
      
      
    }# end run with m along the annotated mc in the stats_table
    
    
    # not needed at this point:
    # prepare the columns we use on X axis to have the proper order:
    #interstingClusterDensity.permouse.data.table$Cond <- factor(interstingClusterDensity.permouse.data.table$Cond,  levels = all.condtions )
    #interstingClusterDensity.permouse.data.table$MC <- factor(interstingClusterDensity.permouse.data.table$MC,  levels = cluster.order )
    stats_table_filtered_yAdjust$MC <- factor(stats_table_filtered_yAdjust$MC,  levels = CD4.clusters )
    
    
    #only process the statistics if you got some in the data in stats_table 
  }else{
    # if the stats_table remained empty since everything was NS, we still need to reset the y-adjust since we absolutely need to protect from drawing stars from the previous run:
    stats_table_filtered_yAdjust <- stats_table # both are empty now
  }
  
  
  frequency.data.table$Ref_Area <- factor(frequency.data.table$Ref_Area,  levels = c("ROI" ,   "stroma",  "tumor") )
  
  
  g<- ggplot(subset(frequency.data.table, MC %in% CD4.clusters)   , 
             aes(x = Ref_Area, 
                 y = density_mm2, 
                 color=Ref_Area
             )
  )+
    facet_wrap2(.~MC,
                strip = strip,
                labeller = labeller(source = facet.labs),
                nrow = 5, scales="free_y"
    )+
    geom_violin(trim=T)+
    geom_quasirandom(alpha=1,size=2) +
    scale_colour_manual( values = tissue_palette, limits = force  )+
    
    # significance brackets are cut off, expand axis:
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
    
    

    labs(title="Local Cluster Densities on ROIs" , 
         subtitle="Split by tumor masks", 
         caption=test.stat.lyrics, 
         y="Cluster density in respective area [cells/mm^2]", 
         x=""
    )+
    
    
    
    theme(panel.background = element_rect(fill = "grey90", colour = "grey55", linewidth = 0.5), # change 'colour' to black for informative axis
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=8),
          #axis.text.x = element_text(angle = 45, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          legend.position="none",
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )
  # only add brackets if there are significances above threshold
  if(nrow(stats_table_filtered_yAdjust) >0){
    g <- g+
      geom_signif(
        data = stats_table_filtered_yAdjust,
        aes(xmin = as.factor(xmin), 
            xmax = as.factor(xmax), 
            annotations = p.signif,#p.format,#p.signif, 
            y_position = y),
        textsize = 6, 
        tip_length=0.03,
        color= "black",
        vjust = 0.5,
        manual = TRUE
      )
  }#only process the statistics if you got some in the data
  
  plot(g)
  ggsave("CD4_densities_byTissue.pdf", g, width = 10, height = 6, dpi=600)
  
  
   
  
  
 
  
  
  
  
  
  
  
  
  
  
  ## ------------------------ PD-L1 expression by clusters -------------
  
 subset(gated.cell.dat, Annotated_metacluster %in% c(CD4.clusters,CD8.clusters,322) ) [, lapply(.SD, "mean", na.rm = TRUE), by='Phenograph_metacluster',  
                                      .SDcols = PD.L1_Nd150_t2]
  
  
  
  ggplot( subset(gated.cell.dat, Annotated_metacluster %in% c(CD4.clusters,CD8.clusters,322) )  , 
          aes( x=as.factor(Annotated_metacluster), y= PD.L1_Nd150_t2 , color = PD.L1_Nd150_t2 )  ) +
    
    
    
    #geom_violin(trim=T)+
    geom_quasirandom(size=0.5,alpha=0.8) + #
    geom_boxplot(outlier.shape = NA, lwd=0.8,alpha=0.7) + 
    #  scale_color_manual(values = cond_palette)+
    
   # scale_color_manual( values = subset.Phenograph_metacluster_palette    )+
    
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),         labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides="l")   +
    
    
    
    theme_classic()+
    
    labs(title=  "PD.L1 expression"      , 
         #subtitle= paste0("Using ", nrow(temp.df), " gated non-lymphnode cells in G2.") ,
         #caption="Created by M.Barone", 
         y="Channel signal [asinh transf.]" 
        # x= sub("_.*", "", m )  
    )+
    
    theme(panel.background = element_rect(fill = "grey93", colour = "grey55", size = 0.5), # change 'colour' to black for informative axis
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=10),
          axis.text.x = element_text(angle = 45, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )
  
  
  
  
  
  
  
}#end protect running old code
   }#end stage 6: do additional CyTOF-y stuff 

########################## STAGE 7 ########################## 
# STAGE 7: Use EDM matrices to calculte cluster-wise distances to cellular signals:
if(whichSTAGE == 7){  
  
  

  # now I want to see where the non-lymphnode cells are found on the UMAPs:
  setwd(OutputDirectory)
  
  
  
if(run.old.code==1){
  for(s in 1:amountofGATEsubsets  ){
    
    temp.df <- subset(masked.gated.cell.dat, GATEsubset == s )
    # watch out, the make.color.plot engine does not accept a full-range Phenograph_metacluster_palette if you just plot a subset of it.
    # And if you plot the heatmap with the full range palette, they get all listed in the legend. So lets create a temporary Palette subset:
    subset.Phenograph_metacluster_palette <- Phenograph_metacluster_palette[ which( 100*(s) < names(Phenograph_metacluster_palette) & names(Phenograph_metacluster_palette) < 100*(s+1)  ) ]
    
    subtitlestring <- paste0("Subset: ", GATEsubsetNames[s], " (",  nrow(temp.df), " cells plotted)" )
    
  
    # prepare the anchor for the labels:
    #centroidX = tapply( temp.df[["UMAP1"]] , temp.df[["Collapsed_metacluster"]], median)
    #centroidY = tapply(temp.df[["UMAP2"]], temp.df[["Collapsed_metacluster"]], median)
    #centroidCol = tapply(temp.df[["Collapsed_metacluster"]], temp.df[["Collapsed_metacluster"]], median)
    #centroidsDf <- data.frame(centroidX, centroidY, centroidCol)
    

    
  g<-  ggplot(data = temp.df, aes(x = UMAP1, y = UMAP2, colour = as.factor(Collapsed_metacluster))) + 
     # geom_point( alpha = 1., pch=".")+#size = 0.9,
      geom_point(size = 1.2, alpha = 0.05)+
      scale_colour_manual( values = subset.Phenograph_metacluster_palette )+ # push the polychrome palette in
      guides(colour = guide_legend(override.aes = list(size=4, alpha=1),  #we also gonna override the dot size how its depicted in the legend
                                   byrow=TRUE) ) +
      geom_point(data= subset(temp.df, inside.lymphnode==0), size = 1.9, alpha = 1.0)+
      labs(title= paste0("UMAP - non-lymphnode cells highlighted")  , 
           subtitle=subtitlestring#, 
           #caption="Created by M.Barone"#, 
           # y="ROI y", 
           # x="ROI x"
      )+
      theme(panel.background = element_rect(fill = "white", 
                                            colour = "black", size = 0.5), 
            axis.title.x = element_text(color = "Black", size = 15), 
            axis.title.y = element_text(color = "Black", size = 15), 
            axis.text.x = element_text(color = "Black",  size = 13), 
            axis.text.y = element_text(color = "Black",  size = 13), 
            panel.border = element_rect(colour = "black",  fill = NA, size = 2), 
            plot.title = element_text(color = "Black", size = 16), 
            aspect.ratio = 1,
            legend.direction = "vertical", 
            legend.position = "right", 
            legend.text = element_text(size = 18), legend.title = element_blank()
      ) #+
  #  geom_label_repel(data = centroidsDf, 
  #                   hjust = "right", # this only takes effect initially and is lost for labels that are pulled...
  #                   force = 30, # repulsion between overlapping text labels. default 1
  #                   force_pull = 0.8, # attraction betweenlabel and datapoint. default 1
  #                   # nudge_x = 0.1*(Xmax-Xmin), 
  #                   #nudge_y = 0.05*(Ymax-Ymin), 
  #                   xlim = c(-Inf, Inf), #c(Xmin, Xmax), # either plots them also outside the graph or restrains the labels to be within the datarange. this assures that no label is cut off
  #                   ylim = c(-Inf, Inf), #c(Ymin,Ymax), # either plots them also outside the graph or restrains the labels to be within the datarange. this assures that no label is cut off
   #                  box.padding = 0.01, # additional padding around each text label 0.25 default
  #                   label.padding = 0.20, #0.25 default
  #                   point.padding = 0, # additional padding around each point
  #                   
  #                   min.segment.length = 0, # draw all line segments
  #                   segment.curvature = -0.1, #pos: more righthand, 0 straight, neg increase left-hand
  #                   segment.ncp = 3, # control points per curve
  #                   segment.angle = 20,
   #                  
   #                  aes(x = centroidX, y = centroidY, label = centroidCol), #, alpha = 0.5
   #                  fill = "white",
   ##                  col = "black", fontface = "bold",size = 4,
  #                   verbose = TRUE) 
    
    
  
  
    plot(g)
  
  filename <- paste0(GATEsubsetNames[s],"_no-lymphnodeCells.jpg")
  ggsave(filename,g, width = 20*plot.ratio, height = 20, units = "cm",dpi = 600)
  
  

  
  g<- ggplot(data = subset(temp.df, inside.lymphnode==0), aes(x = UMAP1, y = UMAP2, colour = distance )) + 
    # geom_point( alpha = 1., pch=".")+#size = 0.9,
    geom_point(size = 0.8, alpha = 0.5)+
    scico::scale_color_scico(palette = "roma", direction=1)+ # bilbao roma
    #scale_colour_manual( values = subset.Phenograph_metacluster_palette )+ # push the polychrome palette in
    #guides(colour = guide_legend(override.aes = list(size=4, alpha=1),  #we also gonna override the dot size how its depicted in the legend
    #                             byrow=TRUE) ) +
    geom_point(data= subset(temp.df, distance==0), size = 1.9, alpha = 1.0)+
    labs(title= paste0("UMAP - no-lymphnode cells distance to islets")  , 
         subtitle=subtitlestring#, 
         #caption="Created by M.Barone"#, 
         # y="ROI y", 
         # x="ROI x"
    )+
    theme(panel.background = element_rect(fill = "white", 
                                          colour = "black", size = 0.5), 
          axis.title.x = element_text(color = "Black", size = 15), 
          axis.title.y = element_text(color = "Black", size = 15), 
          axis.text.x = element_text(color = "Black",  size = 13), 
          axis.text.y = element_text(color = "Black",  size = 13), 
          panel.border = element_rect(colour = "black",  fill = NA, size = 2), 
          plot.title = element_text(color = "Black", size = 16), 
          aspect.ratio = 1,
          legend.direction = "vertical", 
          legend.position = "right", 
          legend.text = element_text(size = 18), legend.title = element_blank()
    )+
    geom_label_repel(data = centroidsDf, 
                     hjust = "right", # this only takes effect initially and is lost for labels that are pulled...
                     force = 30, # repulsion between overlapping text labels. default 1
                     force_pull = 0.8, # attraction betweenlabel and datapoint. default 1
                     # nudge_x = 0.1*(Xmax-Xmin), 
                     #nudge_y = 0.05*(Ymax-Ymin), 
                     xlim = c(-Inf, Inf), #c(Xmin, Xmax), # either plots them also outside the graph or restrains the labels to be within the datarange. this assures that no label is cut off
                     ylim = c(-Inf, Inf), #c(Ymin,Ymax), # either plots them also outside the graph or restrains the labels to be within the datarange. this assures that no label is cut off
                     box.padding = 0.01, # additional padding around each text label 0.25 default
                     label.padding = 0.20, #0.25 default
                     point.padding = 0, # additional padding around each point
                     
                     min.segment.length = 0, # draw all line segments
                     segment.curvature = -0.1, #pos: more righthand, 0 straight, neg increase left-hand
                     segment.ncp = 3, # control points per curve
                     segment.angle = 20,
                     
                     aes(x = centroidX, y = centroidY, label = centroidCol), #, alpha = 0.5
                     fill = "white",
                     col = "black", fontface = "bold",size = 4,
                     verbose = TRUE) 
    
    
    plot(g)
    
    filename <- paste0(GATEsubsetNames[s],"_no-lymphnodeCells_distancetoIslets.jpg")
    ggsave(filename,g, width = 20*plot.ratio, height = 20, units = "cm",dpi = 600)
    
    
    
    
    
    
    
    
    
    
    
  }#end run with s along the GATEsubsets in masked.gated.cell.dat
  
}# end protect running old code  
  
  
  
  
  
  
  
  
  
 ##---------- PLOT distance to MASK (new code as well. works now) ------------------ 
  
  # this is messy, but the code (and the object cluster.distance.dat) was written way before lymph node flags were invented, so we just keep the pre-defined df and subset it in from masked.gated.cell.dat:

  cluster.distance.dat <- subset(masked.gated.cell.dat, inside.lymphnode==0)
  # we gotta tweak the x axis with the conditions:
  
  cluster.distance.dat$Diet_week  <- ordered(cluster.distance.dat$Diet_week,
                                                     levels = c("ND_1",   "ND_12", "ND_14",  "ND_24",  "ND_100", 
                                                                "WD_1",  "WD_12",  "WD_24"   
                                                                )) 
  # again, we take the second pos ctr out of the plottting
  cluster.distance.dat <-  subset(cluster.distance.dat, Diet_week %!in% "ND_14")
  
  

  
  
  setwd(OutputDataDirectory)
  fwrite( cluster.distance.dat , 'distance.gated.cell.dat.csv')
  
  setwd(OutputDirectory)
  # same game as before: make.color.plot engine does not accept a full-range Phenograph_metacluster_palette if you just plot a subset of it.
  # we gotta create a character vector again:  
  s<-2
  subset.Phenograph_metacluster_palette <- subset(annotated.Phenograph_metacluster_palette, Collapsed_metacluster > 100*(s) &  Collapsed_metacluster < 100*(s+1)   ) 
  temp.pal <- subset.Phenograph_metacluster_palette[,color] 
  names(temp.pal ) <-subset.Phenograph_metacluster_palette$Annotated_metacluster
  subset.Phenograph_metacluster_palette <- temp.pal
  rm(temp.pal)
  
  subset.Phenograph_metacluster_palette <-  subset.Phenograph_metacluster_palette[!duplicated(subset.Phenograph_metacluster_palette)] 
  # ^ this can now be pushed into the plotting engine:
  
  # define the offset of the brackets as you add more and more of them
  bracketoffset <- 0.5* max(subset(cluster.distance.dat, Annotated_metacluster>200 &GATEsubset == 2 )$islet.distance, na.rm=T) #0.91
    # lets kick off all NA in islet distance:
    clean.cluster.distance.dat <-cluster.distance.dat %>% drop_na("islet.distance")
  #Annotated_metacluster>200 &
  g<-
    ggplot(data=subset(cluster.distance.dat,   GATEsubset == 2 ) ,  aes(x=Diet_week, y=islet.distance, colour=Annotated_metacluster  ))+ # colour=factor(Collapsed.metacluster)
    
    
    geom_violin(trim=T)+
    geom_quasirandom(alpha=0.3,size=0.4) +
    scale_colour_manual( values = subset.Phenograph_metacluster_palette, limits = force  )+
    ylim(0,2800)+
    
    
    facet_wrap(~Annotated_metacluster)+ #, scales="free_y"
    theme_classic()+
    
    labs(#title=paste0("Size distribution of gated cells in ",d, " and ",w," weeks" )  , 
      title=paste0("Shortest distance of a non-lymph node cell to next islet" )  , 
      subtitle= paste0("Mask based on combined insulin/GLP/NKx6.1 signals, totally ", nrow(subset(cluster.distance.dat, Annotated_metacluster>200 &GATEsubset == 2 )), " gated cells in subset 2")  , 
      #caption="Created by M.Barone", 
      y= "Distance to islet mask [px]", 
      x=""
    )+
    theme(#axis.title.x = element_blank(), 
      legend.title=element_text(size=13), # turn off with element_blank(),
      #legend.text=element_text(size=11),
      axis.text=element_text(size=13),#,axis.text.x=element_text(size=10),
      axis.title.y=element_text(size=12), axis.title.x=element_text(size=15),
      legend.position="none", 
      # legend.key.width=unit(1,"line")
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    )+
    #scale_y_continuous(expand = c(0.2,0))+
    geom_signif(comparisons = list(c("WD_1", "WD_12"),  c("WD_12", "WD_24"),  
                                       c("ND_1", "ND_12"), c("ND_12", "ND_24"), c("ND_24", "ND_100")
                                   ), 
                #map_signif_level=TRUE, 
                map_signif_level=sigFunc,
                test=wilcox.test, #t.test,
                y_position = bracketoffset, 
                color='black',  
                textsize = 3,
                vjust=sigStar.verticaloffset)
    
  
  
  bracketoffset <- bracketoffset* 1.2
  g <- g +
    geom_signif(comparisons = list(c("ND_1", "WD_1")), 
                #map_signif_level=TRUE, 
                map_signif_level=sigFunc,
                test=wilcox.test, #t.test,
                y_position = bracketoffset, 
                color='black',  
                textsize = 3,
                vjust=sigStar.verticaloffset) #+
  
  bracketoffset <- bracketoffset* 1.2
  g <- g +
    
    geom_signif(comparisons = list(
      c("ND_12", "WD_12")), 
                #map_signif_level=TRUE, 
                map_signif_level=sigFunc,
                test=wilcox.test, #t.test, wilcox.test
                y_position = bracketoffset, 
                color='black',  
                textsize = 3,
                vjust=sigStar.verticaloffset) #+
  
  #c("ND_1", "WD_1")
  bracketoffset <- bracketoffset* 1.18
  g <- g +
    
    geom_signif(comparisons = list(c("ND_24", "WD_24"), c("ND_1", "ND_24")), 
                #map_signif_level=TRUE,
                map_signif_level=sigFunc,
                test=wilcox.test, #t.test,
                y_position = bracketoffset, 
                color='black',  
                textsize = 3,
                vjust=sigStar.verticaloffset) #+
  
  bracketoffset <- bracketoffset* 1.18
  g <- g +
    
    geom_signif(comparisons = list(c("ND_1", "ND_12"), c("ND_12", "ND_100"),c("WD_1", "WD_24")),
                #map_signif_level=TRUE,
                map_signif_level=sigFunc,
                test=wilcox.test, #t.test,
                y_position = bracketoffset, 
                color='black',  
                textsize = 2,
                vjust=sigStar.verticaloffset) 
  
  
 
  
  
 
  
  
  
  plot(g)  
  
  
  
  
  

 filename <- paste0("Distance_Clusters_noLymphnodeCells.png")
  ggsave(filename,g, width = 30, height = 20, units = "cm",dpi = 1200)
  
  
  
# table(subset(cluster.distance.dat, Diet_week %in% c("ND_1") & distance==0   )$ROI )
# BeCeFa_20201109_OT11_s0_p7_r5_a5_ac BeCeFa_20210312_OT13_s0_p3_r1_a1_ac BeCeFa_20211008_OT21_s0_p6_r1_a1_ac 
# 1                                   1                                   2 
# BeCeFa_20211015_OT22_s0_p7_r1_a1_ac BeCeFa_20211029_OT24_s0_p6_r1_a1_ac 
# 1                                   4 
# -> 24a1 has many islets in the ROI
 
#  table(subset(cluster.distance.dat, Diet_week %in% c("WD_24") & distance==0   )$ROI )   
#  BeCeFa_20210611_OT18_s0_p5_r5_a5_ac BeCeFa_20210618_OT20_s0_p5_r4_a4_ac BeCeFa_20211022_OT23_s0_p6_r6_a6_ac 
#  1                                   1                                  24 
#  BeCeFa_20211029_OT24_s0_p6_r6_a6_ac 
#  2 
# -> 23a6 provides 24 cells!!! but is has also a ton of islets in the ROI! 
  
  
  
  
  
  
  
  
## ---------------  NEW CODE ------------------  
  
  
# Plot the gated cells that are within the islet masks:
amount.of.cells.in.islets <- nrow(subset(masked.gated.cell.dat, inside.islet==1  ))
amount.of.cells.in.rims <- nrow(subset(masked.gated.cell.dat, inside.rim.islet==1  ))
#cluster.INislets.dat <- subset(cluster.distance.dat, distance==0 )
  

ClustersINisletsfrequency.data.table <- data.frame()


for(i in unique(all.annot.MC$Annotated_metacluster) ){
  
  temp.df <- subset(masked.gated.cell.dat, Annotated_metacluster== i & (inside.islet==1 | inside.rim.islet==1 )   )

  if( nrow(temp.df) >0){
    
    # process inside islet:
    temp.freq <- as.data.frame(table(  subset(temp.df, inside.islet==1 )$Diet_week ))
    
    if(nrow(temp.freq)>0){
    
    names(temp.freq) <- c("Cond", "Count")
    
    # we make two frequencies:
    # 1: whats the distribution among these Condtions?
    temp.freq <- cbind(temp.freq, Condition.freq.percent = round(100 * temp.freq$Count/ nrow( subset(temp.df, inside.islet==1 ) )    ,3))  
    # 2: whats the occurance of that cell relative to all cells in islets?
    temp.freq <- cbind(temp.freq, Cluster.freq.percent = round(100 * temp.freq$Count/amount.of.cells.in.islets,3)  )
    temp.freq <- cbind(temp.freq, Cluster = i)
    temp.freq <- cbind(temp.freq, Location = "islet")
    
    }else{
    # in case there is no cell of that cluster found inside the islet, we still have to create the temp.freq template:  
      temp.freq <- data.frame(Cond="any",
                                 Count=0, 
                                 Condition.freq.percent=0,
                                 Cluster.freq.percent=0,
                                 Cluster = i,
                                 Location= "islet"      )
      
    }
    
    # process islet rim:
    temp.freq2 <- as.data.frame(table(  subset(temp.df, inside.rim.islet==1 )$Diet_week ))
    
    if(nrow(temp.freq2)>0){
    names(temp.freq2) <- c("Cond", "Count")
    
    # we make two frequencies:
    # 1: whats the distribution among these Condtions?
    temp.freq2 <- cbind(temp.freq2, Condition.freq.percent = round(100 * temp.freq2$Count/ nrow( subset(temp.df, inside.rim.islet==1 ) )    ,3))  
    # 2: whats the occurance of that cell relative to all cells in islets?
    temp.freq2 <- cbind(temp.freq2, Cluster.freq.percent = round(100 * temp.freq2$Count/amount.of.cells.in.rims,3)  )
    temp.freq2 <- cbind(temp.freq2, Cluster = i)
    temp.freq2 <- cbind(temp.freq2, Location = "rim")
    
    }else{
      # in case there is no cell of that cluster found inside the islet, we still have to create the temp.freq template:  
      temp.freq2 <- data.frame(Cond="any",
                              Count=0, 
                              Condition.freq.percent=0,
                              Cluster.freq.percent=0,
                              Cluster = i,
                              Location= "rim"      )
      
    }
    
    
    temp.freq <- rbind(temp.freq,temp.freq2)
    rm(temp.freq2)
    
    
    
    
    ClustersINisletsfrequency.data.table <- rbind(ClustersINisletsfrequency.data.table,temp.freq)
    
  }
  
}

# next, we get rid of the the ones that had either rim or islet missing and this absence was tagged as "any" in the condition field:

ClustersINisletsfrequency.data.table <- subset(ClustersINisletsfrequency.data.table, Cond %!in% "any" )
# we also get rid of ND14:
ClustersINisletsfrequency.data.table <- subset(ClustersINisletsfrequency.data.table, Cond %!in% "ND_14" )
# and we get rid of the G1 clusters:
ClustersINisletsfrequency.data.table <- subset(ClustersINisletsfrequency.data.table, Cluster > 200 )


# and now we add more columns:
ClustersINisletsfrequency.data.table <- ClustersINisletsfrequency.data.table %>%
  mutate(Diet = case_when(
    startsWith(as.character(Cond), "ND") ~ "ND",
    startsWith(as.character(Cond), "WD") ~ "WD"
  ))

ClustersINisletsfrequency.data.table <-  ClustersINisletsfrequency.data.table %>%
  mutate(Week = case_when(
    endsWith(as.character(Cond), "_1") ~ 1,
    endsWith(as.character(Cond), "_12") ~ 12,
   # endsWith(as.character(Cond), "_14") ~ 14,
    endsWith(as.character(Cond), "_24") ~ 24,
    endsWith(as.character(Cond), "_100") ~ 100
    
  ))
 
# this is not needed at the moment 
#ClustersINisletsfrequency.data.table$Cluster <- factor(ClustersINisletsfrequency.data.table$Cluster,      # Reordering group factor levels
#                                                       levels = all.annot.MC$Annotated_metacluster ) 



#all.conditions.having.cells.in.islets <- unique(ClustersINisletsfrequency.data.table$Cond)
#isletcells.condions <- c()
#for(i in all.conditions.having.cells.in.islets){
#isletcells.condions <- c(isletcells.condions,  sum(subset(ClustersINisletsfrequency.data.table, Cond == i           )$Count) )
#  }
#names(isletcells.condions) <- all.conditions.having.cells.in.islets
#isletcells.condions[match( ClustersINisletsfrequency.data.table$Condition, names( isletcells.condions)   )]

#not run #ClustersINisletsfrequency.data.table[, "Cluster.freq.Condition"] <- ClustersINisletsfrequency.data.table$Count / isletcells.condions[match( ClustersINisletsfrequency.data.table$Cond names( isletcells.condions)   )]











g <- ggplot(ClustersINisletsfrequency.data.table, aes(x=factor(Week), y = Condition.freq.percent, fill=Diet )) + 
  facet_wrap(Cluster~Location )+
  geom_bar(stat="identity") +
  
  labs(title="Distribution of a cluster among conditions found within the islets/rim " , 
       subtitle=paste0(amount.of.cells.in.islets, " cells within islet mask, ",amount.of.cells.in.rims," cells within ",islet.rim.width,"px rim"), 
       #caption="Created by M.Barone", 
       y="Distribution of a cluster among condition [%]", 
       x=""
  )+
  
  theme(panel.background = element_rect(fill = "grey90", colour = "grey55", size = 0.5), # change 'colour' to black for informative axis
        axis.title.x=element_text(color="grey15", size=11),
        axis.title.y=element_text(color="grey15", size=11),
        axis.text=element_text(size=10),
        #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
        legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
        legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
        #legend.title=element_blank(),
        plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
        plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
  )

plot(g)

filename <- paste0("ClusterDistribution_islet_and",islet.rim.width,"px rim.png")
ggsave(filename,g, width = 30, height = 30, units = "cm",dpi = 1200)





g <- ggplot(ClustersINisletsfrequency.data.table, aes(x=factor(Week), y = Cluster.freq.percent, fill=Diet )) + 
  facet_wrap(Cluster~Location )+
  geom_bar(stat="identity") +
  
  labs(title="Cluster frequencies rel. to all cells found within the islets/rim" , 
       subtitle=paste0(amount.of.cells.in.islets, " cells within islet mask, ",amount.of.cells.in.rims," cells within ",islet.rim.width,"px rim"), 
       #caption="Created by M.Barone", 
       y="Frequency of a cluster [%]", 
       x=""
  )+
  
  theme(panel.background = element_rect(fill = "grey90", colour = "grey55", size = 0.5), # change 'colour' to black for informative axis
        axis.title.x=element_text(color="grey15", size=11),
        axis.title.y=element_text(color="grey15", size=11),
        axis.text=element_text(size=10),
        #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
        legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
        legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
        #legend.title=element_blank(),
        plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
        plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
  )

plot(g)

filename <- paste0("ClusterOccurence_islet_and",islet.rim.width,"px rim.png")
ggsave(filename,g, width = 30, height = 30, units = "cm",dpi = 1200)









#----------------------- DENSITIES ISLETS and RIM --------------------------------------
#----- calculate cluster densities on mice from scratch.

ClusterDensity.permouse.data.table <- data.frame()
#i <- unique(masked.gated.cell.dat$Mouse_ID)[3]
for(i in unique(masked.gated.cell.dat$Mouse_ID) ){
  
  temp.df <- subset(masked.gated.cell.dat, Mouse_ID %in% i  ) 
  
  temp.totalSegmentationArea <- sum( subset(all.cells,  Mouse_ID %in% i & inside.lymphnode==0 )$Area) 
  temp.isletSegmentationArea <- sum( subset(all.cells,  Mouse_ID %in% i & inside.lymphnode==0 & inside.islet == 1)$Area)
  temp.rimSegmentationArea <- sum( subset(all.cells,  Mouse_ID %in% i & inside.lymphnode==0 & inside.rim.islet == 1)$Area)
  
  temp.freq <-  as.data.frame(table(temp.df$Annotated_metacluster))
  names(temp.freq) <- c("MC", "Cells")
  
  
  # at that point, MC is factor:
  temp.freq$MC <- as.character(temp.freq$MC)
  
  # we gonna put more lines in here, namely all clusters that are absent:
  
  temp.absent <- data.frame( "MC" = setdiff(all.annot.MC$Annotated_metacluster, temp.freq$MC),
                             "Cells"= rep(0, length(   setdiff(all.annot.MC$Annotated_metacluster, temp.freq$MC)  )      ) 
  )
  temp.freq <- rbind(temp.freq,temp.absent)
  
  

  
# nrow( subset(temp.df, Annotated_metacluster %in% temp.freq$MC[22] & inside.islet == 1)  )
# temp.df$Annotated_metacluster == temp.freq$MC[22]
#sapply(temp.freq$MC[22], function(x) sum(temp.df$Annotated_metacluster == x & temp.df$inside.islet == 1 , na.rm=T )    ) 
  
# now we gonna scan through all clusters in temp.freq and filter out how many of the "Cells" are within islet and within rim:
# watch out, there is a chance that some mice have ROIs with nasty islet masks, that caused NAs in the flags, you need to take care of that when sapplying:
temp.freq <- cbind(temp.freq, inside.islet =  sapply(temp.freq$MC, function(x) sum(temp.df$Annotated_metacluster == x & temp.df$inside.islet == 1, na.rm=T )))  
temp.freq <- cbind(temp.freq, inside.rim.islet =  sapply(temp.freq$MC, function(x) sum(temp.df$Annotated_metacluster == x & temp.df$inside.rim.islet == 1 , na.rm=T ))) 

  
  # now we perform some normalizations:
# 1) all cells found on that mouse to the total segementation area of all ROIs excluding lymph node area:
  #temp.freq <- cbind(temp.freq, totClusterDensityMOUSE = temp.freq$Cells/temp.totalSegmentationArea     )
  temp.freq <- cbind(temp.freq, totClusterCellsMOUSE_per100sqrpx = temp.freq$Cells/temp.totalSegmentationArea * 10000    )
  
# 2) cells found inside the islet mask. we also write out the amount of islets per ROI, this would be an idea to normalize here...
  #temp.freq <- cbind(temp.freq, isletClusterDensityMOUSE = temp.freq$inside.islet/temp.isletSegmentationArea     )
  temp.freq <- cbind(temp.freq, isletClusterCellsMOUSE_per100sqrpx = temp.freq$inside.islet/temp.isletSegmentationArea * 10000    )    
 

  # 3) cells found in the rim area: if they pile up there, this density is elevated rel. to the totClusterDensityMOUSE
  #temp.freq <- cbind(temp.freq, rimClusterDensityMOUSE = temp.freq$inside.rim.islet/temp.rimSegmentationArea     )
  temp.freq <- cbind(temp.freq, rimClusterCellsMOUSE_per100sqrpx = temp.freq$inside.rim.islet/temp.rimSegmentationArea * 10000    )
  # so lets see if there is a fold increase in that area:
  temp.freq <- cbind(temp.freq, foldincrease.rim.vs.total = temp.freq$rimClusterCellsMOUSE_per100sqrpx/temp.freq$totClusterCellsMOUSE_per100sqrpx    )
  
  
  temp.freq <- cbind(temp.freq, Mouse_ID = unique(temp.df$Mouse_ID)   )
  temp.freq <- cbind(temp.freq, Cond = unique(temp.df$Diet_week) )
  
  ClusterDensity.permouse.data.table <- rbind(ClusterDensity.permouse.data.table,temp.freq)    
  rm(temp.freq,temp.absent)
  
}

# now we get rid of ND14 and sort cond
ClusterDensity.permouse.data.table$Cond <- factor(ClusterDensity.permouse.data.table$Cond,
                                                  levels = all.condtions #c(all.condtions[1],all.condtions[3:5],all.condtions[2],all.condtions[6:8] )         #c(all.condtions         ) # 
)


# Birgit wants to get rid of the unneccessary second control ND14 since its just distrubing to find this control in the plots:

ClusterDensity.permouse.data.table <-  subset(ClusterDensity.permouse.data.table, Cond %!in% "ND_14")


# Birgit has another idea. we plot them temporally and make linear fits:
# for that we need to split up cond into two more columns:
ClusterDensity.permouse.data.table <- ClusterDensity.permouse.data.table %>%
  mutate(Diet = case_when(
    startsWith(as.character(Cond), "ND") ~ "ND",
    startsWith(as.character(Cond), "WD") ~ "WD"
  ))

ClusterDensity.permouse.data.table <-  ClusterDensity.permouse.data.table %>%
  mutate(Week = case_when(
    endsWith(as.character(Cond), "_1") ~ 1,
    endsWith(as.character(Cond), "_12") ~ 12,
    #endsWith(as.character(Cond), "_14") ~ 14,
    endsWith(as.character(Cond), "_24") ~ 24,
    endsWith(as.character(Cond), "_100") ~ 100
    
  ))



#--------------- PLOT Cluster density in Islet and (all ROIs of each mouse summed up) [cells/total px or ROI] ---------------

# watch out, the area.totals is the dimension of the raster image, and not the total cell segmentation area
# so dont use that!
#area.totals <- fread('area.totals.csv')

if(lymphnode.extension>0){
  subtitlestring <- paste0("Using cells within ",islet.rim.width,"px rim around istlet and lymph node mask extended by ",lymphnode.extension, "px")
}else{
  subtitlestring <- paste0("Using cells within ",islet.rim.width,"px rim around istlet. Lymph node mask not extended")
}

# we got these here:    totClusterCellsMOUSE_per100sqrpx  isletClusterCellsMOUSE_per100sqrpx    rimClusterCellsMOUSE_per100sqrpx 

# run through all normalizations in a loop and plot:
normalizations.to.do <- c("totClusterCellsMOUSE_per100sqrpx",  
                          "isletClusterCellsMOUSE_per100sqrpx", 
                          "rimClusterCellsMOUSE_per100sqrpx",
                          "foldincrease.rim.vs.total")

titlestrings <- c(  "G2 cluster density over ROI. Mouse-wise summed up ROIs",
                   "G2 cluster density within islet. Mouse-wise summed up ROIs",
                   "G2 cluster density in rim area of islet. Mouse-wise summed up ROIs",
                   "G2 cluster density increase in rim vs. total non-lymphnode area. Mouse-wise summed up ROIs"
                   )

yaxisstrings <- c( "Cluster density [cells/10000 px^2 non-lymphnode area]",
                   "Cluster density [cells/10000 px^2 islet area]",
                   "Cluster density [cells/10000 px^2 rim area]",
                   "Fold increase"
  
  
)


filestrings <- c("Clusterdesity_Total_",
                 "Clusterdesity_Islet_",
                 "Clusterdesity_Rim_",
                 "Clusterdesity_foldincreaseinRim_"
                 )

for(x in 1:4){
  
  i <- normalizations.to.do[x]
  titlestring <- titlestrings[x]
  filestring <- filestrings[x]
  yaxisstring <- yaxisstrings[x]
  
  bracketoffset <- 0.4*max( subset(ClusterDensity.permouse.data.table, MC>200)[, i], na.rm=T )
  print(bracketoffset)
  

# watch out, we print only the G2 here, the other is not gated for CD45 and doestn make any sense 
g <- ggplot( subset(ClusterDensity.permouse.data.table, MC>200), aes(x=as.factor(Cond), y = subset(ClusterDensity.permouse.data.table, MC>200)[,i] , colour=as.factor(Cond) )) + # the bar wants fill=Cond
  facet_wrap(~MC , scales="free_y")+
  # geom_bar(stat="identity") +
  geom_quasirandom(size=1.5) + #alpha=0.8
  scale_color_manual(values = cond_palette)+
  
  labs(title=titlestring , 
       subtitle=subtitlestring, 
       #caption="Created by M.Barone", 
       y=yaxisstring, 
       x=""
  )+
  
  theme(panel.background = element_rect(fill = "grey90", colour = "grey55", size = 0.5), # change 'colour' to black for informative axis
        axis.title.x=element_text(color="grey15", size=11),
        axis.title.y=element_text(color="grey15", size=11),
        axis.text=element_text(size=8),
        axis.text.x = element_text(angle = 45, hjust=1),
        #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
        legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
        legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
        #legend.title=element_blank(),
        legend.position="none",
        plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
        plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
  )+
  
  scale_y_continuous(expand = c(0.2,0))+
  geom_signif(comparisons = list(c("ND_1", "WD_1")), 
              #map_signif_level=TRUE, 
              map_signif_level=sigFunc,
              test=wilcox.test, #t.test,
              y_position = bracketoffset, 
              color='black',  
              textsize = 3,
              vjust=sigStar.verticaloffset) #+

bracketoffset <- bracketoffset* 1.2
g <- g +
  
  geom_signif(comparisons = list(c("ND_12", "WD_12")), 
              #map_signif_level=TRUE, 
              map_signif_level=sigFunc,
              test=wilcox.test, #t.test, wilcox.test
              y_position = bracketoffset, 
              color='black',  
              textsize = 3,
              vjust=sigStar.verticaloffset) #+


bracketoffset <- bracketoffset* 1.18
g <- g +
  
  geom_signif(comparisons = list(c("ND_24", "WD_24")), 
              #map_signif_level=TRUE,
              map_signif_level=sigFunc,
              test=wilcox.test, #t.test,
              y_position = bracketoffset, 
              color='black',  
              textsize = 3,
              vjust=sigStar.verticaloffset) #+

bracketoffset <- bracketoffset* 1.18
g <- g +
  
  geom_signif(comparisons = list(c("ND_1", "ND_12")), 
              #map_signif_level=TRUE,
              map_signif_level=sigFunc,
              test=wilcox.test, #t.test,
              y_position = bracketoffset, 
              color='black',  
              textsize = 2,
              vjust=sigStar.verticaloffset) +
  geom_signif(comparisons = list(c("ND_12", "ND_24")), 
              #map_signif_level=TRUE,
              map_signif_level=sigFunc,
              test=wilcox.test, #t.test,
              y_position = bracketoffset, 
              color='black',  
              textsize = 2,
              vjust=sigStar.verticaloffset)+
  geom_signif(comparisons = list(c("ND_24", "ND_100")), 
              #map_signif_level=TRUE,
              map_signif_level=sigFunc,
              test=wilcox.test, #t.test,
              y_position = bracketoffset, 
              color='black',  
              textsize = 2,
              vjust=sigStar.verticaloffset)+
  geom_signif(comparisons = list(c("WD_1", "WD_24")), 
              #map_signif_level=F,
              map_signif_level=sigFunc,
              test=wilcox.test, #t.test,
              y_position = bracketoffset, 
              color='black',  
              textsize = 2,
              vjust=sigStar.verticaloffset)

bracketoffset <- bracketoffset* 1.18
g <- g +
  
  geom_signif(comparisons = list(c("WD_1", "WD_12")), 
              #map_signif_level=TRUE,
              map_signif_level=sigFunc,
              test=wilcox.test, #t.test,
              y_position = bracketoffset, 
              color='black',  
              textsize = 2,
              vjust=sigStar.verticaloffset) +
  geom_signif(comparisons = list(c("WD_12", "WD_24")), 
              #map_signif_level=TRUE,
              map_signif_level=sigFunc,
              test=wilcox.test, #t.test,
              y_position = bracketoffset, 
              color='black',  
              textsize = 2,
              vjust=sigStar.verticaloffset)+
  geom_signif(comparisons = list(c("ND_1", "ND_100")), 
              #map_signif_level=TRUE,
              map_signif_level=sigFunc,
              test=wilcox.test, #t.test,
              y_position = bracketoffset, 
              color='black',  
              textsize = 2,
              vjust=sigStar.verticaloffset)





plot(g)




filename <- paste0(filestring,"ftest",lymphnode.extension,"pxLymphmaskextend_newcalculated.png")
ggsave(filename,g, width = 30, height = 25, units = "cm",dpi = 800)




g <- ggplot( subset(ClusterDensity.permouse.data.table, MC>200 ), aes(x=Week, y = subset(ClusterDensity.permouse.data.table, MC>200)[,i], group=Diet , colour=Diet )) + # the bar wants fill=Cond
  #ggplot(subset(ClusterDensity.permouse.data.table, MC>200 & Week<100), aes(x=Week, y = ClusterCellsMOUSE_100sqrpx, group=Diet, colour=Diet ))+ #as.factor(Cond)
  geom_point()+
  facet_wrap(~MC , scales="free_y")+
  scale_color_manual(values = diet_palette)+
  geom_smooth(method='lm')+
  labs(title=titlestring , 
       subtitle=subtitlestring, 
       #caption="Created by M.Barone", 
       y=yaxisstring, 
       x=""
  )+
  
  theme(panel.background = element_rect(fill = "grey90", colour = "grey55", size = 0.5), # change 'colour' to black for informative axis
        axis.title.x=element_text(color="grey15", size=11),
        axis.title.y=element_text(color="grey15", size=11),
        axis.text=element_text(size=8),
        axis.text.x = element_text(angle = 45, hjust=1),
        #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
        legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
        legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
        #legend.title=element_blank(),
        legend.position="none",
        plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
        plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
  )
filename <- paste0(filestring,"_timecourse_",lymphnode.extension,"pxLymphmaskextend_newcalculated.png")
ggsave(filename,g, width = 30, height = 25, units = "cm",dpi = 800)


#watch out its not ordered, and and ND100 on second pos.
}#end run along all normalizations for densities











# For the cluster density change over time, we need to fit each cluster separately, split by diet:

ClusterDensity.fit.linear <- data.frame()

cluster.densities.to.fit <- unique(subset(ClusterDensity.permouse.data.table, MC>200 )$MC)

for(f in cluster.densities.to.fit ){
for(d in diet.vector){  
temp.df <-  subset(ClusterDensity.permouse.data.table, MC==f & Diet==d)
  
#rimClusterCellsMOUSE_per100sqrpx totClusterCellsMOUSE_per100sqrpx
fit.linear = lm(rimClusterCellsMOUSE_per100sqrpx~Week, data = temp.df)  
 summary(fit.linear) 
 
#only extract the coefficients if the fit actually worked
 if( is.character(fit.linear) == FALSE ){
   
   
   intercept <- coef(summary(fit.linear))[1, 1]
   intercept.SE <- coef(summary(fit.linear))[1, 2]
   slope <- coef(summary(fit.linear))[2, 1]
   slope.SE <- coef(summary(fit.linear))[2, 2]
   slope.Pr <- coef(summary(fit.linear))[2, 4]
   df <- df.residual(fit.linear)
   
   ClusterDensity.fit.linear <- rbind(   ClusterDensity.fit.linear, c(unique(temp.df$MC),unique(temp.df$Diet),   slope,slope.SE,slope.Pr,df  ))
   
 }
 

 
  
  g <- ggplot(temp.df, aes(x=Week, y = totClusterCellsMOUSE_per100sqrpx, group=Diet , colour=rimClusterCellsMOUSE_per100sqrpx )) + # the bar wants fill=Cond
    #ggplot(subset(ClusterDensity.permouse.data.table, MC>200 & Week<100), aes(x=Week, y = ClusterCellsMOUSE_100sqrpx, group=Diet, colour=Diet ))+ #as.factor(Cond)
    geom_point()+
    geom_abline(intercept = intercept, slope = slope)+
    
   # scale_color_manual(values = diet_palette)+
 
   labs(title=paste0(f," in Diet ",diet.vector) , 
        subtitle=paste0("Slope ",round(slope,5), ", p-val on slope ",round(slope.Pr,5)), 
  #       #caption="Created by M.Barone", 
         y="Cells per mouse per 10000px", 
  #       x=""
    )+
    
    theme(panel.background = element_rect(fill = "grey90", colour = "grey55", size = 0.5), # change 'colour' to black for informative axis
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=8),
          axis.text.x = element_text(angle = 45, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          legend.position="right",
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )
  
 
  
  plot(g)
  
}# end run with d along both diets  
}#end run with f along all clusters to fit the densities

# this is prone for error. watch out that when you change the rbind loop, change the names here as well:
names(ClusterDensity.fit.linear) <- c("Annotated_cluster", "Diet", "slope", "slope.SE", "slope.Pr", "DOF"   )
# and then the stuff is all character:
ClusterDensity.fit.linear$slope <- as.numeric(ClusterDensity.fit.linear$slope)
ClusterDensity.fit.linear$slope.SE <- as.numeric(ClusterDensity.fit.linear$slope.SE)
ClusterDensity.fit.linear$slope.Pr <- as.numeric(ClusterDensity.fit.linear$slope.Pr)
ClusterDensity.fit.linear$DOF <- as.numeric(ClusterDensity.fit.linear$DOF)

# now we gonna run through once again and check for delta.slope. 
# this works only if every fit worked!!!!
slope.lyrics <-  paste("\n\n")

for(f in cluster.densities.to.fit ){
  
  temp.df <-  subset(ClusterDensity.fit.linear, Annotated_cluster==f)    

    
  
  
DDG = temp.df$slope[1]-  temp.df$slope[2]
#print(paste("Slope diff -> ", round(DDG,2), " cell.dens/week"))

SE.DDG = sqrt(temp.df$slope.SE[1]^2 + temp.df$slope.SE[2]^2)
#print(paste("SE.slope.diff -> ", round(SE.DDG,3), " cell.dens/week" ))

FG = (temp.df$slope.SE[1]^2 + temp.df$slope.SE[2]^2)^2/(temp.df$slope.SE[1]^4/(temp.df$DOF[1]-1.0) + temp.df$slope.SE[2]^4/(temp.df$DOF[1]-1.0))
#print(paste("FG -> ", round(FG) ))

t.quantile = qt(0.975,round(FG))
#print(paste("t-Wert -> ", abs(DDG/sqrt(SE.DDG))))
#print(paste("t-Quantile -> ", t.quantile))

#print(paste("95% CI slope.diff -> ", round( DDG-t.quantile*SE.DDG, 3) , " - ", round( DDG + t.quantile*SE.DDG, 3), " cell.dens/week" ))


#cat("\n Monta-Carlo Simulation Version 1\n")
#maxit =1000000
#DDG.random = sort(rnorm(maxit,temp.df$slope[1],temp.df$slope.SE[1]) - rnorm(maxit,temp.df$slope[1],temp.df$slope.SE[2]))

#lower = which.max(pmin(DDG.random,0)) - 1

#SE.DDG = sd(DDG.random)
#print(paste("SE.DDG -> ", round(SE.DDG,3), " kJ/mol") )
#print(paste("95% CI DDG -> ", round( DDG - 1.96*SE.DDG, 3) , " bis ", round( DDG + 1.96*SE.DDG, 3), " kJ/mol" ))
#print(paste("lower -> ", round(100.0*lower/maxit,3), " ; upper -> ", round(100.0*(maxit-lower)/maxit,3) ))
#hist(DDG.random)

slope.lyrics <-  paste0(slope.lyrics, unique(temp.df$Annotated_cluster)," 95% CI slope.diff -> ", round( DDG-t.quantile*SE.DDG, 7) , " - ", round( DDG + t.quantile*SE.DDG, 3), " cell.dens/week\n")


}#end run with f along all clusters to fit the densities

cat(slope.lyrics)









# birgit wants to correlate activ. eff with the F4.80neg by mouse, and also the F4.80low, to make sure that they are dragged together into the rim and the overall ROI:

temp.df <- subset(ClusterDensity.permouse.data.table, MC %in% c("activ. effector-like CD8+", "F4.80- macrophages", "F4.80lo macrophages"  )& Week==24)

pairwise.wilcox.test(temp.df$totClusterCellsMOUSE_per100sqrpx, temp.df$MC, p.adj="bonferroni")
















s<-2
subset.Phenograph_metacluster_palette <- subset(annotated.Phenograph_metacluster_palette, Collapsed_metacluster > 100*(s) &  Collapsed_metacluster < 100*(s+1)   ) 
temp.pal <- subset.Phenograph_metacluster_palette[,color] 
names(temp.pal ) <-subset.Phenograph_metacluster_palette$Annotated_metacluster
subset.Phenograph_metacluster_palette <- temp.pal
rm(temp.pal)

subset.Phenograph_metacluster_palette <-  subset.Phenograph_metacluster_palette[!duplicated(subset.Phenograph_metacluster_palette)] 
# ^ this can now be pushed into the plotting engine:


g <- 
ggplot(subset(ClusterDensity.permouse.data.table, MC>200 ) , aes(x=factor(Cond), y = inside.islet, fill=MC )) + 
  geom_bar(stat="identity") +
  scale_fill_manual( values = subset.Phenograph_metacluster_palette  ,
                    limits = force 
  )+
  
  labs(title="Absolute count of every cluster found within the islets" , 
       subtitle=paste0(sum(ClusterDensity.permouse.data.table$inside.islet), " gated cells within islet masks"), 
       #caption="Created by M.Barone", 
       y="Summed up cells per cluster", 
       x=""
  )+
  
  #i wanna count them one by one
  scale_y_continuous(breaks = seq(0, 100, by = 1)) +
  
  theme(panel.background = element_rect(fill = "grey90", colour = "grey55", size = 0.5), # change 'colour' to black for informative axis
        axis.title.x=element_text(color="grey15", size=11),
        axis.title.y=element_text(color="grey15", size=11),
        axis.text=element_text(size=10),
        #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
        legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
        legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
        #legend.title=element_blank(),
        plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
        plot.subtitle = element_text(color = "Black", size = 12, hjust = 0),
        panel.grid.minor = element_blank()
  )


plot(g)
filename <- "ClusterCounts_withinIslet_byConditionAndWeek.png"
ggsave(filename,g, width = 30, height = 30, units = "cm",dpi = 1200)








g <- 
  ggplot(subset(ClusterDensity.permouse.data.table, MC>200 ) , aes(x=factor(Cond), y = inside.rim.islet, fill=MC )) + 
  geom_bar(stat="identity") +
  scale_fill_manual( values = subset.Phenograph_metacluster_palette  ,
                     limits = force 
  )+
  
  labs(title="Absolute count of every cluster found within the islet rim" , 
       subtitle=paste0(sum(ClusterDensity.permouse.data.table$inside.rim.islet), " gated cells within ",islet.rim.width,"px rim around istlet rims"), 
       #caption="Created by M.Barone", 
       y="Summed up cells per cluster", 
       x=""
  )+
  
  #i wanna count them one by one
  scale_y_continuous(breaks = seq(0, 600, by = 20)) +
  
  theme(panel.background = element_rect(fill = "grey90", colour = "grey55", size = 0.5), # change 'colour' to black for informative axis
        axis.title.x=element_text(color="grey15", size=11),
        axis.title.y=element_text(color="grey15", size=11),
        axis.text=element_text(size=10),
        #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
        legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
        legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
        #legend.title=element_blank(),
        plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
        plot.subtitle = element_text(color = "Black", size = 12, hjust = 0),
        panel.grid.minor = element_blank()
  )


plot(g)
filename <- paste0("ClusterCounts_within",islet.rim.width,"px rim_byConditionAndWeek.png")
ggsave(filename,g, width = 30, height = 30, units = "cm",dpi = 1200)









# we also want to plot channel siganal raw or transf split by macs or condtions:

# if you re-load the masked.gated.cell.dat. the order is gone:

masked.gated.cell.dat$Diet_week  <- ordered(masked.gated.cell.dat$Diet_week,
                                           levels = c("ND_1",   "ND_12", "ND_14",  "ND_24",  "ND_100", 
                                                      "WD_1",  "WD_12",  "WD_24"   
                                           )) 
# again, we take the second pos ctr out of the plottting
masked.gated.cell.dat <-  subset(masked.gated.cell.dat, Diet_week %!in% "ND_14")




### print channel sigals of CD11b and CD69:
#someshs.channels <- c("CD69_Sm154_t2_z", "CD11b_Nd150_t2_z")

setwd(OutputDirectory)
# you could take the marker.order here for the z-scored, but I wanna see the raw signal, so:
for(m in cellular.cols ){
  
  bracketoffset <- 0.5 * max(subset( masked.gated.cell.dat, GATEsubset==2)[[m]])
  
  g <- 
    ggplot( subset( masked.gated.cell.dat, GATEsubset==2)  , aes( x=as.factor(Diet_week), y= subset( masked.gated.cell.dat, GATEsubset==2)[[m]] , color = Diet_week )  ) +
    
    
    #geom_violin(trim=T)+
    geom_quasirandom(size=0.5,alpha=0.8) + #
    geom_boxplot(outlier.shape = NA, lwd=0.8,alpha=0.7) + 
    scale_color_manual(values = cond_palette[c(1:2,4:8)])+
    
    # geom_area( alpha = 0.15 ,
    #             stat = "bin",
    #             binwidth = function(x) 2 * IQR(x) / (length(x)^(1/3)), 
    #             size = 0.5    )+
    
    # scico::scale_fill_scico(palette = "roma", direction=-1)+ # bilbao
    #scale_color_manual(values = scico(length(all_images), begin = 0.12, palette = "roma"))+
    #scale_fill_manual(values = scico(length(all_images), begin = 0.12, palette = "roma"))+
    #facet_wrap(~Diet, scales = "free_y")+
    
  #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),  labels = trans_format("log10", math_format(10^.x))) +
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),         labels = trans_format("log10", math_format(10^.x))) +
  #annotation_logticks(sides="l")   +
  
  theme_classic()+
    
    labs(title= paste0(sub("_.*", "", m ), " channel signal (f.test)"  )    , 
         subtitle= paste0("Using ", nrow(subset( masked.gated.cell.dat, GATEsubset==2)), " gated non-lymphnode cells in G2, cells split by diet.") ,
         #caption="Created by M.Barone", 
         y="Channel signal [raw]", 
         x= sub("_.*", "", m )  
    )+
    
    theme(panel.background = element_rect(fill = "grey93", colour = "grey55", size = 0.5), # change 'colour' to black for informative axis
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=10),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )+
    
    geom_signif(comparisons = list(c("ND_1", "WD_1")), 
                #map_signif_level=TRUE, 
                map_signif_level=sigFunc,
                test=wilcox.test, #t.test, wilcox.test
                y_position = bracketoffset, 
                color='black',  
                textsize = 3,
                vjust=sigStar.verticaloffset) #+
  
  bracketoffset <- bracketoffset*1.1
  g <- g +
    
    geom_signif(comparisons = list(c("ND_12", "WD_12")), 
                #map_signif_level=TRUE, 
                map_signif_level=sigFunc,
                test=wilcox.test, #t.test,
                y_position = bracketoffset, 
                color='black',  
                textsize = 3,
                vjust=sigStar.verticaloffset) #+
  
  
  bracketoffset <- bracketoffset*1.1
  g <- g +
    
    geom_signif(comparisons = list(c("ND_24", "WD_24")), 
                #map_signif_level=TRUE,
                map_signif_level=sigFunc,
                test=wilcox.test, #t.test,
                y_position = bracketoffset, 
                color='black',  
                textsize = 3,
                vjust=sigStar.verticaloffset) #+
  
  
  bracketoffset <- bracketoffset*1.2
  g <- g +
    
    geom_signif(comparisons = list(c("ND_1", "ND_12"),
                                   c("ND_12", "ND_24"),
                                   c("ND_24", "ND_100"),
                                   c("WD_1", "WD_12"),
                                   c("WD_12", "WD_24")
                                   ), 
                #map_signif_level=TRUE,
                map_signif_level=sigFunc,
                test=wilcox.test, #t.test,
                y_position = bracketoffset, 
                color='black',  
                textsize = 3,
                vjust=sigStar.verticaloffset) #+
  
  
  plot(g)
  
  filename <- paste0("Channel_", sub("_.*", "", m ) , "_byCondition.png")
  ggsave(filename,g, width = 30, height = 20, units = "cm",dpi = 600)
  
}



# for the clustes, we only wanna see different expression levels in our four:

# we are only interested in spatial plots of 4 clusters:
temp.df <- subset(masked.gated.cell.dat, GATEsubset==2 & Annotated_metacluster %in% c("activ. effector-like CD8+",
                                                  "F4.80lo macrophages",
                                                  "F4.80- macrophages",
                                                  "F4.80hi macrophages"))


# and now we only bring out the 4 colors:
subset.Phenograph_metacluster_palette <- subset(annotated.Phenograph_metacluster_palette, Annotated_metacluster%in% c("activ. effector-like CD8+",
                                                                                                                      "F4.80lo macrophages",
                                                                                                                      "F4.80- macrophages",
                                                                                                                      "F4.80hi macrophages"))[,color] 
names(subset.Phenograph_metacluster_palette ) <- c("activ. effector-like CD8+",
                                                   "F4.80lo macrophages",
                                                   "F4.80- macrophages",
                                                   "F4.80hi macrophages")




#install.packages('nortest')
#library(nortest)
#ad.test(  subset( temp.df, Annotated_metacluster==c("mDCs") )[[m]])

for(m in cellular.cols ){
  bracketoffset <- 0.9*max(temp.df[[m]])
  
  g <- 
    ggplot( temp.df  , aes( x=as.factor(Annotated_metacluster), y= temp.df[[m]] , color = Annotated_metacluster )  ) +
  
  
    
    #geom_violin(trim=T)+
    geom_quasirandom(size=0.5,alpha=0.8) + #
    geom_boxplot(outlier.shape = NA, lwd=0.8,alpha=0.7) + 
    #  scale_color_manual(values = cond_palette)+
    
    scale_color_manual( values = subset.Phenograph_metacluster_palette 
                        
    )+
    facet_wrap(~Diet_week)+
    # geom_area( alpha = 0.15 ,
    #             stat = "bin",
    #             binwidth = function(x) 2 * IQR(x) / (length(x)^(1/3)), 
    #             size = 0.5    )+
    
    # scico::scale_fill_scico(palette = "roma", direction=-1)+ # bilbao
    #scale_color_manual(values = scico(length(all_images), begin = 0.12, palette = "roma"))+
    #scale_fill_manual(values = scico(length(all_images), begin = 0.12, palette = "roma"))+
    #facet_wrap(~Diet, scales = "free_y")+
  
  #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),  labels = trans_format("log10", math_format(10^.x))) +
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),         labels = trans_format("log10", math_format(10^.x))) +
  #annotation_logticks(sides="l")   +
  
  theme_classic()+
    
    labs(title= paste0(sub("_.*", "", m ), " channel signal (f.test)"  )    , 
         subtitle= paste0("Using ", nrow(temp.df), " gated non-lymphnode cells in G2.") ,
         #caption="Created by M.Barone", 
         y="Channel signal [transf., z-score]", 
         x= sub("_.*", "", m )  
    )+
    
    theme(panel.background = element_rect(fill = "grey93", colour = "grey55", size = 0.5), # change 'colour' to black for informative axis
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=10),
          axis.text.x = element_text(angle = 45, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )+
    
 #   "activ. effector-like CD8+",

    
    geom_signif(comparisons = list(c("activ. effector-like CD8+", "F4.80- macrophages"),
                                   c("F4.80- macrophages", "F4.80hi macrophages"),
                                   c("F4.80hi macrophages", "F4.80lo macrophages")
                                   
                                   
                                   ), 
                #map_signif_level=TRUE, 
                map_signif_level=sigFunc,
                test=wilcox.test, #t.test, wilcox.test
                y_position = bracketoffset, 
                color='black',  
                textsize = 3,
                vjust=sigStar.verticaloffset) 
  
  
  bracketoffset <- bracketoffset*1.1
  g <- g +
    geom_signif(comparisons = list(c("activ. effector-like CD8+", "F4.80lo macrophages")), 
                #map_signif_level=TRUE, 
                map_signif_level=sigFunc,
                test=wilcox.test, #t.test, wilcox.test
                y_position = bracketoffset, 
                color='black',  
                textsize = 3,
                vjust=sigStar.verticaloffset)
  
  
 # plot(g)
  
#  filename <- paste0("Channel_", sub("_.*", "", m ) , "_byMacs_splitbycondition.png")
#  ggsave(filename,g, width = 30, height = 20, units = "cm",dpi = 600)
  
  
  
  bracketoffset <- 0.6*max(temp.df[[m]])
  g <- 
    ggplot( temp.df  , aes( x=Diet_week, y= temp.df[[m]] , color = Annotated_metacluster )  ) +
    
    
    
    #geom_violin(trim=T)+
    geom_quasirandom(size=0.5,alpha=0.8) + #
    geom_boxplot(outlier.shape = NA, lwd=0.8,alpha=0.7) + 
    #  scale_color_manual(values = cond_palette)+
    
    scale_color_manual( values = subset.Phenograph_metacluster_palette 
                        
    )+
    facet_wrap(~as.factor(Annotated_metacluster))+

  
  theme_classic()+
    
    labs(title= paste0(sub("_.*", "", m ), " channel signal (f.test)"  )    , 
         subtitle= paste0("Using ", nrow(temp.df), " gated non-lymphnode cells in G2.") ,
         #caption="Created by M.Barone", 
         y="Channel signal [transf., z-score]", 
         x= sub("_.*", "", m )  
    )+
    
    theme(panel.background = element_rect(fill = "grey93", colour = "grey55", size = 0.5), # change 'colour' to black for informative axis
          axis.title.x=element_text(color="grey15", size=11),
          axis.title.y=element_text(color="grey15", size=11),
          axis.text=element_text(size=10),
          axis.text.x = element_text(angle = 45, hjust=1),
          #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
          legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank(),
          plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
          plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
    )+
    
    #   "activ. effector-like CD8+",
    
    
    geom_signif(comparisons = list(c("ND_1", "ND_12"),
                                   c("ND_12", "ND_24"),
                                   c("ND_24", "ND_100"),
                                   c("WD_1", "WD_12"),
                                   c("WD_12", "WD_24"),
                                   c("WD_12", "WD_24")    ), 
    #map_signif_level=TRUE, 
    map_signif_level=sigFunc,
    test=wilcox.test, #t.test, wilcox.test
    y_position = bracketoffset, 
    color='black',  
    textsize = 3,
    vjust=sigStar.verticaloffset) 
  
  
  bracketoffset <- bracketoffset*1.1
  g <- g +
    geom_signif(comparisons = list(c("ND_1", "ND_24"),c("ND_24", "ND_100"),
                                   c("WD_1", "WD_24")
                                   ), 
                #map_signif_level=TRUE, 
                map_signif_level=sigFunc,
                test=wilcox.test, #t.test,
                y_position = bracketoffset, 
                color='black',  
                textsize = 3,
                vjust=sigStar.verticaloffset) #+
  
  bracketoffset <- bracketoffset*1.2
  g <- g +
    geom_signif(comparisons = list(c("ND_1", "WD_1")), 
    #map_signif_level=TRUE, 
    map_signif_level=sigFunc,
    test=wilcox.test, #t.test,
    y_position = bracketoffset, 
    color='black',  
    textsize = 3,
    vjust=sigStar.verticaloffset) #+
  
  
  bracketoffset <- bracketoffset*1.1
  g <- g +
    geom_signif(comparisons = list(c("ND_12", "WD_12")), 
                #map_signif_level=TRUE, 
                map_signif_level=sigFunc,
                test=wilcox.test, #t.test,
                y_position = bracketoffset, 
                color='black',  
                textsize = 3,
                vjust=sigStar.verticaloffset) #+
  
  
  bracketoffset <- bracketoffset*1.1
  g <- g +
    geom_signif(comparisons = list(c("ND_24", "WD_24")), 
                #map_signif_level=TRUE, 
                map_signif_level=sigFunc,
                test=wilcox.test, #t.test,
                y_position = bracketoffset, 
                color='black',  
                textsize = 3,
                vjust=sigStar.verticaloffset) #+
 
  
  
  
  
  
  plot(g)
  
  filename <- paste0("Channel_", sub("_.*", "", m ) , "_byCondition_splitbyMacs.png")
  ggsave(filename,g, width = 30, height = 20, units = "cm",dpi = 600)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}
# end plot channel siganl / expression levels 














## --------------- END NEW CODE ------------------



# beside the cluster density in [cells / totalSegementationAreaonROI] Birgit wants to know Cluster presence: cell / iselts on ROI

if( run.old.code==1){

ClustersINisletsPerIslet.data.table <- data.frame()
for(i in plot.ROIs ){
  
  
  temp.df <- subset(cluster.distance.dat, ROI== i & distance==0  )
  temp.amount.of.islets <- fiji.measurement[ match(i, fiji.measurement$ExperimentName   ), ]$Count
  
  # we do not have to protect from temp.amount.of.islets being 0, since we only process ROIs with cells within that mask objects:
  if(nrow(temp.df)>0){
    


  
  temp.freq <-  as.data.frame(table(temp.df$Collapsed_metacluster))
  names(temp.freq) <- c("MC", "Cells")
  
  # at that point, MC is factor:
  temp.freq$MC <- as.numeric(as.character(temp.freq$MC))
  
  
  temp.freq <- cbind(temp.freq, ROI = i)
  temp.freq <- cbind(temp.freq, Condition = unique(temp.df$Diet_week)   )
  
  # now we want to know how many cells per islet are on this ROI. 
  # yes, you would intuively normalize to the total iselt area, creating the density we had before, but here we want to know
  # if we get more cells into the islets...
  temp.freq <- cbind(temp.freq, ClusterDensityPERCENT = 100* temp.freq$Cells/  temp.amount.of.islets   )
  
  
  ClustersINisletsPerIslet.data.table <- rbind(ClustersINisletsPerIslet.data.table,temp.freq)
  
  }#protect from processing ROIs without any cell inside the islet mask 
}#end run through all ROIs to calculate the Cluster presence 


ClustersINisletsPerIslet.data.table <- ClustersINisletsPerIslet.data.table %>%
  mutate(Diet = case_when(
    startsWith(as.character(Condition), "ND") ~ "ND",
    startsWith(as.character(Condition), "WD") ~ "WD"
  ))

ClustersINisletsPerIslet.data.table <-  ClustersINisletsPerIslet.data.table %>%
  mutate(Week = case_when(
    endsWith(as.character(Condition), "_1") ~ 1,
    endsWith(as.character(Condition), "_12") ~ 12,
    endsWith(as.character(Condition), "_14") ~ 14,
    endsWith(as.character(Condition), "_24") ~ 24,
    endsWith(as.character(Condition), "_100") ~ 100
    
  ))

#ClustersINisletsPerIslet.data.table$MC <- factor(ClustersINisletsPerIslet.data.table$MC,      # Reordering group factor levels
#                                                       levels = all.MC ) 





bracketoffset <- 0.1*max(subset(ClustersINisletsPerIslet.data.table, MC>200)$ClusterDensityPERCENT )
# watch out, we print only the G2 here, the other is not gated for CD45 and doestn make any sense 
g <- ggplot(subset(ClustersINisletsPerIslet.data.table, MC>200), aes(x=as.factor(Condition), y = ClusterDensityPERCENT, colour=as.factor(Condition) )) + # the bar wants fill=Cond
  facet_wrap(~MC , scales="free_y")+
  # geom_bar(stat="identity") +
  geom_quasirandom(size=1.5) + #alpha=0.8
  scale_color_manual(values = cond_palette)+
  
  labs(title="G2 cluster frequency per islet on ROI" , 
       subtitle=paste0(nrow(gated.cell.dat), " gated cells in pancreas"), 
       #caption="Created by M.Barone", 
       y="Cluster density [% cells/islets found on ROI]", 
       x=""
  )+
  
  theme(panel.background = element_rect(fill = "grey90", colour = "grey55", size = 0.5), # change 'colour' to black for informative axis
        axis.title.x=element_text(color="grey15", size=11),
        axis.title.y=element_text(color="grey15", size=11),
        axis.text=element_text(size=8),
        axis.text.x = element_text(angle = 45, hjust=1),
        #legend.text=element_text(size=12), # large = 30 # small = 8 # taken out since we define that in the legend layouts just above here. 
        legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
        legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
        #legend.title=element_blank(),
        legend.position="none",
        plot.title = element_text(color="Black", size=14, hjust=0), # size 70 for large, # 18 for small
        plot.subtitle = element_text(color = "Black", size = 12, hjust = 0)
  )+
  
  scale_y_continuous(expand = c(0.2,0))+
  geom_signif(comparisons = list(c("ND_1", "WD_1")), 
              #map_signif_level=TRUE, 
              map_signif_level=sigFunc,
              test=wilcox.test, #t.test,
              y_position = bracketoffset, 
              color='black',  
              textsize = 3,
              vjust=sigStar.verticaloffset) #+

bracketoffset <- bracketoffset* 1.05
g <- g +
  
  geom_signif(comparisons = list(c("ND_12", "WD_12")), 
              #map_signif_level=TRUE, 
              map_signif_level=sigFunc,
              test=wilcox.test, #t.test,
              y_position = bracketoffset, 
              color='black',  
              textsize = 3,
              vjust=sigStar.verticaloffset) #+


bracketoffset <- bracketoffset* 1.05
g <- g +
  
  geom_signif(comparisons = list(c("ND_24", "WD_24")), 
              #map_signif_level=TRUE,
              map_signif_level=sigFunc,
              test=wilcox.test, #t.test,
              y_position = bracketoffset, 
              color='black',  
              textsize = 3,
              vjust=sigStar.verticaloffset) #+


plot(g)

filename <- "ClusterFrequency_withinIslet_byNumberofIsletsonROI.png"
ggsave(filename,g, width = 30, height = 20, units = "cm",dpi = 1200)
}#end protect running old code


}#end stage 7: run proximity index engine with EDM matrices 


########################## STAGE 8 ########################## 
# STAGE 8: Use Adrian Hucks slightly modified run.spatial.analysis function for reg.dat build up
if(whichSTAGE == 8){
  
  celltocell.distance <-function (dat, sample.col, pop.col, annot.cols = NULL, region.col = NULL, 
                                  area.table = NULL, adj.dist = 20, x.col = "x", y.col = "y", 
                                  distribution = TRUE, composition = TRUE, counts = TRUE, counts.per.area = TRUE, 
                                  distance = TRUE, adjacency = TRUE, func = "mean", numCores=4)
  {
    
    
    require(foreach)
    require(doParallel)
    
    
    pops <- unique(dat[[pop.col]])
    samples <- unique(dat[[sample.col]])
    regions <- unique(dat[[region.col]])
    setorderv(dat, sample.col)
    setorderv(dat, pop.col)
    all.counts <- list()
    
    registerDoParallel(numCores)
    foreach (i = samples, .combine=rbind) %dopar% {
   
      message("Processing ", i, " -- (", which(samples == i), " / ", length(samples), ")")
      samp.dat <- dat[dat[[sample.col]] == i, ]
      samp.front <- samp.dat[1, c(sample.col, annot.cols), with = FALSE]
      
      
      
      message(" -- Calculating distance")
      
      
      names.dist <- c(pop.col, region.col, "CellFrom", "CellTo", "Distance")
      
      
      adj.res.list <- data.table()
      regions.sub <- unique(samp.dat[[region.col]])
      
      comb.dist.final <- DataFrame(matrix(nrow = 0, ncol = 5))
      
      for (u in regions.sub) {
        
        message("     ", paste0("Region: ", u))
        reg.temp <- do.filter(samp.dat, region.col, u)
        reg.pops <- unique(reg.temp[[pop.col]])
        message("got till here")
        combs <- gtools::permutations(n = length(reg.pops), r = 2, v = reg.pops, repeats.allowed = TRUE)
        
        #   comb.adj.list <- list()
        comb.dist.list <- list()
        for (a in c(1:nrow(combs))) {
          cell.from <- combs[a, 1]
          cell.to <- combs[a, 2]
          message(paste0(cell.from, " to ", cell.to))
          
          pop.temp <- do.filter(reg.temp, pop.col, c(cell.from,
                                                     cell.to))
          cell.types <- pop.temp[[pop.col]]
          pop.temp <- pop.temp[, c(x.col, y.col), with = FALSE]
          A <- pop.temp[which(cell.types == cell.from), ]
          B <- pop.temp[which(cell.types == cell.to), ]
          res <- proxy::dist(as.data.frame(A), as.data.frame(B))
          if (cell.from == cell.to) {
            diag(res) <- NA
          }
          calc.dist <- sum(rowMeans(res, na.rm = TRUE),   na.rm = TRUE)/nrow(res)
          #comb.adj.list[[nme]] <- sum(res < adj.dist, na.rm = TRUE)/nrow(res)
          
          dist.row <- (c(i, u, cell.from, cell.to, calc.dist))
          comb.dist.list[[a]] <- dist.row
          
        }
        comb.dist.list.region <- do.call(rbind, comb.dist.list)
        colnames(comb.dist.list.region) <- names.dist
        comb.dist.final <- rbind(comb.dist.final, comb.dist.list.region)
        comb.dist.final$Distance <- as.numeric(comb.dist.final$Distance)
      }
      comb.dist.final
    }
    
  }
  
  reg.dat2 <- celltocell.distance(dat = data.table(subset(gated.cell.dat, GATEsubset == 2 & ROI %in% plot.ROIs[c(1:1)])), 
                                       sample.col = 'ROI', 
                                       pop.col = "Annotated_metacluster", 
                                       annot.cols = "UICC_grouped", #"Patient_TT", # this defaulted to: group.col and by doing so only split WD and ND. Lets see if using Diet_Week does the job.
                                       region.col = "tissue.region" , # region.col <- our spatial.dat is already subset and contains only pancreas here.
                                       adj.dist = 20,
                                       numCores=5)
  
  reg.dat2 <- celltocell.distance(dat = cell.dat, 
                                       sample.col = 'ROI', 
                                       pop.col = "Annotated_metacluster", 
                                       annot.cols = "UICC_Staging", #UICC_grouped", #"Patient_TT", # this defaulted to: group.col and by doing so only split WD and ND. Lets see if using Diet_Week does the job.
                                       region.col = "Tissue" , # region.col <- our spatial.dat is already subset and contains only pancreas here.
                                       adj.dist = 20,
                                       numCores=5)
  
  
  
  calc.close.contacts  <- function (dat, 
                                    sample.col, 
                                    pop.col, 
                                    annot.cols = NULL, 
                                    region.col = NULL, 
                                   area.table = NULL, 
                                   adj.dist = 20, 
                                   x.col = "x", 
                                   y.col = "y", 
                                   distribution = TRUE, composition = TRUE, counts = TRUE, counts.per.area = TRUE, 
                                   distance = TRUE, adjacency = TRUE, func = "mean", numCores=4)
  {
    # this is based on the re-written function of Adrian Huck, sent to me 20230706
    
    require(foreach)
    require(doParallel)
    require(progress)
    
    
    pops <- unique(dat[[pop.col]])
    samples <- unique(dat[[sample.col]])
    regions <- unique(dat[[region.col]])
    setorderv(dat, sample.col)
    setorderv(dat, pop.col)
    
    registerDoParallel(numCores)
    

    
    
    
    # foreach (i = samples, .combine=rbind) %dopar% {
     for (i in samples) {
      
      
      message("Processing ", i, " -- (", which(samples == i),      " / ", length(samples), ")")
      samp.dat <- dat[dat[[sample.col]] == i, ]
      samp.front <- samp.dat[1, c(sample.col, annot.cols), with = FALSE] #
      
      
      
     # message(" -- Calculating adjacency")
      
      # watch out, here was:
     # names.adj <- c(pop.col, region.col, "CellFrom", "CellTo", "Adjacency")
      # but we push in "i" afterwards, which is sample.col
      names.adj <- c(sample.col, region.col, "CellFrom", "CellTo", "Adjacency")
      
      regions.sub <- unique(samp.dat[[region.col]])
      
      #comb.adj.final <- DataFrame(matrix(nrow = 0, ncol = length(names.adj)))
      comb.adj.final <- data.frame(DataFrame(matrix(nrow = 0, ncol = length(names.adj))))
      
      for (u in regions.sub) {
        

        
       # message("     ", paste0("Region: ", u))
        reg.temp <- do.filter(samp.dat, region.col, u)
        reg.pops <- unique(reg.temp[[pop.col]])
        combs <- gtools::permutations(n = length(reg.pops), r = 2, v = reg.pops, repeats.allowed = TRUE)
        

        
        pb <- progress_bar$new(format = paste0("[:bar] :percent [Calc. close contacts in ",u," | :eta]"),
                               total = nrow(combs), # we run through the ROIs with three regions
                               show_after=0, #allows to call it right way 
                               current = "|",    # Current bar character
                               clear = FALSE) # make it persist
        pb$tick(0) # call in the progress bar without any progress, just
        
        comb.adj.list <- list()
        for (a in c(1:nrow(combs))) {
          cell.from <- combs[a, 1]
          cell.to <- combs[a, 2]
         # message(paste0(cell.from, " to ", cell.to))
          
          pop.temp <- do.filter(reg.temp, pop.col, c(cell.from,
                                                     cell.to))
          cell.types <- pop.temp[[pop.col]]
          pop.temp <- pop.temp[, c(x.col, y.col), with = FALSE]
          A <- pop.temp[which(cell.types == cell.from), ]
          B <- pop.temp[which(cell.types == cell.to), ]
          res <- proxy::dist(as.data.frame(A), as.data.frame(B))
          if (cell.from == cell.to) {  diag(res) <- NA       }
          
          calc.adj <- sum(res < adj.dist, na.rm = TRUE)/nrow(res)
          
          adj.row <- (c(i, u, cell.from, cell.to, calc.adj))
          comb.adj.list[[a]] <- adj.row
          

          
          pb$tick() 
         
        }
        
        # now we compress the list into an array
        comb.adj.list.region <- do.call(rbind, comb.adj.list)
        colnames(comb.adj.list.region) <- names.adj
        message("got here 1")
        
        # the second region fails HERE! after running the first region through checkpoints 1 and 2, the second region gets to 1 and then pops the error:
        #Error in if (class(from) == "list") { : the condition has length > 1
        # its because  class(comb.adj.final) is #[1] "DFrame" attr(,"package")[1] "S4Vectors"
        comb.adj.final <- rbind(comb.adj.final, comb.adj.list.region)
        comb.adj.final$Adjacency <- as.numeric(comb.adj.final$Adjacency)
        comb.adj.final <- as.data.frame(comb.adj.final)
        message("got here 2")
      }
      message("got here 3")
      comb.adj.final
      #stopImplicitCluster()
      
  #   
    }
    #stopImplicitCluster()
    
    return(comb.adj.final)
  }
  
  
  fast.calc.close.contacts  <- function (dat, 
                                    sample.col, 
                                    pop.col, 
                                    annot.cols = NULL, 
                                    region.col = NULL, 
                                    area.table = NULL, 
                                    adj.dist = 20, 
                                    x.col = "x", 
                                    y.col = "y", 
                                    distribution = TRUE, composition = TRUE, counts = TRUE, counts.per.area = TRUE, 
                                    distance = TRUE, adjacency = TRUE, func = "mean", numCores=4)
  {
    # this is based on the re-written function of Adrian Huck, sent to me 20230706
    
    require(foreach)
    require(doParallel)
    require(progress)
    
    
    pops <- unique(dat[[pop.col]])
    samples <- unique(dat[[sample.col]])
    regions <- unique(dat[[region.col]])
   # dat <- as.data.table(dat)
   # print(class(dat))
    setorderv(dat, sample.col)
    setorderv(dat, pop.col)
    
    registerDoParallel(numCores)
    
    
     #dat <- data.table(dat)
     # print(class(dat))
    
     foreach (i = samples, .combine=rbind) %dopar% {
    #for (i in samples) {
      
      
      message("Processing ", i, " -- (", which(samples == i),      " / ", length(samples), ")")
      #samp.dat <- dat[dat[[sample.col]] == i, ]
      samp.dat <- subset(dat, sample.col==i )
      class(samp.dat)
      samp.front <- samp.dat[1, c(sample.col, annot.cols), with = FALSE] #
      
      
      
       message(" -- Calculating adjacency")
      
      # watch out, here was:
      # names.adj <- c(pop.col, region.col, "CellFrom", "CellTo", "Adjacency")
      # but we push in "i" afterwards, which is sample.col
      names.adj <- c(sample.col, region.col, "CellFrom", "CellTo", "Adjacency")
      
      regions.sub <- unique(samp.dat[[region.col]])
      
      #comb.adj.final <- DataFrame(matrix(nrow = 0, ncol = length(names.adj)))
      comb.adj.final <- DataFrame(matrix(nrow = 0, ncol = length(names.adj)))
      
      for (u in regions.sub) {
        
        
        
         message("     ", paste0("Region: ", u))
        reg.temp <- do.filter(samp.dat, region.col, u)
        reg.pops <- unique(reg.temp[[pop.col]])
        combs <- gtools::permutations(n = length(reg.pops), r = 2, v = reg.pops, repeats.allowed = TRUE)
        
        
        
      #  pb <- progress_bar$new(format = paste0("[:bar] :percent [Calc. close contacts in ",u," | :eta]"),
      #                         total = nrow(combs), # we run through the ROIs with three regions
      #                         show_after=0, #allows to call it right way 
      #                         current = "|",    # Current bar character
      #                         clear = FALSE) # make it persist
      #  pb$tick(0) # call in the progress bar without any progress, just
        
        comb.adj.list <- list()
        for (a in c(1:nrow(combs))) {
          cell.from <- combs[a, 1]
          cell.to <- combs[a, 2]
          message(paste0(cell.from, " to ", cell.to))
          
          pop.temp <- do.filter(reg.temp, pop.col, c(cell.from,
                                                     cell.to))
          cell.types <- pop.temp[[pop.col]]
          pop.temp <- pop.temp[, c(x.col, y.col), with = FALSE]
          A <- pop.temp[which(cell.types == cell.from), ]
          B <- pop.temp[which(cell.types == cell.to), ]
          res <- proxy::dist(as.data.frame(A), as.data.frame(B))
          if (cell.from == cell.to) {  diag(res) <- NA       }
          
          calc.adj <- sum(res < adj.dist, na.rm = TRUE)/nrow(res)
          
          adj.row <- (c(i, u, cell.from, cell.to, calc.adj))
          comb.adj.list[[a]] <- adj.row
          
          
          
       #   pb$tick() 
          
        }
        
        # now we compress the list into an array
        comb.adj.list.region <- do.call(rbind, comb.adj.list)
        colnames(comb.adj.list.region) <- names.adj
    
        
        # the second region fails HERE! after running the first region through checkpoints 1 and 2, the second region gets to 1 and then pops the error:
        #Error in if (class(from) == "list") { : the condition has length > 1
        # its because  class(comb.adj.final) is #[1] "DFrame" attr(,"package")[1] "S4Vectors"
        comb.adj.final <- rbind(comb.adj.final, comb.adj.list.region)
        comb.adj.final$Adjacency <- as.numeric(comb.adj.final$Adjacency)
        #comb.adj.final <- as.data.frame(comb.adj.final)
   
      }
  
      comb.adj.final
      #stopImplicitCluster()
      
      #   
    }
    #stopImplicitCluster()
    
    return(comb.adj.final)
  }
  
  
  reg.dat2 <- fast.calc.close.contacts(dat = subset(gated.cell.dat, GATEsubset == s & ROI %in% plot.ROIs[c(1:1)]), 
                                 sample.col = 'ROI', 
                                 pop.col = "Annotated_metacluster", 
                                 annot.cols = "UICC_grouped", #"Patient_TT", # this defaulted to: group.col and by doing so only split WD and ND. Lets see if using Diet_Week does the job.
                                 region.col = "tissue.region" , # region.col <- our spatial.dat is already subset and contains only pancreas here.
                                 adj.dist = 20,
                                 numCores=5)
  
  
  # lets see how that baby purrs:
  
  # the IM project of Laia used:
#  roi.col <- 'ROI'
#  sample.col <- "Patient_TT"# "Patient_TT" # #'Sample'. 
#  group.col <- "UICC_grouped" #'Group'
#  batch.col <- 'Batch'
#  
#  pop.col <- "Annotated_metacluster" # watch out, Laia had 'Annotated_metacluster' # watch out, the script defaults the string with a space. Since I worked with the data.frame, this is by now Annotated.metacluster with a dot!
#  region.col <- "Tissue"# 'Annotated region' # <- again, this is gonna fail since we dont have the regions defined yet.
  
  # BeCeFa used:
  
#  roi.col <- 'ROI'
#  sample.col <- "Diet_week"# "Sample_week" # #'Sample'. #
#  group.col <- "Diet" #'Group'
#  batch.col <- 'Batch'
  
#  pop.col <- 'Annotated_metacluster' # watch out, the script defaults the string with a space. Since I worked with the data.frame, this is by now Annotated.metacluster with a dot!
#  region.col <- "Tissue"# 'Annotated region' # <- again, this is gonna fail since we dont have the regions defined yet.
  
  
  s<-2 # is probably set, just in case
  #celltocell.adjacency celltocell.distance
  reg.dat <- calc.close.contacts(dat = subset(gated.cell.dat, GATEsubset == s ), 
                                     sample.col = 'ROI', 
                                     pop.col = "Annotated_metacluster", 
                                     annot.cols = "UICC_grouped", #"Patient_TT", # this defaulted to: group.col and by doing so only split WD and ND. Lets see if using Diet_Week does the job.
                                     region.col = "tissue.region" , # region.col <- our spatial.dat is already subset and contains only pancreas here.
                                  adj.dist = 20,
                                  numCores=5)
  
  
}

