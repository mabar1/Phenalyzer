run.experimental.code<-0

########################## STAGE 0 SETTINGS ##########################

# helper functions are stored in the utils file
# define the version you want to use here:
# use it for project-specific stuff
# e.g: "0" will load utils_v0.R 
utils_version <- "0" 

# Do we start from scratch or re-load an old run?
#-------------------------------------------------
# 0: Load all data saved away from previous run. 
# Caution1: for 0 changes in the metadata file wont be incorporated since you skip STAGE 0
# Caution2: Since you essentially skip STAGE 0-3 entirely, you need to supply a good amount of variables manually here: "COPY-PASTE VARIABLES IF YOU DO NOT GENERATE CELL.DAT ANEW"
# 1: recommended. Read cell.dat anew, connect metadata file. Doesn't take long and is the clean way.
rebuild.cell.data.anew <- 1 # 1 recommended 

# for rebuild.cell.data.anew=0: Reload an old run, so get backup data frames (such as gated.cell.dat) from this date:
day.to.reload.data <- 20231018 # date from where to load data: use 18th Oct to include the trashbin gate and proper cluster numbers in G2

# Get the files from this add-masks integration run:
#-------------------------------------------------
# 20230113 comes from a global segementation model with all 4 TMAs
# 20240524 is Steinbock but now with CK18 channel but mask.tiff in cell.dat somehow
# 20240525 finally clean and good
date.of.spatial.dat <- 20240525 # do.fast.extract did that one



# Get the metadata 
#-------------------------------------------------
# Point to Metadata stored as excel sheet with an absolute path. This way you know exaclty what file was loaded.
# for gate thresholds, is is crucial that they have the format:
# CD123.highpass, CD122.lowpass
# furthermore you need area.highpass and area.lowpass columns in square micrometers
# this select cells area.highpass < segmentation areas < .lowpass area penalties
metadatapath <- c("Small TMA_testing_metadata.xlsx") # cross-checked by Juliette, added the gate thresholds 20230313

# Data transformation and normalization:
#-------------------------------------------------
# Transformation will create columns with a "_t5", short for asinh transformed with cofactor 5 (ry different factors and check if heatmap makes sense)
# If this is turned on, the normalization will happen on these columns and not the raw data anymore!
# Min-max normalize transformed signals? Additional to _t5, they will get "_m0m1" showing you the boundaries you have defined
# You can also Z-score the data as well. Additional to _t5, they will get "_z" for z-score data.
# Transform:
transform.rawdata = TRUE
asinh.cofactor = 2 # we need to test them at some point
# Min-Max normalize
minmax.norm = FALSE
min.norm = 0
max.norm = 1
# Z-score data
zscore.norm = TRUE

# which cellular readout is processed by the upper transformations?
cellular.segment.compartment <- "Cell" # or "Cell", "Nucleus", "Membrane" 
cellular.segment.readout <- "Mean" # or "Median", "Min", "Max" 


# Spillover correction via excel sheet:
#-------------------------------------------------
correct.spill.over <- 0 
# 20211021_BeCeFa_spillmat.xlsx is from the panel design
# 20200617_Messung6_sm.xlsx is from Franziska given on 20220601
# the spillover matrix should have the donors in first column, and this column name is "acq.chnl". 
# Whether you supply it in [0-1] or in percent [0-100] does not matter, the script will take care of that. 
# Empty fields do not need to be filled with zeros, we fill them up for you.
#spillover.matrix.path <- c("D:/IMI/IMC/Spillover-Matrix/20200617_Messung6_sm_spillma.xlsx") 



# Tissue types supplied via metadata. Positive Control core name
#-------------------------------------------------
tissue.type.col <- "Tissue.type"
postive.control.Tissue.name <- "spleen" # name of the tissue as defined in the metadata file. Watch out, these cells are deleted later on.

# depreciated:
# Positive Control coreBatch normalization of markers present on pos-ctr. 
#-------------------------------------------------
# if you want to do spillover and batch normalization, run stage 0 the first time with correct.spill.over 0 as well as scaling factors off
# after you checked currently selected "Cellular cols", turn spillover and scaling on for the second+ runs (no worry, I will run it only once)
# IMPORTANT: if you perform batch normalization, leave calculate.scaling.factors on 1 for the upcoming stages. Do NOT turn that back to 0!
calculate.scaling.factors <- 0 # do batch normalization: 1, off: 0
# Batchnormalization: Markers to normalize at the given percentile across all TMA. 
#-------------------------------------------------
# Data will be rescaled onto the n-th percentile of one marker
# The one TMA whose set percentile signal is closest to the median of all TMA percentiles becomes the anchor sample: 
# The anchor sample signal is then used to calculate the individual scaling factors for the remaining TMAs 
# only the following percentiles are plotted and are available: 
# 0.6,0.8,0.85,0.90,0.95,0.99
batch.normalization.matrix <- c(
  
    "CD3_Nd143", 0.80,
    "CD4_Dy161", 0.80,
    "CD8_Nd145", 0.80,
    "CD127_Eu153", 0.80,
    "FoxP3_Tm169", 0.80, 
    "EOMES_Tb159", 0.80,
    "CRTH2_Gd158", 0.80,
    "CD56_Yb172", 0.80,
    "Ki67_Er170", 0.95,
    "CD44_Eu151", 0.99,
    "CD69_Sm154", 0.80,
    "CD16_Sm147", 0.80,
    "CD19_Sm149", 0.95,
    "F4.80_Dy162", 0.95,
    "MHCII_Ho165", 0.90,
    "CD11c_Yb171", 0.95,
    "CD11b_Nd150", 0.85,
    "Ly6G_Gd160", 0.80,
    "B220_Nd142", 0.80,
    "CD138_Nd144", 0.80,
    "CD31_Nd146", 0.85,
    "CD45_Er168", 0.95
)




# Define which markers are transformed and normalized 
#-------------------------------------------------

# since v18-ish we do no longer use column numbers to pull the markers since this is very fragile.
# instead we supply the markers and let the code build up the column name by itself. 
# for this, we need to know how the tsv gets fixed in the code after import. Usually it adds two dots for the separator. change here in case qupath should change that behaviour:
# how are the markers separated from Median, Mean and so on:
qupath.separator <- ".." 

# Important: exclude any marker channels that were not recorded consistently on all TMAs!

# The first time STAGE 0 is executed, the code will show all available columns in cell.dat, among which the raw data column of markers
# You choose which raw signal columns should be processed further for the upcoming clustering. Exclude any

Markers2Process <- c(
#  # major lineage T cells:
"CD3e",
"CD4",
"CD8",
"TRG",
#  # t cell activation status
"HLA.DR",
"CD103",
"CD44",
"CD38.AKZP0110", 
#  "FOXP3",
#  # major lineage B cells: 
#  #"CD19",
"CD20",
#  # immune
"CD45",
"Ki67",
#  # other immune cells
"CD11c", # DC
"CD11b", # NK
"CD68", # macrophages
#  # other structural
"CD147", # MMP inducer
"Podoplanin", # endo and mesothelial cells
"Pan.Cytokeratin", # epithelial cells
"E.cadherin", # epithelial cells
# "Collagen.IV", # basal lamina
"CD31",
#  #cancer?
"CD117", # c-KIT
  "b.Catenin1"
)


# Set the order of the markers to show up in the heatmap 
#-------------------------------------------------
# this needs manual adaptation. All markers that go into clustering will be plotted in the heatmap. 
# you need them all in this list and give them some order that helps looking at the heatmap:


global.marker.order <-Markers2Process
#global.marker.order <- c(
#  "CD45",
#  "CD4",
#  "Ki67",
#  "DAPI"
#)



########################## STAGE 1 SETTINGS ##########################


# Define key columns supplied via metadata file
#-------------------------------------------------
# here we just plug the col names from the sample.meta: see above: # so from now on we need to all Sample as Sample_week Group as Diet Batch as Batch
roi.col <- 'ROI'
#sample.col we had Sample_week, but this does not allow me to split afterwards "Sample_week"
# so you need a tag here that contains the week and sample, some compressed tag
sample.col <- "Tissue" # best is to create a combination of two parameters that allows to split the data set into as many categories as possible
group.col <- "Tissue" #'Group' # here we just sep for Diets to get the diets on the same plot, then use sample.col to sep the weeks from each other.
batch.col <- 'Image'





# Define marker combinations that you need to look at when finding your gating routine
#-------------------------------------------------
# While plotting the scatter plot, the code will use different coloring:
#  1) if no gate is present for any of the two markers, the plot comes with a soothing blue
#  2) if a gate is found, the cells are color-coded according whether they passed the gate or not
#  3) since this is R and we can just easily do that: supply a third marker and let the scatter plot be colored according to this raw signal
#     The color range of this third marker can be adapted. Do so to dampen the influence of outlier signals. 
#     Set the the lower and upper limit of the gradient symmetrically at this percentile:
percentile.colorlimits <- 0.95 # set the upper percentile here, the lower will be symmetrical to that


scatter.plots <- list( 
  c("CD45", "CD4"  , NA), # <- gate in B
  c("CD45", "Ki67", NA), # <- gate in T
  c("CD45", "CD4", "DAPI") # <- gate in other immune
)  





# Quality Control and plots for setting up gating thresholds:
#-------------------------------------------------
# 1: Check positions of gates and segmentation area penalties?
plot.TMAwise.gates <- 0
# 2: Check the distribution of cell segementation area across TMAs and tissues?
# metadata needs area.highpass and area.lowpass columns in square micrometers. typical cells are up to 300μm^2
plot.segementation.areas <- 1 


# 3. Check the gating steps. How does the dataset change while running through the hierarical gating? 
# So this is the "With this subset of cells we apply high/lowpass gate xyz"
plot.scatterplots.gatingsteps <- 0
# when plotting the gating step scatter plots, you can either pool all marker in comparison into one file (one file per "step")
# or you can plot every "marker" comparison separately, separately for every of the gating steps
# watch out, choosing "marker" over "step" is even more time-consuming!
# plotting step can get very messy if you have many subgate populations and many scatter plots with gates to plot
scatterplots.gatingsteps.packaging = "marker"#-wise 




# Set up gating routine (first part)
# -> the second is found in the code at "SET UP GATING ROUTINE")
#-------------------------------------------------

#lets define some names for the gated subsets that we can use to store the files with unique names as we run through the GATEsubset:
# be aware that we gate all bounced back cells in an additional gate, the trash bin gate so to say. 
# But turning on include.trash.bin.gate you can plot and cluster these cells 
#GATEsubsetNames <- c( "G1 CD19+", 
#                      "G2 CD19- AND CD45+ AND CD3+", 
#                      "G3 CD19- AND CD45+ AND CD3-",
#                      "G4 CD19-, CD3-, AND CD45-"#, <- thats the trash.bin.gate name
#                      #"G5 empty"# <- thats the trash.bin.gate name
#)



#Birgit's clustering strategy        
#1 = CD19+ (B cells)
#Markers for clustering: AICL, KLRF1, T.Bet,Granzyme B,KLRG1,CD25,CD107a,PD1,Tim3,CD69,CD38,CD103,HLADR,CD11c,PDL1,Ki67 (16)

#2 = CD19-CD3+ (other immune cells)
#Markers for clustering: AICL, KLRF1, CD4,CD8,CD127,Foxp3,TCRdelta,T.Bet,Granzyme B,KLRG1,CD25,CD107a,CTLA4,PD1,Tim3,CD69,CD38,CD103,CD31,
#                        HLADR,CD45,CD68,CD11c,PDL1,Ki67

#3 = CD19- CD3- CD45+ (other immune cells)
#Markers for clustering: CD127,CD107a,Tim3,CD69,CD38,CD103,CD31,HLADR,CK18,CD324,PDL1,CA19.9,Ki67,aSMA 

#4 G4 CD19-, CD3-, AND CD45- (non immune cells)
# CD127 107a tim3 69 38 103 31 hladr CK18 324 pdl1 CA19.9 ki67 aSMA 
# since she wants to spot PD1+ Tcells, this was added to the G2 and G3

#v19 change all fragile marker numbers to marker names that are arraged in user-defined order that will be used to plot the HM 
# additionally, the marker is selected for clustering by adding the subset number(s) in GATEsubsetNames
markers.to.cluster.per.subset <- list(
  
  "G1: Jurkarts",
  c("CD45",
    # immune
    "CD4",
    "Ki67"
  )
)


if(run.experimental.code==1){
markers.to.cluster.per.subset <- list(
  
  "G1: B Cells",
  c("CD20",
    # immune
    "CD45",
    "CD11b", # some B cells in peritoneal cavity express that
    "CD147", # MMP inducer
    "Ki67"
    ),
  
  "G2 T Cells",
  c("CD3e",
    "CD4",
    "CD8",
    #"TRG", # Adrian does not trust this one
    # t cell activation status
    "HLA.DR",
    "CD103",
    "Ki67",
    "CD44",
    #"FOXP3",
    "CD38.AKZP0110",
    #"Podoplanin", # endo and mesothelial cells
    #"Pan.Cytokeratin", # epithelial cells
    #"E.cadherin", # epithelial cells
    #"CD147", # MMP inducer
    #"Collagen.IV", # basal lamina. Adrian doesnt trust this one
    "CD31"),
  
  "G3 other immune",
  c("CD11c", # DC
    "CD11b", # innate marker such as CD14+CD11b+ moncytes, tiss.res macs, neutros, eosinos conv DCs, NKs
    "CD68", # macrophages
    # other structural
    "Podoplanin", # endo and mesothelial cells
    "Pan.Cytokeratin", # epithelial cells
    "E.cadherin", # epithelial cells
   # "Collagen.IV", # basal lamina. Adrian doesnt trust this one
    "CD147", # MMP inducer
    "CD31",
    #cancer?
    "b.Catenin1"),
  
  "G4 other cells",
  c(# other structural
    "Podoplanin", # endo and mesothelial cells
    "Pan.Cytokeratin", # epithelial cells
    "E.cadherin", # epithelial cells
    "CD147", # MMP inducer
    #"Collagen.IV", # basal lamina. Adrian doesnt trust this one
    "CD31",
    
    #cancer?
    "b.Catenin1")
)
}
# In the code of stage 1 "SET UP GATING ROUTINE" you need to define the gating routine. 
# End the hierarchical gating always by collecting cells that were rejected from every gate into a "flowthrough" subset
# Especially if you set up a multi-step gating, you are probably not interested in the cells that were not gated:
# You can therefore just drop these cells in the last subset before processing, plotting, and being merged with the clustered and annotated subsets.
exclude.trash.bin.gate <- 0 # if you need the last subset, turn off. In that case it is of course not a trash bin subset.

# if you want the above defined set of markers ordered in the same way in every heatmap of each subset, 
# use global marker order defined in stage 0 settings and turn on here: 
use.global.marker.order <- "yes" # yes recommended



# Define which markers are used for clustering each gated subset individually
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Important: also define a marker set for the last subset, whether you will analyze it or not.
# also: keep the -0 in: if you add/remove columns in the metadata file, these numbers need to be corrected to pull the proper columns out!

# The first time STAGE 1 is executed, the code will show all available columns in gated.cell.dat
# You choose from this list the markers for every subset. 
# This allows you to even run the code with differently transformed/normalized data
# Be aware that this is fragile, and adding/removing one column in the metadata file will cause this list the be offset by one

# for the first attempt, we rely on a mail from birgit:
#Birgit's clustering strategy        
#1 = CD19+ (B cells)
#Markers for clustering: AICL, KLRF1, T.Bet,Granzyme B,KLRG1,CD25,CD107a,PD1,Tim3,CD69,CD38,CD103,HLADR,CD11c,PDL1,Ki67 (16)

#2 = CD19-CD3+ (other immune cells)
#Markers for clustering: AICL, KLRF1, CD4,CD8,CD127,Foxp3,TCRdelta,T.Bet,Granzyme B,KLRG1,CD25,CD107a,CTLA4,PD1,Tim3,CD69,CD38,CD103,CD31,
#                        HLADR,CD45,CD68,CD11c,PDL1,Ki67

#3 = CD19- CD3- CD45+ (other immune cells)
#Markers for clustering: CD127,CD107a,Tim3,CD69,CD38,CD103,CD31,HLADR,CK18,CD324,PDL1,CA19.9,Ki67,aSMA 

#4 G4 CD19-, CD3-, AND CD45- (non immune cells)
# CD127 107a tim3 69 38 103 31 hladr CK18 324 pdl1 CA19.9 ki67 aSMA 
# since she wants to spot PD1+ Tcells, this was added to the G2 and G3

# watch out, as of 9th of Oct we gate only for the CD45+ G2 and G3, so CD45 needs to go from these. its 149
# watch out, as of 17th of Oct we add CD4-multisegment line in (its called CD8.lowpass) - and everything is shifted by 1.

#GATEsubset.Markers2Cluster <- list(
#  #"n-th gate:",   c( ..., ...), 
# "Gate Subset 1:",   c(128,146,139,145,136,143,149,140,155,158,151,135,148,131,150,129)+1, # x-th entry in names(gated.cell.dat) vector
# "Gate Subset 2:",   c(128,146,137,127,132,157,139,145,136,143,149,130,140,155,158,151,135,134,148,154,156,131,150,129,160)+1, # watch out this contained CD3 which is gate: 110+30
# "Gate Subset 3:",   c(133,150,156,159,152,136,135,149,160,139,151,127,130,162,157,132)+0, # CD68, CD11c was missing
# "Gate Subset 4:",   c(133,150,156,159,152,136,135,149,160,139,151,127,130,162)+0#,
# #"Trashbin gate:",   c(110) # trash bin gate also gets an existing column number

#)










########################## STAGE 2 SETTINGS ##########################

# Algorithms to run (this is no option anymore...)
# - - - - - - - - - - - - - - - - - - - - - - - - -
run.FlowSOM <- 0 # depreciated. not run
run.Phenograph <- 1
run.tSNE <- 0 # depreciated. not run
run.UMAP <- 1


phenograph.k=30 #k applies to Phenograph and FlowSOM (default 30)


# Cluster engine options

# cytof: Jin Miao Chen's cytof engine
# Rphenoannoy: Stuchly lab's Rphenoannoy engine
# Leiden: TomKellyGenetics Leiden engine
cluster.engine="Rphenoannoy" #"cytof" #choose engine Jin Miao Chen: cytof, stuchly lab: Rphenoannoy
# Rphenoannoy is faster because its a multicore approach https://github.com/stuchly/Rphenoannoy


# When gating into subsets to cluster, should we 
# keep the global z-scores calculated across all subsets or
# re-normalize the subset (the way Somesh does it in scRNA-seq)
marker.zscoring = "global" # "global" "subset"


# major lineage markers might need more weight for clustering and heatmap spread

#For subset 2, give more weight to CD4 and CD8 markers
#if (s == 2) {
#  zscored.exp.mat[, "CD4_Dy161_t2_z"] <- zscored.exp.mat[, "CD4_Dy161_t2_z"] * 2
#  zscored.exp.mat[, "CD8_Nd145_t2_z"] <- zscored.exp.mat[, "CD8_Nd145_t2_z"] * 2
#}


# After clustering a first time, delete cluster numbers and re-run
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# you can use clustering to get rid of rim effects, mis-gated populations, or wrongly segmented cells.
# for that the script has different behavior when running STAGE 2 the first, or the second+ times.
# after the first run, UMAPs and Heatmap are stored in the STAGE 2 First pass folder. 
# select the clusters you want to delete by putting them into the field just down here. 
# The delete step will be executed exactly ones in the second pass!
# the syntax line by line in the list: c("subset number","cluster number as written in the First Pass folder")
delete.clustersinsubsets  <- list(
                                 # c(1, 21), # subset G1, cluster 21
                                 # c(2, 19), # subset G2, cluster 19
                                 # c(2, 20)  # subset G2, cluster 20
                                  )
# Re-run a second time and these cluster numbers get deleted, and Clusters, Heatmap, and UMAP newly calculated.
# upon a third re-run, further deleting is blocked for obvious reasons


# If you run STAGE 2 a second time, whatever is written in the upper list is deleted for good.
# To make sure we do not delete unwanted clusters (please check the output folders before re-running STAGE 2!!)
# Important: Stage 2 has a ticker that wont allow you to delete clusters after running it two times. 
confirm.cluster.delete.after.first.run.STAGE2 <- "no" # yes no # confirm here that you checked them before re-running STAGE 2:




# Also analyse positive control with these algorithms?
# - - - - - - - - - - - - - - - - - - - - - - - - -
cell.analyis.spleen <- 0 # not recommended if you have a lot of ROIs with dense tissue...




    




# not needed in Laias data:
# for the analysis of the pos- control
marker.order.spleen <- c(
  
  "CD3_Nd143_asinh",
  "CD4_Dy161_asinh",
  "CD8_Nd145_asinh",
  "CD127_Eu153_asinh",
  "FoxP3_Tm169_asinh",
  "EOMES_Tb159_asinh",
  "CRTH2_Gd158_asinh",
  "CD56_Yb172_asinh",
  "Ki67_Er170_asinh",
  "CD44_Eu151_asinh",
  "CD69_Sm154_asinh",
  "CD16_Sm147_asinh",
  "F4.80_Dy162_asinh",
  "CD11c_Yb171_asinh",
  "CD11b_Nd150_asinh",
  "Ly6G_Gd160_asinh",
  "B220_Nd142_asinh",
  "CD138_Nd144_asinh",
  "CD31_Nd146_asinh",
  "CD45_Er168_asinh",
  "GLP.1R_Yb174_asinh",
  "NKx6.1_Er166_asinh"
  
)





# ignore:
run.old.code <- 0 # protect the script from running old code that is currently not used (and almost certainly crash the script)
run.experimental.code <- 0 # this is a temporal fix to run code that would probably crash the script in the current state
plot.complexheatmap <- 0 # this is a temporal fix to not run complexheatmap code that is far from running








########################## STAGE 3 SETTINGS ##########################

# if you collapse clusters, the numbers get renumbered. 
# To know which clusters should be collapsed and to which extent we need to collapse branching points in the dentrogram,
# we can collapse first a few, and then keep on collapsing further and further to different levels.
# as you collapse into the first level, clusters will have new numbers. Define if you then sequentially use these new numbers (and work with the new heatmap),
# or if you want to define the collapsed clusters according to the initial heatmap (using the numbers from the original run)

# 0: according to original cluster numbers  
# it comes with the advantage that you can define which original clusters are gonna be pooled in every step, mixing wildly and trying out differnt things
# 1: sequential: use numbers of the beforehand collapse step
# renumber in sequence means that the renumbered cluster numbers are recycled and you pool based on the pooled clusters beforehand
# Watch out. The cluster numbers are changing after every step, make sure you check that every step. This takes effect only AFTER the first collapse round
ClusterCollapsing <- 1 # 1 recommended: only sequential allows for automatic cluster collapse.


# for every subset, you have to decide how many clusters you want to collapse.
# this process happens automatically and the script sequentially collapses the two most similar clusters - one collapse of two clusters at a time
# Every collapse step is numbered: 
#level 0: no collapse,
#level n: Collapse the two most similar clusters and renumber the remaining clusters. Repeat n-1 times.
# For the first run, set it to "max": the script will run a full collapse into one cluster per subset. 
# Choose then the collapse level for each subset: 
# collapse.levels <- c(0, "max", "max") # maximally collapse to two clusters. Check which collapse level fits and re-run with according number supplied.
collapse.levels <- c(0, 0, 0) # <- do not collapse anything, just pool the subsets back together
#collapse.levels <- c(x, y) # collapse closest pair of clusters in subset 1 x times, subset 2 y times iteratively

# Benchmarking the homogeneity of the subsets: we use an accuracy measure to return how homogeneity is lost during collapsing:
#  RMSD; sqrt 1/n sum (x-m)2
#  MAE; 1/n sum abs x-m
clusterhomogeneitymeasure <- "MAE" # better use MAE, which is less sensitive to outlier signals

# Adrian Huck hat ein anderes paket: fviz_nbclust {factoextra}. kmeans, ist also vielleicht suboptimal.


# After collapsing, or not by keeping collapse.levels at 0, the script will annotate the clusters
# other than the original code, you need to explicitly annotate a cluster as "other"!!!
# any un-annotated cluster will keep the three-digit number by the collapse-routine
# if you checked the annotation list down in stage 3, turn annotation on:
use.my.annotation.list <- FALSE # setting this to FALSE will skip the annotation and will just use the three-digit numbers after re-uniting the subsets

# After the annotation, check how each cluster is distributed within the markercomparision we used to set gate thresholds
# there is a chance that some clusters end up at the edge of the threshold and are therefore wrongly put in that subset: 
# e.g a cluster that is "false CD3 gated" would now show up just at the border of the CD3 gate.
# Clusters might end up in certain areas in the scatter plot, that is how clustering works, but you don't want to see a cluster sticking at the gate threshold:
# you might need to address such a cluster by either adjusting the gate thresholds, or removing it entirely already in stage 2
plot.cluster_wise.scatterplots <- 0

# Cluster annotations list.
# Watch out: You will annotate the clusters according to the numbers found in the STAGE2 First(Second) Pass (depending on whether you deleted some)
# Or depending on the numbers you find at the according collapse.level. This number is NOT the one you put in the list, but
# you need to add 100 to subset 1, 200 to subset 2. The code will re-unite the subsets first, giving cluster 3 of subset 2 the number 203. 
# This is the annotated number:
cluster.annots <- list()















# if you wish to see a heatmap&UMAP calculated over the entire dataset (after we re-united the subsets with three-digit number)
# Keep on 0, not sure why we did all the fuss with the gates if you wanna see the entire dataset, but be my guest.
plot.global.heatmap.all.subsets <- 0 # I just cannot remember why this is actually an option in the code



########################## STAGE 4 SETTINGS ##########################

# before we go on, we need to mask lymph nodes in the ROIs, since these provide a ton of cells to any spatial analysis.

# For now, we only consider masks that behave well, meaning that the mask itself covers maximally a certain percentage of the whole ROI
# Any mask that covers more than that on the ROI is ignored without warning!
# This threshold holds for both the lymphnode binary mask, as well as the objects to which the EDM matrices are calculated
max.mask.coverage.percent <-  50 #[%] 25 


mask.lymphnodes <- 0 # yes: 1, no: 0
# Excluding the lymph node binary mask:
#lymphnode.binarymask.path <- c("D:/IMI/IMC/Western-diet/Analysis/Lymph Node ImageJ Mask/Lymphnode_EDM&BM_imageJ_output") # thats the new imageJ script 1.2.2
# You can either use a simple Binary Mask (creating inside.lymphnode flag) or you can use also here EDMs, that will give additional to inside.lymphnode a distance
lymphnode.mask.type <- "EDM" # use either EDM or BM. imageJ produces both for you...
# using EDM has the advantage that you can now extend the mask since you are probably missing a clear lymph node signal and the masks are not covering
# the t cell area that nicely 
# in this case, the script lets you extend the mask, and it will flag cells with a distance of x px also as inside.lymphnode:
lymphnode.extension <- 80 #60[px] extend the mask by these pixels to flag these cells as inside.lymphnode as well

# plot the flags on every ROI, please check these to make sure your extension is applicable!
plot.lymphnodeflags.on.ROIs <- 0




# only then we calculate distances to cellular masks, such as the islets in the ROIs:
run.EDM.calculation <- 1 # yes: 1, no: 0

# the batch process of imageJ stores jpgs, tiffs and the EDM matrix as .txt file in a folder. 
# supply here the absolute path of it
# v1.0 script produces masks based on the signal as we know it:
#EDMpath <- c("D:/IMI/IMC/Western-diet/Analysis/Islet ImageJ EDM/EDM_imageJ_output from v1.0 not extended original masks")

# in the discussion with Birgit 20220422, also seeing the glucagon signal in Juliettes panel, we decide to extend our mask before
# calculating EDM. These here are in v1.1.0 and are extended by 9 iterations on count 1px:
#EDMpath <- c("D:/IMI/IMC/Western-diet/Analysis/Islet ImageJ EDM/EDM_imageJ_output_5IterationsDilate")
EDMpath <- c("D:/IMI/IMC/IM - Matthaus Felsenstein/Analysis/Tumor ImageJ Mask/CK18_EDM&BM_imageJ_output")
EDM.name <- "TE" # used to name the column in gated.cell.dat: inside.*EDM.name*
# now, cross-check with the .jpgs in this folder the following: 
# find an empty mask, which is all-blue and open the EDM.txt file. 
# imageJ puts a number for empty (no mask present), which was in my case 46341 (besides some of them having 0) 
# Adapt if necessary since we gonna skip masks that got this number in:
imageJ.emptymask <- 46341

islet.rim.width <- 40 #[px] define a peri-mask rim of this size in pixel: inside.rim.*EDM.name* . turn off by setting it to 0























########################## STAGE 5 SETTINGS ##########################

# recommended to be run at least once:
# take the dimension reduction engine (UMAP, most probably) and split the dataset by ROI and by TMA
# use this to check if one ablation provides one cluster solely and needs to be taken out accordingly:
plot.DimRed.by.ROI <- 0
plot.DimRed.by.TMA <- 0

#plot the clusters on every ROI. You can choose to plot all subsets together on the ROI, or plot every gated subset separately
plot.all.subset.clusters.on.ROIs <- 0
plot.each.subset.clusters.on.ROIs <- 0 




########################## STAGE 6 SETTINGS ##########################

check.signaldrift.on.TMA <- 0 # check the corrected raw signals by TMA. not needed IMO since they could vary from condition to condion per se, so why check signal drift?

# Birigt plots:
#CD4.clusters <- c(205, 212, 207, 209, 202, 201) # normal threshold gate at CD3
#CD4.clusters <- c(210, 216, 203, 206, 202, 215,209,207,201) # polygon threshold gate at CD3
CD4.clusters <- c(210, 203) # polygon threshold gate at CD3, but now only the AICL+
#CD8.clusters <- c(203,212,211,214,204,213,205) # polygon threshold gate at CD3
CD8.clusters <- c(212,211,214) # polygon threshold gate at CD3, but now only the KLRF1+

# watch out! you gotta have to adapt the p.val filtering if you change these numbers in STAGE 6 code



########################## STAGE 7 SETTINGS  ##########################






