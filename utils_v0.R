# Advanced spatial 2 - modified functions.R changed to: Phenalyzer_utils.R on 20250328 since PCF is incomaptible with the old code

# Lets make a wrapper function to take care of all transformation and normalization in one go:
do.data.normalization <- function (dat, 
                                   use.cols, 
                                   do.transform = FALSE,
                                   cofactor = 5,
                                   do.minmax = FALSE,
                                   do.zscore = FALSE,
                                   new.min = 0, 
                                   new.max = 1, 
                                   append.name = "_rescaled") 
{  require("data.table")
  
  # Lets make sure the chose columns are numeric and make sense to process:
  value <- dat[, use.cols, with = FALSE]
  

  if (isFALSE(all(sapply(value, is.numeric)))) {
    message("It appears that one column in your dataset is non numeric")
    print(sapply(value, is.numeric))
    stop("Transformation/Normalization engine stopped")
  }
  
  
  # Transform first:
  if(do.transform==TRUE){
    #this is the light version of do.sinh
    message("Asinh transformation started...")
    value <- asinh(value/cofactor)
    
    if (length(use.cols) > 1) {
      names(value) <- paste0(names(value), "_t",cofactor)
    }
    if (length(use.cols) == 1) {
      names(value) <- paste0(use.cols, "_t",cofactor)
    }
    
    dat <- cbind(dat, value)
    
    # now if we did transform, the upcoming functions need to work with these newly generated columns
    # Hence, we gonna update the column names we supplied to the function:
    use.cols <- paste0(use.cols, "_t",cofactor)
    

   
    
  }# end transform
  
  
  # Normalize data to 0 to 1:
  if(do.minmax==TRUE){ 
    message("Min-max normalization started...")
    norm.fun <- function(x) {
      (x - min(x))/(max(x) - min(x)) * (new.max - new.min) + 
        new.min
    }
    value <- dat[, use.cols, with = FALSE]
    res <- as.data.table(lapply(value, norm.fun))
    names(res) <- paste0(names(res), "_m",new.min,"m",new.max)
    dat <- cbind(dat, res)

    # We deliberatly do NOT update use.cols now so that z-score also pulls the asinh transformed columns.
    # if you want, for any reason, do min-max AND z-score, you need to update use.cols right here.
    
  }# end min-max normalization
  
  
  # z-score data 
  if(do.zscore==TRUE){ 
    message("Z-score normalization started...")
    value <- apply( dat[, use.cols, with = FALSE], scale, MARGIN = 2)
    res <- as.data.table(value)
    names(res) <- paste0(names(res), "_z")

    dat <- cbind(dat, res)
 
  }# end z-score normalization 
  
  return(dat)
}   # end do.data.normalization function def





# lets make a wrapper function for the spill-over matrix correction of RAW data:
do.spill.over.corr <- function( dat = cell.dat,
                                spillover.matrix.path = spillover.matrix.path,
                                donor.col.name = "acq.chnl"
                                )
{
  
  # Because the code would scale the scaled spilled-in signals over and over again,
  # the function will create an object spill.over.ran only after it ran successfully, and will check its existence before starting.
  # This object does not exist yet, and it will persist after the first run, blocking a second spill-over correction:
  
  if( !exists("spill.over.ran")   ){
    
  # this is experimental and an attempt to get the spillover matrix incorporated, best would be the matrix that comes from the core unit:
  #https://bodenmillergroup.github.io/IMCDataAnalysis/spillover-correction.html
  
  message("\n ")
  
  spillover.mat <- read_excel(spillover.matrix.path)
  
  spillover.mat <- as.data.frame(spillover.mat)
  
  spillover.mat[is.na(spillover.mat)] <- 0

  # we now get rid of that percent counting if the max is 100. Instead we just divide by 1 here:
  spillover.mat[,-1] <- spillover.mat[,-1] / max(spillover.mat[,-1])
  # ...and we set the diagonal to 0 as well, we do not want to affect the own channel of course:
  spillover.mat[spillover.mat=="1"]<-0
  
    
  # now we will extract the isotope names:
  acquired.channels <- pull(spillover.mat[donor.col.name]) #spillover.mat["acq.chnl"]
  spotted.channels <- names(spillover.mat)[-1] 
  
  # and match the channel names onto spillover.mat:
  spillover.mat[donor.col.name] <- names(dat)[match ( acquired.channels, gsub(".*_","",names(dat)) )]
  names(spillover.mat) <- c(donor.col.name, names(dat)[match ( spotted.channels, gsub(".*_","",names(dat)) )])
  
  # finally clean out isotopes that were not part of the panel:
  spillover.mat <- spillover.mat[!is.na(names(spillover.mat))]
  spillover.mat <- spillover.mat %>% drop_na(acq.chnl)
  
  # and update the channel name vector that are part of the spillover matrix:
  acquired.channels <- pull(spillover.mat[donor.col.name]) #spillover.mat["acq.chnl"]
  spotted.channels <- names(spillover.mat)[-1] 

  
  
  

  
  
  
  
  
  pb <- progress_bar$new(format = "[:bar] :percent [Calculating total spillover | :eta]",
                         total = length(acquired.channels) , #*length(spotted.channels)
                         show_after=0, #allows to call it right way 
                         current = "|",    # Current bar character
                         clear = FALSE) # make it persist
  pb$tick(0) # call in the progress bar without any progress, just to show it 
  
  
  
# we will now run along all aquired channels, extract the cell segment signals of that channel and
# then create a dataframe with the relative spill-in from that channel based on the spillover.mat
# we then sum up all spilled over signals for each channel for all cells in a data.frame called total.signal.corr
  
# since we dont know the dimesions of that data.frame, lets create and fill it the first time with values i=1 
# afterwards just sum the existing object up:  
  total.signal.corr <- data.frame()
  i<-0
  for(a in acquired.channels){
    
    
    
    
    
    #and get its relative contribution into all other channels out:
    temp.spill.contribution <-  spillover.mat[  match ( a, spillover.mat$acq.chnl ) ,  ]
    temp.spill.out.channel <- dat[,a, with=F]
    
  # protect from working with a channel that is only party recorded in the dataset and would then yield in NA when calculating the spill-out in another channel:
    if( !any(is.na(temp.spill.out.channel))  ){
      

      if(i==0){
        total.signal.corr <- data.frame(mapply(`*`, temp.spill.contribution[-1],temp.spill.out.channel,SIMPLIFY=FALSE))
        i<-1
        pb$tick()
      }else{
        total.signal.corr <- total.signal.corr + data.frame(mapply(`*`, temp.spill.contribution[-1],temp.spill.out.channel,SIMPLIFY=FALSE))
        pb$tick()
      }
      
   }else{pb$tick()} # end protect from trying to calculate the spill in in other channels if spill out channel is not recorded on all ROIs
      
   
  } # end run with a along acquired channels                             
  

  # now there are some ideas as how to substract values in one df from the same column names in another
  #https://stackoverflow.com/questions/18708395/subtract-values-in-one-dataframe-from-another
  # but they need a unique column with unique entries by row, sth we do not have here: dat contains much more cols! 
  # due to missing ideas, lets do it the brutal way by substracting the spilled in signal column-by-column:
  
  pb <- progress_bar$new(format = "[:bar] :percent [Correcting total spillover | :eta]",
                         total = length(spotted.channels) , #*length(spotted.channels)
                         show_after=0, #allows to call it right way 
                         current = "|",    # Current bar character
                         clear = FALSE) # make it persist
  pb$tick(0) # call in the progress bar without any progress, just to show it 
  
  for(s in spotted.channels){

    dat[, s[1]  ] <-  dat[,s, with=F]  -  total.signal.corr[,s]
    pb$tick()
    
  }
  
  
  
  
  
  return(dat)
  
  
  }else{
    # if spill.over routine ran, the object spill.over.ran exists and we skip it:
    message(paste0("Spill-over correction ran already on ",data, " and calculation blocked!")   )  
  }
  
} # end function def spill over correction





# the old heatmap wrapper is gonna be depreciated soon. we use complex heatmap for this. Here are the functions we need:
fh = function(x) fastcluster::hclust(dist(x))
# use the faster function for dendrogram sorting:
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))


#heatmap wrapper saved the png by default, but somehow this doesnt work today (yday worked but I dunno what function was loaded.) 
# so we orverride the function here since the error comes back unused argument (file.name = paste0("somepath")):
# the rest is updated to the function found in v1.0.0; dendrograms.sort routine, which, funnily is turned off by default anyway..
make.pheatmap <- function (dat, 
                           sample.col, 
                           plot.cols, 
                           annot.cols = NULL, 
                           feature.annots = NULL, 
                           annotation_colors = NULL, 
                           plot.title = paste0(sample.col,  " heatmap"), 
                           # path = NULL, # this would give you the option to store it somewhere else. Its all commented out since we control that via outputdirs
                           file.name = NA, # if file.name is not provided, print rather than save
                           transpose = FALSE, 
                           normalise = TRUE, 
                           is.fold = FALSE, 
                           fold.range = NULL, 
                           dendrograms = "both", 
                           dendrograms.sort = FALSE, # if TRUE, will SORT dendrograms
                           cutree_rows = 1, 
                           cutree_cols = 1, 
                           row.sep = c(), 
                           col.sep = c(), 
                           cell.size = NA, # was 15 
                           standard.colours = "BuPu", 
                           fold.colours = "Spectre") 
{
  if (!is.element("pheatmap", installed.packages()[, 1])) 
    stop("pheatmap is required but not installed")
  if (!is.element("RColorBrewer", installed.packages()[, 1])) 
    stop("RColorBrewer is required but not installed")
  if (!is.element("scales", installed.packages()[, 1])) 
    stop("scales is required but not installed")
  require(pheatmap)
  require(RColorBrewer)
  require(scales)
  if (standard.colours == "BuPu") {
    colour.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                                 "BuPu"))(31))
  }
  if (standard.colours == "RdYlBu") {
    colour.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                                 "RdYlBu"))(31))
    colour.palette <- rev(colour.palette)
  }
  if (standard.colours == "rev(RdBu)") {
    colour.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                                 "RdBu"))(31))
    colour.palette <- rev(colour.palette)
  }
  if (standard.colours == "Blues") {
    colour.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                                 "Blues"))(31))
  }
  if (standard.colours == "Reds") {
    colour.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                                 "Reds"))(31))
  }
  if (standard.colours == "Greys") {
    colour.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                                 "Greys"))(31))
  }
  if (standard.colours == "YlGnBu") {
    colour.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                                 "YlGnBu"))(31))
  }
  if (standard.colours == "viridis") {
    colour.palette <- colorRampPalette(c((scales::viridis_pal(option = "viridis"))(50)))
    colour.palette <- colour.palette(31)
  }
  if (standard.colours == "spectral") {
    spectral.list <- colorRampPalette(RColorBrewer::brewer.pal(11, 
                                                               "Spectral"))(50)
    spectral.list <- rev(spectral.list)
    colour.palette <- colorRampPalette(c(spectral.list))
    colour.palette <- colour.palette(31)
  }
  if (standard.colours == "magma") {
    colour.palette <- colorRampPalette(c((scales::viridis_pal(option = "magma"))(50)))
    colour.palette <- colour.palette(31)
  }
  if (standard.colours == "inferno") {
    colour.palette <- colorRampPalette(c((scales::viridis_pal(option = "inferno"))(50)))
    colour.palette <- colour.palette(31)
  }
  if (fold.colours == "Spectre") {
    fold.palette <- colorRampPalette(rev(c("#ffeda0", "#fed976", 
                                           "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026", 
                                           "#800026", "black", "#023858", "#045a8d", "#0570b0", 
                                           "#3690c0", "#74a9cf", "#a6bddb", "#d0d1e6", "#ece7f2")))
    fold.palette <- fold.palette(31)
  }
  if (fold.colours == "BuPu") {
    fold.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                               "BuPu"))(31))
  }
  if (fold.colours == "RdYlBu") {
    fold.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                               "RdYlBu"))(31))
    fold.palette <- rev(fold.palette)
  }
  if (fold.colours == "rev(RdBu)") {
    fold.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                               "RdBu"))(31))
    fold.palette <- rev(fold.palette)
  }
  if (fold.colours == "Blues") {
    fold.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                               "Blues"))(31))
  }
  if (fold.colours == "Reds") {
    fold.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                               "Reds"))(31))
  }
  if (fold.colours == "Greys") {
    fold.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                               "Greys"))(31))
  }
  if (fold.colours == "YlGnBu") {
    fold.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                               "YlGnBu"))(31))
  }
  if (fold.colours == "viridis") {
    fold.palette <- colorRampPalette(c((scales::viridis_pal(option = "viridis"))(50)))
    fold.palette <- fold.palette(31)
  }
  if (fold.colours == "spectral") {
    spectral.list <- colorRampPalette(RColorBrewer::brewer.pal(11, 
                                                               "Spectral"))(50)
    spectral.list <- rev(spectral.list)
    fold.palette <- colorRampPalette(c(spectral.list))
    fold.palette <- fold.palette(31)
  }
  if (fold.colours == "magma") {
    fold.palette <- colorRampPalette(c((scales::viridis_pal(option = "magma"))(50)))
    fold.palette <- fold.palette(31)
  }
  if (fold.colours == "inferno") {
    fold.palette <- colorRampPalette(c((scales::viridis_pal(option = "inferno"))(50)))
    fold.palette <- fold.palette(31)
  }
  dat <- as.data.frame(dat)
  heatmap.data <- dat
  rownames(heatmap.data) <- t(dat[sample.col])
  heatmap.data
  if (is.null(annot.cols) == FALSE) {
    annot <- heatmap.data[annot.cols]
    heatmap.data <- heatmap.data[plot.cols]
    heatmap.data
  }
  ### Transpose (ONLY IF REQUIRED) -- the longest set (clusters or parameters) on x-axis -- 
  # by default MARKERS are columns, CLUSTERS are rows -- transpose to flip these defaults
  if (is.null(annot.cols) == TRUE) {
    annot <- NULL
    heatmap.data <- heatmap.data[plot.cols]
    heatmap.data
  }
  if (transpose == TRUE) {
    heatmap.data.t <- as.data.frame(t(heatmap.data))
    heatmap.data <- heatmap.data.t
  }
  ### NORMALISE BY COLUMN (i.e. each column/parameter has a max of 1 and a minimum of 0) # This is optional, but allows for better comparison between markers
  if (normalise == TRUE) {
    if (is.fold == FALSE) {
      row.nam <- row.names(heatmap.data)
      col.nam <- names(heatmap.data)
      norm.fun <- function(x) {
        (x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - 
                                      min(x, na.rm = TRUE))
      }
      heatmap.data.norm <- as.data.frame(lapply(heatmap.data, 
                                                norm.fun))
      names(heatmap.data.norm) <- col.nam
      max(heatmap.data.norm)
      heatmap.data.norm <- as.matrix(heatmap.data.norm)
      heatmap.data <- heatmap.data.norm
      rownames(heatmap.data) <- row.nam
    }
  }#end normalize=T
  
  heatmap.data <- as.matrix(heatmap.data)
  
  ### Set up clustering
  if (dendrograms == "none") {
    row.clustering <- FALSE
    col.clustering <- FALSE
  }
  if (dendrograms != "none") {
    # set the custom distance and clustering functions, per your example
    hclustfunc <- function(x) hclust(x, method = "complete")
    distfunc <- function(x) dist(x, method = "euclidean")
    
    # perform clustering on rows and columns
    if (dendrograms == "both") {
      row.clustering <- TRUE
      col.clustering <- TRUE
      
      #this is the new block as by 0.5.4->1.0.0
      if(isTRUE(dendrograms.sort)){
        
        row.clustering <- hclustfunc(distfunc(heatmap.data))
        col.clustering <- hclustfunc(distfunc(t(heatmap.data)))
        
        require(dendsort)
        sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
        
        row.clustering <- sort_hclust(row.clustering)
        col.clustering <- sort_hclust(col.clustering)
      }
      #end new block
      
    }
    if (dendrograms == "column") {
      row.clustering <- FALSE
      col.clustering <- TRUE
      
      
      #this is the new block as by 0.5.4->1.0.0
      if(isTRUE(dendrograms.sort)){
        
        col.clustering <- hclustfunc(distfunc(t(heatmap.data)))
        
        require(dendsort)
        sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
        
        col.clustering <- sort_hclust(col.clustering)
      }
      #end new block
      
      
    }
    if (dendrograms == "row") {
      row.clustering <- TRUE
      col.clustering <- FALSE
      
      
      #this is the new block as by 0.5.4->1.0.0
      if(isTRUE(dendrograms.sort)){
        
        row.clustering <- hclustfunc(distfunc(t(heatmap.data)))
        
        require(dendsort)
        sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
        
        row.clustering <- sort_hclust(row.clustering)
      }
      #end new block
    }
  }
  if (is.fold == TRUE) {
    map.colour <- fold.palette
    sym.key <- FALSE
    sym.breaks <- TRUE
    if (is.null(fold.range)) {
      fld.max <- max(heatmap.data, na.rm = TRUE)
      fld.min <- min(heatmap.data, na.rm = TRUE)
      if (fld.max == -fld.min) {
        fold.max.range <- fld.max
        fold.min.range <- fld.min
      }
      if (fld.max > -fld.min) {
        fold.max.range <- fld.max
        fold.min.range <- -fld.max
      }
      if (fld.max < -fld.min) {
        fold.max.range <- -fld.min
        fold.min.range <- fld.min
      }
    }
    if (!is.null(fold.range)) {
      fold.max.range <- fold.range[1]
      fold.min.range <- fold.range[2]
    }
    my.breaks <- seq(fold.min.range, fold.max.range, length.out = 32)
  }
  if (is.fold == FALSE) {
    map.colour <- colour.palette
    sym.key <- FALSE
    sym.breaks <- FALSE
    heatmap.data
    my.max <- function(x) ifelse(!all(is.na(x)), max(x, 
                                                     na.rm = T), NA)
    my.min <- function(x) ifelse(!all(is.na(x)), min(x, 
                                                     na.rm = T), NA)
    my.breaks <- seq(my.min(heatmap.data), my.max(heatmap.data), 
                     length.out = 32)
  }
  scale.set <- "none"
  title.text <- plot.title
  
  
  
  # this here is a chunk of code deleted in v1.0, but I liked it and so I bring it back:
  # for this to work, you need to supply path and you need to uncomment path= in the function input
  # Specify directory heatmap will be saved
  # if(is.null(path)){  flnm <- file.name }
  #  if(!is.null(path)){ flnm <- paste0(path, '/', file.name)  }
  
  # this is a bit stupid, but I want it pushed into storage and printed. so we need to pipe it to pheatmap twice:
  
  pheatmap::pheatmap(mat = as.matrix(heatmap.data), main = title.text, 
                     cellwidth = cell.size, cellheight = cell.size, cluster_rows = row.clustering, 
                     cluster_cols = col.clustering, breaks = my.breaks, cutree_rows = cutree_rows, 
                     cutree_cols = cutree_cols, gaps_row = row.sep, gaps_col = col.sep, 
                     annotation_row = annot, annotation_col = feature.annots, 
                     annotation_colors = annotation_colors, color = map.colour ,
                     filename = file.name # I brought this back in to control the filename when calling the funtion
  )
  
  
  
  
} # end function def make.pheatmap




make.colour.plot <- function (dat, x.axis, y.axis, 
                              point.alpha,  # I added an alpha channel for the geom_point for col.type=factor coloring plots
                              col.axis = NULL, col.type = "continuous", 
                              add.label = FALSE, hex = FALSE, hex.bins = 30, 
                              colours = "spectral", 
                              polychromepalette = NULL, # this is the go-to way to color, you need to provide a table with colors via polychrome
                              col.min.threshold = 0.01, col.max.threshold = 0.995, 
                              align.xy.by = dat, 
                              align.col.by = dat, regression.line = NULL, 
                              titlestring = col.axis, # that used to be title, but I changed the code to allow title and subtitle:
                              subtitlestring = NULL,
                              filename = NULL, dot.size = 1, plot.width = 9, plot.height = 7, 
                              nudge_x = 0.5, # <- depreciated, as its calculated dynamically from the X-axis range
                              nudge_y = 0.5, # <- depreciated...
                              square = TRUE, 
                              legend.loc = "right", 
                              legend.text.size = 18, # was 18 default
                              save.to.disk = TRUE, path = getwd(), blank.axis = FALSE,
                              axis.title.size = 15, #font size of the axis titles. new, was 28
                              axis.text.size = 13, #font size of the axis numbers. new, was 24
                              title.size = 16 #font size of the title, new, was 32
) 
{
  if (!is.element("ggplot2", installed.packages()[, 1])) 
    stop("ggplot2 is required but not installed")
  if (!is.element("scales", installed.packages()[, 1])) 
    stop("scales is required but not installed")
  if (!is.element("colorRamps", installed.packages()[, 1])) 
    stop("colorRamps is required but not installed")
  if (!is.element("ggthemes", installed.packages()[, 1])) 
    stop("ggthemes is required but not installed")
  if (!is.element("RColorBrewer", installed.packages()[, 1])) 
    stop("RColorBrewer is required but not installed")
  require(ggplot2)
  require(scales)
  require(colorRamps)
  require(ggthemes)
  require(RColorBrewer)
  if (hex == TRUE) {
    if (is.null(col.axis)) {
      message("Note: hex bins do not currently work for density plots, only for colour plots when col.axis is specified and can be plotted as a continuous numeric variable")
    }
    if (!is.null(col.axis)) {
      if (!is.numeric(dat[[col.axis]])) {
        stop("Sorry, hex bins only work when col.type is specified, and can be plotted as a continuous numeric variable")
      }
    }
  }
  if (!is.null(col.axis)) {
    if (col.type == "continuous") {
      if (!is.numeric(dat[[col.axis]])) {
        message("Non-numeric values detected in col.axis -- using col.type = 'factor'")
        col.type <- "factor"
      }
    }
    if (col.type == "factor") {
      if (length(unique(as.factor(dat[[col.axis]]))) > 200) {
        message("Over 200 factors detected, using continuous scale instead of a factor scale")
        col.type <- "continuous"
      }
    }
  }
  if (colours == "jet") {
    colour.scheme <- colorRampPalette(c("#00007F", "blue", 
                                        "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", 
                                        "red", "#7F0000"))
  }
  if (colours == "spectral") {
    spectral.list <- colorRampPalette(brewer.pal(11, "Spectral"))(50)
    spectral.list <- rev(spectral.list)
    colour.scheme <- colorRampPalette(c(spectral.list))
  }
  if (colours == "viridis") {
    colour.scheme <- colorRampPalette(c(viridis_pal(option = "viridis")(50)))
  }
  if (colours == "inferno") {
    colour.scheme <- colorRampPalette(c(viridis_pal(option = "inferno")(50)))
  }
  if (colours == "magma") {
    colour.scheme <- colorRampPalette(c(viridis_pal(option = "magma")(50)))
  }
  if (colours == "BuPu") {
    colour.list <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                              "BuPu"))(31))
    colour.scheme <- colorRampPalette(c(colour.list))
  }
  if (colours == "turbo") {
    colour.scheme <- colorRampPalette(c(viridis_pal(option = "turbo")(50)))
  }
  if (colours == "mako") {
    colour.scheme <- colorRampPalette(c(viridis_pal(option = "mako")(50)))
  }
  if (colours == "rocket") {
    colour.scheme <- colorRampPalette(c(viridis_pal(option = "rocket")(50)))
  }
  if (is.null(align.xy.by)) {
    Xmax <- max(dat[[x.axis]])
    Xmin <- min(dat[[x.axis]])
  }
  else {
    Xmax <- max(align.xy.by[[x.axis]])
    Xmin <- min(align.xy.by[[x.axis]])
  }
  if (is.null(align.xy.by)) {
    Ymax <- max(dat[[y.axis]])
    Ymin <- min(dat[[y.axis]])
  }
  else {
    Ymax <- max(align.xy.by[[y.axis]])
    Ymin <- min(align.xy.by[[y.axis]])
  }
  if (!is.null(col.axis)) {
    if (col.type == "continuous") {
      if (is.null(align.col.by)) {
        ColrMin <- quantile(dat[[col.axis]], probs = c(col.min.threshold), 
                            na.rm = TRUE)
        ColrMax <- quantile(dat[[col.axis]], probs = c(col.max.threshold), 
                            na.rm = TRUE)
      }
      else {
        ColrMin <- quantile(align.col.by[[col.axis]], 
                            probs = c(col.min.threshold), na.rm = TRUE)
        ColrMax <- quantile(align.col.by[[col.axis]], 
                            probs = c(col.max.threshold), na.rm = TRUE)
      }
    }
    if (col.type == "factor") {
      if (is.null(align.col.by)) {
        colRange <- unique(dat[[col.axis]])
        colRange <- colRange[order(colRange)]
        colRange <- as.character(colRange)
      }
      else {
        colRange <- unique(align.col.by[[col.axis]])
        colRange <- colRange[order(colRange)]
        colRange <- as.character(colRange)
      }
    }
  }
  if (!is.null(col.axis)) {
    if (col.type == "continuous") {
      p <- ggplot(data = dat, aes(x = .data[[x.axis]], 
                                  y = .data[[y.axis]], colour = .data[[col.axis]]))
      if (hex == TRUE) {
        p <- p + stat_summary_hex(aes(z = dat[[col.axis]]), 
                                  fun = "mean", bins = hex.bins)
        p <- p + scale_fill_gradientn(colours = c(colour.scheme(50)), 
                                      limits = c(ColrMin, ColrMax), oob = squish)
      }
      else {
        p <- p + geom_point(size = dot.size)
        p <- p + scale_colour_gradientn(colours = colour.scheme(50), 
                                        limits = c(ColrMin, ColrMax), oob = squish, 
                                        na.value = "grey50")
      }
    }
    else if (col.type == "factor") {
      # these are tSNEs and UMAPs. lets not overdraw the points here:
      
      #this specral thing is annoying since you cannot see the clusters and the colors get re-assigned everytime we print another subset.
      #so lets predefine a palette and then use it in here:
      if (colours == "polychrome") {
        p <- ggplot(data = dat, aes(x = .data[[x.axis]], 
                                    y = .data[[y.axis]], colour = as.factor(.data[[col.axis]]))) + 
          geom_point(size = dot.size, alpha = point.alpha)+
          scale_colour_manual( values = polychromepalette )+ # push the polychrome palette in
          guides(colour = guide_legend(override.aes = list(size=4, alpha=1),  #we also gonna override the dot size how its depicted in the legend
                                       byrow=TRUE) ) # and I want the elements listed row-by-row
        
      } else {
        
        
        
        p <- ggplot(data = dat, aes(x = .data[[x.axis]], 
                                    y = .data[[y.axis]], colour = as.factor(.data[[col.axis]]))) + 
          geom_point(size = dot.size, alpha = point.alpha) + 
          lims(colour = colRange)+
          guides(colour = guide_legend(override.aes = list(size=4) ) ) #we also gonna override the dot size how its depicted in the legend
        
      }#end default factor plot with col.axis set
      
    }
    
    
    
  }
  if (is.null(col.axis)) {
    p <- ggplot(data = dat, aes(x = .data[[x.axis]], y = .data[[y.axis]])) + 
      ggpointdensity::geom_pointdensity(size = dot.size)
    if (colours == "viridis" || colours == "magma" || colours == 
        "inferno") {
      p <- p + viridis::scale_colour_viridis(option = colours)
    }
    else if (colours == "jet") {
      p <- p + ggplot2::scale_colour_gradientn(colours = c("#00007F", 
                                                           "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", 
                                                           "#FF7F00", "red", "#7F0000"))
    }
    else if (colours == "spectral") {
      p <- p + ggplot2::scale_colour_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, 
                                                                                                       "Spectral"))(50)))
    }
    else if (colours == "BuPu") {
      colour.list <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                                "BuPu"))(31))
      p <- p + ggplot2::scale_colour_gradientn(colours = colour.list)
    }
  }
  
  if (!is.null(regression.line)) {
    p <- p + geom_smooth(method = regression.line)
  }
  if (is.null(title)) {
    title <- "Density"
  }
  p <- p +  labs(title= titlestring , 
                 subtitle=subtitlestring#, 
                 #caption="Created by M.Barone"#, 
                 # y="ROI y", 
                 # x="ROI x"
  )
  
  
  p <- p + scale_x_continuous(breaks = scales::pretty_breaks(n = 8), 
                              name = x.axis, limits = c(Xmin, Xmax))
  p <- p + scale_y_continuous(breaks = scales::pretty_breaks(n = 8), 
                              name = y.axis, limits = c(Ymin, Ymax))
  if (col.type == "continuous") {
    p <- p + theme(panel.background = element_rect(fill = "white", 
                                                   colour = "black", size = 0.5), 
                   axis.title.x = element_text(color = "Black", size = 28), 
                   axis.title.y = element_text(color = "Black", size = 28), axis.text.x = element_text(color = "Black", size = 24), 
                   axis.text.y = element_text(color = "Black",  size = 24), 
                   panel.border = element_rect(colour = "black",  fill = NA, size = 2), 
                   plot.title = element_text(color = "Black",  face = "bold", size = 32, hjust = 0)
    )
  }
  if (col.type == "factor") {
    p <- p + theme(panel.background = element_rect(fill = "white", 
                                                   colour = "black", size = 0.5), 
                   axis.title.x = element_text(color = "Black", size = axis.title.size), 
                   axis.title.y = element_text(color = "Black", size = axis.title.size), 
                   axis.text.x = element_text(color = "Black",  size = axis.text.size), 
                   axis.text.y = element_text(color = "Black",  size = axis.text.size), 
                   panel.border = element_rect(colour = "black",  fill = NA, size = 2), 
                   plot.title = element_text(color = "Black", size = title.size) # was also ,face = "bold", hjust = 0
                   #plot.subtitle = element_text(color = "Black", size = subtitle.size)
    )
  }
  if (square == TRUE) {
    p <- p + theme(aspect.ratio = 1)
  }
  if (legend.loc %in% c("top", "bottom")) {
    p <- p + theme(legend.direction = "horizontal", legend.position = legend.loc, 
                   legend.text = element_text(size = legend.text.size), legend.title = element_blank())
  }
  if (legend.loc %in% c("left", "right")) {
    p <- p + theme(legend.direction = "vertical", legend.position = legend.loc, 
                   legend.text = element_text(size = legend.text.size), legend.title = element_blank())
  }
  if (col.type == "factor") {
    if (add.label == TRUE) {
      if (is.numeric(dat[[col.axis]])) {
        centroidX = tapply(dat[[x.axis]], dat[[col.axis]], 
                           median)
        centroidY = tapply(dat[[y.axis]], dat[[col.axis]], 
                           median)
        centroidCol = tapply(dat[[col.axis]], dat[[col.axis]], 
                             median)
        centroidsDf <- data.frame(centroidX, centroidY, 
                                  centroidCol)
      }
      if (!is.numeric(dat[[col.axis]])) {
        labels <- sort(unique(dat[[col.axis]]))
        centroidsDf <- data.frame(centroidX = tapply(dat[[x.axis]], 
                                                     dat[[col.axis]], median), centroidY = tapply(dat[[y.axis]], 
                                                                                                  dat[[col.axis]], median), centroidCol = labels)
      }
      #https://ggrepel.slowkow.com/reference/geom_text_repel.html
      #https://ggrepel.slowkow.com/articles/examples.html
      p <- p + geom_label_repel(data = centroidsDf, 
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
                                
                                max.overlaps = Inf, # default 10
                                
                                aes(x = centroidX, y = centroidY, label = centroidCol), #, alpha = 0.5
                                fill = "white",
                                col = "black", fontface = "bold",size = 5,
                                verbose = TRUE) 
      
      
      
      # Put geom_point() of the labels after geom_label_repel, so that its point is on top layer
      p <- p + geom_point(data = centroidsDf, aes(x = centroidX, 
                                                  y = centroidY), col = "black", size = 2 , alpha = 0.3)
      
      p <- p +coord_cartesian(clip = "off") #this ensures that the label is not clipped off
      
      #p <- p + guides(alpha = "none")
      

    }
  }
  if (blank.axis == TRUE) {
    p <- p + theme(axis.line = element_blank(), axis.text.x = element_blank(), 
                   axis.text.y = element_blank(), axis.ticks = element_blank(), 
                   axis.title.x = element_blank(), axis.title.y = element_blank(), 
                   panel.grid.major = element_blank(), panel.background = element_blank(), 
                   panel.border = element_blank(), panel.grid.minor = element_blank(), 
                   plot.background = element_blank(), )
  }
  if (save.to.disk == TRUE) {
    if (!is.null(col.axis)) {
      if (col.type == "continuous") {
        lb <- "Colour"
      }
      if (col.type == "factor") {
        lb <- "Factor"
      }
    }
    if (is.null(col.axis)) {
      lb <- "Density plot"
    }
    if (is.null(filename)) {
      filename <- paste0(lb, " plot - ", title, " - plotted on ", 
                         x.axis, " by ", y.axis, ".png")
    }
    ggsave(filename = filename, plot = p, path = path, width = plot.width, 
           height = plot.height, limitsize = FALSE)
  }
  else {
    # watch out, there is a cryptic error message if the Plots tab in RStudio is too small to send the plot to the viewer
    # so protect from doing that:

    if (min(dev.size("in")) > 2.2   ) {  # Minimum width & height of 4 inches
      print(p)  
    } else {
      message("Plot window is too small. Resize the RStudio plot window to see the plot.")
    }
    
  
  }
  return(p)
} # end function def make.colour.plot








make.colour.plot.adapted <- function (dat, x.axis, y.axis, 
                              col.axis = NULL, col.type = "continuous", 
                              add.label = FALSE, hex = FALSE, hex.bins = 30, 
                              colours = "spectral", 
                              polychromepalette = NULL, # this is the go-to way to color, you need to provide a table with colors via polychrome
                              NA.color = "#F1F1F1",
                              min.line.length.tolabel = 0,
                              
                              col.min.threshold = 0.01, col.max.threshold = 0.995, 
                              point.alpha,  # I added an alpha channel for the geom_point for col.type=factor coloring plots
                              align.xy.by = dat, 
                              align.col.by = dat, regression.line = NULL, 
                              titlestring = col.axis, # that used to be title, but I changed the code to allow title and subtitle:
                              subtitlestring = NULL,
                              filename = NULL, dot.size = 1, plot.width = 9, plot.height = 7, 
                              nudge_x = 0.5, # <- depreciated, as its calculated dynamically from the X-axis range
                              nudge_y = 0.5, # <- depreciated...
                              draw.repel.labels="yes",
                              repel.be.verbose = FALSE,
                              repel.label.size = 4, # repel labels in UMAP  
                              label.label.repelforce = 30,
                              label.datapoint.pullforce = 0.8,
                              square = TRUE, 
                              legend.loc = "right", 
                              legend.text.size = 18, # was 18 default
                              save.to.disk = TRUE, path = getwd(), blank.axis = FALSE,
                              axis.title.size = 15, #font size of the axis titles. new, was 28
                              axis.text.size = 13, #font size of the axis numbers. new, was 24
                              title.size = 16 #font size of the title, new, was 32
) 
{
  if (!is.element("ggplot2", installed.packages()[, 1])) 
    stop("ggplot2 is required but not installed")
  if (!is.element("scales", installed.packages()[, 1])) 
    stop("scales is required but not installed")
  if (!is.element("colorRamps", installed.packages()[, 1])) 
    stop("colorRamps is required but not installed")
  if (!is.element("ggthemes", installed.packages()[, 1])) 
    stop("ggthemes is required but not installed")
  if (!is.element("RColorBrewer", installed.packages()[, 1])) 
    stop("RColorBrewer is required but not installed")
  require(ggplot2)
  require(scales)
  require(colorRamps)
  require(ggthemes)
  require(RColorBrewer)
  if (hex == TRUE) {
    if (is.null(col.axis)) {
      message("Note: hex bins do not currently work for density plots, only for colour plots when col.axis is specified and can be plotted as a continuous numeric variable")
    }
    if (!is.null(col.axis)) {
      if (!is.numeric(dat[[col.axis]])) {
        stop("Sorry, hex bins only work when col.type is specified, and can be plotted as a continuous numeric variable")
      }
    }
  }
  if (!is.null(col.axis)) {
    if (col.type == "continuous") {
      if (!is.numeric(dat[[col.axis]])) {
        message("Non-numeric values detected in col.axis -- using col.type = 'factor'")
        col.type <- "factor"
      }
    }
    if (col.type == "factor") {
      if (length(unique(as.factor(dat[[col.axis]]))) > 200) {
        message("Over 200 factors detected, using continuous scale instead of a factor scale")
        col.type <- "continuous"
      }
    }
  }
  if (colours == "jet") {
    colour.scheme <- colorRampPalette(c("#00007F", "blue", 
                                        "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", 
                                        "red", "#7F0000"))
  }
  if (colours == "spectral") {
    spectral.list <- colorRampPalette(brewer.pal(11, "Spectral"))(50)
    spectral.list <- rev(spectral.list)
    colour.scheme <- colorRampPalette(c(spectral.list))
  }
  if (colours == "viridis") {
    colour.scheme <- colorRampPalette(c(viridis_pal(option = "viridis")(50)))
  }
  if (colours == "inferno") {
    colour.scheme <- colorRampPalette(c(viridis_pal(option = "inferno")(50)))
  }
  if (colours == "magma") {
    colour.scheme <- colorRampPalette(c(viridis_pal(option = "magma")(50)))
  }
  if (colours == "BuPu") {
    colour.list <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                              "BuPu"))(31))
    colour.scheme <- colorRampPalette(c(colour.list))
  }
  if (colours == "turbo") {
    colour.scheme <- colorRampPalette(c(viridis_pal(option = "turbo")(50)))
  }
  if (colours == "mako") {
    colour.scheme <- colorRampPalette(c(viridis_pal(option = "mako")(50)))
  }
  if (colours == "rocket") {
    colour.scheme <- colorRampPalette(c(viridis_pal(option = "rocket")(50)))
  }
  if (is.null(align.xy.by)) {
    Xmax <- max(dat[[x.axis]])
    Xmin <- min(dat[[x.axis]])
  }
  else {
    Xmax <- max(align.xy.by[[x.axis]])
    Xmin <- min(align.xy.by[[x.axis]])
  }
  if (is.null(align.xy.by)) {
    Ymax <- max(dat[[y.axis]])
    Ymin <- min(dat[[y.axis]])
  }
  else {
    Ymax <- max(align.xy.by[[y.axis]])
    Ymin <- min(align.xy.by[[y.axis]])
  }
  if (!is.null(col.axis)) {
    if (col.type == "continuous") {
      if (is.null(align.col.by)) {
        ColrMin <- quantile(dat[[col.axis]], probs = c(col.min.threshold), 
                            na.rm = TRUE)
        ColrMax <- quantile(dat[[col.axis]], probs = c(col.max.threshold), 
                            na.rm = TRUE)
      }
      else {
        ColrMin <- quantile(align.col.by[[col.axis]], 
                            probs = c(col.min.threshold), na.rm = TRUE)
        ColrMax <- quantile(align.col.by[[col.axis]], 
                            probs = c(col.max.threshold), na.rm = TRUE)
      }
    }
    if (col.type == "factor") {
      if (is.null(align.col.by)) {
        colRange <- unique(dat[[col.axis]])
        colRange <- colRange[order(colRange)]
        colRange <- as.character(colRange)
      }
      else {
        colRange <- unique(align.col.by[[col.axis]])
        colRange <- colRange[order(colRange)]
        colRange <- as.character(colRange)
      }
    }
  }
  if (!is.null(col.axis)) {
    if (col.type == "continuous") {
      p <- ggplot(data = dat, aes(x = .data[[x.axis]], 
                                  y = .data[[y.axis]], colour = .data[[col.axis]]))
      if (hex == TRUE) {
        p <- p + stat_summary_hex(aes(z = dat[[col.axis]]), 
                                  fun = "mean", bins = hex.bins)
        p <- p + scale_fill_gradientn(colours = c(colour.scheme(50)), 
                                      limits = c(ColrMin, ColrMax), oob = squish)
      }
      else {
        p <- p + geom_point(size = dot.size)
        p <- p + scale_colour_gradientn(colours = colour.scheme(50), 
                                        limits = c(ColrMin, ColrMax), oob = squish, 
                                        na.value = "grey50")
      }
    }
    else if (col.type == "factor") {
      # these are tSNEs and UMAPs. lets not overdraw the points here:
      
      #this specral thing is annoying since you cannot see the clusters and the colors get re-assigned everytime we print another subset.
      #so lets predefine a palette and then use it in here:
      if (colours == "polychrome") {
        if (legend.loc %in% c("none")) {
          p <- ggplot(data = dat, aes(x = .data[[x.axis]], 
                                      y = .data[[y.axis]], colour = as.factor(.data[[col.axis]]))) + 
            geom_point(size = dot.size, alpha = point.alpha)+
            scale_colour_manual( values = polychromepalette,
                                 na.value = NA.color) 
          
        }else{
          p <- ggplot(data = dat, aes(x = .data[[x.axis]], 
                                      y = .data[[y.axis]], colour = as.factor(.data[[col.axis]]))) + 
            geom_point(size = dot.size, alpha = point.alpha)+
            scale_colour_manual( values = polychromepalette ,
                                 na.value = NA.color)+ # push the polychrome palette in
            guides(colour = guide_legend(override.aes = list(size=4, alpha=1),  #we also gonna override the dot size how its depicted in the legend
                                         byrow=TRUE) ) # and I want the elements listed row-by-row
        }
      } else {
        
        
        
        p <- ggplot(data = dat, aes(x = .data[[x.axis]], 
                                    y = .data[[y.axis]], colour = as.factor(.data[[col.axis]]))) + 
          geom_point(size = dot.size, alpha = point.alpha) + 
          lims(colour = colRange)+
          guides(colour = guide_legend(override.aes = list(size=4) ) ) #we also gonna override the dot size how its depicted in the legend
        
      }#end default factor plot with col.axis set
      
    }
    
    
    
  }
  if (is.null(col.axis)) {
    p <- ggplot(data = dat, aes(x = .data[[x.axis]], y = .data[[y.axis]])) + 
      ggpointdensity::geom_pointdensity(size = dot.size)
    if (colours == "viridis" || colours == "magma" || colours == 
        "inferno") {
      p <- p + viridis::scale_colour_viridis(option = colours)
    }
    else if (colours == "jet") {
      p <- p + ggplot2::scale_colour_gradientn(colours = c("#00007F", 
                                                           "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", 
                                                           "#FF7F00", "red", "#7F0000"))
    }
    else if (colours == "spectral") {
      p <- p + ggplot2::scale_colour_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, 
                                                                                                       "Spectral"))(50)))
    }
    else if (colours == "BuPu") {
      colour.list <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                                "BuPu"))(31))
      p <- p + ggplot2::scale_colour_gradientn(colours = colour.list)
    }
  }
  
  if (!is.null(regression.line)) {
    p <- p + geom_smooth(method = regression.line)
  }
  if (is.null(title)) {
    title <- "Density"
  }
  p <- p +  labs(title= titlestring , 
                 subtitle=subtitlestring#, 
                 #caption="Created by M.Barone"#, 
                 # y="ROI y", 
                 # x="ROI x"
  )
  
  
  p <- p + scale_x_continuous(breaks = scales::pretty_breaks(n = 8), 
                              name = x.axis, limits = c(Xmin, Xmax))
  p <- p + scale_y_continuous(breaks = scales::pretty_breaks(n = 8), 
                              name = y.axis, limits = c(Ymin, Ymax))
  if (col.type == "continuous") {
    p <- p + theme(panel.background = element_rect(fill = "white", 
                                                   colour = "black", size = 0.5), 
                   axis.title.x = element_text(color = "Black", size = 28), 
                   axis.title.y = element_text(color = "Black", size = 28), axis.text.x = element_text(color = "Black", size = 24), 
                   axis.text.y = element_text(color = "Black",  size = 24), 
                   panel.border = element_rect(colour = "black",  fill = NA, size = 2), 
                   plot.title = element_text(color = "Black",  face = "bold", size = 32, hjust = 0)
    )
  }
  if (col.type == "factor") {
    p <- p + theme(panel.background = element_rect(fill = "white", 
                                                   colour = "black", size = 0.5), 
                   axis.title.x = element_text(color = "Black", size = axis.title.size), 
                   axis.title.y = element_text(color = "Black", size = axis.title.size), 
                   axis.text.x = element_text(color = "Black",  size = axis.text.size), 
                   axis.text.y = element_text(color = "Black",  size = axis.text.size), 
                   panel.border = element_rect(colour = "black",  fill = NA, size = 2), 
                   plot.title = element_text(color = "Black", size = title.size) # was also ,face = "bold", hjust = 0
                   #plot.subtitle = element_text(color = "Black", size = subtitle.size)
    )
  }
  if (square == TRUE) {
    p <- p + theme(aspect.ratio = 1)
  }
  if (legend.loc %in% c("top", "bottom")) {
    p <- p + theme(legend.direction = "horizontal", legend.position = legend.loc, 
                   legend.text = element_text(size = legend.text.size), legend.title = element_blank())
  }
  if (legend.loc %in% c("left", "right")) {
    p <- p + theme(legend.direction = "vertical", legend.position = legend.loc, 
                   legend.text = element_text(size = legend.text.size), legend.title = element_blank())
  }
  if (legend.loc %in% c("none")) {
    p <- p + theme(legend.direction = "vertical", legend.position = legend.loc, 
                   legend.text = element_text(size = legend.text.size), legend.title = element_blank())
  }
  if (col.type == "factor") {
    if (add.label == TRUE) {
      if (is.numeric(dat[[col.axis]])) {
        centroidX = tapply(dat[[x.axis]], dat[[col.axis]], 
                           median)
        centroidY = tapply(dat[[y.axis]], dat[[col.axis]], 
                           median)
        centroidCol = tapply(dat[[col.axis]], dat[[col.axis]], 
                             median)
        centroidsDf <- data.frame(centroidX, centroidY, 
                                  centroidCol)
      }
      if (!is.numeric(dat[[col.axis]])) {
        labels <- sort(unique(dat[[col.axis]]))
        centroidsDf <- data.frame(centroidX = tapply(dat[[x.axis]], 
                                                     dat[[col.axis]], median), centroidY = tapply(dat[[y.axis]], 
                                                                                                  dat[[col.axis]], median), centroidCol = labels)
      }
      
      
      
      
      if(draw.repel.labels=="yes"){
      #https://ggrepel.slowkow.com/reference/geom_text_repel.html
      #https://ggrepel.slowkow.com/articles/examples.html
      p <- p + geom_label_repel(data = centroidsDf, 
                                hjust = "right", # this only takes effect initially and is lost for labels that are pulled...
                                force = label.label.repelforce, # repulsion between overlapping text labels. default 1
                                force_pull = label.datapoint.pullforce, # attraction betweenlabel and datapoint. default 1
                                # nudge_x = 0.1*(Xmax-Xmin), 
                                #nudge_y = 0.05*(Ymax-Ymin), 
                                xlim = c(NA, NA), # repel away from the edges of the plot, #c(-Inf, Inf), # either plots them also outside the graph or restrains the labels to be within the datarange. this assures that no label is cut off
                                ylim = c(NA, NA), # repel away from the edges of the plot, c(-Inf, Inf), # either plots them also outside the graph or restrains the labels to be within the datarange. this assures that no label is cut off
                                box.padding = 0.01, # additional padding around each text label 0.25 default
                                label.padding = 0.20, #0.25 default
                                point.padding = 0, # additional padding around each point
                                
                                min.segment.length = min.line.length.tolabel, 
                                segment.curvature = -0.1, #pos: more righthand, 0 straight, neg increase left-hand
                                segment.ncp = 3, # control points per curve
                                segment.angle = 20,
                                max.overlaps = Inf,
                                aes(x = centroidX, y = centroidY, label = centroidCol), #, alpha = 0.5
                                fill = "white",
                                col = "black", fontface = "bold", size = repel.label.size,
                                verbose = repel.be.verbose) 
      
      }
      
      
      
      # Put geom_point() of the labels after geom_label_repel, so that its point is on top layer
      if (colours == "polychrome") {
        if(draw.repel.labels=="yes"){
        #in case we supply the palette, we can even match the label dot color. for that we need to supply another column to centroidsDf:
        # centroidsDf <- do.add.cols(centroidsDf, 'centroidCol', all.annot.MC, 'Values')
        
        p <- p + geom_point(data = centroidsDf, aes(x = centroidX, 
                                                    y = centroidY,
                                                    col = as.factor(centroidCol)), size = 4.5 , alpha = 1)+
          # scale_colour_manual( values = polychromepalette ) +
          geom_point(data = centroidsDf, aes(x = centroidX, 
                                             y = centroidY), shape=3, col = "black", size = 1 , alpha = 0.5)
        
        
        }
        
        
      } else {
        
        
        p <- p + geom_point(data = centroidsDf, aes(x = centroidX, 
                                                    y = centroidY), col = "black", size = 2 , alpha = 0.3)
        
      }
      
      
      
      
      
      
      p <- p +coord_cartesian(clip = "off") #this ensures that the label is not clipped off
      
      #p <- p + guides(alpha = "none")
      
      
    }
  }
  if (blank.axis == TRUE) {
    p <- p + theme(axis.line = element_blank(), axis.text.x = element_blank(), 
                   axis.text.y = element_blank(), axis.ticks = element_blank(), 
                   axis.title.x = element_blank(), axis.title.y = element_blank(), 
                   panel.grid.major = element_blank(), panel.background = element_blank(), 
                   panel.border = element_blank(), panel.grid.minor = element_blank(), 
                   plot.background = element_blank(), )
  }
  if (save.to.disk == TRUE) {
    if (!is.null(col.axis)) {
      if (col.type == "continuous") {
        lb <- "Colour"
      }
      if (col.type == "factor") {
        lb <- "Factor"
      }
    }
    if (is.null(col.axis)) {
      lb <- "Density plot"
    }
    if (is.null(filename)) {
      filename <- paste0(lb, " plot - ", title, " - plotted on ", 
                         x.axis, " by ", y.axis, ".png")
    }
    ggsave(filename = filename, plot = p, path = path, width = plot.width, 
           height = plot.height, limitsize = FALSE)
  }
  else {
    # watch out, there is a cryptic error message if the Plots tab in RStudio is too small to send the plot to the viewer
    # so protect from doing that:
    
    if (min(dev.size("in")) > 2.2   ) {  # Minimum width & height of 4 inches
      print(p)  
    } else {
      message("Plot window is too small. Resize the RStudio plot window to see the plot.")
    }
  }
  return(p)
} # end function def make.colour.plot.adapted
















# Adapt some functions. These are based on v1.0.0 functions, but the adaptation is coming from when I tweaked them in 0.5.4. 
# the setup of make.spatial.plot is still the same:
#this is just a copy-paste of Advaced spatial 1 v2, but now we call it .adapted, since I changed the annotated cells and use global color palette...

make.spatial.plot.adapted <- function(dat, # spatial data object
                                      image.roi, # name of ROI
                                      
                                      image.channel, # name of channel # <- this is the marker signal to plot 
                                      
                                      image.min.threshold = 0.00, # this is the marker signal low percentile. this is by default set to black
                                      image.med.threshold = 0.50, # this is the marker signal mid percentile. Kinda makes sense to put that on the median
                                      markersignal.mid.col = "#2B739B", # this color is picked on the direct line of markersignal.high.col towards black
                                      image.max.threshold = 1.00, # this is the marker signal high percentile which is colored by:
                                      markersignal.high.col = "#56B4E9", # colorblind friendly skyblue alternative would be dark blue 0072B2 or "white"
                                      
                                      ## Options for adding cell outlines
                                      mask.outlines = NULL, # character -- the outlines in dat object
                                      
                                      ## Options for adding cellular data
                                      cell.dat = NULL, # can be character (if it's data within dat) or a data.table
                                      cell.col = NULL, # column for colouration
                                      
                                      ## Other settings (with defaults)
                                      image.y.flip = TRUE,
                                      image.mask.size = 0.1,
                                      image.mask.colour = "#CC79A7", # colorblind friendly reddish purple
                                      
                                      
                                      
                                      image.blank = FALSE,
                                      
                                      cell.x = "x",
                                      cell.y = "y",
                                      cell.col.type = "numeric",
                                      amountofdiscretcols = 50, # this is for no-numeric colors, defaults to 50 further down the code
                                      cell.colours = "spectral",
                                      polychromepalette = NULL, # this is the go-to way to color, you need to provide a table with colors via polychrome
                                      cell.col.min.threshold = 0.01,
                                      cell.col.max.threshold = 0.995,
                                      
                                      #legend location:
                                      legend.loc = "right",

                                      # title = paste0(image.roi), # depreciated as title. not quite sure why they used title as way to save the file name, thats bad. instead:
                                      titlestring = NULL, subtitlestring = NULL, # I added more control over how to call the ROIs...
                                      dot.size = 1,
                                      dot.alpha = 1,
                                      align.xy.by = cell.dat, # choose a data frame to set absolute limits for X/Y/colour
                                      align.col.by = cell.dat,
                                      save.to.disk = TRUE,
                                      file.extension = ".png",
                                      path = getwd(),
                                      plot.width = 30,
                                      plot.height = 20,
                                      blank.axis = FALSE)
{
  
  ### TESTING
  # library(raster)
  # library(data.table)
  # library(tiff)
  # library(ggplot2)
  #
  # dat = dat
  #
  # dat$meta.data
  #
  # roi = "20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac"
  # roi.marker = "CD20_Dy161"
  
  # cell.dat <- dat$cell.dat.means.filtered
  # cell.dat <- cell.dat[cell.dat[["ImageName"]] == "20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac_ilastik_s2_Probabilities_mask.tiff",]
  # cell.dat = cell.dat
  # cell.x = "X"
  # cell.y = "Y"
  # cell.colour = 'CD20'
  # 
  # add.outlines = TRUE
  # flip.y.axis = TRUE
  
  ### Check that necessary packages are installed
  if(!is.element('Spectre', installed.packages()[,1])) stop('Spectre is required but not installed')
  if(!is.element('ggplot2', installed.packages()[,1])) stop('ggplot2 is required but not installed')
  if(!is.element('scales', installed.packages()[,1])) stop('scales is required but not installed')
  if(!is.element('colorRamps', installed.packages()[,1])) stop('colorRamps is required but not installed')
  if(!is.element('ggthemes', installed.packages()[,1])) stop('ggthemes is required but not installed')
  if(!is.element('RColorBrewer', installed.packages()[,1])) stop('RColorBrewer is required but not installed')
  if(!is.element('raster', installed.packages()[,1])) stop('raster is required but not installed')
  if(!is.element('rgeos', installed.packages()[,1])) stop('rgeos is required but not installed')
  
  ### Require packages
  require(ggplot2)
  require(scales)
  require(colorRamps)
  require(ggthemes)
  require(RColorBrewer)
  require(raster)
  require(rgeos)
  
  ### Compatability conversions
  
  roi <- image.roi
  roi.marker <- image.channel
  
  #cell.dat
  # if('data.table' %in% class(cell.dat)){
  #   cell.dat <- cell.dat[cell.dat[['ROI']] == image.roi,]
  # }
  
  cell.colour <- cell.col
  
  add.outlines <- image.outlines <- mask.outlines
  flip.y.axis <- image.y.flip
  
  cell.colour.type <- cell.col.type
  
  raster.mask.size <- image.mask.size
  raster.mask.colour <- image.mask.colour
  raster.min.threshold <- image.min.threshold
  raster.max.threshold <- image.max.threshold
  
  col.min.threshold <- cell.col.min.threshold
  col.max.threshold <- cell.col.max.threshold
  
  colours <- cell.colours
  
  ### Colour setup
  
  # Jet
  if(colours == "jet"){
    colour.scheme <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  }
  
  # Spectral
  if(colours == "spectral"){
    spectral.list <- colorRampPalette(brewer.pal(11,"Spectral"))(amountofdiscretcols) # these were all set to 50 cols, which is overkill..
    spectral.list <- rev(spectral.list)
    colour.scheme <- colorRampPalette(c(spectral.list))
  }
  
  # Viridis
  if(colours == "viridis"){
    colour.scheme <- colorRampPalette(c(viridis_pal(option = "viridis")(amountofdiscretcols)))
  }
  
  # Inferno
  if(colours == "inferno"){
    colour.scheme <- colorRampPalette(c(viridis_pal(option = "inferno")(amountofdiscretcols)))
  }
  
  #Magma
  if(colours == "magma"){
    colour.scheme <- colorRampPalette(c(viridis_pal(option = "magma")(amountofdiscretcols)))
  }
  
  ### cell.dat setup
  
  if(!is.null(cell.dat)){
    
    if(is.character(cell.dat) == TRUE){
      temp <- dat[[roi]]@DATA[[cell.dat]]
      cell.dat <- temp
    }
    
    
    if(cell.colour.type == "numeric"){
      
      # Dot point colouration
      if(is.null(align.col.by) == TRUE){
        ColrMin <- quantile(cell.dat[[cell.colour]], probs = c(col.min.threshold))
        ColrMax <- quantile(cell.dat[[cell.colour]], probs = c(col.max.threshold))
      }
      
      if(is.null(align.col.by) == FALSE){
        ColrMin <- quantile(align.col.by[[cell.colour]], probs = c(col.min.threshold))
        ColrMax <- quantile(align.col.by[[cell.colour]], probs = c(col.max.threshold))
      }
      
    }
    
  }
  
  
  ### Preparat the raster data
  
  ## Image prep
  
  raster.image <- dat[[roi]]@RASTERS[[roi.marker]]
  
  tiff.p <- rasterToPoints(raster.image)
  tiff.df <- data.frame(tiff.p)
  raster.label <- names(tiff.df)[3]
  colnames(tiff.df) <- c("x_axis", "y_axis", raster.label)
  
  ## Create cell outlines
  if(!is.null(mask.outlines)){
    outline <- dat[[roi]]@MASKS[[mask.outlines]]$outlines
    centroids <- dat[[roi]]@MASKS[[mask.outlines]]$centroids
    
    centroid.xmin <- centroids@bbox[1]
    centroid.xmax <- centroids@bbox[3]
    
    centroid.ymin <- centroids@bbox[2]
    centroid.ymax <- centroids@bbox[4]
  }
  
  ## Flip y-axis values
  
  # if(flip.y.axis == TRUE){
  #   dat <- invert.y.axis(dat, y.axis)
  # }
  
  ## Normalise XY for cell centroids
  plot.normalize <- function(dat, min, max){
    return(((dat- min(dat)) / (max(dat)-min(dat))) * (max - min) + min)
  }
  
  # if(!is.null(cell.dat)){
  #   # X AXIS
  #   cell.dat[[cell.x]] <- plot.normalize(cell.dat[[cell.x]], min = centroid.xmin, max = centroid.xmax)
  #
  #   # Y AXIS
  #   cell.dat[[cell.y]] <- plot.normalize(cell.dat[[cell.y]], min = centroid.ymin, max = centroid.ymax)
  #
  # }
  
  ## Raster colour limits
  
  RastMin <- quantile(tiff.df[[3]], probs = c(raster.min.threshold))
  RastMax <- quantile(tiff.df[[3]], probs = c(raster.max.threshold))
  
  
  # ok, Im not fond of renaming an input variable and renaming again, thats some bad scripting here. so we do it our way:
  
  marker.percentile <- quantile(tiff.df[[3]], probs = c(raster.min.threshold,image.med.threshold,image.max.threshold))
  
  marker.lower.threshold <- marker.percentile[1]
  marker.median <- marker.percentile[2]
  marker.upper.threshold <- marker.percentile[3]
  
  ###############################################
  ### Add a check to see if centroids line up ###
  ###############################################
  
  
  
  ### Generate and show coloured plot
  
  # this plots the marker tiff channel signal:
  
  if(image.blank == FALSE){
    p <- ggplot(data=tiff.df, aes(x=tiff.df[[1]], y=tiff.df[[2]])) +
      
      ## Plot the raster (IMC image)
      geom_raster(aes(fill=tiff.df[[3]])) +
      #depreciated:
      #scale_fill_gradient(raster.label,
      #                    low = "black", #"black"
      #                    high = markersignal.high.col, # defaults to a colorblind friendly blue with grey 56B4E9 
      #                    limits=c(RastMin,RastMax),
      #                    oob=squish)
      scale_fill_gradient2(raster.label,
                           low="black", # "#0047a3", #8C30A4FF
                           mid= markersignal.mid.col, #"#339c3c",  #00A087FF
                           high=markersignal.high.col,# "#f77225", #F39B7FFF
                           midpoint = marker.median, #space='Lab' other values are depreciated
                           limits=c(marker.lower.threshold,marker.upper.threshold), 
                           #na.value = "#C00DC6FF" # as you see, we squish but plot the outside-percentiles pink
                           oob=squish
      )
    
    
  }
  
  if(image.blank == TRUE){
    p <- ggplot(data=tiff.df, aes(x=tiff.df[[1]], y=tiff.df[[2]])) +
      
      ## Plot the raster (IMC image)
      geom_raster(aes(fill=tiff.df[[3]])) +
      scale_fill_gradient(raster.label,
                          low = "black", #"black"
                          high = markersignal.high.col, 
                          limits=c(RastMin,RastMax),
                          oob=squish)
  }
  
  
  
  ### Plot the cell mask boundaries
  
  if(!is.null(image.outlines)){
    p <- p + geom_path(aes(x = long, y = lat, group = group),
                       data = outline,
                       size = raster.mask.size,
                       col = raster.mask.colour)
  }
  
  ## Plot the cellular data
  
  if(!is.null(cell.dat)){
    if(cell.colour.type == "numeric"){
      p <- p + geom_point(data=cell.dat,
                          aes(x=cell.dat[[cell.x]], y=cell.dat[[cell.y]], color = cell.dat[[cell.colour]]),  #as.numeric(as.character(col))
                          size = dot.size, #dot.size
                          alpha = dot.alpha # shape = 1
      ) +
        
        scale_color_gradientn(colours = colour.scheme(50),
                              limits = (c(ColrMin,ColrMax)),
                              oob=squish,
                              name = cell.colour)
    }
    
    # this here plots the annotated clusters:
    if(cell.colour.type != "numeric"){
      #this covers the cell.color.type, now for the special case of polychrome
      
      #this specral thing is annoying since you cannot see the clusters and the colors get re-assigned everytime we print another subset.
      #so lets predefine a palette and then use it in here:
      if (colours == "polychrome") {
      
        p <- p + geom_point(data=cell.dat,
                            aes(x=cell.dat[[cell.x]], y=cell.dat[[cell.y]], color = as.factor(cell.dat[[cell.colour]])),
                            size = dot.size, 
                            alpha = dot.alpha)+
          # push the polychrome palette in
          # now there is a bug in ggplot, which would plot the whole palette here, like 37 clusters or so. 
          #you need to force the limits so that only the current subset ends in the legend
          scale_colour_manual( values = polychromepalette,
                               limits = force 
                               )+ 
          guides(colour = guide_legend(override.aes = list(size=4, alpha=1),  #we also gonna override the dot size how its depicted in the legend
                                       byrow=TRUE) ) # and I want the elements listed row-by-row
      
      } else {
      #old code:
      p <- p + geom_point(data=cell.dat,
                          aes(x=cell.dat[[cell.x]], y=cell.dat[[cell.y]], color = as.factor(cell.dat[[cell.colour]])),  #as.numeric(as.character(col))
                          size = dot.size, #dot.size 
                          alpha = dot.alpha # shape = 1
      ) +
        scale_colour_discrete(name = cell.colour)
      
      }#end default factor plot with col.axis set
      
    }#end cell.color.type not numeric
    
  }
  
  ## Setup some themes
  p <- p + theme_bw() +
    coord_equal() +
    xlab(cell.x)+
    ylab(cell.y)+
    #ggtitle(title) # depreciated, we use titlestring and subtitlestring:
    labs(title= titlestring , 
         subtitle=paste0(subtitlestring, "        (Marker: ", image.channel, ")" ), # subtitlestring, 
         #caption="Created by M.Barone"#, 
         y="ROI y", 
         x="ROI x"
    )
  
  
  # define where to go with the legend(s)
  if (legend.loc %in% c("top", "bottom")) {
    p <- p + theme(legend.direction = "horizontal", 
                   legend.position = legend.loc, 
                   legend.text = element_text(size = 12), 
                   legend.title = element_blank() # title is out, we plot the marker in the subtitlelyric
                   )
    }
  if (legend.loc %in% c("left", "right")) {
    p <- p + theme(legend.direction = "vertical", 
                   legend.position = legend.loc, 
                   legend.text = element_text(size = 12), 
                   legend.title = element_blank() # title is out, we plot the marker in the subtitlelyric
                   )
  }
  
  
  ## More themes
  p <- p + theme(panel.background = element_rect(fill = "grey5", colour = "grey15", size = 0.5), # change 'colour' to black for informative axis
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
  
  if(flip.y.axis == TRUE){
    
    p <- p + scale_y_reverse()
    
  }
  
  ### Save ggplot to disk if desired
  if(save.to.disk == TRUE){
    ggsave(filename = paste0(image.roi, " ",image.channel , file.extension), #paste0(image.roi, "_ROI_", roi.marker, "_marker_", cell.colour,".png"), 
           plot = p,
           path = path,
           width = plot.width,
           height = plot.height,
           limitsize = FALSE)
  }
  
  print(p)
  
}#end make.spatial.plot adaptation






#CD45_Er168, ~GLP.1R_Yb174

# this function requires three more variables that get produced as the engine runs trhough the markercomparision in a m-loop:
# it will strip off the _*isotope* and create a marker string of the AB:
#xMarker <- sub("_[^_]+$", "", markercomparision[[m]][1] )
#yMarker <- sub("_[^_]+$", "", markercomparision[[m]][2] )
#colorMarker <- sub("_[^_]+$", "", markercomparision[[m]][3] )  
# whith this strings, we match back in the meta or the multisegmented lines to plot additional stuff


plot.scatter.gate <- function(indata, 
                             inx, 
                             iny, 
                             incol=NA, 
                             dot.size = 1.0,
                             dot.alpha = 0.2,
                             gate=NULL,
                             titlestring=paste0(i) # this used to be a function to run through all TMA with i. hence we use that as default title here.
){
  
  #p=NULL
  p <- ggplot(indata, 
              aes_string(x = inx, 
                         y = iny, 
                         colour = indata[[incol]]),  #indata[,Insulin_blabal] worked
              environment = environment()
  ) #  ggplot needs to search within the function env and this needs explicit pass
  
  
  
  
  
  
  
  
  if (!is.na(incol)) {
    p <- p + geom_point(alpha=dot.alpha,
                        size=dot.size) +
      scale_color_gradient2(low="#0047a3", #8C30A4FF
                            mid="#339c3c",  #00A087FF
                            high="#f77225", #F39B7FFF
                            midpoint = median, #space='Lab' other values are depreciated
                            limits=c(lower.threshold,upper.threshold), 
                            #limits=c(quartiles[1],quartiles[3]), #oob=squish, # logical: squish outofbound vals into limits
                            na.value = "#C00DC6FF", name= sub("_[^_]+$", "", incol )  # the sub just removes everything after the underscore
      )
    
  }else{
    
    # now we check if xMarker or yMarker are part of any outside.*AB*.highpass or lowpass col and draw in and out with according colors:
    # the problem here is that if you have markercombinations that both carry flags, you dont know if to plot boolean AND or simply OR. #
    # We go with default AND here:
    
    if ( inx %in% gate_affected.markers & iny %in% gate_affected.markers ) {
      
      
      #we need to check if they are part of the multisegmented lines. in that case we need to color them accordingly and plot the multisegmented line.
      
      
      # this is the ugly condition where we need to plot both gates and dont know if they are related at all and if with AND or OR?
      # so we dont separate colors:
     p <- p +  geom_point(color="#030280",
                         alpha=dot.alpha,
                         size=dot.size) 
     
     
     
     
     # end inx and iny are part of gate_affected.markers
    }else if ( inx %in% gate_affected.markers ) {
      
      # If the x-axis marker carries a gate, use the flags to draw different color
      # we do not care if its a highpass or lowpass but just brute-force our way through both options:
      
      p <- p +  geom_point(data= subset(indata, indata[[paste0("outside.",xMarker,".highpass")]]==1  ),
                           color="#EB8563", # that used to be a violet 8D60AC
                           alpha=dot.alpha,
                           size=dot.size) +
        
                 geom_point(data=subset(indata, indata[[paste0("outside.",xMarker,".highpass")]]==0   ),
                           color="#308441", # that used to be dark blue 030280
                           alpha=dot.alpha,
                           size=dot.size) +
        
        geom_point(data= subset(indata, indata[[paste0("outside.",xMarker,".lowpass")]]==1  ),
                   color="#EB8563",
                   alpha=dot.alpha,
                   size=dot.size) +
        
        geom_point(data=subset(indata, indata[[paste0("outside.",xMarker,".lowpass")]]==0   ),
                   color="#308441",
                   alpha=dot.alpha,
                   size=dot.size) 
        
   
      
    }else if ( iny %in% gate_affected.markers ) {
      
      # If the y-axis marker carries a gate, use the flags to draw different color
      # we do not care if its a highpass or lowpass but just brute-force our way through both options:
      
      p <- p +  geom_point(data= subset(indata, indata[[paste0("outside.",yMarker,".highpass")]]==1  ),
                           color="#EB8563",
                           alpha=dot.alpha,
                           size=dot.size) +
        
        geom_point(data=subset(indata, indata[[paste0("outside.",yMarker,".highpass")]]==0   ),
                   color="#308441",
                   alpha=dot.alpha,
                   size=dot.size) +
        
        geom_point(data= subset(indata, indata[[paste0("outside.",yMarker,".lowpass")]]==1  ),
                   color="#EB8563",
                   alpha=dot.alpha,
                   size=dot.size) +
        
        geom_point(data=subset(indata, indata[[paste0("outside.",yMarker,".lowpass")]]==0   ),
                   color="#308441",
                   alpha=dot.alpha,
                   size=dot.size)
        
        
        

    } else {
      
      # whats left here are the markers that are not part of the flag columns, including the multi-segmented gates. 
      
      p <- p +  geom_point(color="#030280",
                           alpha=dot.alpha,
                           size=dot.size) 
      
      
    }
    
    
    

  }# end else bracket to take care of the geom_point plotting for when there is no incol marker given to draw a gradient 
  
  
  
  
  # after the points are drawn, lets put the gate lines on top:
  
  # for the multisegemented lines this is still messy and needs a generalized approach.
  # if the combination is part of the multisegmented lines, we take care of that first, and only then check back in the metadata:
  
 # if( exists("CD45segmentlinePOINTS") ){}
  
  if ( xMarker %in% "CD45" & yMarker %in% "GLP.1R" ) {

    # so we scale the line on the fly by just using the scaling factors of the fine-tune dataset. its therefore important that the Batch is stored in "i"...:
    scaled.CD45segmentlinePOINTS <- CD45segmentlinePOINTS %>%
      mutate(CD45_Er168 = CD45_Er168 * finetune.CD45segmentline[  match (i, finetune.CD45segmentline$Batch.ID), names(CD45segmentlinePOINTS)[1]] ) %>%
      mutate(GLP.1R_Yb174 = GLP.1R_Yb174 * finetune.CD45segmentline[  match (i, finetune.CD45segmentline$Batch.ID), names(CD45segmentlinePOINTS)[2]] )
    
    
    #scaled.CD45segmentlinePOINTS <- CD45segmentlinePOINTS 
    
    # p <- p + geom_path(data = scaled.CD45segmentlinePOINTS,  size = 1.6, colour = "grey")+
    #    geom_path(data = scaled.CD45segmentlinePOINTS, color = "grey",size=1.1, alpha=0.55)
    
    
    p <- p + geom_path(data = scaled.CD45segmentlinePOINTS,  size = 1.1, colour = "#A100AB", alpha=0.4)
    
  }else{
    
    # we first check if inx and iny are part of the gate_affected.marker, and if so, plot the according line: xMarker %in%  
    # of course we need to set gates TMA wise. a ROI-wise gating is not allowed since the "unique" command in the geom_hline will crash the script here:   
    # as of v7, this is changed. We pull the proper gate column via threshold.gates in one go and do not have to care if its an Hgate or Lgate:
    
    if ( inx %in% gate_affected.markers & inx %!in% names(CD45segmentlinePOINTS)  ) {
      # you need to protect here from plotting an inexistent threshold gate if inx was part of the multisegmented gating routine:
      
      p <- p +  geom_vline(xintercept = as.numeric(unique( indata[[ threshold.gates[str_detect(threshold.gates, xMarker)]    ]]    )) ,  size = 1.1, colour = "#A100AB", alpha=0.4)
    }
    
    if ( iny %in% gate_affected.markers ) {
      p <- p +  geom_hline(yintercept = as.numeric(unique( indata[[ threshold.gates[str_detect(threshold.gates, yMarker)]    ]]    )) ,  size = 1.1, colour = "#A100AB", alpha=0.4)
    }
    
  }#end else bracket to paint line thresholds instead of the multisegmented line
  
  
  
  
  
  
  
  
  
  
  

  
  
  
  
  
  p <- p +  
    scale_x_log10(breaks = breaks, 
                  minor_breaks = minor_breaks,
                  labels = trans_format("log10", math_format(10^.x))) +
    
    scale_y_log10(breaks = breaks, 
                  minor_breaks = minor_breaks,
                  labels = trans_format("log10", math_format(10^.x))) +
    
    
    #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
    #              minor_breaks = trans_breaks("log10", math_format(10^.x)),
    #              labels = trans_format("log10", math_format(10^.x))) +
    
    
   # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #                labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks()   +
    
    theme_minimal()  
  # guides(colour = guide_legend(nrow = 8)) + # watch out, this induces binning into the continous color scale
  
  
  if (!is.na(incol)) {
    # here we prepare the title string. if incol is given, mention it in the title:
    
    # I think by now gate is depreciated since the function searched automatically if either marker is part of a gate.
    #  if (!is.null(gate)) {
    #    p <- p +  labs(title=paste0(i, " (gate set at ",gate,")")  , 
    #      subtitle=paste0(nrow(indata), " cells, ",sub("_[^_]+$", "", incol )," colored from ",(1-percentile.colorlimits)*100,"-", percentile.colorlimits*100,"th percentile") 
    #caption="Created by M.Barone", 
    # y="Maximal Velocity [AU/sec]", 
    #  x="Condition"
    #  )
    #  }
    
    p <- p +  labs(title=titlestring , 
                   subtitle=paste0(nrow(indata), " pancreatic cells, ",sub("_[^_]+$", "", incol )," colored ",(1-percentile.colorlimits)*100,"-", percentile.colorlimits*100,"th percentile") 
                   #caption="Created by M.Barone", 
                   # y="Maximal Velocity [AU/sec]", 
                   #  x="Condition"
    )
    
    
    
  }else{
    
    p <- p +  labs(title=titlestring  , 
                   subtitle=paste0(nrow(indata), " pancreatic cells") 
                   #caption="Created by M.Barone", 
                   # y="Maximal Velocity [AU/sec]", 
                   #  x="Condition"
    )
    
  }# end else bracket for when incol was NA and there is no need to plot this in the subtitle string 
  
  
  p <- p +  coord_fixed() + # this forces equidistant axis scaling. turn off aspect.ratio, of course...
    
    ylab( yMarker ) + 
    xlab( xMarker ) +  
    theme(#axis.title.x = element_blank(), 
      legend.title=element_text(size=13), # turn off with element_blank(),
      legend.text=element_text(size=11),
      axis.text=element_text(size=13),#,axis.text.x=element_text(size=10),
      axis.title.y=element_text(size=12), axis.title.x=element_text(size=12),
      # legend.position="top", 
      legend.key.width=unit(1,"line"),
      #aspect.ratio=1 #force ggplot on squared plotting # turn off if you force fixed coords of course
    ) 
  
  
  print(p)
  
  
} # end function def  plot.scatter.gate




plot.scatter.gate.v2 <- function(indata, 
                                 inx, 
                                 iny, 
                                 incol=NA, 
                                 dot.size = 1.0,
                                 dot.alpha = 0.2,
                                 gate.flag = NA,
                                 gate.line.alpha=0.4,
                                 gate.line.size=1.1, 
                                 gate.line.color = "#A100AB",
                                 passed.gate.color = "#308441",
                                 rejected.gate.color = "#EB8563",
                                 gate=NULL,
                                 titlestring=paste0(i) # this used to be a function to run through all TMA with i. hence we use that as default title here.
){
  
  #p=NULL
  p <- ggplot(indata, 
              aes_string(x = inx, 
                         y = iny, 
                         colour = indata[[incol]]), 
              environment = environment()
  ) #  ggplot needs to search within the function env and this needs explicit pass
  
  
  if (!is.na(incol)) {
    p <- p + geom_point(alpha=dot.alpha,
                        size=dot.size) +
      scale_color_gradient2(low="#0047a3", #8C30A4FF
                            mid="#339c3c",  #00A087FF
                            high="#f77225", #F39B7FFF
                            midpoint = median, #space='Lab' other values are depreciated
                            limits=c(lower.threshold,upper.threshold), 
                            #limits=c(quartiles[1],quartiles[3]), #oob=squish, # logical: squish outofbound vals into limits
                            na.value = "#C00DC6FF", name= sub("_[^_]+$", "", incol )  # the sub just removes everything after the underscore
      )
    
  }
  
  
  # we set already some variable that help us now to order the plotting routine:
  
  #
  # ...........  the current gate flags are not part of the scatter plot: ...........
  if(gate.printing.needed == 0){
    #  the markers are not part of the flag columns, no action needed
    p <- p +  geom_point(color="#030280",
                         alpha=dot.alpha,
                         size=dot.size) 
    
  }
  
  
  # ...........  the current gate flags are part of the scatter plot and are thresholds from the metadata file: ...........
  if(gate.printing.needed == 1 & is.metadata.threshold == 1 ){
    
    # next, we differenctiate between high and lowpass. This was not done in the first plotting routine and lead to some weird plotting artifacts:
    
    
    
    
    p <- p +  geom_point(data= subset(indata, indata[[ gate.flag ]]==1  ),
                         color=rejected.gate.color, 
                         alpha=dot.alpha,
                         size=dot.size) +
      
      geom_point(data=subset(indata, indata[[ gate.flag ]]==0   ),
                 color= passed.gate.color,
                 alpha=dot.alpha,
                 size=dot.size) 
    
    
    # the only thing left is to check if our threshold gate needs to be drawn vertically or horizonatlly:
    if(  grepl(temp.gated.marker, inx) ){ 
      p <- p +  geom_vline(xintercept = as.numeric(unique( indata[[ sub("^.*?\\.", "", gate.flag)     ]]    )) ,  
                           size = gate.line.size, 
                           colour = gate.line.color, 
                           alpha=gate.line.alpha)
    }else{
      p <- p +  geom_hline(yintercept = as.numeric(unique( indata[[ sub("^.*?\\.", "", gate.flag)    ]]    )) ,  
                           size = gate.line.size, 
                           colour = gate.line.color, 
                           alpha=gate.line.alpha) 
    }
    
    
    
    
    
    
    
  }
  
  
  # first, lets take care of the intersting case where inx is part of gate(s)
  
  # now we check if xMarker or yMarker are part of any outside.*AB*.highpass or lowpass col and draw in and out with according colors:
  # the problem here is that if you have markercombinations that both carry flags, you dont know if to plot boolean AND or simply OR. #
  # We go with default AND here:
  
  
  
  
  #we need to check if they are part of the multisegmented lines. in that case we need to color them accordingly and plot the multisegmented line.
  if(exists("polygon.gates") ){
    
    for(g in names(polygon.gates)){
      
      if ( inx %in%  names(polygon.gates[[g]])[1] & iny %in% names(polygon.gates[[g]])[2] ) {
        
        
        p <- p +  geom_point(data= subset(indata, indata[[  paste0("outside.", g)  ]]==1  ),
                             color="#EB8563", 
                             alpha=dot.alpha,
                             size=dot.size) +
          
          geom_point(data=subset(indata, indata[[ paste0("outside.", g)  ]]==0   ),
                     color="#308441", 
                     alpha=dot.alpha,
                     size=dot.size)
        
        # for plotting the line, we could have fine-tuned it TMA-wise:
        
        if(exists("finetune.segmentline")){
          
          # so we scale the line on the fly by just using the scaling factors of the fine-tune dataset. its therefore important that the Batch is stored in "i"...:
          #scaled.segmentlinePOINTS <- polygon.gates[[g]] %>%
          #mutate(
          # polygon.gates[[g]][1] = polygon.gates[[g]][1] * finetune.segmentline[  match (i, finetune.segmentline$Batch.ID), names(polygon.gates[[g]])[1] ] ,
          # polygon.gates[[g]][2] = polygon.gates[[g]][2] * finetune.segmentline[  match (i, finetune.segmentline$Batch.ID), names(polygon.gates[[g]])[2] ] 
          # )
          
          
          #  p <- p + geom_path(data = scaled.segmentlinePOINTS,  size = 1.1, colour = "#A100AB", alpha=0.4)
        }else{
          p <- p + geom_path(data = polygon.gates[[g]],  size = 1.1, colour = "#A100AB", alpha=0.4) 
        }
        
      }# end when you found a matching polygon that contains the marker in correct orientation
    }#end run along the polygon.gate list
    
    
    # end both inx is in gate_affected.markers and polygon.gates exists 
  }
  
  
  
  p <- p +  
    scale_x_log10(breaks = breaks, 
                  minor_breaks = minor_breaks,
                  labels = trans_format("log10", math_format(10^.x))) +
    
    scale_y_log10(breaks = breaks, 
                  minor_breaks = minor_breaks,
                  labels = trans_format("log10", math_format(10^.x))) +
    
    
    #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
    #              minor_breaks = trans_breaks("log10", math_format(10^.x)),
    #              labels = trans_format("log10", math_format(10^.x))) +
    
    
    # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
    #                labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks()   +
    
    theme_minimal()  
  # guides(colour = guide_legend(nrow = 8)) + # watch out, this induces binning into the continous color scale
  
  
  if (!is.na(incol)) {
    # here we prepare the title string. if incol is given, mention it in the title:
    
    p <- p +  labs(title=titlestring , 
                   subtitle=paste0("Plotted ", nrow(indata), " cells, ",sub("_[^_]+$", "", incol )," colored ",(1-percentile.colorlimits)*100,"-", percentile.colorlimits*100,"th percentile") 
                   #caption="Created by M.Barone", 
                   # y="Maximal Velocity [AU/sec]", 
                   #  x="Condition"
    )
    
    
    
  }else{
    
    p <- p +  labs(title=titlestring  , 
                   subtitle=paste0("Plotted ", nrow(indata), " cells") 
                   #caption="Created by M.Barone", 
                   # y="Maximal Velocity [AU/sec]", 
                   #  x="Condition"
    )
    
  }# end else bracket for when incol was NA and there is no need to plot this in the subtitle string 
  
  
  p <- p +  coord_fixed() + # this forces equidistant axis scaling. turn off aspect.ratio, of course...
    
    ylab( yMarker ) + 
    xlab( xMarker ) +  
    theme(#axis.title.x = element_blank(), 
      legend.title=element_text(size=13), # turn off with element_blank(),
      legend.text=element_text(size=11),
      axis.text=element_text(size=13),#,axis.text.x=element_text(size=10),
      axis.title.y=element_text(size=12), axis.title.x=element_text(size=12),
      # legend.position="top", 
      legend.key.width=unit(1,"line"),
      #aspect.ratio=1 #force ggplot on squared plotting # turn off if you force fixed coords of course
    ) 
  
  
  print(p)
  
  
} # end function def  plot.scatter.gate.v2




plot.scatter.density<- function(indata, 
                                inx, 
                                iny, 
                                titlestring,
                                dot.size = NULL
                                ){
  
 
 # indata <- subset(indata, indata[[inx]]>0 & indata[[iny]]>0 )  
  
  
  p <- ggplot(indata, 
              aes_string(x = inx, 
                         y = iny  
              ), 
              environment = environment() )  #  ggplot needs to search within the function env and this needs explicit pass
  
  


  
  # we first check if inx and iny are part of the gate_affected.marker, and if so, plot the according line: xMarker %in%  
  # the intercept is pushed into line by just pasting _Hgate back onto the marker antibody, choose this column and make unique.
  # of course we need to set gates TMA wise. a ROI-wise gating is not allowed since unique will crash the script here:
  
  # if the combination is part of the multisegmented lines, we take care of that first, and only then check back in the metadata:
  
  if ( xMarker %in% "CD45" & yMarker %in% "GLP.1R" ) {
    
    scaled.CD45segmentlinePOINTS <- CD45segmentlinePOINTS %>%
      mutate(CD45_Er168 = CD45_Er168 * finetune.CD45segmentline[  match (i, finetune.CD45segmentline$Batch.ID), names(CD45segmentlinePOINTS)[1]] ) %>%
      mutate(GLP.1R_Yb174 = GLP.1R_Yb174 * finetune.CD45segmentline[  match (i, finetune.CD45segmentline$Batch.ID), names(CD45segmentlinePOINTS)[2]] )
    
   # scaled.CD45segmentlinePOINTS <- CD45segmentlinePOINTS 
    
    p <- p + geom_path(data = scaled.CD45segmentlinePOINTS,  size = 1.1, colour = "#A100AB", alpha=0.4)
    
  }else{
   
    # we first check if inx and iny are part of the gate_affected.marker, and if so, plot the according line: xMarker %in%  
    # of course we need to set gates TMA wise. a ROI-wise gating is not allowed since the "unique" command in the geom_hline will crash the script here:   
    # as of v7, this is changed. We pull the proper gate column via threshold.gates in one go and do not have to care if its an Hgate or Lgate:
    
    if ( inx %in% gate_affected.markers & inx %!in% names(CD45segmentlinePOINTS)  ) {
      # you need to protect here from plotting an inexistent threshold gate if inx was part of the multisegmented gating routine:
      
      p <- p +  geom_vline(xintercept = as.numeric(unique( indata[[ threshold.gates[str_detect(threshold.gates, xMarker)]    ]]    )) ,  size = 1.1, colour = "#A100AB", alpha=0.4)
    }
    
    if ( iny %in% gate_affected.markers ) {
      p <- p +  geom_hline(yintercept = as.numeric(unique( indata[[ threshold.gates[str_detect(threshold.gates, yMarker)]    ]]    )) ,  size = 1.1, colour = "#A100AB", alpha=0.4)
    }
    
  }#end else bracket to paint line thresholds instead of the multisegmented line
  
  
  
  
  # by default, the density comes with painting dots, and not circles to speed up the plotting:
  if (is.null(dot.size)) { 
    p <- p + geom_point(alpha = 0.1,colour = "#030280", pch='.')
  }else{
    p <- p + geom_point(alpha=0.1,colour = "#030280", size=dot.size) 
  }
  


    p <- p +  stat_density_2d(aes( fill = ..level..),        #fill = ..level..), 
                    alpha = 0.2, 
                    contour = TRUE,
                    contour_var = "density",
                    n = 120, # default 100, computationally intense, do not exceed that value too much..
                    geom = "polygon", #geom: polygon or raster (equidistant lines or heatmap-ish picture)
                    colour="#99CBFE")+ 
    
    
    scico::scale_fill_scico(palette = "roma", direction=-1)+ # bilbao
    
      
      # the multisegment line is defined farther than data and will be plotted fully, which makes not much sense. So:
      scale_x_log10(breaks = breaks, 
                    minor_breaks = minor_breaks,
                    labels = trans_format("log10", math_format(10^.x))) +
      
      scale_y_log10(breaks = breaks, 
                    minor_breaks = minor_breaks,
                    labels = trans_format("log10", math_format(10^.x))) +
                        #limits = c( log10(min( indata[[inx]] ))  ,  log10(max( indata[[inx]]  )) ) 
                  
      

      
      
      
    annotation_logticks()   +
    
    theme_minimal() + 
    
    labs(title= titlestring , 
         subtitle=paste0(nrow(indata), " cells plotted")#, 
         #caption="Created by M.Barone"#, 
         #       # y="Maximal Velocity [AU/sec]", 
         #       #  x="Condition"
    )+
    
    coord_fixed() + # this forces equidistant axis scaling. turn off aspect.ratio, of course...
    
    ylab( yMarker ) + 
    xlab( xMarker ) +  
      

      
    theme(# this was on, but lets see what happens if title and subtitle is left to ggplot: plot.title = element_text(size=20, face="bold"),plot.subtitle = element_text(size=16),
      #axis.title.x = element_blank(), 
      legend.title=element_text(size=13), # turn off with element_blank(),
      legend.text=element_text(size=11),
      axis.text=element_text(size=13),#,axis.text.x=element_text(size=10),
      axis.title.y=element_text(size=12), axis.title.x=element_text(size=12),
      # legend.position="bottom", 
      legend.key.width=unit(1,"line")#,
      #aspect.ratio=1 #force ggplot on squared plotting
    )
  
    print(p)
  
} # end function def  plot.scatter.density






plot.scatter.density.v2 <- function(indata, 
                                    inx, 
                                    iny, 
                                    incol=NA, 
                                    dot.size = NULL,
                                    dot.alpha = 0.1,
                                    gate.flag = NA,
                                    gate.line.alpha=0.4,
                                    gate.line.size=1.1, 
                                    gate.line.color = "#A100AB",
                                    passed.gate.color = "#308441",
                                    rejected.gate.color = "#F6CD7A",
                                    gate=NULL,
                                    titlestring=paste0(i), # this used to be a function to run through all TMA with i. hence we use that as default title here.
                                    subtitlestring=paste0("Plotted ", nrow(indata), " cells") 
){
  
  #p=NULL
  p <- ggplot(indata, 
              aes_string(x = inx, 
                         y = iny, 
                         colour = indata[[incol]]), 
              environment = environment()
  ) #  ggplot needs to search within the function env and this needs explicit pass
  
  
  if (!is.na(incol)) {
    
    
    # by default, the density comes with painting dots, and not circles to speed up the plotting:
    if (is.null(dot.size)) { 
      p <- p + geom_point(alpha = 0.5, pch='.')+
        
   #     scale_color_gradient2(low="#0047a3", #8C30A4FF
  #                            mid="#339c3c",  #00A087FF
  #                            high="#f77225", #F39B7FFF
  #                            midpoint = median, #space='Lab' other values are depreciated
   #                           limits=c(lower.threshold,upper.threshold), 
  #                            #limits=c(quartiles[1],quartiles[3]), #oob=squish, # logical: squish outofbound vals into limits
  #                            na.value = "#C00DC6FF", name= sub("_[^_]+$", "", incol )  # the sub just removes everything after the underscore
  #      )
                              
       scico::scale_color_scico(palette = "roma", direction=-1, 
                                                       midpoint = median,
                                                       limits=c(lower.threshold,upper.threshold),
                                                       na.value = "#8C0172"                              ) # bilbao roma                      
                              
      
    }else{
      p <- p + geom_point(alpha=0.1, size=dot.size)+
     #   scale_color_gradient2(low="#0047a3", #8C30A4FF
    #                          mid="#339c3c",  #00A087FF
    #                          high="#f77225", #F39B7FFF
    #                          midpoint = median, #space='Lab' other values are depreciated
     #                         limits=c(lower.threshold,upper.threshold), 
     #                         #limits=c(quartiles[1],quartiles[3]), #oob=squish, # logical: squish outofbound vals into limits
     #                         na.value = "#C00DC6FF", name= sub("_[^_]+$", "", incol )  # the sub just removes everything after the underscore
     #   )
      scico::scale_color_scico(palette = "roma", direction=-1, 
                               midpoint = median,
                               limits=c(lower.threshold,upper.threshold),
                               na.value = "#8C0172"      ) # bilbao roma   
    }
    
    
    
    
  }# end set color gradient over third marker experssion, rather than the flags.
  
  
  # we set already some variable that help us now to order the plotting routine:
  
  #
  # ...........  the current gate flags are not part of the scatter plot: ...........
  if(gate.printing.needed == 0){
    
    # if no gates are to be printed, we just need to set the dot colors:
    
    #  the markers are not part of the flag columns, no action needed
    # by default, the density comes with painting dots, and not circles to speed up the plotting:
    if (is.null(dot.size)) { 
      p <- p + geom_point(alpha = 0.1,colour = "#030280", pch='.')
    }else{
      p <- p + geom_point(alpha=0.1,colour = "#030280", size=dot.size) 
    }
    
  } # no gates need printing
  
  
  # ...........  the current gate flags are part of the scatter plot and are thresholds from the metadata file: ...........
  if(gate.printing.needed == 1 ){
    
    # next, we differenctiate between high and lowpass. This was not done in the first plotting routine and lead to some weird plotting artifacts:
    
    # now, if we have a color marker, the dot colors are set already. we just have to take care of a missing color marker. 
    # in that case we of course dont plot blue dots, but plot the cells by their gate line status in or out:
    if (is.na(incol)){
    # by default, the density comes with painting dots, and not circles to speed up the plotting:
    if (is.null(dot.size) ) { 
      p <- p +  geom_point(data= subset(indata, indata[[ gate.flag ]]==1  ),
                           color=rejected.gate.color, 
                           alpha=dot.alpha,
                           pch='.') +
        
        geom_point(data=subset(indata, indata[[ gate.flag ]]==0   ),
                   color= passed.gate.color,
                   alpha=dot.alpha,
                   pch='.') 
      
      
    }else{
      
      p <- p +  geom_point(data= subset(indata, indata[[ gate.flag ]]==1  ),
                           color=rejected.gate.color, 
                           alpha=dot.alpha,
                           size=dot.size) +
        
        geom_point(data=subset(indata, indata[[ gate.flag ]]==0   ),
                   color= passed.gate.color,
                   alpha=dot.alpha,
                   size=dot.size) 
      
    }
  }#set colors if incol is NA
    
    
    
    
    
    # so we are within gate.printing.needed yes.
    
    # if its a metadata threshold, this flag is also set to 1
    if( is.metadata.threshold == 1 ){
      
      # the only thing left is to check if our threshold gate needs to be drawn vertically or horizonatlly:
      if(  grepl(temp.gated.marker, inx) ){ 
        p <- p +  geom_vline(xintercept = as.numeric(unique( indata[[ sub("^.*?\\.", "", gate.flag)     ]]    )) ,  
                             size = gate.line.size, 
                             colour = gate.line.color, 
                             alpha=gate.line.alpha)
      }else{
        p <- p +  geom_hline(yintercept = as.numeric(unique( indata[[ sub("^.*?\\.", "", gate.flag)    ]]    )) ,  
                             size = gate.line.size, 
                             colour = gate.line.color, 
                             alpha=gate.line.alpha) 
      }
      
      
    }else{
      # gate printing 1 but metadata threshold 0 leaves only the polygons 
      
      
      #we need to check if they are part of the multisegmented lines. in that case we need to color them accordingly and plot the multisegmented line.
      if(exists("polygon.gates") ){
        
        for(g in names(polygon.gates)){
          
          if ( inx %in%  names(polygon.gates[[g]])[1] & iny %in% names(polygon.gates[[g]])[2] ) {
            
            #  message("Plotting polygon")
            #   p <- p +  geom_point(data= subset(indata, indata[[  paste0("outside.", g)  ]]==1  ),
            #                         color="#EB8563", 
            #                         alpha=dot.alpha,
            #                         size=dot.size) +
            #      
            #      geom_point(data=subset(indata, indata[[ paste0("outside.", g)  ]]==0   ),
            #                 color="#308441", 
            #                 alpha=dot.alpha,
            #                 size=dot.size)
            
            # for plotting the line, we could have fine-tuned it TMA-wise:
            
            if( exists("finetune.segmentline") ){
              
              # watch out, you need to make sure that if you dont bring a long a i for all.TMAs, you are probably plottling the entire dataset,
              # so protect from running into a crash when the function wants to match the TMA:
              
              if(!is.na(i)){
                
                scaled.segmentlinePOINTS <- polygon.gates[[g]] %>%
                  mutate( 
                    polygon.gates[[g]][1] * finetune.segmentline[  match (i, finetune.segmentline$Batch.ID), names(polygon.gates[[g]])[1] ] ,
                    polygon.gates[[g]][2] * finetune.segmentline[  match (i, finetune.segmentline$Batch.ID), names(polygon.gates[[g]])[2] ] 
                  )
                
                
                p <- p + geom_path(data = scaled.segmentlinePOINTS,  
                                   size = gate.line.size, 
                                   colour = gate.line.color, 
                                   alpha=gate.line.alpha)
                
              }else{
                
                # so, i was set to NA, meaning we need to plot all lines subsequently by running through all TMAs:
                
                for(j in all.TMAs){
                  
                  scaled.segmentlinePOINTS <- polygon.gates[[g]] %>%
                    mutate( 
                      polygon.gates[[g]][1] * finetune.segmentline[  match (j, finetune.segmentline$Batch.ID), names(polygon.gates[[g]])[1] ] ,
                      polygon.gates[[g]][2] * finetune.segmentline[  match (j, finetune.segmentline$Batch.ID), names(polygon.gates[[g]])[2] ] 
                    )
                  
                  #    message(paste0("i set to NA, plotting all polygons: ",j))
                  #    print(scaled.segmentlinePOINTS)
                  
                  p <- p + geom_path(data = scaled.segmentlinePOINTS,  
                                     size = gate.line.size, 
                                     colour = gate.line.color, 
                                     alpha=gate.line.alpha)
                  
                }
                
                
                
              }
              
              
              # so we scale the line on the fly by just using the scaling factors of the fine-tune dataset. its therefore important that the Batch is stored in "i"...:
              scaled.segmentlinePOINTS <- polygon.gates[[g]] %>%
                mutate( 
                  polygon.gates[[g]][1] * finetune.segmentline[  match (i, finetune.segmentline$Batch.ID), names(polygon.gates[[g]])[1] ] ,
                  polygon.gates[[g]][2] * finetune.segmentline[  match (i, finetune.segmentline$Batch.ID), names(polygon.gates[[g]])[2] ] 
                )
              
              
              p <- p + geom_path(data = scaled.segmentlinePOINTS,  
                                 size = gate.line.size, 
                                 colour = gate.line.color, 
                                 alpha=gate.line.alpha)
              
              
              # end finetune.segmentline exists, so plot them accordingly 
            }else{
              #all TMAs use the same gateline, so no need to finetune.
              p <- p + geom_path(data = polygon.gates[[g]],  size = gate.line.size, 
                                 colour = gate.line.color, 
                                 alpha=gate.line.alpha) 
            }
            
            
            
            
          }# end when you found a matching polygon that contains the marker in correct orientation
        }#end run along the polygon.gate list
        
        
        # end both inx is in gate_affected.markers and polygon.gates exists 
      }# end polygon.gates actually does exist. lucky us
    }# end plotting polygons 
    
    
  } # end plot gates and color the dots with in and out of flags.
  
  
  
  p <- p +  stat_density_2d(aes( fill = ..level..),        #fill = ..level..), 
                            alpha = 0.2, 
                            contour = TRUE,
                            contour_var = "density",
                            n = 120, # default 100 grid points in either dimension. computationally intense, do not exceed that value too much..
                            geom = "polygon", #geom: polygon or raster (equidistant lines or heatmap-ish picture)
                            colour="#99CBFE")+ 
    
    
    scico::scale_fill_scico(palette = "roma", direction=-1) # bilbao
  
  
  # first, lets take care of the intersting case where inx is part of gate(s)
  
  # now we check if xMarker or yMarker are part of any outside.*AB*.highpass or lowpass col and draw in and out with according colors:
  # the problem here is that if you have markercombinations that both carry flags, you dont know if to plot boolean AND or simply OR. #
  # We go with default AND here:

  
  
  
  
  p <- p +  
    scale_x_log10(#breaks = trans_breaks("log10", function(x) 10^x),
                  #labels = trans_format("log10", math_format(10^.x))
      ) +
    
    scale_y_log10(#breaks = trans_breaks("log10", function(x) 10^x),
                  #labels = trans_format("log10", math_format(10^.x))
      ) +
    
    
    #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
    #              minor_breaks = trans_breaks("log10", math_format(10^.x)),
    #              labels = trans_format("log10", math_format(10^.x))) +
    
    
    # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
    #                labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "trbl")+
    theme_bw()
    
    #theme_minimal()  
  # guides(colour = guide_legend(nrow = 8)) + # watch out, this induces binning into the continous color scale
  
  
  if (!is.na(incol)) {
    # here we prepare the title string. if incol is given, mention it in the title:
    
    p <- p +  labs(title=titlestring , 
                   subtitle=paste0("Plotted ", nrow(indata), " cells, ",sub("_[^_]+$", "", incol )," colored ",(1-percentile.colorlimits)*100,"-", percentile.colorlimits*100,"th percentile") 
                   #caption="Created by M.Barone", 
                   # y="Maximal Velocity [AU/sec]", 
                   #  x="Condition"
    )
    
    
    
  }else{
    
    p <- p +  labs(title=titlestring  , 
                   subtitle=subtitlestring
                   #caption="Created by M.Barone", 
                   # y="Maximal Velocity [AU/sec]", 
                   #  x="Condition"
    )
    
  }# end else bracket for when incol was NA and there is no need to plot this in the subtitle string 
  
  
  p <- p +  coord_fixed() + # this forces equidistant axis scaling. turn off aspect.ratio, of course...
    
    ylab( yMarker ) + 
    xlab( xMarker ) +  
    theme(#axis.title.x = element_blank(), 
      legend.title=element_text(size=13), # turn off with element_blank(),
      legend.text=element_text(size=11),
      axis.text=element_text(size=13),#,axis.text.x=element_text(size=10),
      axis.title.y=element_text(size=12), axis.title.x=element_text(size=12),
      # legend.position="top", 
      legend.key.width=unit(1,"line"),
      #aspect.ratio=1 #force ggplot on squared plotting # turn off if you force fixed coords of course
    ) 
  
  
  print(p)
  
  
} # end function def  plot.scatter.density.v2












# ............. Batch alignement percentages are globally defined and called like that within the function ......................
# if you decide to change these, you need to adjust them inside the function as well!!
# these are the percentiles we will explore when we attempt to align the batches.
# we will choose one of these per channel and drop that value into the batch.normalization.matrix.
batch.aligning.percentiles = c(0.6,0.8,0.85,0.90,0.95,0.99)
percentile_colors <- c('p60'='#a361c7','p80'='#5ba962','p85'='#c75a87','p90'='#ab973d','p95'='#648ace','p99'='#cb6342') 









# ............. GATTING via OUTLIER Flags 1 or active data point 0 flag: ......................
# The advantage of flagging outliers as 1 comes from my crstallography background, where reflections NOT participating in the model calculations are flagged with 1. 
# These belonged to a Rfree set. active data was always 0. 
#The second advantage is that if a datapoint passes x gates successfully, the sum of its flags remains 0, no matter how many gates you applied. So its easy to fish out such AND-connected gated cells if needed.

# just to clear out: the lower.penalty sets a lower limit under whose threshold the signal gets flagged as outlier. datapoints that remain active get a 0, outliers a 1.
# since we focus on outliers, the according flag column will be called outside.*AB*.highpass. So with the lower penalty we make a highpass gate. And flag the outliers 


# we define a stupid function that checks if i is smaller than j, i being the number to be checked against the gate/threshold j
# if this condition is met, we gonna flag the outlier with 1, if not we dont flip the flag


function.lower.penalty <- function(i,j){
  
  if( i < j ){
    #if gate condition is met, flip the gate flag:
    n <- 1
  }else{n<-0}
  return(n)
}


# and accordingly a second funtion to flag if i is above the penalty j
function.upper.penalty <- function(i,j){
  
  if( i > j ){
    #if gate condition is met, flip the gate flag:
    n <- 1
  }else{n<-0}
  return(n)
}


# for the area we create a function that checks if i is between the low j and high k threshold
#this one is different because we dont flag if i is in between, but flip the flag if either threshold is violated
function.range.penalty <- function(i,j,k){
  
  if( j <= i & i <= k ){
    #Only keep the cell unflagged, if i is both, bigger than j and the same time smaller than k
    n <- 0
  }else{n<-1}
  return(n)
}


# we further define a funtion to exclude strings in subset. This is because I wanna have spleen out and I dunno if we gonna have other 
# strings than pancreas further down the line. So lets make it the complicated way to exclude the string that will remain constant for sure: spleen
'%!in%' <- function(x,y)!('%in%'(x,y))









# this is the function to return the values of the multisegmented line looking along y!

multisegmentGATEline <- function(x.var, y.var){
  
  stopifnot(length(x.var) == length(y.var), sum(duplicated(x.var)) == 0)
  
  p <- order(x.var)
  x.var <- x.var[p]
  y.var <- y.var[p]
  
  k <- diff(y.var) / diff(x.var)
  l <- -1 * k * head(x.var, -1) + head(y.var, -1)
  
  function(x){
    
    ind <- findInterval(x, x.var)
    if(!all(between(ind, 1, length(x.var) - 1))) stop("\n\nGate threshold line is not defined over entire data range!")
    
    x * k[ind] + l[ind]
  }
}




# cluster collapse functions:




RMSD = function(x){
  sqrt( 1/length(x) * sum(  ( x - mean(x) )^2       )   )
}

MAE = function(x){
  1/length(x) * sum(   abs( x - mean(x) )  )
}



calc.homogeneity = function(m){
  
  # we first extract the currently present cluster numbers:
  all.clusters <-  sort(unique(m[[ ncol(m) ]] ))
  s<-0
  
  for(i in all.clusters){
    
    #message("Subsetting into cluster ",i)
    signals.in.cluster <- subset(m, V2==i)
    
    #now we drop out the last column:
    #signals.in.cluster <- select(signals.in.cluster,-last_col() )
    signals.in.cluster <- signals.in.cluster[, 1:( ncol(signals.in.cluster)-1  ) ]
    
    
    # there is certainly more to do with channel-wise RMSD calculations. but for now Im just interested in a global read-out of homegeneity of each cluster
    if (clusterhomogeneitymeasure == "MAE") {  cluster.sum <- sum( mapply( MAE , signals.in.cluster  ))   }
    if (clusterhomogeneitymeasure == "RMSD") {  cluster.sum <- sum( mapply( RMSD ,signals.in.cluster  ))  }
    
    
    
    
    #message(paste0("Cluster ", i, " RMSD: ", round( cluster.sum,3)  ))
    
    # we calculate RMSD of every marker, then sum up and add the RMSD-cluster-sum to the RMSD-subset-sum
    s <- s +  cluster.sum
    
    
  } # end run along all clusters to calculate each homogeneity 
  
  # Im pretty sure we need to correct for the amount of clusters we just ran through, since decreasing the cluster amount massively decreases the collected RMDS over them:
  s <- s / length(all.clusters)
  
  return(s) 
  
}# end function def calc.homogeneity


sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))


# proximity index functions:
extract.distance <- function(x=x,
                             y=y,
                             mask=mask){
  
  n =   mapply( function(x,y) mask[y,x] ,y=y, x=x   )
  return(n)
}

lookup.binary.mask <- function(x=x,
                               y=y#,
                              # mask=lymphnode.mask
                               ){
  ifelse(lymphnode.mask[y,x]==0, n = 0, n=1)
  return(n)
}



# for plotting log scale, we define some major and minor breaks:
breaks <- 10^(-10:10)
minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))


# for plotting and statistics, i just want to see the starts, no NS

# we wanna skip the comparisions that dont have significance:
sigFunc = function(x){
  if(x < 0.0001){"****"} 
  else if(x < 0.001){"***"} 
  else if(x < 0.01){"**"}
  else if(x < 0.05){"*"}
  else{"NS"}}

# just to bring the stars closer to the horizontal bar:
sigStar.verticaloffset <- 0.1







make.autograph.adapted <- function (dat, x.axis, y.axis, colour.by = x.axis, y.axis.label = y.axis, 
                                    grp.order = NULL, 
                                    colours = NULL, 
                                    polychromepalette = NULL, # this is the go-to way to color, you need to provide a table with colors via polychrome
                                    my_comparisons = NULL, 
                                    comparison.bar.tip.height = 0.04,
                                    inter.bracket.spacing = 0.1,
                                    Variance_test = NULL, 
                                    test= NULL,
                                    title = paste0(y.axis), 
                                    subtitle = NULL, 
                                    filename = paste0(y.axis, ".pdf"), 
                                    dot.alpha=1,
                                    violin.plot = FALSE,
                                    box.plot = FALSE, 
                                    errorbar.plot =FALSE,
                                    scale = "lin", dot.size = 3, width = 5, height = 5, max.y = 1.4, 
                                    path = getwd()) 
{
  if (!is.element("Spectre", installed.packages()[, 1])) 
    stop("Spectre is required but not installed")
  if (!is.element("ggplot2", installed.packages()[, 1])) 
    stop("ggplot2 is required but not installed")
  if (!is.element("data.table", installed.packages()[, 1])) 
    stop("data.table is required but not installed")
  if (!is.element("ggpubr", installed.packages()[, 1])) 
    stop("ggpubr is required but not installed")
  
  if (!is.element("ggbeeswarm", installed.packages()[, 1])) 
    stop("ggbeeswarm is required but not installed")
  
  
  require(ggplot2)
  require(data.table)
  require(ggpubr)
  require(ggbeeswarm)
  
  if (!is.null(colours)) {
    if (length(unique(dat[[colour.by]])) != length(colours)) {
      stop("The length of factors you want to colour the plot by does not match the number of colours you have provided.")
    }
  }
  message(paste0("AutoGraph for `", y.axis.label, " - ", y.axis, 
                 "` started"))
  #message("AutoGraph - setup started")
  spectral.list <- colorRampPalette(brewer.pal(11, "Spectral"))(50)
  spectral.list <- rev(spectral.list)
  colour.scheme <- colorRampPalette(c(spectral.list))
  dat <- data.table::as.data.table(dat)
  max_y_value <- max(dat[, y.axis, with = FALSE], na.rm = TRUE)
  max_y_value_p40 <- max_y_value * max.y
  min_y_value <- min(dat[, y.axis, with = FALSE], na.rm = TRUE)
  bottom_y <- min_y_value
  Xaxis <- dat[[x.axis]] <- as.factor(dat[[x.axis]])
  
  
  if (is.null(grp.order)) {
    if (is.null(colours)) { colours <- colour.scheme(length(unique(dat[[x.axis]])))       }
    #if (colours == "polychrome") {  colours <-  polychromepalette   }
  }
  if (!is.null(grp.order)) {
    Xaxis <- dat[[x.axis]] <- factor(dat[[x.axis]], levels = grp.order)
    if (is.null(colours)) { colours <- colour.scheme(length(as.vector(grp.order))) }
    # if (colours == "polychrome") {   colours <-  polychromepalette }
    
    
  }
  # message("AutoGraph - setup complete")
  #message("AutoGraph - plotting started")
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = Xaxis, 
                                         y = dat[[y.axis]],   
                                        # fill = as.factor(dat[[colour.by]]), # is this needed?
                                         colour = as.factor(dat[[colour.by]]))    
  )+
    geom_quasirandom(size=dot.size, alpha=dot.alpha) +
    guides(colour = guide_legend(override.aes = list(size=4, alpha=1)  ) )+ #we gonna override the dot size how its depicted in the legend  byrow=TRUE)  # and I want the elements listed row-by-row
    
   # scale_fill_manual(name = colour.by, values = colours)+ # is this needed?
    scale_color_manual(name = colour.by, values = colours)+
    
    labs(title=title  , 
         subtitle=subtitle,
         #caption="Created by M.Barone", 
         x = paste0(x.axis), 
         y = y.axis.label
    )+
    theme_classic(base_size = 30)+
    theme(legend.position = "right", 
          legend.text = ggplot2::element_text(colour = "black",
                                              size = 10), #, angle = 0, hjust = 0, vjust = 0.5, face = "bold" all were bold...
          legend.title = element_blank(), 
          axis.text.x = ggplot2::element_text(colour = "black", 
                                              size = 12, angle = 45, hjust = 1, vjust = 1), 
          axis.text.y = ggplot2::element_text(colour = "black", 
                                              size = 12, angle = 0, hjust = 1, vjust = 0), 
          axis.title.x = ggplot2::element_text(color="grey15", size=11, angle = 0, hjust = 0.5, vjust = 0), 
          axis.title.y = ggplot2::element_text(colour = "black", 
                                               size = 12, angle = 90, hjust = 0.5, vjust = 1), 
          plot.title = ggplot2::element_text(lineheight = 0.8,, hjust = 0, size = 12), 
          plot.subtitle = ggplot2::element_text(size = 10, color = "black"), #,face = "italic"
          axis.line = ggplot2::element_line(colour = "black",   size = 1)
    )
  
  
  
  
  
  
  if(box.plot == TRUE){ p <- p + geom_boxplot(outlier.shape = NA, lwd=0.4,alpha=0.1,show.legend = FALSE)  }
  if (violin.plot == TRUE) { p <- p + geom_violin(ggplot2::aes(fill = as.factor(dat[[colour.by]]), 
                                                               colour = as.factor(dat[[colour.by]])), trim = FALSE, 
                                                  show.legend = FALSE, alpha = 0.1)
  }
  
  if(errorbar.plot ==  TRUE  ){
    p <- p + ggplot2::stat_summary(fun.max = function(i) mean(i) + 
                                     sd(i)/sqrt(length(i)), fun.min = function(i) mean(i) - 
                                     sd(i)/sqrt(length(i)), geom = "errorbar", width = 0.5, 
                                   size = 1)
    p <- p + ggplot2::stat_summary(fun = mean, fun.min = mean, 
                                   fun.max = mean, geom = "crossbar", width = 0.7, size = 0.5)
    
    
  }
  
  
  
  
  if (!is.null(my_comparisons)) {
    if (!is.null(test)) {
      p <- p + stat_compare_means(comparisons = my_comparisons, 
                                  method = test,
                                  tip.length=comparison.bar.tip.height,
                                  step.increase = inter.bracket.spacing
      )
    }
  }
  if (!is.null(Variance_test)) {
    p <- p + stat_compare_means(method = Variance_test, 
                                label.y = max_y_value_p40, size = 4)
  }
  
  
  if (scale == "lin") {
    p <- p + ggplot2::scale_y_continuous(limits = c(0, max_y_value_p40))
  }
  if (scale == "sci") {
    p <- p + ggplot2::scale_y_continuous(labels = scales::scientific, 
                                         limits = c(0, max_y_value_p40))
  }
  p
  ggplot2::ggsave(plot = p, filename = paste0(filename), width = width, 
                  height = height, path = path)
  # message(paste0("AutoGraph for `", y.axis.label, " - ", y.axis, 
  #                "` saved to disk"))
  return(p)
} # end function def make.autograph.adapted

