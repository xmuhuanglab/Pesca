mySpatialDimPlot <- function(metadata, # the meta.data in the seurat objectBDB76B
                          row, # the column which contains the x coordinate information
                          col, # the column which contains the y coordinate information
                          color, # the column contain cluster column
                          pointsize = 2.5, # the spot size for visualization
                          colors = NULL # the color for the cluster
                          ){
  library(ggplot2)
  color_ref <- c("#F5DEB3", "#AFEEEE", "#0000FF", "#40E0D0", "#00CED1", "#A52A2A", "#FFD700", "#FF7F50", "#6495ED", "#00FA9A", "#9370DB", "#00FFFF", "#FFA07A", "#F08080", "#008000", "#FF69B4", "#EEE8AA", "#FF69B4", "#B22222", "#C71585", "#006400", "#556B2F", "#F5F5DC", "#9400D3", "#FFFF00", "#ADFF2F", "#DAA520", "#F0E68C", "#FFEBCD", "#8A2BE2", "#BDB76B", "#32CD32", "#FFC0CB", "#F5DEB3", "#F0FFF0", "#FFE4B5", "#FF4500", "#1E90FF", "#FAEBD7")
  
  if(is.null(colors)){
    colors <- color_ref[1:length(unique(metadata[,color]))]
  }
  temp <- metadata[,c(row, col, color)]
  colnames(temp) <- c("a", "b", "c")

  p <- ggplot(temp, aes(a, b, color = c)) + geom_point(size = pointsize) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = 'white'), plot.background = element_rect(fill = "white")) + theme(legend.title = element_blank(), legend.key = element_rect(fill = 'white'), legend.text = element_text(size = 10), legend.key.size = unit(1, 'cm') )+scale_color_manual(values = colors)
  return(p)
}



mySpatialFeaturePlot <- function(object, # seurat object
                                 feature, # the features used to plot 
                                 pointsize = 2.5, # spot size
                                 scale = FALSE, # whether scale the data
                                 ncol = 2, # the number of column for plots when the number of features is more than two
                                 assay, # the assay data used to plot 
                                 slot = "data" # the slot in the assay 
                                 ){
  print(feature)
  data <- object@meta.data
  data <- data[,c("row", "col")]
  data$cell <- rownames(data)
  rownames(data) <- NULL
  
  assaydata <- GetAssayData(object = object, assay = assay, slot = slot)
  finaldata <- data.frame()
  
  # scale the data
  if(scale){
    for (i in feature) {
      temp <- data
      temp$feature <- i
      feamat <- as.data.frame(assaydata[i,])
      temp$expre <- apply(feamat, 2, function(x){
        (x-min(x))/(max(x)-min(x))
      })
      finaldata <- rbind(finaldata, temp)
    }
  }else{
    for (i in feature) {
      temp <- data
      temp$feature <- i
      feamat <- as.data.frame(assaydata[i,])
      temp$expre <- feamat$`assaydata[i, ]`
      finaldata <- rbind(finaldata,temp)
    }
  }
  finaldata$expre <- as.numeric(finaldata$expre)
  
  # plot data
  if(length(feature) == 1){
    p <- ggplot(finaldata, aes(row, col, color = expre)) + geom_point(size = pointsize) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = 'white'), plot.background = element_rect(fill = "white")) + theme(legend.title = element_blank(), legend.key = element_rect(fill = 'white'), legend.text = element_text(size = 10), legend.key.size = unit(1, 'cm') ) + scale_color_gradientn(colours = c(colorRampPalette(c("#5b51a3", "#79c9a4", "#f2faac", "#fdb465", "#a4104d"))(90)), limits = c(min(finaldata$expre), max(finaldata$expre)))
    return(p)
  }else{
    p <- ggplot(finaldata, aes(row, col, color = expre)) + geom_point(size = pointsize) + facet_grid(~feature) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = 'white'), plot.background = element_rect(fill = "white"))+theme(legend.title = element_blank(), legend.key = element_rect(fill = 'white'), legend.text = element_text(size = 10),legend.key.size=unit(1, 'cm') ) + scale_color_gradientn(colours = c(colorRampPalette(c("#5b51a3", "#79c9a4", "#f2faac", "#fdb465", "#a4104d"))(90)), limits = c(min(finaldata$expre), max(finaldata$expre)));p
    return(p)
  }
}




