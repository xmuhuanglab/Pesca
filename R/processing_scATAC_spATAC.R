# This script is to call peak and generated Seurat object for spATAC and scATAC seq data. For the scATAC-seq data, the gene score is also calculated using Signac.


# 1. Call peaks for scATAC-seq data
#######################################################################
Callpeak_scATAC <- function(cell_type, # cell_type related to the separate fragment file names used to call peak
                            fragpath, # path storing fragment file (tabix)
                            macs2_path # path storing  macs2 software
){
  # call peak for scATAC data separately and filtered
  if(dir.exists("./Processing_scATAC/")){
    print("Directory exists")
  }else{
    dir.create("./Processing_scATAC/", recursive = T)
  }
  
  allpeaks <- list()
  for(i in cell_type){
    print("Callpeaks: 1/2.call peak for scATAC-seq data")
    path <- grep(i, fragpath, value = T, ignore.case = T)
    peaks_all <- CallPeaks(path, macs2.path = macs2_path)
    peaks_all <- keepStandardChromosomes(peaks_all, pruning.mode = "coarse")
    peaks_all <- subsetByOverlaps(x = peaks_all, ranges = blacklist_mm10, invert = TRUE)
    
    peaks <- as.data.frame(peaks_all)
    peaks <- peaks[which(peaks$neg_log10qvalue_summit > 2 & peaks$fold_change > 2),]
    
    peaks <- makeGRangesFromDataFrame(peaks)
    allpeaks[[i]] <- peaks
  }
  saveRDS(allpeaks, file = "./Processing_scATAC/1.scATAC_allpeaks.rds")
  print("Finished calling peak for scATAC-seq data!")
  
  # merge all peaks
  print("Callpeaks: 2/2.merge peaks for scATAC-seq data")
  extracted_data <- list()
  for (i in seq_along(allpeaks)) {
    extracted_data[[i]] <- allpeaks[[i]]
  }
  combined_iranges <- do.call(c, extracted_data)
  
  peaks <- reduce(x = combined_iranges, drop.empty.ranges = T)
  peaks_width <- width(peaks)
  
  saveRDS(peaks, file = paste0("./Processing_scATAC/2.scATAC_mergedpeaks.rds"))
  print(paste0("Finished merging peaks: total ", length(peaks), " peaks!"))
  return(peaks)
}


# 2.Generate Matrix and create SeuratObject for scATAC-seq data
#######################################################################
scATAC_processing <- function(fragpath = NULL, # path storing fragment file (tabix)
                              peaks = NULL, # merged peaks in scATAC-seq data
                              cell_type = NULL, # vector containing the organ names
                              sc_mat = NULL,
                              mincell = 100,
                              minfeature = 1000,
                              cell_anno = NULL, # a vector contain the cell identify
                              genome = "mm10", # mm10 or hg19 or hg38
                              sep = c("-", "-"), # separate the peak
                              genome_anno = NULL,
                              savedata = TRUE,
                              gene.activity = TRUE
){
  ###################################
  if(is.null(fragpath) & is.null(peaks)){  # using the peak_by_cell matrix of scATAC-seq data to create seurat object
    message("Using peak-by-cell matrix to create seurat object, please note that gene activity will not be calculated because of the lack of fragment data")
    chrom_assay <- CreateChromatinAssay(
      counts = sc_mat,
      sep = sep,
      genome = genome,
      min.cells = mincell,
      min.features = minfeature,
    )
    sc_atac <- CreateSeuratObject(counts = chrom_assay,assay = "peaks")
    print(paste0("Number of total cells:", ncol(sc_atac)))
    print(paste0("Number of total peaks:", nrow(sc_atac)))
    
    if(!is.null(cell_anno)){
      cells <- intersect(Cells(sc_atac), cell_anno$cell)
      if(length(cells) == 0){print("Could not add annotation, please check the cell barcode in cell_anno and sc_mat")}else{
        sc_atac <- sc_atac[,cells]
        rownames(cell_anno) <- cell_anno$cell
        sc_atac$annotation = cell_anno[Cells(sc_atac),"anno"]
      }
    }
    
    if(savedata){
      message("saving data: scATAC_object.rds")
      saveRDS(sc_atac, file = paste0("./Processing_scATAC/5.scATAC_object.rds"))
    }
    return(sc_atac)
    
  }else if(!is.null(fragpath) & !is.null(peaks)){ # using fragment data of scATAC-seq data to create seurat object
    ###################################
    
    message("Using fragment data to create seurat object")
    allmatrix <- list()
    cellanno <- list()
    fragObjects <- list()
    
    for (f in cell_type) {
      fragment <- data.table::fread(grep(f, fragpath, value = T, ignore.case = T))
      cells <- unique(as.vector(fragment$V4))
      temp <- CreateFragmentObject(path = grep(f, fragpath, value = T, ignore.case = T), cells = cells, validate.fragments = F)
      mat <- FeatureMatrix(fragments = temp, features = peaks)
      cellanno[[f]] <- colnames(mat)
      allmatrix[[f]] <- mat
      fragObjects[[f]] <- temp
      print(paste0("Finished matrix generation:", f, "!"))
    }
    rm(fragment)
    sc_mat <- do.call(cbind, allmatrix)
    sc_mat <- as(sc_mat, "sparseMatrix")
    
    if(savedata){
      message("save data: scATAC_matrix.rds")
      saveRDS(sc_mat, file = paste0("./Processing_scATAC/3.scATAC_matrix.rds"))
      message("save data: scATAC_cell_annotation.rds")
      saveRDS(cellanno, file = paste0("./Processing_scATAC/4.scATAC_cell_annotation.rds"))
    }
    
    # create Seurat object
    print("Processing scATAC: 1/2.dimension reduction")
    chrom_assay <- CreateChromatinAssay(
      counts = sc_mat,
      sep = sep,
      genome = genome,
      fragments = fragObjects,
      min.cells = mincell,
      min.features = minfeature,
    )
    
    sc_atac <- CreateSeuratObject(
      counts = chrom_assay,
      assay = "peaks"
    )
    
    print(paste0("Number of total cells:", ncol(sc_atac)))
    print(paste0("Number of total peaks:", nrow(sc_atac)))
    
    if(genome == "mm10"){
      annotations <- genome_anno$v79
    }else if(
      genome == "hg19"
    ){
      annotations <- genome_anno$v75
    }else{
      annotations <- genome_anno$v86
    }
    seqlevelsStyle(annotations) <- 'UCSC'
    Annotation(sc_atac) <- annotations
    
    ### add annotation information
    sc_atac$annotation = NA
    sc_atac$cells = rownames(sc_atac@meta.data)
    
    for (o in names(cellanno)) {
      sc_atac@meta.data[which(sc_atac@meta.data$cells %in% cellanno[[o]]), "annotation"] = o
    }
    
    if(gene.activity){
      print("Processing scATAC: 2/2.calculate gene activity")
      # add gene activity
      ga <- GeneActivity(sc_atac, extend.upstream = 2000, extend.downstream = 2000)
      sc_atac[['gene_activity']] <- CreateAssayObject(counts = ga[,intersect(colnames(ga), colnames(sc_atac))])
      sc_atac <- NormalizeData(
        object = sc_atac,
        assay = 'gene_activity',
        normalization.method = 'LogNormalize',
        scale.factor = median(sc_atac$nCount_gene_activity)
      )
    }
    if(savedata){
      message("saving data: scATAC_object.rds")
      saveRDS(sc_atac, file = paste0("./Processing_scATAC/5.scATAC_object.rds"))
    }
    return(sc_atac)
  } 
  }



# 3. Call peak and create Seurat Object for spATAC data
#######################################################################
spATAC_processing <- function(fragpath_spatac = NULL, # path storing ATAC fragment file
                              peaks = NULL, # merged peaks set generated from corresponding scATAC-seq data
                              spotIDs = NULL, # spots used to generate matrix
                              spot_anno = NULL, # a data frame containing the identity of spot and spot barcode
                              sp_mat = NULL,
                              genome = "mm10", # genome annotation 
                              minspot = 5, # the minimum number of spots which peaks are accessible in 
                              minfeature = 10, # the minimum number of peaks that a spot contains
                              savedata = FALSE
){
  if(dir.exists("./Processing_spATAC/")){
    print("Directory exists")
  }else{
    dir.create("./Processing_spATAC/", recursive = T)
  }
  ###################################
  if(!is.null(sp_mat)){
    message("Using peak-by-spot matrix to create seurat object.")
    if(is.null(spotIDs)){
      chrom_assay <- CreateChromatinAssay(
        counts = sp_mat,
        sep = sep,
        genome = genome,
        min.cells = minspot,
        min.features = minfeature,
      )
    }else{
      chrom_assay <- CreateChromatinAssay(
        counts = sp_mat[, spotIDs],
        sep = sep,
        genome = genome,
        min.cells = minspot,
        min.features = minfeature,
      )
    }
    sp_atac <- CreateSeuratObject(counts = chrom_assay,assay = "peaks")
    print(paste0("Number of total spots:", ncol(sp_atac)))
    print(paste0("Number of total peaks:", nrow(sp_atac)))
    
    if(!is.null(spot_anno)){
      if(length(spot_anno) != ncol(sp_mat)){message("The length of cell identify is not equal to cell number.")} else{
        spots <- intersect(Cells(sp_atac), spot_anno$spot)
        if(length(spots) == 0){print("Could not add annotation, please check the spot barcode in spot_anno and sp_mat")}else{
          sp_atac <- subset(sp_atac, cells = spots)
          rownames(spot_anno) <- spot_anno$spot
          sp_atac$annotation = spot_anno[Cells(sp_atac),"anno"]
        }
        }
    }
    if(savedata){
      message("saving data: spATAC_object.rds")
      saveRDS(sc_atac, file = paste0("./Processing_spATAC/1.spATAC_object.rds"))
    }
    return(sp_atac)
  }
  
  ###################################
  if(is.null(spotIDs)){
    fragment <- data.table::fread(fragpath_spatac)
    spotIDs <- unique(fragment$V4)
    print("Processing spATAC: 1/3.using all spots data to call peaks")
  }
  fragObject <- CreateFragmentObject(path = fragpath_spatac, cells = spotIDs)
  spatac_mat <- FeatureMatrix(fragments = fragObject, features = peaks)
  spatac_mat <- as(spatac_mat, "sparseMatrix")
  
  print("Processing spATAC: 2/3.creatate SeuratObject")
  chrom_assay <- CreateChromatinAssay(
    counts = spatac_mat,
    sep = c("-", "-"),
    genome = genome,
    fragments = fragObject,
    min.cells = minspot,
    min.features = minfeature
  )
  sp_atac <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks"
  )
  
  print(paste0("Number of total cells:", ncol(sp_atac), "in spatial ATAC data"))
  print(paste0("Number of total cells:", ncol(sp_atac), "in spatial ATAC data"))
  
  ### add annotation information
  if(!is.null(spot_anno)){
      spots <- intersect(Cells(sp_atac), spot_anno$spot)
      if(length(spots) == 0){print("Could not add annotation, please check the spot barcode in spot_anno and fragment file")}else{
        sp_atac <- sp_atac[,spots]
        rownames(spot_anno) <- spot_anno$spot
        sp_atac$annotation = spot_anno[Cells(sp_atac),"anno"]
      }
  }

  if(savedata){
    message("saving data: spATAC_object.rds")
    saveRDS(sp_atac, file = paste0("./Processing_spATAC/1.spATAC_object.rds"))
  }
  return(sp_atac)
}

