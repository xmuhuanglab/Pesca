# This script is to call peak and generated Seurat object for spATAC and scATAC seq data. For the scATAC-seq data, the gene score is also calculated using Signac.


# 1. Call peaks for scATAC-seq data
#######################################################################
Callpeak_scATAC <- function(organs, # organs related to the separate fragment file names used to call peak
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
  for(i in organs){
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
scATAC_processing <- function(fragpath, # path storing fragment file (tabix)
                              peaks, # merged peaks in scATAC-seq data
                              organs, # vector containing the organ names
                              mincell = 100,
                              minfeature = 1000,
                              genome = "mm10", # mm10 or hg19 or hg38
                              sep = c("-", "-"), # separate the peak
                              annotation
){
  allmatrix <- list()
  cellanno <- list()
  fragObjects <- list()
  
  for (f in organs) {
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
  
  saveRDS(sc_mat, file = paste0("./Processing_scATAC/3.scATAC_matrix.rds"))
  saveRDS(cellanno, file = paste0("./Processing_scATAC/4.scATAC_cell_annotation.rds"))
  
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
    annotations <- annotation$v79
  }else if(
    genome == "hg19"
  ){
    annotations <- annotation$v75
  }else{
    annotations <- annotation$v86
  }
  seqlevelsStyle(annotations) <- 'UCSC'
  Annotation(sc_atac) <- annotations
  
  ### add annotation information
  sc_atac$annotation = NA
  sc_atac$cells = rownames(sc_atac@meta.data)
 
  for (o in names(cellanno)) {
    sc_atac@meta.data[which(sc_atac@meta.data$cells %in% cellanno[[o]]), "annotation"] = o
  }

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
  saveRDS(sc_atac, file = paste0("./Processing_scATAC/5.scATAC_seuratobject.rds"))
  return(sc_atac)
                                                               }


# 3. Call peak and create Seurat Object for spATAC data
#######################################################################
spATAC_processing <- function(fragpath_spatac, # path storing ATAC fragment file
                              peaks, # merged peaks set generated from corresponding scATAC-seq data
                              spotIDs = NULL, # spots used to generate matrix
                              spots_anno = NULL, # list containing spots annotation
                              genome = "mm10", # genome annotation 
                              minspot = 5, # the minimum number of spots which peaks are accessible in 
                              minfeature = 10, # the minimum number of peaks that a spot contains
                              annotation
){
  if(dir.exists("./Processing_spATAC/")){
    print("Directory exists")
  }else{
    dir.create("./Processing_spATAC/", recursive = T)
  }
  
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
  
  if(genome == "mm10"){
    annotations <- annotation$v79
  }else if(
    genome == "hg19"
  ){
    annotations <- annotation$v75
  }else{
    annotations <- annotation$v86
  }
  seqlevelsStyle(annotations) <- 'UCSC'
  Annotation(sp_atac) <- annotations
  
  ### add annotation information
  sp_atac$annotation = NA
  sp_atac$cells = rownames(sp_atac@meta.data)
  
  if(is.null(spots_anno)){
    sp_atac$annotation <- NA
  }else{
    for (o in names(spots_anno)) {
      sp_atac@meta.data[which(sp_atac@meta.data$cells %in% spots_anno[[o]]), "annotation"] <-  o
    }
  }
  
  print("Processing spATAC: 3/3.calculate gene activity")
  ga <- GeneActivity(sp_atac, extend.upstream = 2000, extend.downstream = 2000)
  sp_atac[['gene_activity']] <- CreateAssayObject(counts = ga[, intersect(colnames(ga), colnames(sp_atac))])
  sp_atac <- NormalizeData(
    object = sp_atac,
    assay = 'gene_activity',
    normalization.method = 'LogNormalize',
    scale.factor = median(sp_atac$nCount_gene_activity)
  )
  saveRDS(sp_atac, file = paste0("./Processing_spATAC/1.SeuratObject_spATAC.rds"))
  return(sp_atac)
}

