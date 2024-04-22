# This script is to enhance the spatial ATAC data using scATAC-seq and spatial ATAC data. The neighbor identifying process in enhancement could be performed on the cell-type pairs level (more accurate) or the all-data level

# enhance the spatial ATAC data using scATAC-seq data and spatial ATAC data
#######################################################################
Enhancement <- function(scATAC, # scATAC seurat object
                        spATAC, # spATAC seurat object
                        method = "rlsi", # integration method, lsi or cca
                        processing_scATAC = TRUE, # whether scATAC-seq data will be processed
                        processing_spATAC = TRUE, # whether spATAC data will be processed
                        spATAC.assay = "peaks", # the assay data used for integration in spATAC data
                        scATAC.assay = "peaks", # the assay data used for integration in scATAC-seq data
                        initial.neighbor = 100, # the initial neighbors identified for spots, default 100
                        cell_type, # cell types data for enhancement
                        dims = 2:30, # the dimensions when perform the integration
                        k = 10, # the number of final single-cell neighbors for each spot
                        cell_type_pair = TRUE, # whether the neighbor identification restricted to cell types. Annotation should be represent in scATAC and spATAC data when it is true.
                        k.weight = 20, # k.weight when transferring data
                        lambda = 0.9, # lambda parameter when performing Harmony to integrate
                        savedata = TRUE, 
                        saveplot = TRUE
){
  if(dir.exists("./Enhancement/Plot")){
    print("Directory exists")
  }else{
    dir.create("./Enhancement/Plot", recursive = T)
  }
  
  if(cell_type_pair == TRUE){
    print("Neighbor identification will be restrict to the identical cell type")
  }else{
    print("Neighbor identification will not be restrict to the identical cell type")
  }
  
  print("Enhance spatial ATAC profile on: ")
  print(cell_type)
  scATAC <- subset(scATAC, annotation %in% cell_type)
  spATAC <- subset(spATAC, annotation %in% cell_type)

  
  spATAC$annotation <- paste0("spatial_", spATAC$annotation)
  scATAC$annotation <- paste0("scATAC_", scATAC$annotation)
  spATAC$dataset <- "spATAC"
  scATAC$dataset <- "scATAC"
  DefaultAssay(scATAC) <- spATAC.assay
  DefaultAssay(spATAC)<- scATAC.assay  
  
  # processing scATAC-seq data 
 if(processing_scATAC){
   print("Enhance spATAC: 1/8.scATAC processing")
   sc_anno <- c()
   for(i in unique(cell_type)){
     sc_anno <- append(sc_anno, grep(i, unique(scATAC$annotation), ignore.case = T, value = T))
   }
   scATAC <- subset(scATAC, annotation%in% sc_anno)
   scATAC <- RunTFIDF(scATAC)
   scATAC <- FindTopFeatures(scATAC, min.cutoff = 'q5')
   scATAC <- RunSVD(scATAC)
   scATAC <- RunUMAP(object = scATAC, reduction = 'lsi', dims = 2:20)
 }else{print("scATAC-seq data will be not processed.")}
  
  # processing spATAC data
 if(processing_spATAC){
   print("Enhance spATAC: 2/8.spATAC processing")
   sp_anno <- c()
   for(i in unique(cell_type)){
     sp_anno <- append(sp_anno, grep(i, unique(spATAC$annotation), ignore.case = T, value = T))
   }
   spATAC <- subset(spATAC, annotation%in% sp_anno)
   spATAC <- RunTFIDF(spATAC)
   spATAC <- FindTopFeatures(spATAC, min.cutoff = 'q5')
   spATAC <- RunSVD(spATAC)
   spATAC <- RunUMAP(object = spATAC, reduction = 'lsi', dims = 2:30)
 }else{print("spATC data will be not processed.")}
  
  # combine datasets
  combined <- merge(spATAC, scATAC)
  
  # integrate two datasets
  print("Enhance spATAC: 3/8.Integrate two datasets")
  integration.anchors <- FindIntegrationAnchors(
    object.list = list(scATAC, spATAC),
    anchor.features = intersect(rownames(spATAC), rownames(scATAC)),
    reduction = method,
    assay = c(spATAC.assay, scATAC.assay),
    dims = 2:20
  )
  
  integrated <- IntegrateData(anchorset = integration.anchors)
  integrated <- RunSVD(integrated)
  integrated <- RunUMAP(object = integrated, reduction = 'lsi', dims = 2:20)
  
  if(saveplot){
    p1 <- DimPlot(integrated, group.by = "dataset", pt.size = 0.5)
    ggsave(p1, filename = "./Enhancement/Plot/F1.UMAP_integrated(dataset).pdf", height = 6, width = 7)
    p2  <-  DimPlot(integrated, group.by = "annotation", pt.size = 0.5)
    ggsave(p2, filename = "./Enhancement/Plot/F2.UMAP_integrated(annotation).pdf", height = 6, width = 7)
  }
  
  # run Harmony to integrate again
  print("Enhance spATAC: 4/8.Second-round integration")
  integrated <- suppressWarnings(harmony::RunHarmony(
    object = integrated,
    group.by.vars = 'dataset',
    reduction = 'lsi',
    plot_convergence = F,
    theta = 3,
    lambda = lambda,
    assay.use = "integrated",
    verbose = FALSE,
    max.iter.harmony = 60,
    project.dim = F
  ))
  
  integrated <- RunUMAP(object = integrated, reduction = "harmony", dims = 2:20)
  
  if(saveplot){
    p1 <- DimPlot(integrated, group.by = "dataset", pt.size = 0.5)
    ggsave(p1, filename = "./Enhancement/Plot/F3.UMAP_integrated_harmonized(dataset).pdf", height = 6, width = 7)
    p2 <- DimPlot(integrated, group.by = "annotation", pt.size = 0.5)
    ggsave(p2, filename = "./Enhancement/Plot/F4.UMAP_integrated_harmonized(annotation).pdf", height = 6, width = 7)
  }
  
  if(savedata){
    message("saving data: Integrated_Harmonized_data.rds")
    saveRDS(integrated, file="./Enhancement/1.Integrated_Harmonized_data.rds")
  }
  
  # find neighbors
  if(cell_type_pair == TRUE){
    print("Enhance spATAC: 5/8.Find neighbors (cell type level)")
    final_neigbors <- list()
    
    for (i in cell_type) {
      print(paste0("Find neighbors:", i))
      anno1 <- grep(i, unique(integrated$annotation), value = T, ignore.case = T)
      temp <- subset(integrated, annotation %in% anno1)
      temp <- Seurat::FindNeighbors(temp, 
                                    k.param = initial.neighbor,
                                    reduction = "harmony",
                                    dims = dims, 
                                    return.neighbor = T, 
                                    assay = "integrated")
      neigbors <- suppressWarnings(sapply(colnames(subset(temp, subset = dataset == "scATAC", invert = TRUE)), function(i){
        Seurat::TopNeighbors(temp@neighbors$integrated.nn, i, n = initial.neighbor)
      }))
      final_neigbors[[i]] <- neigbors
      
    }}else{
      print("Enhance spATAC: 5/8.Find neighbors (all-data level)")
      integrated <- Seurat::FindNeighbors(integrated, 
                                          k.param = initial.neighbor,
                                          reduction = "harmony",
                                          dims = dims, 
                                          return.neighbor = T, 
                                          assay = "integrated")
      neigbors <- suppressWarnings(sapply(colnames(subset(integrated, subset = dataset == "scATAC", invert = TRUE)), function(i){
        Seurat::TopNeighbors(integrated@neighbors$integrated.nn, i, n = initial.neighbor)
      }))
      final_neigbors <- neigbors
    }
  
  final_neigbor <- as.data.frame(final_neigbors)
  
  if(savedata){
    message("saving data: All_initial_neighbors.rds")
    saveRDS(final_neigbor,file = "./Enhancement/2.All_initial_neighbors.rds") 
  }

  # filtered neighbors
  print("Enhance spATAC: 6/8.Filtered neighbors")
  integrated_scATAC <- suppressWarnings(subset(integrated, subset = dataset == "scATAC"))
  integrated_spATAC <- suppressWarnings(subset(integrated, subset = dataset == "spATAC"))
  int_cell_scATAC <- colnames(integrated_scATAC)
  int_cell_spATAC <- colnames(integrated_spATAC)
  
  neigbor <- final_neigbor[1:k,]
  for (i in colnames(final_neigbor)) {
    scATAC <- intersect(final_neigbor[,i],int_cell_scATAC)[1:k]
    if(length(scATAC) == 0){
      print("erro: 0 scATAC found")
    }
    neigbor[,i] <- scATAC
  }

  if(savedata){
    message("saving data: Filtered_scNeighbors.rds")
    saveRDS(neigbor, file = "./Enhancement/3.Filtered_scNeighbors.rds")
  }
  
  # find neighbors pairs
  print("Enhance spATAC: 7/8.Find neighbors pairs")
  colnames(neigbor) <- sapply(colnames(neigbor), function(x){
    unlist(strsplit(x, "[.]"))[2]
  }) %>% unlist()
  
  data <- data.frame()
  for ( i in names(neigbor)) {
    cells <- neigbor[[i]]
    indx1 <- which(int_cell_scATAC %in% cells)
    indx2 <- which(int_cell_spATAC %in% paste0(i,"-1"))
    if(length(indx1) == 0){
      indx1 <- "none"
      print("warning: 0 scATAC found")
    }
    tmp <- data.frame(cell1 = indx1, cell2 = indx2, score = 1) 
    data <- rbind(data, tmp)
  }
  
  if(length(which(data$indx1 == "none")) != 0){
    data <- data[-which(data$indx1 == "none"),]
  }
  data$cell1 <- as.numeric(data$cell1)
  data$cell2 <- as.numeric(data$cell2)
  data$score <- as.numeric(data$score)
  
  
  if(savedata){
    message("saving data: Neighbors_pairs.rds")
    saveRDS(data, file = "./Enhancement/4.Neighbors_pairs.rds")
  }
  
  # transfer data
  print("Enhance spATAC: 8/8.Transfer data")
  transfer.anchors <- FindTransferAnchors(
    reference = integrated_scATAC,
    query = integrated_spATAC,
    reference.reduction = "lsi",
    reduction = "lsiproject",
    dims = 2:30
  )
  transfer.anchors@anchors <- data
  enhanced_spATAC <- TransferData(
    anchorset = transfer.anchors,
    refdata = GetAssayData(integrated_scATAC, assay = scATAC.assay, slot = "counts"),
    weight.reduction = integrated_spATAC[["harmony"]],
    dims = 2:30, 
    k.weight = k.weight
  )
  return(enhanced_spATAC)
  if(savedata){
    message(paste0("saving data: ","Enhanced_spATAC_assay_k", k, ".rds"))
    saveRDS(enhanced_spATAC, file = paste0("./Enhancement/5.Enhanced_spATAC_assay_k", k, ".rds"))
  }
  return(enhanced_spATAC)
  print("Finished enhancement!!!")
  
}
