# This script is to predict the spatial ATAC data using scATAC-seq and spatial RNA data. The neighbor identifying process in prediction could be performed on the cell-type pairs level (more accurate) or the all-data level

# predict the spatial ATAC data using scATAC-seq data and spatial RNA data
#######################################################################
Prediction <- function(scATAC, # scATAC seurat object
                       spRNA, # spRNA seurat object
                       processing_scATAC = TRUE, # whether scATAC-seq data will be processed
                       processing_spRNA = TRUE, # whether spRNA data will be processed
                       method = "rpca", # integration method, rpca or cca
                       k.weight = 50, # # k.weight when transferring data
                       spRNA.assay = "spatial", # the name of count assay in spRNA data
                       scATAC.assay = NULL, # the name of count assay in spRNA data
                       impute.assay = "peaks",
                       cell_type, # cell type data used for enhancement
                       initial.neighbor = 100, # the initial neighbors identified for spots
                       dims = 1:30, # the dimensions when perform the integration
                       k, # the number of final single-cell neighbors for each spot
                       cell_type_pair = TRUE, 
                       lambda = 0.9, # lambda parameter when performing Harmony to integrate
                       savedata = TRUE,
                       saveplot = TRUE
){
  if(dir.exists("./Prediction/Plot")){
    print("Directory exists")
  }else{
    dir.create("./Prediction/Plot/", recursive = T)
  }
  spRNA$dataset <- "spRNA"
  scATAC$dataset <- "scATAC"
  DefaultAssay(scATAC) <- scATAC.assay
  DefaultAssay(spRNA) <- spRNA.assay  
  
  print("Predict spatial ATAC profile on: ")
  print(cell_type)
  scATAC <- subset(scATAC, annotation %in% cell_type)
  spATAC <- subset(spATAC, annotation %in% cell_type)
  
  # processing scATAC 
  if(processing_scATAC){
    print("Predict spATAC: 1/7.scATAC processing")
    sc_anno <- c()
    for(i in unique(cell_type)){
      sc_anno <- append(sc_anno, grep(i, unique(scATAC$annotation), ignore.case = T, value = T))
    }
    scATAC <- NormalizeData(object = scATAC)
    scATAC <- FindVariableFeatures(object = scATAC)
    scATAC <- ScaleData(object = scATAC)
    scATAC <- RunPCA(object = scATAC)
    scATAC <- RunUMAP(object = scATAC, dims = 1:30)
  }else{print("scATAC-seq will be not processed.")}
  
  # processing spRNA 
  if(processing_spRNA){
    print("Predict spATAC: 2/7.spATAC processing")
    sp_anno <- c()
    for(i in unique(cell_type)){
      sp_anno <- append(sp_anno, grep(i, unique(spRNA$annotation), ignore.case = T, value = T))
    }
    spRNA <- subset(spRNA, annotation %in% sp_anno)
    spRNA <- NormalizeData(object = spRNA)
    spRNA <- FindVariableFeatures(object = spRNA)
    spRNA <- ScaleData(object = spRNA)
    spRNA <- RunPCA(object = spRNA)
    spRNA <- RunUMAP(object = spRNA, dims = 1:30)
  }else{print("spATAC data will be not processed.")}
  
  # combine datasets
  message("Combine two datasets")
  combined <- merge(spRNA, scATAC)
  
  # integrate two datasets
  print("Predict spATAC: 3/7.Integrate two datasets")
  integration.anchors <- FindIntegrationAnchors(
    object.list = list(scATAC, spRNA),
    anchor.features = intersect(rownames(spRNA),rownames(scATAC)),
    reduction = method,
    assay = c(scATAC.assay, spRNA.assay),
    dims = 1:30
  )
  integrated <- IntegrateData(anchorset = integration.anchors, dims = dims)
  integrated <- FindVariableFeatures(object = integrated)
  integrated <- ScaleData(object = integrated)
  integrated <- RunPCA(object = integrated)
  integrated <- RunUMAP(object = integrated, dims = 1:30)
  
  if(saveplot){
    p1 <- DimPlot(integrated, group.by = "dataset", pt.size = 0.8)
    ggsave(p1, filename = "./Prediction/Plot/F1.UMAP_integrated(dataset).pdf", height = 6, width = 7)
    p2  <-  DimPlot(integrated, group.by = "annotation", pt.size = 0.8)
    ggsave(p2, filename = "./Prediction/Plot/F2.UMAP_integrated(annotation).pdf", height = 6, width = 7)
  }
  
  integrated <- suppressWarnings(harmony::RunHarmony(
    object = integrated,
    group.by.vars = 'dataset',
    reduction = 'pca',
    plot_convergence = F,
    theta = 3,
    lambda = lambda,
    assay.use = "integrated",
    verbose = FALSE,
    max.iter.harmony = 20,
    project.dim = F
  ))
  
  integrated <- RunUMAP(object = integrated, reduction = "harmony", dims = 1:30)
  
  if(saveplot){
    p1 <- DimPlot(integrated, group.by = "dataset", pt.size = 0.5)
    ggsave(p1, filename = "./Prediction/Plot/F3.UMAP_integrated_harmonized(dataset).pdf", height = 6, width = 7)
    p2 <- DimPlot(integrated, group.by = "annotation", pt.size = 0.5)
    ggsave(p2, filename = "./Prediction/Plot/F4.UMAP_integrated_harmonized(annotation).pdf", height = 6, width = 7)
  }
  if(saveplot){
    message("saving data: Integrated_Harmonized_data.rds")
    saveRDS(integrated, file = "./Prediction/1.Integrated_Harmonized_data.rds")
  }
  # return(integrated)
  
  # find neighbors
  if(cell_type_pair == TRUE){
    print("Predict spATAC: 4/7.Find neighbors (cell-type pairs level)")
    final_neigbors <- list()
    for (i in cell_type) {
      print(paste0("Find neighbors:",i))
      anno <- grep(i,unique(integrated$annotation),value = T,ignore.case = T)
      temp <- subset(integrated,annotation %in% anno)
      temp <- Seurat::FindNeighbors(temp, 
                                    k.param = initial.neighbor,
                                    reduction = "harmony",
                                    dims = 1:30, 
                                    return.neighbor = T, 
                                    assay = "integrated")
      
      neigbors <- suppressWarnings(sapply(colnames(subset(temp, subset = dataset == "scATAC", invert = TRUE)), function(i){
        Seurat::TopNeighbors(temp@neighbors$integrated.nn, i, n = initial.neighbor)
      }))
      final_neigbors[[i]] <- neigbors
    }
  }else{
    final_neigbors <- list()
    for (i in cell_type) {
      print("Predict spATAC: 4/7.Find neighbors (all-data level)")
      integrated <- Seurat::FindNeighbors(integrated, 
                                    k.param = initial.neighbor,
                                    reduction = "harmony",
                                    dims = 1:30, 
                                    return.neighbor = T, 
                                    assay = "integrated")
      neigbors <- suppressWarnings(sapply(colnames(subset(integrated, subset = dataset == "scATAC", invert = TRUE)), function(i){
        Seurat::TopNeighbors(integrated@neighbors$integrated.nn, i, n = initial.neighbor)
      }))
      final_neigbors <- neigbors
    }
  }
  
  final_neigbor <- as.data.frame(final_neigbors)
  
  if(savedata){
    message("saving data: Integrated_Harmonized_data.rds")
    saveRDS(final_neigbor,file = "./Prediction/2.All_initial_neighbors.rds")
    
  }
  # filtered neighbors
  print("Predict spATAC: 5/7.Filtered neighbors")
  integrated_snATAC <- suppressWarnings(subset(integrated, subset = dataset == "scATAC"))
  integrated_spRNA <- suppressWarnings(subset(integrated, subset = dataset == "spRNA"))
  int_cell_snATAC <- colnames(integrated_snATAC)
  int_cell_spRNA <- colnames(integrated_spRNA)
  
  neigbor <- final_neigbor[1:k,]
  for (i in colnames(final_neigbor)) {
    snatac <- intersect(final_neigbor[,i], int_cell_snATAC)[1:k]
    if(length(snatac) == 0){
      print("erro: 0 snATAC found")
    }
    neigbor[,i] <- snatac
  }

  if(savedata){
    message("saving data: Filtered_scNeighbors.rds")
    saveRDS(neigbor, file = "./Prediction/3.Filtered_scNeighbors.rds")
  }
  
  # find neighbors' pairs
  print("Predict spATAC: 6/7.Find neighbors pairs")
  colnames(neigbor) <- sapply(colnames(neigbor), function(x){
    unlist(strsplit(x, "[.]"))[2]
  }) %>% unlist()
  
  data <- data.frame()
  i <-  names(neigbor)[1]
  for ( i in names(neigbor)) {
    cells <- neigbor[[i]]
    indx1 <- which(int_cell_snATAC %in% cells)
    indx2 <- which(int_cell_spRNA %in% i)
    if(length(indx1) == 0){
      indx1 <- "none"
      print("warning: 0 snATAC found")
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
    saveRDS(data, file = "./Prediction/4.Neighbors_pairs.rds")
  }
  
  # transfer data
  print("Predict spATAC: 7/7.Transfer data")
  transfer.anchors <- FindTransferAnchors(
    reference = integrated_snATAC,
    query = integrated_spRNA,
    reference.reduction = "pca",
    reduction = "rpca",
    dims = 1:30
  )
  
  transfer.anchors@anchors <- data
  predicted_spATAC <- TransferData(
    anchorset = transfer.anchors,
    refdata = GetAssayData(integrated_snATAC, assay = impute.assay, slot = "counts"),
    weight.reduction = integrated_spRNA[["harmony"]],
    dims = 1:30,
    k.weight = k.weight
  )
  
  if(savedata){
    message(paste0("saving data: ", "predicted_spATAC_assay_k", k, "_.rds"))
    saveRDS(predicted_spATAC, file = paste0("./Prediction/5.predicted_spATAC_assay_k", k, ".rds"))
    print("Finished prediction!!!")
  }
  return(predicted_spATAC)
}


