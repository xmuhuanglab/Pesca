# This script is to enhance the spatial ATAC data using scATAC-seq and spatial ATAC data. The neighbor identifying process in enhancement could be performed on the cell-type pairs level (more accurate) or the all-data level

# enhance the spatial ATAC data using scATAC-seq data and spatial ATAC data
#######################################################################
Enhancement <- function(scATAC, # scATAC seurat object
                        spATAC, # spATAC seurat object
                        organs, # organs/tissues data for enhancement
                        k.neighbor, # the initial neighbors identified for spots
                        dims = 2:30, # the dimensions when perform the integration
                        k, # the number of final single-cell neighbors for each spot
                        cell_type_pair = TRUE, 
                        k.weight = 20, # k.weight when integrating two datasets
                        lambda = 0.9 # lambda parameter when performing Harmony to integrate
){
  if(dir.exists("./Enhancemnet/Plot")){
    print("Directory exists")
  }else{
    dir.create("./Enhancemnet/Plot", recursive = T)
  }
  
  if(cell_type_pair == TRUE){
    print("You choose to find the neighbors for spots on cell-type pairs level")
  }else{
    print("You choose to find the neighbors for spots on all-data pairs level")
  }
  
  spATAC$annotation <- paste0("spatial_", spATAC$annotation)
  scATAC$annotation <- paste0("sc_", scATAC$annotation)
  spATAC$dataset <- "spATAC"
  scATAC$dataset <- "scATAC"
  DefaultAssay(scATAC) <- "peaks"
  DefaultAssay(spATAC)<- "peaks"  
  
  # processing scATAC-seq data 
  print("Enhance spATAC: 1/4.scATAC processing")
  sc_anno <- c()
  for(i in unique(organs)){
    sc_anno <- append(sc_anno, grep(i, unique(scATAC$annotation), ignore.case = T, value = T))
  }
  scATAC <- subset(scATAC, annotation%in% sc_anno)
  scATAC <- RunTFIDF(scATAC)
  scATAC <- FindTopFeatures(scATAC, min.cutoff = 'q5')
  scATAC <- RunSVD(scATAC)
  scATAC <- RunUMAP(object = scATAC, reduction = 'lsi', dims = 2:20)
  
  # processing spATAC data
  print("Enhance spATAC: 2/4.spATAC processing")
  sp_anno <- c()
  for(i in unique(organs)){
    sp_anno <- append(sp_anno, grep(i, unique(spATAC$annotation), ignore.case = T, value = T))
  }
  spATAC <- subset(spATAC, annotation%in% sp_anno)
  spATAC <- RunTFIDF(spATAC)
  spATAC <- FindTopFeatures(spATAC, min.cutoff = 'q5')
  spATAC <- RunSVD(spATAC)
  spATAC <- RunUMAP(object = spATAC, reduction = 'lsi', dims = 2:30)
  
  # combine datasets
  print("3/3.Combine two datasets")
  combined <- merge(spATAC, scATAC)
  
  # integrate two datasets
  print("Enhance spATAC: 3/4.Integrate two datasets")
  integration.anchors <- FindIntegrationAnchors(
    object.list = list(scATAC, spATAC),
    anchor.features = intersect(rownames(spATAC), rownames(scATAC)),
    reduction = "rlsi",
    assay = c("peaks", "peaks"),
    dims = 2:20
  )
  
  integrated <- IntegrateData(anchorset = integration.anchors, k.weight = k.weight)
  integrated <- RunSVD(integrated)
  integrated <- RunUMAP(object = integrated, reduction = 'lsi', dims = 2:20)
  
  p <- DimPlot(object = integrated, label = TRUE, group.by = "annotation") + NoLegend()
  ggsave(p, filename = "./Enhancement/Plot/F1.UMAP_integrated.pdf", height = 6, width = 7)
  
  # run Harmony to integrate again
  print("Enhance spATAC: 4/4.Harmonize two datasets")
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
  
  integrated <- RunUMAP(integrated, dims = 2:20, reduction = 'harmony')
  
  p1 <- DimPlot(integrated, group.by ="dataset", pt.size = 0.5)
  ggsave(p1, filename = "./Enhancement/Plot/F2.UMAP_integrated_harmonized(dataset).pdf", height = 6, width = 7)
  p2 <- DimPlot(integrated, group.by = "annotation", pt.size = 0.5)
  ggsave(p2, filename = "./Enhancement/Plot/F2.UMAP_integrated_harmonized(annotation).pdf", height = 6, width = 7)
  
  saveRDS(integrated, file="./Enhancement/1.Integrated_Harmonized_data.rds")
  
  # find neighbors
  if(cell_type_pair == TRUE){
    print("Enhance spATAC: 1/4.Find neighbors (cell-type pairs level)")
    final_neigbors <- list()
    
    for (i in organs) {
      print(paste0("Find neighbors:", i))
      anno1 <- grep(i, unique(integrated$annotation), value = T, ignore.case = T)
      temp <- subset(integrated, annotation %in% anno1)
      temp <- Seurat::FindNeighbors(temp, 
                                    k.param = k.neighbor,
                                    reduction = "harmony",
                                    dims = dims, 
                                    return.neighbor = T, 
                                    assay = "integrated")
      neigbors <- suppressWarnings(sapply(colnames(subset(temp, subset = dataset == "scATAC", invert = TRUE)), function(i){
        Seurat::TopNeighbors(temp@neighbors$integrated.nn, i, n = k.neighbor)
      }))
      final_neigbors[[i]] <- neigbors
      
    }}else{
      print("Enhance spATAC: 1/4.Find neighbors (all-data level)")
      integrated <- Seurat::FindNeighbors(integrated, 
                                          k.param = k.neighbor,
                                          reduction = "harmony",
                                          dims = dims, 
                                          return.neighbor = T, 
                                          assay = "integrated")
      neigbors <- suppressWarnings(sapply(colnames(subset(integrated, subset = dataset == "scATAC", invert = TRUE)), function(i){
        Seurat::TopNeighbors(integrated@neighbors$integrated.nn, i, n = k.neighbor)
      }))
      final_neigbors <- neigbors
    }
  
  final_neigbor <- as.data.frame(final_neigbors)
  saveRDS(final_neigbor,file = "./Enhancement/2.All_initial_neighbors.rds") 
  
  # filtered neighbors
  print("2/4.Filtered neighbors")
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
  saveRDS(neigbor, file = "./Enhancement/3.Filtered_scNeighbors.rds")
  
  # find neighbors pairs
  print("3/4.Find neighbors pairs")
  colnames(neigbor) <- sapply(colnames(neigbor), function(x){
    unlist(strsplit(x, "[.]"))[2]
  }) %>% unlist()
  
  data <- data.frame()
  for ( i in names(neigbor)) {
    cells <- neigbor[[i]]
    indx1 <- which(int_cell_scATAC %in% cells)
    indx2 <- which(int_cell_spATAC %in% paste0(i,"-1"))
    if(length(indx1) == 0){
      indx1 <- "mei"
      print("erro: 0 scATAC found")
    }
    tmp <- data.frame(cell1 =indx1, cell2 = indx2, score = 1) 
    data <- rbind(data, tmp)
  }
  saveRDS(data, file = "./Enhancement/4.Neighbors_pairs.rds")
  
  # transfer data
  print("4/4.Transfer data")
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
    refdata = GetAssayData(integrated_scATAC, assay = "peaks", slot = "counts"),
    weight.reduction = integrated_spATAC[["harmony"]],
    dims = 2:30
  )
  saveRDS(enhanced_spATAC, file = paste0("./Enhancement/5.Enhanced_spATAC_assay_k", k, "_.rds"))
  print("Finished enhancement!!!")
}





    

  

  

  
  

  


  
  
  
  
  
  
  
  
