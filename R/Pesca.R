######################################                                                         ######################################
######################################   Generate enhanced and predicted spatial ATAC Assay    ######################################                                                       

#Call peaks for single-cell ATAC data
#######################################################################
callpeak_scATAC=function(organs, # Organ related to the separate fragment filename used to call peak
                         fragpath, # Path to stored fragment file (tabix)
                         macs2_path # Path to stored macs2 software
){
  
  #### call peak for scATAC data separately and filtered
  allpeaks=list()
  for(i in organs){
    print("Callpeak: 1.call peak for scATAC-seq data")
    path=grep(i,fragpath,value = T,ignore.case = T)
    peaks_all <- CallPeaks(path, macs2.path =macs2_path)
    peaks_all <- keepStandardChromosomes(peaks_all, pruning.mode = "coarse")
    peaks_all <- subsetByOverlaps(x = peaks_all, ranges =  blacklist_mm10, invert = TRUE)
    
    peaks=as.data.frame(peaks_all)
    peaks=peaks[which(peaks$neg_log10qvalue_summit>2 & peaks$fold_change>2),]
    
    peaks=makeGRangesFromDataFrame(peaks)
    allpeaks=append(allpeaks,peaks)
    allpeaks[[i]]=peaks
    
  }
  saveRDS(allpeaks,file = "Step0_1.Allpeaks.rds")
  print("Finished calling peak for scATAC-seq data")
  
  ###### merge all peaks
  print("Callpeak: 2.merge peaks for scATAC-seq data")
  
  extracted_data <- list()
  for (i in seq_along(allpeaks)) {
    extracted_data[[i]] <- allpeaks[[i]]
  }
  combined_iranges <- do.call(c, extracted_data)
  
  peaks <- reduce(x = combined_iranges,drop.empty.ranges = T)
  peaks_width <- width(peaks)
  
  saveRDS(peaks,file = paste0("Step0_2.MergedPeaks.rds"))
  print("Finished merging peaks: total ", length(peaks)," peaks")
  return(peaks)
  
  
}

######### 2.Generate Matrix and create SeuratObject for scATAC-seq data
scATAC_processing=function( fragpath, # Path to stored fragment file (tabix)
                            peaks, # Merged peaked in scATAC
                            organs, # vector contain the organ names
                            mincell=100,
                            minfeature=1000,
                            genome="mm10", # mm10,hg19
                            annotation # 
                            ){
  allmatrix=list()
  cellanno=list()
  fragObjects=list()
  
  for (f in organs) {
    fragment=data.table::fread(grep(f,fragpath,value = T,ignore.case = T))
    cells=unique(as.vector(fragment$V4))
    temp=CreateFragmentObject(path = grep(f,fragpath,value = T,ignore.case = T),cells = cells,validate.fragments = F)
    mat=FeatureMatrix(fragments =temp,features = peaks)
    cellanno[[f]]=colnames(mat)
    allmatrix[[f]]=mat
    fragObjects[[f]]=temp
    print("Finished matrix generation:",f)
  }
  rm(fragment)
  sc_mat=do.call(cbind,allmatrix)
  sc_mat=as(sc_mat, "sparseMatrix")
  
  saveRDS(sc_mat,file = paste0("Step0_3.scMatrix.rds"))
  saveRDS(cellanno,file = paste0("Step0_4.Cell_annotation.rds"))
  
  ### create Seurat object
  chrom_assay <- CreateChromatinAssay(
    counts = sc_mat,
    sep = c("-", "-"),
    genome =genome,
    fragments = fragObjects,
    min.cells = mincell,
    min.features = minfeature,
  )
  
  sc_atac <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks"
  )
  
  print("Number of total cells:",ncol(sc_atac))
  print("Number of total peaks:",nrow(sc_atac))
  
  if(genome=="mm10"){
    annotations <- ann_data$v79
  }else if(
    genome=="hg19"
  ){
    annotations <- ann_data$v75
  }else{
    annotations <- ann_data$v86
  }
  seqlevelsStyle(annotations) <- 'UCSC'
  Annotation(sc_atac) <- annotations
  
  ### add annotation information
  sc_atac$annotation=NA
  sc_atac$cells=rownames(sc_atac@meta.data)
  for (o in names(cellanno)) {
    sc_atac@meta.data[which(sc_atac@meta.data$cells %in% cellanno[[o]]),"annotation"]=o
  }
  
  ### add gene activity
  ga <- GeneActivity(sc_atac,extend.upstream = 2000,extend.downstream = 2000)
  sc_atac[['gene_activity']] <- CreateAssayObject(counts = ga[,intersect(colnames(ga),colnames(sc_atac))])
  sc_atac <- NormalizeData(
    object = sc_atac,
    assay = 'gene_activity',
    normalization.method = 'LogNormalize',
    scale.factor = median(sc_atac$nCount_gene_activity)
  )
  saveRDS(sc_atac,file = paste0("Step0_5.SeuratObject_scATAC.rds"))
  return(sc_atac)
}


#Call peaks for spatial ATAC data
#######################################################################
spATAC_processing=function(fragpath_spatac, # The path of spatial ATAC fragment file
                           peaks, # Merged peaks set generated from corresponding scATAC-seq data
                           spotIDs=NULL,# Spots used to generate matrix
                           anno, # List containing spots annotation
                           genome="mm10",# Genome annotation 
                           minspot=5,# The minimum number of spots which peaks are accessible in 
                           minfeature=10# The minimum number of peaks that a spot contains
                           ){
  
  if(is.null(spotIDs)){
    fragment=data.table::fread(fragpath_spatac)
    spotIDs=unique(fragment$V4)
    print("1/3.Using all spots data to call peaks")
  }
  fragObject=CreateFragmentObject(path = fragpath_spatac,cells = spotIDs)
  spatac_mat=FeatureMatrix(fragments =fragObject,features = peaks)
  spatac_mat=as(spatac_mat, "sparseMatrix")
  
  print("2/3.Creatate SeuratObject")
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
  
  print("Number of total cells:",ncol(sp_atac),"in spatial ATAC data")
  print("Number of total peaks:",nrow(sp_atac),"in spatial ATAC data")
  
  if(genome=="mm10"){
    annotations <- ann_data$v79
  }else if(
    genome=="hg19"
  ){
    annotations <- ann_data$v75
  }else{
    annotations <- ann_data$v86
  }
  seqlevelsStyle(annotations) <- 'UCSC'
  Annotation(sp_atac) <- annotations
  
  ### add annotation information
  sp_atac$annotation=NA
  sp_atac$cells=rownames(sp_atac@meta.data)
  for (o in names(anno)) {
    sp_atac@meta.data[which(sp_atac@meta.data$cells %in% anno[[o]]),"annotation"]=o
  }
  
  print("3/3.Calculate spatial ATAC gene activity")
  ga <- GeneActivity(sp_atac,extend.upstream = 2000,extend.downstream = 2000)
  sp_atac[['gene_activity']] <- CreateAssayObject(counts = ga[,intersect(colnames(ga),colnames(sp_atac))])
  sp_atac <- NormalizeData(
    object = sp_atac,
    assay = 'gene_activity',
    normalization.method = 'LogNormalize',
    scale.factor = median(sp_atac$nCount_gene_activity)
  )
  saveRDS(sp_atac,file = paste0("Step1_1.SeuratObject_spATAC.rds"))
  return(sp_atac)
}



# Pesca Enhancement
#######################################################################
Enhancement_processing=function(scATAC, # seuratObject
                    spATAC, # SeuratObject
                    organs, # Organs/tissues data for enhancement
                    k.weight = 20, # k.weight when integrating two datasets
                    lambda = 0.9 # Lambda parameter when using Harmony
                    ){
  if(dir.exists("./Enhancemnet/Plot")){
    print("Directory exists")
    
  }else{
    dir.create("./Enhancemnet/Plot",recursive=T)
  }
  spATAC$annotation=paste0("spatial_",spATAC$annotation)
  scATAC$annotation=paste0("sc_",scATAC$annotation)
  spATAC$dataset <- "spATAC"
  scATAC$dataset <- "scATAC"
  DefaultAssay(scATAC)="peaks"
  DefaultAssay(spATAC)="peaks"  
  
  ######## scATAC processing
  print("1/4.scATAC processing")
  sc_anno=c()
  for(i in unique(organs)){
    sc_anno=append(sc_anno,grep(i,unique(scATAC$annotation),ignore.case = T,value = T))
  }
  scATAC=subset(scATAC, annotation%in% sc_anno)
  scATAC <- RunTFIDF(scATAC)
  scATAC <- FindTopFeatures(scATAC, min.cutoff = 'q5')
  scATAC <- RunSVD(scATAC)
  scATAC <- RunUMAP(object = scATAC, reduction = 'lsi', dims = 2:20)

  ######## spATAC processing
  print("2/4.spATAC processing")
  sp_anno=c()
  for(i in unique(organs)){
    sp_anno=append(sp_anno,grep(i,unique(spATAC$annotation),ignore.case = T,value = T))
  }
  spATAC=subset(spATAC, annotation%in% sp_anno)
  spATAC <- RunTFIDF(spATAC)
  spATAC <- FindTopFeatures(spATAC, min.cutoff = 'q5')
  spATAC <- RunSVD(spATAC)
  spATAC <- RunUMAP(object = spATAC, reduction = 'lsi', dims = 2:30)
  
  ######## Combine datasets
  print("3/3.Combine two datasets")
  combined <- merge(spATAC, scATAC)
 
  ######## Integrate datasets
  print("3/4.Integrate two datasets")
  integration.anchors <- FindIntegrationAnchors(
    object.list = list(scATAC, spATAC),
    anchor.features = intersect(rownames(spATAC),rownames(scATAC)),
    reduction = "rlsi",
    assay = c("peaks","peaks"),
    dims = 2:20
  )
  
  integrated <- IntegrateData(anchorset = integration.anchors,k.weight = k.weight)
  integrated <- RunSVD(integrated)
  integrated <- RunUMAP(object = integrated, reduction = 'lsi', dims = 2:20)

  p=DimPlot(object = integrated, label = TRUE,group.by = "annotation") + NoLegend()
  ggsave(p,filename = "./Enhancement/Plot/F1.UMAP_integrated.pdf",height = 6,width = 7)
  
  ######## Harmony 
  print("4/4.Harmonize two datasets")
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

  p1=DimPlot(integrated,group.by ="dataset",pt.size = 0.5)
  ggsave(p1,filename = "./Enhancement/Plot/F2.UMAP_integrated_harmonized(dataset).pdf",height = 6,width = 7)
  p2 = DimPlot(integrated, group.by = "annotation",pt.size = 0.5)
  ggsave(p2,filename = "./Enhancement/Plot/F2.UMAP_integrated_harmonized(annotation).pdf",height = 6,width = 7)
  
  # saveRDS(integrated,file="./Enhancement/Step2_1.Integrated_Harmonized_data.rds")
  return(integrated)
}



#Pesca prediction processing
#######################################################################
Enhancemment=function( IntegratedData, # Integrated dataset
                                 k=10, # the number of single-cell neighbors for each spot
                                 k.neighbor = 200, # Initial neighbors
                                 dims = 1:30,
                                 organs
                                 ){
  
  ######## Find neighbors 
  print("1/4.Find neighbors")
  final_neigbors=list()
  for (i in organs) {
    print(paste0("Find neighbors:",i))
    anno1=grep(i,unique(integrated$annotation),value = T,ignore.case = T)
    temp=subset(integrated,annotation %in% anno1)
    temp = Seurat::FindNeighbors(temp, 
                                 k.param = k.neighbor,
                                 reduction="harmony",
                                 dims=dims, 
                                 return.neighbor = T, 
                                 assay = "integrated")
    neigbors = suppressWarnings(sapply(colnames(subset(temp, subset = dataset == "scATAC", invert = TRUE)), function(i){
      Seurat::TopNeighbors(temp@neighbors$integrated.nn, i, n = k.neighbor)
    }))
    final_neigbors[[i]]=neigbors
  }
  
  final_neigbor=as.data.frame(final_neigbors)
  saveRDS(final_neigbor,file = "./Enhancement/Step3_1.All_Neighbors.rds")
  
  ######## Filtered neighbors
  print("2/4.Filtered neighbors")
  integrated_scATAC = suppressWarnings(subset(integrated, subset = dataset == "scATAC"))
  integrated_spATAC = suppressWarnings(subset(integrated, subset = dataset == "spATAC"))
  int_cell_scATAC=colnames(integrated_scATAC)
  int_cell_spATAC=colnames(integrated_spATAC)
  
  neigbor=final_neigbor[1:k,]
  for (i in colnames(final_neigbor)) {
    scATAC=intersect(final_neigbor[,i],int_cell_scATAC)[1:k]
    if(length(scATAC)==0){
      print("erro: 0 scATAC found")
    }
    neigbor[,i]=scATAC
  }
  saveRDS(neigbor,file = "./Enhancement/Step3_2.Filtered_scNeighbors.rds")

  ######## Find neighbors pairs
  print("3/4.Find neighbors pairs")
  colnames(neigbor)=sapply(colnames(neigbor), function(x){
    unlist(strsplit(x,"[.]"))[2]
  }) %>% unlist()
  
  data <- data.frame()
  for ( i in names(neigbor)) {
    cells = neigbor[[i]]
    indx1 = which(int_cell_scATAC %in% cells)
    indx2 = which(int_cell_spATAC %in% paste0(i,"-1"))
    if(length(indx1)==0){
      indx1="mei"
      print("erro: 0 scATAC found")
    }
    tmp = data.frame(cell1 =indx1, cell2 = indx2, score = 1) 
    data = rbind(data,tmp)
  }
  saveRDS(data,file = "./Enhancement/Step3_3.Neighbors_Pairs.rds")
  
  ######## Transfer data
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
  saveRDS(enhanced_spATAC,file = paste0("./Enhancement/Step3_4.Assay_enhanced_spATAC_k",k,"_.rds"))
  return(enhanced_spATAC)
}



#Pesca Prediction
#######################################################################
Prediction_processing=function(scATAC,
                               spRNA,
                               k.weight = 50, # k.weight in integration step
                               default.assay="spatial" # the count assay of spRNA 
                               ){
  if(dir.exists("./Prediction/Plot")){
    print("Directory exists")
  }else{
      dir.create("./Prediction/Plot/",recursive=T)
  }
  spRNA$dataset <- "spRNA"
  scATAC$dataset <- "scATAC"
  DefaultAssay(scATAC)="gene_activity"
  DefaultAssay(spRNA)=default.assay  
  
  ######## scATAC processing
  print("2/4.spATAC processing")
  sc_anno=c()
  for(i in unique(organs)){
    sc_anno=append(sc_anno,grep(i,unique(scATAC$annotation),ignore.case = T,value = T))
  }
  scATAC <- NormalizeData(object = scATAC)
  scATAC <- FindVariableFeatures(object = scATAC)
  scATAC <- ScaleData(object = scATAC)
  scATAC <- RunPCA(object = scATAC)
  scATAC <- RunUMAP(object = scATAC,dims = 1:30)
  
  ########  spRNA processing
  print("3/4.spATAC processing")
  sp_anno=c()
  for(i in unique(organs)){
    sp_anno=append(sp_anno,grep(i,unique(spRNA$annotation),ignore.case = T,value = T))
  }
  spRNA=subset(spRNA, annotation%in% sp_anno)
  spRNA <- NormalizeData(object = spRNA)
  spRNA <- FindVariableFeatures(object = spRNA)
  spRNA <- ScaleData(object = spRNA)
  spRNA <- RunPCA(object = spRNA)
  spRNA <- RunUMAP(object = spRNA,dims = 1:30)
  
  ######## Combine datasets
  message("Combine two datasets")
  combined <- merge(spRNA, scATAC)

  ######## Integrate two datasets
  message("Integrate two datasets")
  integration.anchors <- FindIntegrationAnchors(
    object.list = list(scATAC, spRNA),
    anchor.features = intersect(rownames(spRNA),rownames(scATAC)),
    reduction = "cca",
    # reduction="rpca",
    assay = c("gene_activity","spatial"),
    dims = 2:30
  )
  integrated <- IntegrateData(anchorset = integration.anchors,k.weight = k.weight)
  integrated <- FindVariableFeatures(object = integrated)
  integrated <- ScaleData(object = integrated)
  integrated <- RunPCA(object = integrated)
  integrated <- RunUMAP(object = integrated,dims = 1:30)
  
  integrated <- suppressWarnings( harmony::RunHarmony(
    object = integrated,
    group.by.vars = 'dataset',
    reduction = 'pca',
    plot_convergence = F,
    theta = 3,
    lambda = 0.5,
    assay.use = "integrated",
    verbose = FALSE,
    max.iter.harmony = 20,
    project.dim = F
  ))
  
  p=DimPlot(object = integrated, label = TRUE,group.by = "annotation") + NoLegend()
  ggsave(p,filename = "./Prediction/Plot/F1.UMAP_integrated.pdf",height = 6,width = 7)
  
  p1=DimPlot(integrated,group.by ="dataset",pt.size = 1.0)
  ggsave(p1,filename = "./Prediction/Plot/F2.UMAP_integrated_harmonized(dataset).pdf",height = 6,width = 7)
  p2 = DimPlot(integrated, group.by = "annotation",pt.size = 0.8)
  ggsave(p2,filename = "./Prediction/Plot/F2.UMAP_integrated_harmonized(annotation).pdf",height = 6,width = 7)
  
  # saveRDS(integrated,file="./Prediction/Step2_1.Integrated_Harmonized_data.rds")
  return(integrated)
}


Prediction=function(k.neighbor,
                     dims=1:30,
                     organs,
                     k # the number of single-cell neighbors for each spot
                    ){
  
  ######## Find neighbors
  print("1/4.Find neighbors")
  final_neigbors=list()
  for (i in organs) {
    print(paste0("Find neighbors:",i))
    anno=grep(i,unique(integrated$annotation),value = T,ignore.case = T)
    temp=subset(integrated,annotation %in% anno)
    temp = Seurat::FindNeighbors(temp, 
                                 k.param = k.neighbor,
                                 reduction="harmony",
                                 dims=dims, 
                                 return.neighbor = T, 
                                 assay = "integrated")
    
    neigbors = suppressWarnings(sapply(colnames(subset(temp, subset = dataset == "scATAC", invert = TRUE)), function(i){
      Seurat::TopNeighbors(temp@neighbors$integrated.nn, i, n = k.neighbor)
    }))
    final_neigbors[[i]]=neigbors
  }
  final_neigbor=as.data.frame(final_neigbors)
  saveRDS(final_neigbor,file = "./Prediction/Step3_1.All_Neighbors.rds")
  
  ######## Filtered neighbors
  print("2/4.Filtered neighbors")
  integrated_snATAC = suppressWarnings(subset(integrated, subset = dataset == "scATAC"))
  integrated_spRNA = suppressWarnings(subset(integrated, subset = dataset == "spRNA"))
  int_cell_snATAC=colnames(integrated_snATAC)
  int_cell_spRNA=colnames(integrated_spRNA)
  
  neigbor=final_neigbor[1:k,]
  for (i in colnames(final_neigbor)) {
    snatac=intersect(final_neigbor[,i],int_cell_snATAC)[1:k]
    if(length(snatac)==0){
      print("erro: 0 snATAC found")
    }
    neigbor[,i]=snatac
  }
  saveRDS(neigbor,file = "./Prediction/Step3_2.Filtered_scNeighbors.rds")
  
  ######## Find  neighbors' pairs
  print("3/4.Find neighbors pairs")
  colnames(neigbor)=sapply(colnames(neigbor), function(x){
    unlist(strsplit(x,"[.]"))[2]
  }) %>% unlist()
  
  data <- data.frame()
  i = names(neigbor)[1]
  for ( i in names(neigbor)) {
    cells = neigbor[[i]]
    indx1 = which(int_cell_snATAC %in% cells)
    indx2 = which(int_cell_spRNA %in% i)
    if(length(indx1)==0){
      indx1="mei"
      print("erro: 0 snATAC found")
    }
    
    tmp = data.frame(cell1 =indx1, cell2 = indx2, score = 1) 
    data = rbind(data,tmp)
    
  }
  saveRDS(data,file = "./Prediction/Step3_3.Neighbors_Pairs.rds")
  
  ######## Transfer data
  print("4/4.Transfer data")
  transfer.anchors <- FindTransferAnchors(
    reference = integrated_snATAC,
    query = integrated_spRNA,
    reference.reduction = "pca",
    reduction = "rpca",
    dims = 2:30
  )

  transfer.anchors@anchors <- data
  predicted_spATAC <- TransferData(
    anchorset = transfer.anchors,
    refdata = GetAssayData(integrated_snATAC, assay = "peaks", slot = "counts"),
    weight.reduction = integrated_spRNA[["harmony"]],
    dims = 2:30
  )
  saveRDS(predicted_spATAC,file = paste0("./Prediction/Step3_4.Assay_predicted_spATAC_k",k,"_.rds"))
}
