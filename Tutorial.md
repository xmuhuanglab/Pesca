
# Tutorial for Pesca

## Introduction
Pesca is a computational method to predict and enhance spatial profiling of chromatin accessibility using scATAC-seq data. The accuracy and generalizability of the Pesca-enhanced and -predicted data are validated by existing various ground-truth data including spatial co-profiling data, VISTA results, and epigenomic MEIFISH data.

![11](https://github.com/xmuhuanglab/Pesca/assets/95668602/849cc8cf-124c-4e96-a7ee-d80a54791eb8)


## Input data
#### Input data for enhancement
When enhancing the spatial ATAC profile, the peak set in scATA-seq data and spATAC data should be identical. Otherwise, please provide the sequencing fragment files of these datasets to generate the common peaks.<br>
(1) spATAC: peak-by-spot matrix or fragment files and spatial coordinates for each spot; <br>
(2) scATAC-seq: peak-by-cell matrix or fragment files; <br>
(3) cell type annotations for spots and cells (optional); <br>

#### Input data for prediction
When predicting spatial ATAC profile, the fragment files of scATA-seq data should be provided to generate the gene activity score using Signac. The gene activity will be regarded as the bridge to integrate the spatial transcriptome (spRNA) data and scATAC-seq data.<br>
(1) spRNA: gene-by-spot matrix and spatial coordinates for each spot; <br>
(2) scATAC-seq: peak-by-cell matrix and fragment files; <br>
(3) cell type annotations for spots and cells (optional); <br>

As an example, here we only enhanced the spatial ATAC profiles of forebrain and midbrain from E13 mouse, which are from research 10.1038/s41586-023-05795-1. The corresponding scATAC-seq data of E13.5 mouse used for enhancement are from 10.1016/j.celrep.2023.112210.

## Data Processing 
#### Data Processing (input data are fragment files)
If the input data of scATAC-seq and spATAC are fragment files, we should generate the reliable peak set first using scATAC-seq data and create Seurat objects for both datasets separately.

```r
setwd("Workdir")
rm(list=ls())
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(patchwork)
library(GenomicRanges)
library(ggplot2)
library(Matrix)
library(data.table)
library(scran)
library(presto)
library(dplyr)
library(ggplot2)

source("./processing_scATAC_spATAC.R")
source("./predict_spATAC.R")
source("./enhance_spATAC.R")
source("./visualization.R")
```
#### Step1. Call peaks using fragment files of scATAC-seq data
```r
fragpath_scatac <- dir('/path/fragment_data/', "*gz$", full.names = TRUE)
macs2_path <- "/path/macs2"
cell_type <- c("forebrain", "midbrain", "hindbrain", "forelimb", "eye")

peaks <- Callpeak_scATAC(cell_type, # the cell types related to the separate fragment filenames used to call peak
                         fragpath_scatac, # the path to stored fragment file (tabix)
                         macs2_path # the path to stored MACS2 software
)
```

#### Step2. Generate peak-by-cell matrix and create Seurat object
```r
annotation <- readRDS(("./data/Ref1.annotation.rds")) # genome annotation
scATAC <- scATAC_processing(fragpath = fragpath_scatac, # path storing fragment file (tabix)
                            peaks = peaks, # merged peaks in scATAC-seq data
                            cell_type = cell_type, # vector containing the organ names
                            sc_mat = NULL, # a peak-by-cell matrix, only needed when fragment data are not available
                            mincell = 100,
                            minfeature = 1000,
                            cell_anno = NULL, # a vector containing the cell identify
                            genome = "mm10", # mm10 or hg19 or hg38
                            sep = c("-", "-"), 
                            genome_anno = NULL, # genome annotation
                            savedata = TRUE,
                            gene.activity = TRUE)
```

#### Step3. Generate peak-by-spot matrix and create Seurat Object

```r
fragpath_spatac <- dir('/path/fragment_data/', "*gz$", "*gz$", full.names = TRUE)
spot_anno <- readRDS("./data/Ref2.spATAC_spots_annotations.rds")

spATAC_Object <- spATAC_processing(fragpath_spatac = NULL, # path storing spATAC fragment file
                                   peaks = peaks, # merged peaks set generated from corresponding scATAC-seq data
                                   spotIDs = NULL, # spots used to generate matrix
                                   spot_anno = spot_anno, # a vector containing spot identity
                                   sp_mat = NULL, # a peak-by-spot matrix, only needed when fragment data are not available
                                   genome = "mm10", # genome annotation 
                                   minspot = 5, # the minimum number of spots which peaks are accessible in 
                                   minfeature = 10, # the minimum number of peaks that a spot contains
                                   savedata = TRUE)
```

## Pesca enhancement 
Pesca enhances spatial ATAC profile using scATAC-seq data. Firstly, scATAC-seq and spATA data will be integrated using LSI or CCA reduction. Since CCA-based integration may also lead to overcorrection, so we recommend “rlsi” as the first choice. Then Harmony will be used for second-round integration. In the common space, Pesca will identify the k nearest single-cell neighbors for each spatial spot, the number of k depends on the spot size. Notably, if the annotations are given, the neighbor identification will be restricted to the same cell type, assuring the origin similarity of the spot and their neighbors. After that, the chromatin accessibility of single-cell neighbors will be used to enhance the spatial profile where the distances of spatial neighbors and the anchor score will be used to calculate the anchor weight.


```r
scATAC <- readRDS("./Processing_scATAC/5.scATAC_object.rds")
spATAC <- readRDS("./Processing_spATAC/1.spATAC_object.rds")

cell_type <- c("forebrain", "midbrain")
enhanced_spATAC <- Enhancement(scATAC = scATAC, # Seurat object of scATAC-seq
                               spATAC = spATAC , # Seurat object of spATAC data
                               processing_scATAC = TRUE, # whether scATAC-seq data will be processed
                               processing_spATAC = TRUE, # whether spATAC data will be processed
                               method = "rlsi", # integration method, lsi or cca
                               spATAC.assay = "peaks", # the assay data used for integration in spATAC data
                               scATAC.assay = "peaks",  # the assay data used for integration in scATAC-seq data
                               cell_type = cell_type, # cell type data used for enhancement
                               initial.neighbor = 200, # the initial neighbors identified for spots, default 100
                               dims = 2:30, # the dimensions when performing the integration
                               k = 10, # the number of final single-cell neighbors for each spot
                               cell_type_pair = TRUE, # whether the neighbor identification is restricted to cell types. Annotation should be represented in scATAC and spATAC data when it is true.
                               k.weight = 20, # k.weight when transferring data
                               lambda = 0.9, # lambda parameter when performing Harmony to integrate
                               savedata = TRUE,
                               saveplot = TRUE)

spATAC[["enhanced"]] <- enhanced_spATAC

# Add spatial information before plotting
coordinate <- readRDS("./data/spRNA_coordinate.rds")
spATAC$x <- coordinate$x
spATAC$y <- coordinate$y

plot0 <- mySpatialDimPlot(metadata = spATAC@meta.data, 
                          x = "x", # x coordinate
                          y = "y", # y coordinate
                          color = "annotation", # the variable visualized
                          pointsize = 2.5, 
                          colors = NULL 
);print(plot0)

features <- rownames(spATAC)[1:10]
plot1 <- mySpatialFeaturePlot(object = spATAC, # Seurat object
                              feature = features[1], # the features used to plot
                              x = "x", # x coordinate
                              y = "y", # y coordinate
                              pointsize = 2.5, # spot size
                              scale = FALSE, # whether scale the data
                              ncol = 2, # the number of columns for plots when the number of features is more than two
                              assay = "enhanced" , # the assay data used to plot 
                              slot = "data" # the slot in the assay 
);print(plot1)


plot2 <- mySpatialFeaturePlot(object = spATAC, # Seurat object
                              feature = features[1], # the features used to plot
                              x = "x", # x coordinate
                              y = "y", # y coordinate
                              pointsize = 2.5, # spot size
                              scale = FALSE, # whether scale the data
                              ncol = 2, # the number of columns for plots when the number of features is more than two
                              assay = "peaks" , # the assay data used to plot 
                              slot = "data" # the slot in the assay 
);print(plot2)

```

# Pesca prediction
Pesca predicts spatial ATAC profile using scATAC-seq data and spRNA data. Firstly, scATAC-seq (gene activity) and spRNA data will be integrated using RPCA or CCA reduction. Since CCA-based integration may also lead to overcorrection, so we recommend “RPCA” as the first choice. Then Harmony will be used for second-round integration.In the common space, Pesca will identify the k nearest single-cell neighbors for each spatial spot, the number of k depends on the spot size. Notably, if the annotations are given, the neighbor identification will be restricted to the same cell type, assuring the origin similarity of the spot and their neighbors. After that, the chromatin accessibility of single-cell neighbors will be imputed to the corresponding spot, where the distances of spatial neighbors and the anchor score will be used to calculate the anchor weight.

```r
spRNA <- readRDS("/path/spRNA.rds")
scATAC <- readRDS("./Processing_scATAC/5.scATAC_object.rds")

predicted_spATAC <- Prediction(scATAC = scATAC, # scATAC Seurat object
                               spRNA = spRNA, # spRNA Seurat object
                               processing_scATAC = TRUE, # whether scATAC-seq data will be processed
                               processing_spRNA = TRUE, # whether spRNA data will be processed
                               method = "rpca", # integration method, rpca or cca
                               k.weight = 50, # k.weight in integration step
                               spRNA.assay = "spatial", # the assay data used for integration in spRNA data
                               scATAC.assay = "gene_activity", # the assay data used for integration in scATAC-seq data
                               impute.assay = "peaks", # the assay data imputed to spatial spot
                               cell_type = cell_type, # cell types data used for enhancement
                               initial.neighbor = 100, # the initial neighbors identified for spots, default 100
                               dims = 1:30, # the dimensions when performing the integration
                               k = 10, # the number of final single-cell neighbors for each spot
                               cell_type_pair = TRUE, # whether the neighbor identification is restricted to cell types. Annotation should be represented in scATAC and spATAC data when it is true.
                               lambda = 0.9, # lambda parameter when performing Harmony to integrate
                               savedata = TRUE,
                               saveplot = TRUE
)

```



