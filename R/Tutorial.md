
## Tutorial for Pesca
### Data for enhancement
(1) fragment files from spatial-ATAC-seq data; <br>
(2) fragment files from scATAC-seq data; # should contain corresponding tissues/organs data as those in spatial data; <br>
(3) annotations list for spots; # a list containing annotations used to find neighbors in matching spatial and single-cell data <br>

### Data for prediction
(1) gene-by-spot matrix and coordinates matrix from spatial transcriptome data; <br>
(2) fragment files from scATAC-seq data; # should contain corresponding tissues/organs data as in spatial data; <br>
(3) annotations list for spots; # the annotations are used to find neighbors in matching spatial and single-cell data <br>

```r
setwd("Workdir")
rm(list=ls())
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(ggplot2)
library(Matrix)
library(data.table)
library(data.table)
library(parallel)
library(scran)
library(presto)
library(dplyr)
library(ggplot2)
```

## Part1: Data Processing 
```r
source("./processing_scATAC_spATAC.R")
source("./predict_spATAC.R")
source("./enhance_spATAC.R")
source("./enhance_fragment.R")
source("./visualization.R")

organs <- c("forebrain")
fragpath <- dir('/path/fragment_data/', "*gz$", full.names = TRUE)
macs2_path <- "/path/macs2"
```

##### step1: call peaks using fragment file of scATAC-seq data
```r
peaks <- Callpeak_scATAC(organs, #the organs/tissues related to the separate fragment filenames used to call peak
                      fragpath, #the path to stored fragment file (tabix)
                      macs2_path #the path to stored macs2 software
                      )
```

##### step2: generate scATAC peak-by-cell matrix and create SeuratObject
```r
organs <- c("forebrain")
peaks <- readRDS("./Processing_scATAC/1.scATAC_allpeaks.rds")
annotation <- readRDS(("./data/Ref1.annotation.rds"))
scATAC <- scATAC_processing (
                            fragpath = fragpath, 
                            peaks = peaks, 
                            organs = organs,
                            mincell = 100,
                            minfeature = 1000,
                            genome = "mm10",
                            sep = c("-", "-"), 
                            annotation = annotation)
```

##### step3: generate spATAC peak-by-spot matrix and create SeuratObject
```r
fragpath_spatac <- dir('/path/0.spATAC_data_/Fragmentfile', "*gz$", full.names = TRUE)
annotation <- readRDS("./data/Ref1.annotation.rds")
peaks <- readRDS("./Processing_scATAC/2.scATAC_mergedpeaks.rds")
spots_anno <- readRDS("./data/Ref2.spATAC_spots_annotations.rds")
spATAC_Object <- spATAC_processing(fragpath = fragpath_spatac, 
                                   peaks = peaks, 
                                   spots_anno = spots_anno,
                                   annotation = annotation)
```

## Part2: Pesca enhancement 
```r
spATAC <- readRDS("./Processing_spATAC/1.SeuratObject_spATAC.rds")
scATAC <- readRDS("./Processing_scATAC/5.scATAC_seuratobject.rds")
organs <- c("forebrain")
enhanced_spATAC <- Enhancement(scATAC = scATAC, # scATAC seurat object
                                spATAC = spATAC , # spATAC seurat object
                                organs = organs, # organs/tissues data for enhancement
                                k.neighbor = 100, # the initial neighbors identified for spots
                                dims = 2:30, # the dimensions when perform the integration
                                k = 10, # the number of final single-cell neighbors for each spot
                                cell_type_pair = FALSE, 
                                k.weight = 20, # k.weight when integrating two datasets
                                lambda = 0.9 # lambda parameter when performing Harmony to integrate
)

spATAC[["enhanced"]] <- enhanced_spATAC
plot <- mySpatialDimPlot(metadata, 
                             row,
                             col, 
                             color, 
                             pointsize = 2.5, 
                             colors = NULL 
)
plot <- mySpatialFeaturePlot(object, # seurat object
                                 feature, # the features used to plot 
                                 pointsize = 2.5, # spot size
                                 scale = FALSE, # whether scale the data
                                 ncol = 2, # the number of column for plots when the number of features is more than two
                                 assay, # the assay data used to plot 
                                 slot = "data" # the slot in the assay 
)
```

## Part3: Pesca prediction
```r
spRNA <- readRDS("/path/spRNA.rds")
scATAC <- readRDS("./Processing_scATAC/5.scATAC_seuratobject.rds")

predicted_spATAC <- Prediction(scATAC = scATAC, # scATAC seurat object
                               spRNA = spRNA, # spRNA seurat object
                               k.weight = 50, # k.weight in integration step
                               default.assay="spatial", # the name of count assay in spRNA data
                               organs = organs, # organs/tissues data for enhancement
                               k.neighbor = 100, # the initial neighbors identified for spots
                               dims = 2:30, # the dimensions when perform the integration
                               k = 10, # the number of final single-cell neighbors for each spot
                               cell_type_pair = TRUE, 
                               lambda = 0.9 # lambda parameter when performing Harmony to integrate
)
```

### Add spatial information before plotting
```r
coordinate <- readRDS("./data/spRNA_coordinate.rds")
spRNA$x <- coordinate$x
spRNA$y <- coordinate$y
```



