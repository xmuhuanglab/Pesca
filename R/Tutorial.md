
## Tutorial for Pesca
### Data for enhancement
(1) fragment files from spatial-ATAC-seq data; <br>
(2) fragment files from scATAC-seq data; # should contain corresponding tissues/organs data as in spatial data; <br>
(3) annotations list for spots; # a list containing annotations used to find neighbors in matching spatial and single-cell data <br>

### Data for prediction
(1) gene-by-spot matrix and coordinates matrix from spatial transcriptome data; <br>
(2) fragment files from scATAC-seq data; # should contain corresponding tissues/organs data as in spatial data; <br>
(3) annotation list for spots; # the annotations are used to find neighbors in matching spatial and single-cell data <br>

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
load("./Data/Ref1.annotation.rdata")
source("./Pesca.R")
```

## Part1: Data Processing 
```r
organs=c("forebrain")
fragpath = dir('/scATAC_fragment_files_path',"*gz$",full.names=TRUE)
macs2_path="/macs2_path"
```

##### step1: call peaks using fragment file from scATAC-seq data
```r
peaks=callpeak_scATAC(organs, # the organs/tissues whose data are used to call peaks
                      fragpath, # the path storing fragment file (tabix), filename of fragment should contain tissue info, like forebrain_fragments.tsv.gz
                      macs2_path # the path storing macs2 software
                      )
```

##### step2: generate scATAC peak-by-cell matrix and create SeuratObject
```r
scATAC_Object=scATAC_processing(fragpath =fragpath, # fragment path 
                                organs = organs, # the organs/tissues whose data are used to generate matrix
                                peaks=peaks  # merged peak set in step1,  Granges format
                               )
```

##### step3: generate spATAC peak-by-spot matrix and create SeuratObject
```r
fragpath_spatac=dir('/spATAC_fragment_files_path',"*gz$",full.names=TRUE)
anno=readRDS("./Ref2.spATAC_spots_annotations.rds") # spots annotations
peaks=readRDS("./Step0_2.MergedPeaks.rds") # merged peak set obtained in step1

spATAC_Object=spATAC_processing(fragpath=fragpath_spatac, # fragment path 
                                peaks=peaks, # merged peak set in step1, Granges format
                                anno = anno # spot annotation
                               )
```

## Part2: Pesca enhancement 
```r
# Integrate scATAC-seq data and spATAC data
Integrated=Enhancement_processing(scATAC = scATAC_Object,
                                  spATAC = spATAC_Object,
                                  organs = organs)
# Enhance spatial ATAC data
Assay_enhanced=Enhancemment(IntegratedData=Integrated, # integrated datasets
                            k.neighbor = 200, # the number of initial neighbors
                            k=10, # the number of final single-cell neighbors for each spot
                            dims = 1:30, # dims when finding neighbors
                            organs)
```

## Part3: Pesca prediction
```r
# Integrate scATAC-seq data and spATAC data
load("./2.ST.Rdata")
anno=readRDS("./Ref3.spRNA_spots_annotations.rds")
Integrated=Prediction_processing(scATAC = scATAC_Object,
                                 spRNA = spRNA_Object,
                                 organs = organs,
                                 default.assay="spatial", # the count assay of spRNA
                                 anno=anno)
# predict spatial ATAC data
Assay_predicted=Enhancemment(IntegratedData=Integrated, # integrated datasets
                             k.neighbor = 200, # initial number of neighbors
                             k=10, # the number of single-cell neighbors for each spot
                             dims = 1:30, # dims used when finding neighbors
                             organs)

```


