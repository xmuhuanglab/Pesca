# Pesca: a computational method to predict and enhance spatial profiling of chromatin accessibility
## Introduction
Pesca is a computational method to predict and enhance spatial profiling of chromatin accessibility using scATAC-seq data.The accuracy and generalizability of the Pesca-enhanced and -predicted data are validated by existing various ground-truth data including spatial co-profiling data, VISTA results, and epigenomic MEIFISH data.

## Workflow
![image](![Pesca](https://github.com/xmuhuanglab/Pesca/assets/95668602/eb23d0c4-a1fa-4baf-86d0-8613892e7146))
See tutorial for more details.

## Tutorial for Pesca
## Pesca enhance spatial ATAC data
### (1)fragment files from spatial-ATAC-seq;
### (2)fragment files from scATAC-seq; # should contain corresponding tissues/organs data as in spatial data
### (3)annotation list for spots; # the annotations are used to find neighbors in matching spatial and single-cell data

## Pesca predict spatial ATAC data
### (1)gene-by-spot matrix and coordinates matrix from spatial transcriptome data;
### (2)fragment files from scATAC-seq; # should contain corresponding tissues/organs data as in spatial data
### (3)annotation list for spots; # the annotations are used to find neighbors in matching spatial and single-cell data
