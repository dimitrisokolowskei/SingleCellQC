library(Seurat)
library(tidyverse)
library(scDblFinder)
library(BiocParallel)
library(Matrix)
library(SingleCellExperiment)

# Doublet Identification ----
set.seed(12345)
sce <- as.SingleCellExperiment(seurat_obj)
bp <- SnowParam(4, RNGseed=1234)
sce <- scDblFinder(sce, samples="sample", BPPARAM=bp)
merged_seurat <- as.Seurat(sce)




