library(Seurat)
library(tidyverse)
library(celda)

# @seurat_obj Seurat Object
decontX_remove <- function(seurat_obj) {
  decontX_results = celda::decontX(seurat_obj@assays$RNA@counts)
  decontaminated_matrix = decontX_results$decontXcounts
  seurat_obj@assays$RNA@counts = decontaminated_matrix
  return(seurat_obj)
}

deconta <- decontX_remove(seurat_obj)
