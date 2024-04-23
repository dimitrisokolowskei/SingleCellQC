library(Seurat)
library(SeuratObject)
library(tidyverse)
library(svglite)

seurat <- readRDS("pbmc_small.rds")
seurat <- readRDS("Zhou_2020.rds")
View(seurat@meta.data)

# Annotation function ----
seurat <- readRDS("Zhou_2020.rds")
class(seurat)

QCalc <- function(obj, species = c("Human", "Mouse")) {
  
  if(class(obj) != "Seurat")
    stop("Not a seurat object")
  
  obj$log10GenesPerUMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)
  
  if(species == "Human") {
    obj$mitoRatio <- PercentageFeatureSet(object = obj, pattern = "^MT-")
    obj$mitoRatio <- obj@meta.data$mitoRatio / 100
    obj$riboRatio <- PercentageFeatureSet(object = obj, pattern = "^RP[SL]")
    obj$riboRatio <- obj@meta.data$riboRatio / 100
  
  } else if(species == "Mouse") {
    obj$mitoRatio <- PercentageFeatureSet(object = data, pattern = "^mt-")
    obj$mitoRatio <- obj@meta.data$mitoRatio / 100
    obj$riboRatio <- PercentageFeatureSet(object = data, pattern = "^Rpl|^Rps")
    obj$riboRatio <- obj@meta.data$riboRatio / 100
  }
  
  return(obj)
}
x <- QCalc(seurat, species = "Human")
View(x@meta.data)


# Quality Controrl Cut-off ----
QC <- function(obj, gene = 200, mito = 0.25, ribo = 0.20, umi = 400, genes_umi = 0.65, doublet = 0.0, min.cells = 5) {
  
  if(class(obj) != "Seurat") 
    stop("Not a seurat object")
  message("Seurat object found. Processing...")
  
  if(is.numeric(gene) & is.numeric(mito) & is.numeric(ribo) & is.numeric(umi)) {
    filtered_seurat <- subset(x = obj, subset = (nCount_RNA >= umi) & 
                                                (nFeature_RNA >= gene) &
                                                (log10GenesPerUMI > genes_umi) & 
                                                (mitoRatio < mito) &
                                                (riboRatio <= ribo) &
                                                (scDblFinder.score < doublet))
    
    counts <- GetAssayData(filtered_seurat, slot = "counts") # v5 -> count = "counts"
    nonzero <- counts > 0
    keep_genes <- Matrix::rowSums(nonzero) >= min.cells
    filtered_counts <- counts[keep_genes, ]
    filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
    
    return(filtered_seurat)
  }
}

a <- QC(x)  
View(a@meta.data)

# Quality Control Visualization ----
QCvis <- function(obj, feature = NULL, param = NULL, xintercept = NULL, log10 = T) {
  metadata <- obj@meta.data
  
  if (log10) {
    
    ggplot(metadata, aes_string(color=feature, x=param, fill=feature)) + 
      geom_density(alpha = 0.2) + 
      scale_x_log10() + 
      theme_classic() +
      ylab("Cell density") +
      geom_vline(xintercept = 600)
    
  } else {
    ggplot(metadata, aes_string(color=feature, x=param, fill=feature)) + 
      geom_density(alpha = 0.2) + 
      theme_classic() +
      geom_vline(xintercept = 0.020)
  }
    
}

QCvis(x, feature = "Sample", param = "log10GenesPerUMI", log10 = F)
