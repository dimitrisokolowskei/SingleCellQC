library(Seurat)
library(SeuratObject)

seurat <- readRDS("Zhou_2020.rds")
View(seurat@meta.data)




QCalc <- function(obj, species = c("Human", "Mouse")) {
  
  data <- readRDS(obj)
  
  if(class(data) != "Seurat") 
    stop("Not a seurat object")
  
  if(species == "Human") {
    data$log10GenesPerUMI <- log10(data$nFeature_RNA) / log10(data$nCount_RNA)
    data$mitoRatio <- PercentageFeatureSet(object = data, pattern = "^MT-")
    data$mitoRatio <- data@meta.data$mitoRatio / 100
    data$riboRatio <- PercentageFeatureSet(object = data, pattern = "^RP[SL]")
    data$riboRatio <- data@meta.data$riboRatio / 100
  }
  
  data$log10GenesPerUMI <- log10(data$nFeature_RNA) / log10(data$nCount_RNA)
  data$mitoRatio <- PercentageFeatureSet(object = data, pattern = "^Mt-")
  data$mitoRatio <- data@meta.data$mitoRatio / 100
  data$riboRatio <- PercentageFeatureSet(object = data, pattern = "^Rp[sl]")
  data$riboRatio <- data@meta.data$riboRatio / 100
  
  data
}  


x <- QCalc("Zhou_2020.rds", species = "Mouse")
View(x@meta.data)

QC <- function(obj, annot = F) {
  data <- readRDS(obj)
  
  if(class(data) != "Seurat") 
    stop("Not a seurat object")
  
  if(annot) {
    
  }
  
}

  
y <- TRUE

if(y) {
  print("Hi")
}


x <- QualityControl("Zhou_2020.rds")
View(x@meta.data)
