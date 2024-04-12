library(Seurat)
library(SeuratObject)



seurat <- readRDS("Zhou_2020.rds")
View(seurat@meta.data)

class(seurat)

QualityControl <- function(obj) {
  seurat <- readRDS(obj)
  
  if(class(obj) != "SeuratObject") 
    stop("Not a seurat object")
  print("Everything is fine around here...")
}

  



data <- QualityControl("Zhou_2020.rds")

d <- class(seurat)
View(class(seurat))


is.Seurat(seurat)

is.
SeuratObject::Seurat-