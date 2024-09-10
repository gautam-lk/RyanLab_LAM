
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(SeuratObject)



##HD5 Seurat Object Saving Sc_LAMCore
LAM2_Spatial <- readRDS("spatial_data.RDS")



SaveH5Seurat(LAM2_Spatial, filename = "LAM2_Spatial.h5Seurat")
Convert("LAM2_Spatial.h5Seurat", dest = "h5ad") 
