##HD5 Seurat Object for Azimuth analysis

library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(SeuratObject)

##Load Object
LAM2_Spatial <- readRDS("spatial_data.RDS")

SaveH5Seurat(LAM2_Spatial, filename = "LAM2_Spatial.h5Seurat")
Convert("LAM2_Spatial.h5Seurat", dest = "h5ad") 
