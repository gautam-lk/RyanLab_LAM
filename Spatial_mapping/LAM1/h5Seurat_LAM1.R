##HD5 Seurat Object for Azimuth analysis


library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(SeuratObject)


## Load Object
LAM1_Spatial <- readRDS("spatial_data.RDS")



SaveH5Seurat(LAM1_Spatial, filename = "LAM1_Spatial.h5Seurat")
Convert("LAM1_Spatial.h5Seurat", dest = "h5ad") 
