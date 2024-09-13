library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(SeuratObject)


library(dplyr)
library(Seurat)
library(patchwork)
library(spacexr)
library(ggplot2)
library(patchwork)
library(limma) # optional
options(Seurat.object.assay.version = "v3")


# Set path for saving data
getwd()

data_path <- "data/"
fig_path <- "fig/"


# Convenience functions
SaveFigure <- function(plots, name, type = "png", width, height, res){
  if(type == "png") {
    png(paste0(fig_path, name, ".", type),
        width = width, height = height, units = "in", res = 200)
  } else {
    pdf(paste0(fig_path, name, ".", type),
        width = width, height = height)
  }
  print(plots)
  dev.off()
}

SaveObject <- function(object, name){
  saveRDS(object, paste0(data_path, name, ".RDS"))
}

ReadObject <- function(name){
  readRDS(paste0(data_path, name, ".RDS"))
}



# Convenience function ends here



#Read/Load Data
LAM1_Spatial_data <- readRDS(".........LAM1/spatial_data.RDS")


#####Subclustering  ### Cluster7

Idents(LAM1_Spatial_data)

#Step 1: Identify cells from the specific cluster
LAM1_Cluster7 <- WhichCells(LAM1_Spatial_data, ident = 7)

#Step 2: Create a new Seurat object with subsetted cells
LAM1_Cluster7_obj <- subset(LAM1_Spatial_data, cells = LAM1_Cluster7)

#Step 3: Recluster the subsetted Seurat object with increased resolution
#LAM1_Cluster7_obj <- FindNeighbors(LAM1_Cluster7_obj, dims = 1:30)  # You can adjust 'dims' based on your data
#LAM1_Cluster7_obj <- FindClusters(LAMcore_Sc_obj, resolution = 0.5)  # Adjust 'resolution' as needed

DimPlot(LAM1_Cluster7_obj, reduction = 'umap',  label = T)

saveRDS(LAM1_Cluster7_obj, "LAM1_Cluster7_obj.RDS")


####
#spatial_data <- NormalizeData(LAM1_subset)
#spatial_data <- ScaleData(spatial_data)


plot3 <- DimPlot(spatial_data_C7, reduction = "umap", label = FALSE) + NoLegend()
plot4 <- SpatialDimPlot(spatial_data_C7, label = FALSE, label.size = 3) + NoLegend()

UMAP_Spatial_DimPlot_C7 <- plot3 + plot4
SaveFigure(UMAP_Spatial_DimPlot_C7, "UMAP_Spatial_DimPlot_C7.png", width = 12, height = 6, res = 400)
