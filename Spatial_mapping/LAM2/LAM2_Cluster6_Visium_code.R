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
LAM2_Spatial_data <- readRDS(".../LAM2/spatial_data.RDS")


#####Subclustering  ### Cluster6

Idents(LAM2_Spatial_data)

#Step 1: Identify cells from the specific cluster
LAM2_Cluster6 <- WhichCells(LAM2_Spatial_data, ident = 6)

#Step 2: Create a new Seurat object with subsetted cells
LAM2_Cluster6_obj <- subset(LAM2_Spatial_data, cells = LAM2_Cluster6)

#Step 3: Recluster the subsetted Seurat object with increased resolution
#LAM2_Cluster6_obj <- FindNeighbors(LAM2_Cluster6_obj, dims = 1:30)  # You can adjust 'dims' based on your data
#LAM2_Cluster6_obj <- FindClusters(LAMcore_Sc_obj, resolution = 0.5)  # Adjust 'resolution' as needed

DimPlot(LAM2_Cluster6_obj, reduction = 'umap',  label = T)

saveRDS(LAM2_Cluster6_obj, "LAM2_Cluster6_obj.RDS")


####
#spatial_data <- NormalizeData(LAM2_subset)
#spatial_data <- ScaleData(spatial_data)


plot3 <- DimPlot(LAM2_Cluster6_obj, reduction = "umap", label = FALSE) + NoLegend()
plot4 <- SpatialDimPlot(LAM2_Cluster6_obj, label = FALSE, label.size = 3) + NoLegend()

UMAP_Spatial_DimPlot_C6 <- plot3 + plot4
SaveFigure(UMAP_Spatial_DimPlot_C6, "UMAP_Spatial_DimPlot_C6.png", width = 12, height = 6, res = 400)
