library(Seurat)
library(SeuratObject)
library(SeuratData)
library(SeuratDisk)
library(dplyr)
library(tidyverse)
library(Seurat)
library(patchwork)
library(spacexr)
library(ggplot2)
library(patchwork)
library(limma) # optional
options(Seurat.object.assay.version = "v3")

# Set path for saving data
data_path <- "........./LAM2/"
fig_path <- "........../LAM2/"


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





####Spatial Analysis starts here
#### Analysis used the the RDS filed saved for each LAM sections

LAM2data <- readRDS(".../LAM2/spatial_data.RDS")
spatial_data <- NormalizeData(LAM2data)
spatial_data <- ScaleData(spatial_data)


spatial_data <- FindVariableFeatures(spatial_data, selection.method = "vst", nfeatures = 2000)

spatial_data <- RunPCA(spatial_data, features = VariableFeatures(object = spatial_data))
spatial_data <- FindNeighbors(spatial_data, reduction = "pca", dims = 1:30)
spatial_data <- FindClusters(spatial_data, resolution = 1.0, verbose = FALSE)

spatial_data <- RunUMAP(spatial_data, reduction = "pca", dims = 1:30)


plot3 <- DimPlot(spatial_data, reduction = "umap", label = TRUE) + NoLegend()
plot4 <- SpatialDimPlot(spatial_data, label = TRUE, label.size = 3) + NoLegend()

UMAP_Spatial_DimPlot <- plot3 + plot4
SaveFigure(UMAP_Spatial_DimPlot, "UMAP_Spatial_DimPlot.png", width = 12, height = 6, res = 400)

#saveRDS(spatial_data, "spatial_data.RDS")
#make use of Sc dataset clusters here

#HeatMap
DoHeatmap(
  LAM2data,
  features = NULL,
  cells = NULL,
  group.by = "ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = NULL,
  label = FALSE, # Change ro True for LABEL
  size = 5.5,
  hjust = 0,
  vjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)
#HeatMap of selected genes

DoHeatmap(
  LAM2data,
  features = c("COL3A1", "ACTA2", "VEGFR2", "EFNB2", "DLL4","KLF2", "EPHB4", "BRG1", "LYVE1", "PROX1", "D2-40", "VEGFR3", "PDPN", "SOX18", "PECAM1", "CD31",  "VECAD", "VWF", "VEGFC", "VEGFD", 
               "P4HB", "MAS516", "PDGFRA", "FLT4"),
  cells = NULL,
  group.by = "ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = NULL,
  label = TRUE,
  size = 5.5,
  hjust = 0,
  vjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)
#Save

#...Identify spatially variable genes


#######Read Moransi data file from directory, if already saved.

getwd()
#Spatial location of Specific LAM core features/genes
LAM2_moransi <- readRDS('LAM2_moransi.rds')

LAM_feature_spatialFeature_test <- SpatialFeaturePlot(LAM2_moransi, 
features = c("COL3A1", "ACTA2", "VEGFR2", "EFNB2", "DLL4","KLF2", "EPHB4", "BRG1", "LYVE1", "PROX1", "D2-40", "VEGFR3", "PDPN", "SOX18", "PECAM1", "CD31",  "VECAD", "VWF", "VEGFC", "VEGFD", 
             "P4HB", "MAS516", "PDGFRA", "FLT4"), ncol = 4, alpha = c(0.1, 1)) +  
  plot_annotation(
    title = "large_geneset_test",
    subtitle = "LAM2")
SaveFigure(LAM_feature_spatialFeature_test, "large_geneset_test.png", width = 40, height = 40)

#LAM2_Moransi_selected_gene_extended
SpatialFeaturePlot(LAM2_moransi, 
                   features = c("COL3A1", "ACTA2", "VEGFR2", "EFNB2", "DLL4","KLF2", "EPHB4", "BRG1", "LYVE1", "PROX1", "D2-40", "VEGFR3", "PDPN", "SOX18", "PECAM1", "CD31",  "VECAD", "VWF", "VEGFC", "VEGFD", 
                                "P4HB", "MAS516", "PDGFRA", "PDGFRB", "FLT4", "ESR1", "VEGFD", "ACTA2", "RAMP1", "PMEL"), ncol = 4, alpha = c(0.1, 1)) +  
  plot_annotation(
    title = "large_geneset_test",
    subtitle = "LAM2")


DoHeatmap(
  LAM2data,
  features = c("COL3A1", "ACTA2", "VEGFR2", "EFNB2", "DLL4","KLF2", "EPHB4", "BRG1", "LYVE1", "PROX1", "D2-40", "VEGFR3", "PDPN", "SOX18", "PECAM1", "CD31",  "VECAD", "VWF", "VEGFC", "VEGFD", 
               "P4HB", "MAS516", "PDGFRA", "FLT4"),
  cells = NULL,
  group.by = "ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = NULL,
  label = TRUE,
  size = 5.5,
  hjust = 0,
  vjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)

all

DoHeatmap(
  LAM2data,
  features = NULL,
  cells = NULL,
  group.by = "ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = NULL,
  label = FALSE,
  size = 5.5,
  hjust = 0,
  vjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)
#Code ends here..

