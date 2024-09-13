library(Seurat)
library(SeuratObject)
library(SeuratData)
library(SeuratDisk)
library(dplyr)
library(Seurat)
library(patchwork)
library(spacexr)
library(ggplot2)
library(patchwork)
library(limma) # optional
options(Seurat.object.assay.version = "v3")

# Set path for saving data
data_path <- "..../"
fig_path <- "...../"


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

LAM_SPAT1_SPAT2_data <- readRDS("LAM_SPAT1_SPAT2.rds")
int_spatial_data <- NormalizeData(LAM_SPAT1_SPAT2_data)
int_spatial_data <- ScaleData(int_spatial_data)


int_spatial_data <- FindVariableFeatures(int_spatial_data, selection.method = "vst", nfeatures = 2000)

int_spatial_data <- RunPCA(int_spatial_data, features = VariableFeatures(object = int_spatial_data))
int_spatial_data <- FindNeighbors(int_spatial_data, reduction = "pca", dims = 1:30)
int_spatial_data <- FindClusters(int_spatial_data, resolution = 1.0, verbose = FALSE)

int_spatial_data <- RunUMAP(int_spatial_data, reduction = "pca", dims = 1:30)

SaveFigure(UMAP_Spatial_DimPlot_int, "UMAP_Spatial_DimPlot_int.png", width = 6, height = 6, res = 400)

#saveRDS(int_spatial_data, "int_spatial_data.RDS")
#make use of Sc dataset clusters here



#HeatMap for Differential gene expression across clusters in the integrated dataset

SeuratObj_integrated <- readRDS("LAM_SPAT1_SPAT2.rds")

# Step 1: View the clusters in your Seurat object
# Check cluster identities to determine which clusters to compare
Idents(SeuratObj_integrated) <- "seurat_clusters"  # The default cluster identity

# Step 2: Find DEGs across all clusters
# Use FindAllMarkers to compare multiple clusters and find DEGs for each cluster
all_degs <- FindAllMarkers(
  SeuratObj_integrated, 
  only.pos = FALSE,          # TRUE, if you want to find positive markers (upregulated genes)
  min.pct = 0.1,            # Minimum percentage of cells expressing the gene
  logfc.threshold = 0.25    # Minimum log2 fold change threshold
)

# Step 3: Filter for the top DEGs per cluster
# Here, we'll take the top 10 DEGs per cluster based on avg_log2FC
top10_degs <- all_degs %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)  # Top 10 genes by log2 fold change

# Step 3a: Filter for significant DEGs per cluster
# Here, we'll take the all DEGs per cluster based on avg_log2FC
ALL_degs <- all_degs %>% 
  group_by(cluster)  #  genes by log2 fold change

# Step 4: Extract the expression matrix for the top DEGs
# This is to plot in the heatmap
ALL_deg_genes <- ALL_degs$gene
expr_matrix <- as.matrix(GetAssayData(SeuratObj_integrated, slot = "data")[ALL_deg_genes, ])



# Step 4a: Extract the expression matrix for the top DEGs
# This is to plot in the heatmap
top10_deg_genes <- top10_degs$gene
expr_matrix <- as.matrix(GetAssayData(SeuratObj_integrated, slot = "data")[top10_deg_genes, ])


# Step 5: Generate a heatmap using pheatmap for the top DEGs across clusters
pheatmap_result <- pheatmap(
  expr_matrix,
  scale = "cluster",  
  show_rownames = FALSE,
  show_colnames = FALSE,  # Adjust to TRUE if you want to show cell barcodes
  annotation_col = as.data.frame(Idents(SeuratObj_integrated))  # Use clusters as annotation
)

# Step 6: Export DEGs to a CSV file
write.csv(ALL_degs, "IntegratedLAM1_LAM2_DEGs_across_Clusters.csv", row.names = FALSE)

#################################
# Generate a heatmap using DoHeatmap for the DEGs across clusters
all_heatmap <- DoHeatmap(
  SeuratObj_integrated, 
  features = ALL_deg_genes,
  group.by = "seurat_clusters",
  label = TRUE
)

# Show the heatmap
print(all_heatmap)

#################################
# Generate a heatmap using DoHeatmap for the top 10 DEGs across clusters
top10_heatmap <- DoHeatmap(
  SeuratObj_integrated, 
  features = top10_deg_genes,
  group.by = "seurat_clusters",
  label = TRUE
)

# Show the heatmap
print(top10_heatmap)


#genes significant to Cluster4
cluster4_sig <- read.csv("clusters_4_sig.csv")
cluster4_sig_gene = c(cluster4_sig$gene)

# Step 3: Sort genes based on log2FC
# You can sort the genes by log2 fold change if desired
sorted_degs_cluster4_sig <- cluster4_sig %>%
  arrange(desc(avg_log2FC))  # Sort by log2 fold change

# Step 4: Select all genes to plot in the heatmap
# Extract the list of all DEGs (or all genes if desired)
all_genes_cluster4_sig <- unique(sorted_degs_cluster4_sig$gene)  # Get the unique list of all genes

# Generate a heatmap using DoHeatmap for the top Cluster4 sig DEGs across clusters
all_gene_cluster4_sig_heatmap <- DoHeatmap(
  SeuratObj_integrated, 
  features = all_genes_cluster4_sig,
  group.by = "seurat_clusters",
  label = TRUE
)

# Show the heatmap
print(all_gene_cluster4_sig_heatmap)

#25up_25Down_genes significant in Cluster4 across cluster
cluster4_UP_DN_sig <- read.csv("25up_25down_clus4_sig.csv")
cluster4_UP_DN_sig_gene = c(cluster4_UP_DN_sig$gene)

# Generate a heatmap using DoHeatmap for the top Cluster4 sig DEGs across clusters
cluster4_UP_DN_sig_heatmap <- DoHeatmap(
  SeuratObj_integrated, 
  features = cluster4_UP_DN_sig_gene,
  group.by = "seurat_clusters",
  label = TRUE
)

# Show the heatmap
print(cluster4_UP_DN_sig_heatmap)

##
#100Upregulated gene in Cluster4

cluster4_100UP <- read.csv("100_up.csv")
cluster4_100UP_gene = c(cluster4_100UP$gene)

# Generate a heatmap using DoHeatmap for the top Cluster4 sig DEGs across clusters
cluster4_100UP_sig_heatmap <- DoHeatmap(
  SeuratObj_integrated, 
  features = cluster4_100UP_gene,
  group.by = "seurat_clusters",
  label = TRUE
)

# Show the heatmap
print(cluster4_100UP_sig_heatmap)


#100 downregulated gene in Cluster4

cluster4_100DN <- read.csv("100_down.csv")
cluster4_100DN_gene = c(cluster4_100DN$gene)

# Generate a heatmap using DoHeatmap for the top Cluster4 sig DEGs across clusters
cluster4_100DN_sig_heatmap <- DoHeatmap(
  SeuratObj_integrated, 
  features = cluster4_100DN_gene,
  group.by = "seurat_clusters",
  label = TRUE
)

# Show the heatmap
print(cluster4_100DN_sig_heatmap)

