LAM1data <- readRDS("...../LAM1/spatial_data.RDS")

Idents(LAM1data) <- "seurat_clusters"  # The default cluster identity

# Step 2: Find DEGs across all clusters
# Use FindAllMarkers to compare multiple clusters and find DEGs for each cluster
all_degs <- FindAllMarkers(
  LAM1data, 
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

# Step 6: Export DEGs to a CSV file
write.csv(ALL_degs, "all_DEGs_across_clusters_LAM1.csv", row.names = FALSE)

