
getwd()
library(pheatmap)
library(grid)
library(ggplot2)
library(matrixStats)
library(tidyverse)
data <-read.csv('Averaged_Kinase_Vs_SEQ_2.csv')
dim(data)
head(data, scale = 'row')
data1 <- as.matrix(data[,-1])
rownames(data1) <- data[,1]
pheatmap(data1, cluster_rows=FALSE, cluster_cols=FALSE)
pheatmap(log2(data1 +1), scale = 'row',row.names = F, show_rownames = T, treeheight_row = 100, cutree_rows = 3, cutree_cols = )
pheatmap(log2(data1 +1), scale = 'row', row.names = TRUE)
my_plot <- pheatmap(log2(data1 +1), scale = 'row', row.names = F, show_rownames = F, treeheight_row = 50)

#Cut the heatmap to pieces
my_plot <- pheatmap(log2(data1 +1), scale = 'row', row.names = F, show_rownames = T, treeheight_row = 100, cutree_rows = 3, cutree_cols = )


out <- pheatmap(data1, 
                show_rownames=F, cluster_cols=T, cluster_rows=T, scale="row",
                cex=1, clustering_distance_rows="euclidean", cex=1,
                clustering_distance_cols="euclidean", clustering_method="complete", border_color=FALSE,
                annotation_col=metadata,
                annotation_row=metadata_gene)

