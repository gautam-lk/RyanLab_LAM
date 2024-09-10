#https://www.10xgenomics.com/resources/analysis-guides/integrating-10x-visium-and-chromium-data-with-r
#library(Seurat)
#library(SeuratObject)
#library(dplyr)
#library(Seurat)
#library(patchwork)


getwd()
library(Seurat)
library(SeuratObject)
library(dplyr)
library(Seurat)
library(patchwork)
options(Seurat.object.assay.version = "v3")
library(spacexr)
library(xlsx)

LAM1.data <- Read10X_h5("LAM1_filt_feat_bc_matrix.h5")


library(ggplot2) #for saving result plots






#########################


# Load Visium data directly from Space Ranger output directory

# Read in the tissue_positions_list.csv file (path dependent on data download method from step 2)
vis_coords <- read.csv('spatial/tissue_positions.csv', 
                       header=TRUE)


head(vis_coords)
# Col 3 is the x coord
# Col 4 is the y coord

# Multiply x and y coordinates by -1 to flip/mirror
vis_coords[,3] <- vis_coords[,3]*1
vis_coords[,4] <- vis_coords[,4]*-1
head(vis_coords)


# Resave file.
# Note that due to the file name requirement, we will write over the original file with this new one
write.table(vis_coords, 'spatial/tissue_positions.csv', quote=FALSE, row.names=FALSE, sep=',')

# Load Visium data directly from Space Ranger output directory
VisiumData<-read.VisiumSpatialRNA("C:/Users/lgautam/Documents/LAMDATA/07222024/Azimuth/LAM1")

# Create a list of barcodes from the column names of the count matrix
barcodes <- colnames(VisiumData@counts)


# Plot number of UMIs per barcode (spot)
plot_puck_continuous(puck=VisiumData, barcodes=barcodes, plot_val=VisiumData@nUMI, 
                     size=1, ylimit=c(0,round(quantile(VisiumData@nUMI,0.9))), 
                     title='plot of nUMI') 



# Load single cell data from Cell Ranger output directory
Counts<-Read10X_h5("LAM1_filt_feat_bc_matrix.h5")

# Create new Seurat Object from the count input
count_SeuratObject<-CreateSeuratObject(Counts)

# Create a count matrix object and take a look at the matrix
sc_counts <- count_SeuratObject@assays$RNA@counts
sc_counts[1:5,1:5]


# Load cell types file

#____--------------------------------------------------
# Load cell types file
cell_types <- read.csv("azimuth_pred_LAM1_finest.tsv", sep="\t")
head(cell_types)

# Check unique cell types
unique(cell_types$predicted.ann_finest_level)
# Check how many barcodes per cell type annotation
# If needed, update annotation file to meet minimum requirements: spacexr requires at least 25 cells (barcodes) per cell type
table(cell_types$predicted.ann_finest_level)

# Group Cell types (....) with (............)

cell_types[cell_types == "Plasma"] <- "Plasma_CD8_T"
cell_types[cell_types == "CD8_T"] <- "Plasma_CD8_T"


cell_types[cell_types == "AT1"] <- "AT1_AT2"
cell_types[cell_types == "AT2"] <- "AT1_AT2"



cell_types[cell_types == "SM"] <- "SM_AF_PF"
cell_types[cell_types == "AF"] <- "SM_AF_PF"
cell_types[cell_types == "PF"] <- "SM_AF_PF"



cell_types_filter <- cell_types[cell_types$predicted.ann_finest_level.score >= 0.5 & cell_types$mapping.score >=0.5, ]
table(cell_types_filter$predicted.ann_finest_level)




#head(cell_types)
# The cell_types object is a dataframe with 4 columns: cell (barcode), predicted subclass (cell type name), predicted subclass score, mapping score
# Set the cell type name and barcode as object names
cell_types <- setNames(cell_types_filter[[2]], cell_types_filter[[1]])

# Convert to factor data type
cell_types <- as.factor(cell_types) 
head(cell_types)


#There were errors in the cell type 


#Step 7c: Create nUMI object
head(count_SeuratObject@meta.data)

# Format for spacexr
sc_umis <- count_SeuratObject@meta.data[,c(1,2)]
# Set object names as nCount_RNA and barcode
sc_umis <- setNames(sc_umis[[2]], rownames(sc_umis))
head(sc_umis)

#Create single cell reference object
SCreference <- Reference(sc_counts, cell_types, sc_umis)



# Optional: spacexr runs faster with more cores, but you may need to adjust system environment settings
Sys.setenv("OPENBLAS_NUM_THREADS"=2)

# Create and run RCTD algorithm
myRCTD <- create.RCTD(VisiumData, SCreference, max_cores = 2)
myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")


# Create the output directory in your working directory
resultsdir <- "RCTD_output_plots/level_finest/" # An output directory for plots, can be any name
dir.create(resultsdir) 

# Create variables from the myRCTD object to plot results
barcodes <- colnames(myRCTD@spatialRNA@counts) # list of spatial barcodes
weights <- myRCTD@results$weights # Weights for each cell type per barcode

# Normalize per spot weights so cell type probabilities sum to 1 for each spot
norm_weights <- normalize_weights(weights) 
cell_type_names<-colnames(norm_weights) # List of cell types
dim(norm_weights) # matrix of 3343 spots and 7 cell types
# For each spot barcode (row), you can see the normalized weight for each cell type (column)

# Look at cell type normalized weights for 2 spots
subset_df <- as.data.frame(t(as.data.frame(norm_weights[1:2,])))
subset_df$celltypes <- rownames(subset_df); rownames(subset_df) <- NULL
subset_df[order(subset_df$'AAACAGAGCGACTCCT-1', decreasing=T),]

# Plot cell type probabilities (normalized) per spot (red = 1, blue = 0 probability)
# Save each plot as a jpg file
for(i in 1:length(cell_types)){
  plot_puck_continuous(myRCTD@spatialRNA, barcodes, norm_weights[,cell_type_names[i]], title =cell_type_names[i], size=1)
  ggsave(paste(resultsdir, cell_type_names[i],'_weights_5.jpg', sep=''), height=5, width=5, units='in', dpi=300)
}

#_______________________________________________________________________________________________
# mapping.score >=0.25

# Load cell types file
cell_types <- read.csv("azimuth_pred_LAM1_finest.tsv", sep="\t")
head(cell_types)

# Check unique cell types
unique(cell_types$predicted.ann_finest_level)
# Check how many barcodes per cell type annotation
# If needed, update annotation file to meet minimum requirements: spacexr requires at least 25 cells (barcodes) per cell type
table(cell_types$predicted.ann_finest_level)

cell_types[cell_types == "T_Cell"] <- "T_cell_B_cell"
cell_types[cell_types == "B_Cell"] <- "T_cell_B_cell"


cell_types[cell_types == "AT1"] <- "AT1_AT2"
cell_types[cell_types == "AT2"] <- "AT1_AT2"


cell_types[cell_types == "SM"] <- "SM_PF"
cell_types[cell_types == "PF"] <- "SM_PF"

cell_types_filter_2_5 <- cell_types[cell_types$predicted.ann_finest_level.score >= 0.25 & cell_types$mapping.score >=0.25, ]


#head(cell_types)
# The cell_types object is a dataframe with 4 columns: cell (barcode), predicted subclass (cell type name), predicted subclass score, mapping score
# Set the cell type name and barcode as object names
cell_types_2_5 <- setNames(cell_types_filter_2_5[[2]], cell_types_filter_2_5[[1]])

# Convert to factor data type
cell_types_2_5 <- as.factor(cell_types_2_5) 
head(cell_types_2_5)


#There were errors in the cell type 


#Step 7c: Create nUMI object
head(count_SeuratObject@meta.data)

# Format for spacexr
sc_umis <- count_SeuratObject@meta.data[,c(1,2)]
# Set object names as nCount_RNA and barcode
sc_umis <- setNames(sc_umis[[2]], rownames(sc_umis))
head(sc_umis)

#Create single cell reference object
SCreference <- Reference(sc_counts, cell_types_2_5, sc_umis)



# Optional: spacexr runs faster with more cores, but you may need to adjust system environment settings
Sys.setenv("OPENBLAS_NUM_THREADS"=2)

# Create and run RCTD algorithm
myRCTD_25 <- create.RCTD(VisiumData, SCreference, max_cores = 2)
myRCTD_25 <- run.RCTD(myRCTD_25, doublet_mode = "full")


# Create the output directory in your working directory
resultsdir <- "RCTD_output_plots/level_finest/" # An output directory for plots, can be any name
dir.create(resultsdir) 

# Create variables from the myRCTD object to plot results
barcodes <- colnames(myRCTD_25@spatialRNA@counts) # list of spatial barcodes
weights <- myRCTD_25@results$weights # Weights for each cell type per barcode

# Normalize per spot weights so cell type probabilities sum to 1 for each spot
norm_weights <- normalize_weights(weights) 
cell_type_names_2_5<-colnames(norm_weights) # List of cell types
dim(norm_weights) # matrix of 3343 spots and 7 cell types
# For each spot barcode (row), you can see the normalized weight for each cell type (column)

# Look at cell type normalized weights for 2 spots
subset_df <- as.data.frame(t(as.data.frame(norm_weights[1:2,])))
subset_df$celltypes_2_5 <- rownames(subset_df); rownames(subset_df) <- NULL
subset_df[order(subset_df$'AAACAGAGCGACTCCT-1', decreasing=T),]

# Plot cell type probabilities (normalized) per spot (red = 1, blue = 0 probability)
# Save each plot as a jpg file
for(i in 1:length(cell_types_2_5)){
  plot_puck_continuous(myRCTD_2_5@spatialRNA, barcodes, norm_weights[,cell_type_names[i]], title =cell_type_names[i], size=1)
  ggsave(paste(resultsdir, cell_type_names[i],'_weights_2_5.jpg', sep=''), height=5, width=5, units='in', dpi=300)
}
#Code Ends here##