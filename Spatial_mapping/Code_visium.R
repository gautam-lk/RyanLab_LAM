#https://www.10xgenomics.com/resources/analysis-guides/integrating-10x-visium-and-chromium-data-with-r
#library(Seurat)
#library(SeuratObject)
#library(dplyr)
#library(Seurat)
#library(patchwork)
#getwd()
#LAM.data <- Read10X_h5("LAM1_filt_feat_bc_matrix.h5")



# Install R packages (only need to do this once, or if you need to update packages)
install.packages("devtools")
# install.packages("devtools")
options(timeout = 600000000) ### set this to avoid timeout error
devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)

# Load libraries
library(spacexr)
library(Seurat)
library(ggplot2) #for saving result plots

setwd("C:/Users/lgautam/Documents/LAMDATA/New_Folder")

# Read in the tissue_positions_list.csv file (path dependent on data download method from step 2)
vis_coords <- read.csv('spatial/tissue_positions.csv', 
                       header=TRUE)

# Read in the tissue_positions_list.csv file (path dependent on data download method from step 2)
vis_coords <- read.csv('spatial/tissue_positions.csv', 
                       header=TRUE)

head(vis_coords)
# Col 3 is the x coord
# Col 4 is the y coord

# Multiply x and y coordinates by -1 to flip/mirror
vis_coords[,3] <- vis_coords[,3]*-1
vis_coords[,4] <- vis_coords[,4]*-1
head(vis_coords)


# Resave file.
# Note that due to the file name requirement, we will write over the original file with this new one
write.table(vis_coords, 'spatial/tissue_positions.csv', quote=FALSE, row.names=FALSE, sep=',')

# Load Visium data directly from Space Ranger output directory
VisiumData<-read.VisiumSpatialRNA("C:/Users/lgautam/Documents/LAMDATA/New_Folder")

# Create a list of barcodes from the column names of the count matrix
barcodes <- colnames(VisiumData@counts)

# Plot number of UMIs per barcode (spot)
plot_puck_continuous(puck=VisiumData, barcodes=barcodes, plot_val=VisiumData@nUMI, 
                     size=1, ylimit=c(0,round(quantile(VisiumData@nUMI,0.9))), 
                     title='plot of nUMI') 




# Load single cell data from Cell Ranger output directory
Counts<-Read10X_h5("filtered_feature_bc_matrix.h5")
# Or if you downloaded the outputs from the 10x Genomics public datasets page:
#Counts<-Read10X_h5("/path/to/adult_mouse_brain_SC/5k_mouse_brain_CNIK_3pv3_filtered_feature_bc_matrix.h5")

# Create new Seurat Object from the count input
count_SeuratObject<-CreateSeuratObject(Counts)

# Create a count matrix object and take a look at the matrix
sc_counts <- count_SeuratObject@assays$RNA@counts
sc_counts[1:5,1:5]


# Load cell types file
cell_types <- read.csv("azimuth_pred.tsv", sep="\t")
head(cell_types)

# Check unique cell types
unique(cell_types$redicted.celltype_level3)
# Check how many barcodes per cell type annotation
# If needed, update annotation file to meet minimum requirements: spacexr requires at least 25 cells (barcodes) per cell type
table(cell_types$predicted.celltype_level3)

# Group CD4_T (....) with CD8_T (............)
cell_types[cell_types == "CD4_T"] <- "CD4_8_T"
cell_types[cell_types == "CD8_T"] <- "CD4_8_T"

# Group AEC (....) with AM (............)
cell_types[cell_types == "AEC"] <- "AEC_AM"
cell_types[cell_types == "AM"] <- "AEC_AM"


# Group AT1 (....) with AT2 (............)
cell_types[cell_types == "AT1"] <- "AT1_AT2"
cell_types[cell_types == "AT2"] <- "AT1_AT2"

# Group VSMC (....) with SCMF (............)
cell_types[cell_types == "VSMC"] <- "VSMC_SCMF_AF2"
cell_types[cell_types == "SCMF"] <- "VSMC_SCMF_AF2"
cell_types[cell_types == "AF2"] <- "VSMC_SCMF_AF2"



cell_types_filter <- cell_types[cell_types$predicted.celltype_level3.score >= 0.5 & cell_types$mapping.score >=0.5, ]


#head(cell_types)
# The cell_types object is a dataframe with 4 columns: cell (barcode), predicted subclass (cell type name), predicted subclass score, mapping score
# Set the cell type name and barcode as object names
cell_types <- setNames(cell_types_filter[[2]], cell_types_filter[[1]])

# Convert to factor data type
cell_types <- as.factor(cell_types) 
head(cell_types)


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
resultsdir <- "RCTD_output_plots/" # An output directory for plots, can be any name
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
  ggsave(paste(resultsdir, cell_type_names[i],'_weights.jpg', sep=''), height=5, width=5, units='in', dpi=300)
}


