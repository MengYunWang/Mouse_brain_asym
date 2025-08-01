## script to create the cell annotation refernces for hip and auc from Allen institute

# scripts coming from here: https://community.brain-map.org/t/gene-expression-matrix-cvs-is-too-large-to-load-it/566/17
# data coming from here: https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x: metadata.csv and expression_matrix.hdf5

# reference: Yao, Z., Van Velthoven, C. T., Nguyen, T. N., Goldy, J., Sedeno-Cortes, A. E., Baftizadeh, F., ... & Zeng, H. (2021). 
# A taxonomy of transcriptomic cell types across the isocortex and hippocampal formation. Cell, 184(12), 3222-3241.

# Adopted by M.-Y. Wang 03-Sep-2024

setwd("/data/workspaces/lag/workspaces/lg-func-asym/working_data/Mengyun/")

rm(list = ls())

################# Inherit the argument from the command line
args <- commandArgs(trailingOnly=TRUE)
region <- args[1]  # This is the region "AUD" or "HIP"


## Load necessary libraries
# library(rhdf5)      # For reading hdf5; https://www.bioconductor.org/packages/devel/bioc/vignettes/rhdf5/inst/doc/rhdf5.html
library(HDF5Array)  # Alternatively option for reading the data matrix with less memory: https://rdrr.io/github/Bioconductor/HDF5Array/
library(data.table) # For fast reading of csv files
library(Seurat)

################################################# Step 1: Loading and subsetting the data
## Read in sample and gene names for use later
genes   <- h5read("expression_matrix.hdf5","/data/gene")
samples <- h5read("expression_matrix.hdf5","/data/samples")

## Read in the metadata using fread (fast!)
metadata <- fread("metadata.csv")
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$sample_name
# Note that the order of metadata and counts and the number of cells are DIFFERENT 
#   (there are 107 cells with no metadata). This will be critical later!


## Subsample data into AUD and HIP
## -- If you want to analyze these data in R and especially using Seurat, you will
## --   have a better chance of success if you only work on parts of the data at a
## --   time.  We suggest, using information in the metadata (e.g., cell type or
## --   brain region columns) to subset.   
metadata_select <- metadata[metadata[,6] == region, ]
use_samples <- intersect(sample(rownames(metadata_select)), samples)
read_samples <- sort(match(use_samples,samples))


########################################## Step 2: reading the count 

## Strategy #2: Read in only a relevant subset of data using h5read. This is the 
##   method that probably works best in most situations.
counts <- h5read("expression_matrix.hdf5", "/data/counts", index = list(read_samples, NULL))
counts <- t(counts)
subcounts <- as(counts, "dgCMatrix")


## Add gene and sample names to the data matrix
## -- Note: We'll works with the recommended strategy (#2) for the rest of the script
rownames(subcounts) <- as.character(genes)
colnames(subcounts) <- as.character(samples) [read_samples]


######################################## Step3: Read the subsetted data and metadata into Seurat
## -- Note that this **WILL NOT WORK** for the full 1 Million+ cell data set!
## -- I would strongly encourage other methods for analysis when dealing with more than ~100,000 cells at a time
## -- Also recall that the order of data and meta do not match, so reorder here
allen_refence <- CreateSeuratObject(counts=subcounts)          # Put in a Seurat object
met <- as.data.frame(metadata[colnames(subcounts),]) # Format the metadata
allen_refence <- AddMetaData(allen_refence,met)                           # Add the metadata



################################# Step4: create reference to be used in RCTD
library(spacexr)

# Create the reference object with allen mouse brain cortex atlas

allen_refence <- UpdateSeuratObject(allen_refence)
Idents(allen_refence) <- "subclass_label"
sort(table(Idents(allen_refence))) # check the number of cells of each cell types

# remove the following cells because there aren't enough of them for annotation (<25)
if (region == "AUD") {
  allen_refence <- subset(allen_refence, idents = c("CR", "L2/3 IT ENTl", "L2/3 IT PPP", "L6 IT ENTl", "Micro-PVM", "SMC-Peri", "VLMC"), invert = TRUE)
} else {
  allen_refence <- subset(allen_refence, idents = c( "L2/3 IT CTX", "L6b CTX", "L6 CT CTX", "L5 PT CTX", "L6 IT CTX", "L2 IT ENTl", "L2/3 IT ENTl", 
                                                     "CT SUB", "L2/3 IT PPP", "Sst Chodl","L4/5 IT CTX", "L2/3 IT RHP", "VLMC") , invert = TRUE)
}

sort(table(Idents(allen_refence))) # check the number of cells of each cell types


allen_refence <- subset(allen_refence, downsample = 10000) # max=10000
counts <-  GetAssayData(allen_refence, assay = "RNA", layer = "counts")
cluster <- as.factor(allen_refence$subclass_label)
names(cluster) <- colnames(allen_refence)
nUMI <- allen_refence$nCount_RNA
names(nUMI) <- colnames(allen_refence)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
reference <- Reference(counts, cluster, nUMI)

if (region == "AUD") {
  saveRDS(reference, file = "aud.reference.rds")
} else {
  saveRDS(reference, file = "hip.reference.rds")
}

