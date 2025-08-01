## This script is to plot the example hip cell type figure

## created by M.-Y. Wang 17-Jul-2025

# Remove all objects created before to prevent clashing
rm(list = ls())

# Set the working directory to the path where your files are located
# setwd("/Users/joeywang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Mouse_brain/")
# setwd("/Users/wang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Mouse_brain/")
setwd("/data/workspaces/lag/workspaces/lg-func-asym/working_data/Mengyun")


################# Step 0. Inherit the argument from the command line
args <- commandArgs(trailingOnly=TRUE)
file_name <- "harmony"  # This is the filenames "orig", "harmony","integrated_rpca"
mouse_id <- "F687" # see the lookup variable
hemi <- "left" # "left" or "right" 
region <- "CA1"  # CA1 CA3 DG

print(paste("File Name:", file_name))
print(paste("Mouse ID:", mouse_id))
print(paste("Hemi:", hemi))
print(paste("Region:", region))

# load essential libraries
library(Seurat)
library(future)
library(ggplot2)
library(sf)
library(concaveman)
library(dplyr)
library(readxl)
library(tidyr)

#plan("multisession", workers = 10) #
options(future.globals.maxSize = 20000 * 1024^2) # to define the max RAM for each worker

################################# Step 1: Load and update the data
mouse_asym <- readRDS("intermedia_data/cell_type/F687_lhip_cell_type_harmony.rds")

# function to process the coordinates
process_coordi <- function(region) {
  HIP <- read_excel("Spatial_transcriptomics_batch1_AUD_HIP.xlsx", sheet = paste0("Coordinates_", region))
  colnames(HIP) <- HIP[1, ]
  HIP <- HIP[-1:-4, ]
  matching_columns <- grep(mouse_id, colnames(HIP))
  
  if (hemi == "left") {
    coordi <- HIP[, matching_columns[1]:(matching_columns[1] + 1)]
  } else {
    coordi <- HIP[, matching_columns[2]:(matching_columns[2] + 1)]
  }
  
  return (coordi)
}

# get the coordi of the CA1 region
coordi <- process_coordi(region)

coordi <- na.omit(coordi)
colnames(coordi) <- c("X", "Y")
coordi <- as.data.frame(coordi)
coordi$x <- as.numeric(coordi$X)
coordi$y <- as.numeric(coordi$Y)
polygon_df <- coordi  # or directly from clipboard/data frame

# Extract the whole spatial coordinates using Seurat's proper accessor
coords <- GetTissueCoordinates(mouse_asym, image = "left")

# Identify CA1 inside the whole HIP region
inside <- sp::point.in.polygon(
  point.x = coords$x,
  point.y = coords$y,
  pol.x   = polygon_df$X,
  pol.y   = polygon_df$Y
)

# Get cell barcodes inside the polygon
inside_cells <- coords$cell[inside == 1]

# Update predicted cell types
mouse_asym$predicted.celltype <- as.character(mouse_asym$predicted.celltype)

mouse_asym$predicted.celltype[
  colnames(mouse_asym) %in% inside_cells &
    mouse_asym$predicted.celltype == "CA3"
] <- "unknown"

mouse_asym$predicted.celltype <- factor(mouse_asym$predicted.celltype)

################################# Step 2: plot the image
# group cells either by their cell type identity or their niche identity.
cell_color_mapping <- c(
  "Astro"       = "#E41A1C",  # red
  "Endo"        = "#377EB8",  # blue
  "Car3"        = "#4DAF4A",  # green
  "Lamp5"       = "#984EA3",  # purple
  "Oligo"       = "#FF7F00",  # orange
  "Pvalb"       = "#FFFF33",  # yellow
  "Sncg"        = "#A65628",  # brown
  "Sst"         = "#F781BF",  # pink
  "Vip"         = "#999999",  # grey-neutral (just one)
  "CR"          = "#66C2A5",  # teal
  "Micro-PVM"   = "#FC8D62",  # salmon
  "NP SUB"      = "#8DA0CB",  # light blue
  "SMC-Peri"    = "#E78AC3",  # magenta
  "SUB-ProS"    = "#A6D854",  # light green
  "Sst Chodl"   = "#FFD92F",   # gold
  "CA1-ProS" = "#008000",     # Green for CA1
  "CA2-IG-FC" = "#8A2BE2",    # Blue-violet for CA2 (unchanged)
  "CA3" = "#7FFFD4",          # Aquamarine for CA3
  "DG" = "#0000FF",
  "unknown" = "#D3D3D3"
)
  
  # "Astro" = "#EE82EE", "Endo" = "#FF4500", "Lamp5" = "#A52A2A", 
  #                       "Oligo" = "#D2691E", "Pvalb" = "#00008B", "Sncg" = "#9ACD32",
  #                       "Sst" = "#4B0082", "Vip" = "#F08080", "CA1-ProS" =  "#008000",
  #                       "CA2-IG-FC" = "#8A2BE2", "CA3" = "#7FFFD4", "CR" = "#FF0000", 
  #                       "DG" = "#0000FF", "Micro-PVM" = "#FFFF00", "NP SUB" = "#FFA500",
  #                       "SMC-Peri" = "#2E8B57", "SUB-ProS" = "#ADFF2F", "unknown" = "#D3D3D3"

celltype.plot <-
  ImageDimPlot(
    mouse_asym,
    group.by = "predicted.celltype",
    size = 1.2,
    cols = cell_color_mapping,
    dark.background = F
  ) + 
  ggtitle("Cell type") +
  
  theme(
    legend.text  = element_text(family = "Helvetica", size = 8),
    legend.title = element_text(family = "Helvetica", face = "bold", size = 9),
    # optional: make the whole plot use the same family
    text         = element_text(family = "Helvetica")
  )


plot2save <- paste0("output/2_Cell_type/acr_sample/batch1_2/", region, "/Plot_", mouse_id, "_lhip_cell_type_", file_name, "_merged.pdf")

ggsave(plot2save, plot = celltype.plot, device   = cairo_pdf, width = 5, height = 5, units = 'in')



