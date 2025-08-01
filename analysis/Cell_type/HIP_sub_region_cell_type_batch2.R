## This script is to do cell type annotation in CA1 CA3 AND DG

## created by M.-Y. Wang 19-Sep-2024; updated 21-Mar-2025

# Remove all objects created before to prevent clashing
rm(list = ls())


# Set the working directory to the path where your files are located
# setwd("/Users/joeywang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Mouse_brain/")
# setwd("/Users/wang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Mouse_brain/")
setwd("/data/workspaces/lag/workspaces/lg-func-asym/working_data/Mengyun")


################# Inherit the argument from the command line
args <- commandArgs(trailingOnly=TRUE)
file_name <- args[1]  # This is the filenames "orig", "harmony"
mouse_id <- args[2] # see the lookup variable
hemi <- args[3] # "left" or "right" 
region <- args[4] # CA1 CA3 DG

print(paste("File Name:", file_name))
print(paste("Mouse ID:", mouse_id))
print(paste("Hemi:", hemi))
print(paste("Region:", region))

# load essential libraries
library(Seurat)
library(future)
library(ggplot2)
library(sf)
library(tidyr)
library(dplyr)
library(readxl)
# plan("multisession", workers = 10) # 
options(future.globals.maxSize = 50000 * 1024^2) # to define the max RAM for each worker

#################################Step 1: Load the data
data2read <- paste("intermedia_data/mouse_asym_acr_sample_", file_name, ".rds", sep="")
mouse_asym <- readRDS(file=data2read)

#################################Step 2: crop the image
# load the coordinates
HIP <- read_excel(paste0("Spatial_transcriptomics_batch1_batch2_Xenium_extractions_", region,".xlsx"), sheet = paste0("Coordinates_", region) )
# give the sampleID to the column names
colnames(HIP) <- HIP[2,]
# delete other unrelavant rows
HIP <- HIP[-1:-5, ]
# find the matching rows
matching_columns <- grep(mouse_id, colnames(HIP))
# get the x y coordinates
if (hemi=="left") {
  coordi <- HIP[,matching_columns[1]] %>%
    separate(col=mouse_id, into = c("x", "y"), sep = ",")
} else {
  coordi <- HIP[,matching_columns[2]] %>%
    separate(col=mouse_id, into = c("x", "y"), sep = ",")
}
# get rid of nan values and numerical them
coordi <- na.omit(coordi)
coordi <- as.data.frame(coordi)
coordi <- lapply(coordi, as.numeric)

# Create the SpatialPolygons object for the brain region
create_roi_polygon <- function(xy_pts) {
  p <- Polygon(cbind(xy_pts$x, xy_pts$y))
  ps <- Polygons(list(p), 1)
  roi_to_crop <- SpatialPolygons(list(ps))
  return(roi_to_crop)
}

# Create the polygon
HIP_to_crop <- create_roi_polygon(coordi)


lookup <- list("M234" = "b814",
               "M253" = "b814",
               "M071" = "b917",
               "M083" = "b918",
               "M650" = "b969",
               "M638" = "b969",
               "M076" = "b970",
               "M236" = "b970",
               
               "F073" = "b814",
               "F078" = "b918",
               "F087" = "b969",
               "F090" = "b970"
)

fov_id <- lookup[[mouse_id]]
mouse_asym[[hemi]] <- Overlay(mouse_asym[[fov_id]], HIP_to_crop, invert = FALSE)


#################################Step 3: cell annotation with spacexr

library(spacexr)

# load the reference object with allen mouse brain cell type atlas
reference <- readRDS("hip.reference.rds")


# Create query object including coords, counts, and UMI
#Counts, cluster, and spot information is extracted
query.counts <- GetAssayData(mouse_asym,
                             assay = "Xenium",
                             layer = paste0("counts.", mouse_id))[, Cells(mouse_asym[[hemi]])]
coords <- GetTissueCoordinates(mouse_asym[[hemi]], which = "centroids")
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))

# run RCTD
RCTD <- create.RCTD(query, reference, max_cores = 10)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

# put the results back to the object
annotations.df <- RCTD@results$results_df # pixel by cell type
annotations <- annotations.df$first_type # first possible cell type
names(annotations) <- rownames(annotations.df)

mouse_asym$predicted.celltype <- annotations
keep.cells <- Cells(mouse_asym)[!is.na(mouse_asym$predicted.celltype)]
mouse_asym <- subset(mouse_asym, cells = keep.cells)

######################################### Step 4: create niche >> layers
# construct a new assay called niche containing the cell type composition spatially neighboring each cell.
mouse_asym <-  BuildNicheAssay(object = mouse_asym,
                               fov = hemi,
                               group.by = "predicted.celltype",
                               niches.k = 4, # number of niches, can change to 6 or soemthing
                               neighbors.k = 30 #number of neighbors to consider
)
if(hemi=="left") {
  saveRDS(mouse_asym,
          file = paste0("intermedia_data/cell_type/",
                        mouse_id,
                        "_l",
                        region,
                        "_cell_type_",
                        file_name,
                        ".rds"
          )
  )
} else {
  saveRDS(mouse_asym,
          file = paste0("intermedia_data/cell_type/",
                        mouse_id,
                        "_r",
                        region,
                        "_cell_type_",
                        file_name,
                        ".rds"
          )
  )
}
###################################### Step 5. Plot the AC according to the niche
# group cells either by their cell type identity, or their niche identity.
cell_color_mapping <- c("Astro" = "#EE82EE", "Endo" = "#FF4500", "Lamp5" = "#A52A2A", 
                        "Oligo" = "#D2691E", "Pvalb" = "#00008B", "Sncg" = "#9ACD32",
                        "Sst" = "#4B0082", "Vip" = "#F08080", "CA1-ProS" =  "#008000",
                        "CA2-IG-FC" = "#8A2BE2", "CA3" = "#7FFFD4", "CR" = "#FF0000", 
                        "DG" = "#0000FF", "Micro-PVM" = "#FFFF00", "NP SUB" = "#FFA500",
                        "SMC-Peri" = "#2E8B57", "SUB-ProS" = "#ADFF2F")

celltype.plot <-
  ImageDimPlot(
    mouse_asym,
    group.by = "predicted.celltype",
    size = 1,
    cols = cell_color_mapping,
    dark.background = F
  ) + 
  ggtitle("Cell type")

niche.plot <-
  ImageDimPlot(
    mouse_asym,
    group.by = "niches",
    size = 1.2,
    dark.background = F
  ) + 
  scale_fill_manual(values = c("1"="#EB7D5B", "2"= "#6CA2EA","3"= "#B5D33D","4" = "#FED23F")) +
  ggtitle("Niches")

plot_annotation <- celltype.plot | niche.plot

if(hemi=="left"){
  plot2save <- paste0("output/2_Cell_type/acr_sample/batch1_2/", region, "/Plot_", mouse_id, "_l", region, "_cell_type_", file_name, ".png")
} else {
  plot2save <- paste0("output/2_Cell_type/acr_sample/batch1_2/", region, "/Plot_", mouse_id, "_r", region, "_cell_type_", file_name, ".png")
}
ggsave(plot2save, plot = plot_annotation, width = 8, height = 5, units = 'in', dpi = 300)

###################################### save the data
cell_per_layer <- table(mouse_asym$predicted.celltype, mouse_asym$niches) # cell types in each layer

library(openxlsx)
if(hemi=="left"){
  cell_type2save <- paste("output/2_Cell_type/acr_sample/batch1_2/", region, "/", mouse_id, "_", region, "_cell_type_", file_name, ".xlsx", sep = "")
  write.xlsx(cell_per_layer, cell_type2save, sheetName = "left", rowNames=TRUE, colNames=TRUE)
  
} else {
  cell_type2save <- paste("output/2_Cell_type/acr_sample/batch1_2/", region, "/", mouse_id, "_", region, "_cell_type_", file_name, ".xlsx", sep = "")
  wb <- loadWorkbook(cell_type2save)
  addWorksheet(wb, "right")
  writeData(wb, "right", cell_per_layer)
  saveWorkbook(wb, cell_type2save, overwrite=TRUE)
}




