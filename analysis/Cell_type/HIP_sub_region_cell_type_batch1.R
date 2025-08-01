## This script is to do cell type annotation in hippocampus

## created by M.-Y. Wang 19-Sep-2024

# Remove all objects created before to prevent clusing
rm(list = ls())

# Set the working directory to the path where your files are located
# setwd("/Users/joeywang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Mouse_brain/")
# setwd("/Users/wang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Mouse_brain/")
setwd("/data/workspaces/lag/workspaces/lg-func-asym/working_data/Mengyun")


################# Inherit the argument from the command line
args <- commandArgs(trailingOnly=TRUE)
file_name <- args[1]  # This is the filenames "orig", "harmony","integrated_rpca"
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
library(dplyr)
library(readxl)
#plan("multisession", workers = 10) # 
options(future.globals.maxSize = 20000 * 1024^2) # to define the max RAM for each worker

#################################Step 1: Load the data
data2read <- paste0("intermedia_data/mouse_asym_acr_sample_", file_name, ".rds")
mouse_asym <- readRDS(file=data2read)

#################################Step 2: crop the image
# load the coordinates
HIP <- read_excel("Spatial_transcriptomics_batch1_AUD_HIP.xlsx", sheet = paste0("Coordinates_", region))
colnames(HIP) <- HIP[1,]
HIP <- HIP[-1:-4,]
matching_columns <- grep(mouse_id, colnames(HIP))

if (hemi=="left") {
  coordi <- HIP[,matching_columns[1]:(matching_columns[1]+1)]
} else {
  coordi <- HIP[,matching_columns[2]:(matching_columns[2]+1)]
}

coordi <- na.omit(coordi)
colnames(coordi) <- c("x", "y")
coordi <- as.data.frame(coordi)
coordi <- lapply(coordi, as.numeric)

# Create the SpatialPolygons object for the brain region
create_roi_polygon <- function(xy_pts) {
  p <- Polygon(cbind(xy_pts$x, xy_pts$y))
  ps <- Polygons(list(p), 1)
  roi_to_crop <- SpatialPolygons(list(ps))
  return(roi_to_crop)
}

HIP_to_crop <- create_roi_polygon(coordi)


lookup <- list("M669" = "b174",
               "M670" = "b174",
               "M671" = "b174",
               
               "M672" = "b277",
               "M674" = "b277",

               "M676" = "b339",
               "M677" = "b339",
               
               "M678" = "b442",
               "F679" = "b442",
               "M680" = "b442", # mistake here M680 should be F680
               
               "F682" = "b597",
               "F683" = "b597",
               
               "F685" = "b695",
               "F686" = "b695",
               "F687" = "b695",
               
               "F688" = "b768"
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
                               niches.k = 3, # number of niches, can change to 6 or soemthing
                               neighbors.k = 30 #number of neighbors to consider
)
if(hemi=="left") {
  saveRDS(mouse_asym,
          file = paste0("intermedia_data/cell_type/",
                        mouse_id,
                        "_l", region, "_cell_type_",
                        file_name,
                        ".rds"
          )
  )
} else {
  saveRDS(mouse_asym,
          file = paste0("intermedia_data/cell_type/",
                        mouse_id,
                        "_r", region, "_cell_type_",
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
    size = 1.2,
    cols = cell_color_mapping,
    dark.background = F
  ) + 
  ggtitle("Cell type")

niche.plot <-
  ImageDimPlot(
    mouse_asym,
    group.by = "niches",
    size = 1.5,
    dark.background = F
  ) + 
  scale_fill_manual(values = c("1"="#EB7D5B", "2"= "#6CA2EA","3"= "#B5D33D")) +
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
  cell_type2save <- paste0("output/2_Cell_type/acr_sample/batch1_2/", region, "/", mouse_id, "_", region, "_cell_type_", file_name, ".xlsx")
  write.xlsx(cell_per_layer, cell_type2save, sheetName = "left", rowNames=TRUE, colNames=TRUE)
  
} else {
  cell_type2save <- paste0("output/2_Cell_type/acr_sample/batch1_2/", region, "/", mouse_id, "_", region, "_cell_type_", file_name, ".xlsx")
  wb <- loadWorkbook(cell_type2save)
  addWorksheet(wb, "right")
  writeData(wb, "right", cell_per_layer)
  saveWorkbook(wb, cell_type2save, overwrite=TRUE)
}





