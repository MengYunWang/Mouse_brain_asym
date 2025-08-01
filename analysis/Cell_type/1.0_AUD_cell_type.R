## This script is to do cell type annotation in auditory cortex

## created by M.-Y. Wang 19-Sep-2024



# Remove all objects created before to prevent clashing
rm(list = ls())

# Set the working directory to the path where the files are located
# setwd("/Users/joeywang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Mouse_brain/")
# setwd("/Users/wang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Mouse_brain/")
setwd("/data/workspaces/lag/workspaces/lg-func-asym/working_data/Mengyun")


#### Inherit the argument from the command line
args <- commandArgs(trailingOnly=TRUE)
file_name <- args[1]  # This is the filenames "orig", "harmony"
mouse_id <- args[2] #see the lookup variable
hemi <- args[3] # "left" or "right" 

print(paste("File Name:", file_name))
print(paste("Mouse ID:", mouse_id))
print(paste("Hemi:", hemi))

# load essential libraries
library(Seurat)
library(future)
library(ggplot2)
library(sf)
library(dplyr)
library(tidyr)
library(readxl)
# plan("multisession", workers = 16) # 
options(future.globals.maxSize = 50000 * 1024^2) # to define the max RAM for each worker

#################################Step 1: Load the data
data2read <- paste("intermedia_data/mouse_asym_acr_sample_", file_name, ".rds", sep="")
mouse_asym <- readRDS(file=data2read)

#################################Step 2: crop the image
# load the coordinates
AUD <- read_excel("Spatial_transcriptomics_batch1_batch2_Xenium_extractions_AUD.xlsx", sheet = "Coordinates_AUD")
# give the sampleID to the column names
colnames(AUD) <- AUD[2,]
# delete other unrelavant rows
AUD <- AUD[-1:-5, ]
# find the matching rows
matching_columns <- grep(mouse_id, colnames(AUD))
# get the x y coordinates
if (hemi=="left") {
  coordi <- AUD[,matching_columns[1]] %>%
    separate(col=mouse_id, into = c("x", "y"), sep = ",")
} else {
  coordi <- AUD[,matching_columns[2]] %>%
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

AC_to_crop <- create_roi_polygon(coordi)


lookup <- list("M670" = "b174",
               "M671" = "b174",
               "M672" = "b277",
               "M673" = "b277",
               "M234" = "b814",
               "M253" = "b814",
               "M071" = "b917",
               "M083" = "b918",
               "M650" = "b969",
               "M638" = "b969",
               "M076" = "b970",
               "M236" = "b970",
               
               "F679" = "b442",
               "M680" = "b442", # here, to be noticed>>F680 is labelled as M680
               "F682" = "b597",
               "F683" = "b597",
               "F687" = "b695",
               "F688" = "b768",
               "F078" = "b918",
               "F087" = "b969",
               "F090" = "b970"
)

fov_id <- lookup[[mouse_id]]
mouse_asym[[hemi]] <- Overlay(mouse_asym[[fov_id]], AC_to_crop, invert = FALSE)


#################################Step 3: cell annotation with spacexr

library(spacexr)

# load the reference object with allen mouse brain cell type atlas
reference <- readRDS("aud.reference.rds")


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
                               niches.k = 5, # number of niches, can change to 6 or soemthing
                               neighbors.k = 30 #number of neighbors to consider
                               )
if(hemi=="left") {
  saveRDS(mouse_asym,
          file = paste0("intermedia_data/cell_type/",
                        mouse_id,
                        "_lAC_cell_type_",
                        file_name,
                        ".rds"
                        )
          )
} else {
  saveRDS(mouse_asym,
          file = paste0("intermedia_data/cell_type/",
                        mouse_id,
                        "_rAC_cell_type_",
                        file_name,
                        ".rds"
          )
  )
}
###################################### Step 5. Plot the AC according to the niche
# group cells either by their cell type identity, or their niche identity.
cell_color_mapping <- c("Astro" = "#FF0000", "Endo" = "#32CD32", "L2-3 IT CTX" = "#0000FF", "L4-5 IT CTX" = "#EE82EE", "L5 PT CTX" = "#FFFF00", 
                        "L5 IT CTX" = "#FFA500", "L6 CT CTX" = "#006400", "L6 IT CTX" = "#ADFF2F", "L6b CTX" = "#FF00FF", "Lamp5" = "#A52A2A", 
                        "Car3" = "#FF4500", "L5-6 IT TPE-ENT" = "#8A2BE2", "L5-6 NP CTX" = "#7FFFD4", "Oligo" = "#D2691E", "Pvalb" = "#00008B", 
                        "Sncg" = "#9ACD32", "Sst" = "#4B0082", "Vip" = "#F08080", "Sst Chodl" = "#20B2AA")

celltype.plot <-
  ImageDimPlot(
    mouse_asym,
    group.by = "predicted.celltype",
    size = 1,
    cols = cell_color_mapping,
    dark.background = F
  ) + 
  ggtitle("Cell type")+ 
  coord_flip()

# niche.plot <-
#   ImageDimPlot(
#     mouse_asym,
#     group.by = "niches",
#     size = 1.5,
#     dark.background = F
#   ) + 
#   scale_fill_manual(values = c("1"="#442288", "2"= "#6CA2EA","3"= "#B5D33D", "4" = "#FED23F","5"= "#EB7D5B")) +
#   ggtitle("Niches")

#plot_annotation <- celltype.plot | niche.plot
plot_annotation <- celltype.plot

if(hemi=="left"){
  plot2save <- paste0("output/2_Cell_type/acr_sample/batch1_2/AUD/Plot_", mouse_id, "_lAC_cell_type_", file_name, ".png")
} else {
  plot2save <- paste0("output/2_Cell_type/acr_sample/batch1_2/AUD/Plot_", mouse_id, "_rAC_cell_type_", file_name, ".png")
}
ggsave(plot2save, plot = plot_annotation, width = 5, height = 5, units = 'in', dpi = 300)

###################################### save the data
cell_per_layer <- table(mouse_asym$predicted.celltype, mouse_asym$niches) # cell types in each layer

library(openxlsx)
if(hemi=="left"){
  cell_type2save <- paste("output/2_Cell_type/acr_sample/batch1_2/AUD/", mouse_id, "_AC_cell_type_", file_name, ".xlsx", sep = "")
  write.xlsx(cell_per_layer, cell_type2save, sheetName = "left", rowNames=TRUE, colNames=TRUE)
  
} else {
  cell_type2save <- paste("output/2_Cell_type/acr_sample/batch1_2/AUD/", mouse_id, "_AC_cell_type_", file_name, ".xlsx", sep = "")
  wb <- loadWorkbook(cell_type2save)
  addWorksheet(wb, "right")
  writeData(wb, "right", cell_per_layer)
  saveWorkbook(wb, cell_type2save, overwrite=TRUE)
}


# if(hemi=="left"){
#   cell_type2save <- paste("output/2_Cell_type/acr_batch/within_batch_95/MB686_lAC_cell_type_", file_name, ".csv", sep = "")
# } else {
#   cell_type2save <- paste("output/2_Cell_type/acr_batch/within_batch_95/MB686_rAC_cell_type_", file_name, ".csv", sep = "")
# }
# # Save cell_per_layer into CSV file
# write.csv(cell_per_layer, cell_type2save, row.names=TRUE)
