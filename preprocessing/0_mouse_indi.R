# This script is to load the data and split it into samples

# Created by M.-Y. WANG
# 05-Sep-2024; updated 17-March-2025

# clear workspace
rm(list=ls())


args <- commandArgs(trailingOnly=TRUE)
file_id <- args[1]  # This is the SGE_TASK_ID:21474 21477 21439 21442 20697 20695 21068 58314 58317 58318 58569 58570
batch_id <- args[2]  # This is the 174 277 339 442 597 695 768 814 917 918 969 970


# Set the working directory to the path where your files are located
# setwd("/data/clusterfs/lag/users/men/mouse_asym/")
setwd("/data/workspaces/lag/workspaces/lg-func-asym/working_data/Mengyun")

# load essential libraries
library(Seurat)
library(future)
library(dplyr)

options(future.globals.maxSize = 200000 * 1024^2, future.seed = TRUE) 

####-----------------------------------Step1. Load the data  
input_file= Sys.glob(paste0("XETG_data/output-XETG00101__00", file_id, "*"))
print(input_file)
mouse_asym <- LoadXenium(input_file , fov = paste0("b", batch_id)) # fov will be b174 b277 ...

save_file=paste0("intermedia_data/mouse_asym_", batch_id, ".rds", sep="")
saveRDS(mouse_asym,file=save_file)
mouse_asym
###-----------------------------------Step2. split the data into three samples
cordi <- mouse_asym@images[[paste0("b", batch_id)]]@boundaries[["centroids"]]@coords


batch_conditions <- list(
  "174" = c("M669", "M670", "M671"),
  "277" = c("M672", "M673", "M674"),
  "442" = c("M678", "F679", "F680"),
  "597" = c("F681", "F682", "F683"),
  "695" = c("F685", "F686", "F687"),
  "814" = c("M234", "M253", "F073"),
  "969" = c("M650", "M638", "F087"),
  "970" = c("M076", "M236", "F090")
)

if (batch_id %in% names(batch_conditions)) {
  mouse_asym@meta.data$sample_id <- case_when(
    cordi[,2] < 7000 ~ batch_conditions[[batch_id]][1],
    cordi[,2] > 7000 & cordi[,2] < 16000 ~ batch_conditions[[batch_id]][2],
    cordi[,2] > 16000 ~ batch_conditions[[batch_id]][3],
    TRUE ~ "Other"
  )
  mouse_asym$Xenium <- split(mouse_asym$Xenium, f = mouse_asym$sample_id)
  
} else if(batch_id == "339") {
  mouse_asym@meta.data$sample_id <- case_when(
    cordi[,2] < 7000 ~ "M675",
    cordi[,2] > 7000 & cordi[,2] < 14500 ~ "M676",
    cordi[,2] > 15000 ~ "M677",
    TRUE ~ "Other"
  )
  
  mouse_asym$Xenium <- split(mouse_asym$Xenium, f = mouse_asym$sample_id)
  
} else if(batch_id == "768") {
  mouse_asym@meta.data$sample_id <- case_when(
    cordi[,2] < 7000 ~ "F688",
    cordi[,2] > 7000 & cordi[,2] < 15200 ~ "F880",
    cordi[,2] > 15300 ~ "F890",
    TRUE ~ "Other"
  )
  
  mouse_asym$Xenium <- split(mouse_asym$Xenium, f = mouse_asym$sample_id)
  mouse_asym <- subset(mouse_asym, sample_id=="F688", invert = FALSE)
  
} else if(batch_id == "917") {
  mouse_asym@meta.data$sample_id <- case_when(
    cordi[,2] < 7000 ~ "M071",
    cordi[,2] > 7000 & cordi[,2] < 15000 ~ "M083",
    cordi[,2] > 15000 ~ "F078",
    TRUE ~ "Other"
  )
  
  mouse_asym$Xenium <- split(mouse_asym$Xenium, f = mouse_asym$sample_id)
  mouse_asym <- subset(mouse_asym, sample_id=="M071", invert = FALSE)
  
} else if(batch_id == "918") {
  mouse_asym@meta.data$sample_id <- case_when(
    cordi[,2] < 7000 ~ "M071",
    cordi[,2] > 7000 & cordi[,2] < 15000 ~ "M083",
    cordi[,2] > 15000 ~ "F078",
    TRUE ~ "Other"
  )
  
  mouse_asym$Xenium <- split(mouse_asym$Xenium, f = mouse_asym$sample_id)
  mouse_asym <- subset(mouse_asym, sample_id=="M071", invert = TRUE)
}

mouse_asym
save_file_split=paste0("intermedia_data/mouse_asym_", batch_id, "_split.rds")
saveRDS(mouse_asym, file=save_file_split)



