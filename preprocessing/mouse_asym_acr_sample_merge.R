# This script is to do harmonization across the 12 batches

# Created by M.-Y. WANG
# 31-Aug-2024; updated 17-Mar-2025

# clear workspace
rm(list=ls())

# Set the working directory to the path where your files are located
# setwd("/data/clusterfs/lag/users/menwan2/mouse_asym/")
setwd("/data/workspaces/lag/workspaces/lg-func-asym/working_data/Mengyun")

# load essential libraries
library(Seurat)
library(parallel)
library(future)
library(dplyr)

plan("cluster", workers = 15)
options(future.globals.maxSize = 350000 * 1024^2, future.seed = TRUE) 

#######################################
#-------Step1. Load the data and merge 
#######################################
# load the data
mouse_asym_1 <- readRDS(file="intermedia_data/mouse_asym_174_split.rds")
mouse_asym_2 <-readRDS(file="intermedia_data/mouse_asym_277_split.rds")
mouse_asym_3 <-readRDS(file="intermedia_data/mouse_asym_339_split.rds")
mouse_asym_4 <-readRDS(file="intermedia_data/mouse_asym_442_split.rds")
mouse_asym_5 <-readRDS(file="intermedia_data/mouse_asym_597_split.rds")
mouse_asym_6 <-readRDS(file="intermedia_data/mouse_asym_695_split.rds")
mouse_asym_7 <-readRDS(file="intermedia_data/mouse_asym_768_split.rds")
mouse_asym_8 <-readRDS(file="intermedia_data/mouse_asym_814_split.rds")
mouse_asym_9 <-readRDS(file="intermedia_data/mouse_asym_917_split.rds")
mouse_asym_10 <-readRDS(file="intermedia_data/mouse_asym_918_split.rds")
mouse_asym_11 <-readRDS(file="intermedia_data/mouse_asym_969_split.rds")
mouse_asym_12 <-readRDS(file="intermedia_data/mouse_asym_970_split.rds")

# merge the data
mouse_asym <- merge(mouse_asym_1,
                    y = c(mouse_asym_2,mouse_asym_3,mouse_asym_4,mouse_asym_5,
                          mouse_asym_6,mouse_asym_7,mouse_asym_8,mouse_asym_9,
                          mouse_asym_10,mouse_asym_11,mouse_asym_12),
                    add.cell.ids = c("B174", "B277", "B339", "B442", "B597", "B695","B768",
                                     "B814", "B917", "B918", "B969", "B970"),
                    project = "mouse_asymmetry"
)

#remove the individual to release the memory
rm(mouse_asym_1,mouse_asym_2,mouse_asym_3,mouse_asym_4,mouse_asym_5,mouse_asym_6,mouse_asym_7,
   mouse_asym_8,mouse_asym_9,mouse_asym_10,mouse_asym_11,mouse_asym_12)


#save the data
saveRDS(mouse_asym,file="intermedia_data/mouse_asym_acr_sample_merged.rds")

# 
mouse_asym <- JoinLayers(mouse_asym)
mouse_asym[["Xenium"]] <- split(mouse_asym[["Xenium"]], f = mouse_asym$sample_id)

#save the data
saveRDS(mouse_asym,file="intermedia_data/mouse_asym_acr_sample_merged.rds")

# filter the cell without any transcripts
mouse_asym <-  subset(mouse_asym, subset = nCount_Xenium > 0) # nCount > transcripts

#########################################
#----- Step2. do the standard processing
#########################################

#normalizes the data, detects high-variance features
mouse_asym <- mouse_asym %>% 
  SCTransform(assay = "Xenium") %>%
  RunPCA(npcs = 30, features = rownames(mouse_asym)) %>%
  RunUMAP(reduction="pca",dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.3)

# save the data
saveRDS(mouse_asym,file="intermedia_data/mouse_asym_acr_sample_orig.rds")
