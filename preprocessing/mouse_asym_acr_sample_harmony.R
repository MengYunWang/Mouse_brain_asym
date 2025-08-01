# This script is to do harmonization across the 12 samples

# Created by M.-Y. WANG
# 31-Aug-2024

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
options(future.globals.maxSize = 300000 * 1024^2, future.seed = TRUE) 

###############-------------------Step1: Load data
mouse_asym <- readRDS(file="intermedia_data/mouse_asym_acr_sample_orig.rds")

###############-------------------Step2. Integration and do the standard again
#in integrate the data with HarmonyIntegration
mouse_asym_acr_sample_harmony <- mouse_asym %>%
  IntegrateLayers(method = HarmonyIntegration,
                  normalization.method = "SCT",
                  orig.reduction = "pca",
                  new.reduction = "harmony",
                  verbose = TRUE
  ) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.3)

## save the data
saveRDS(mouse_asym_acr_sample_harmony,file="intermedia_data/mouse_asym_acr_sample_harmony.rds")


