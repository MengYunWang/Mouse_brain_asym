# This code is going to do gene expression data curation and save it into data2analysis_*.xlsx

# Created 30-Aug-2024; updated 14-OCT-2024; updated 14-Mar-2025
# Created by M.-Y. WANG 

# Remove all objects created before to prevent clashing
rm(list = ls())

# Set the working directory to the path where your files are located
# setwd("/Users/joeywang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Mouse_brain") # change it to the file directory
# setwd("/Users/wang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Mouse_brain/")
setwd("/data/workspaces/lag/workspaces/lg-func-asym/working_data/Mengyun")


library(readxl)
library(openxlsx)
library(dplyr)
library(purrr)
args <- commandArgs(trailingOnly=TRUE)
roi <- args[1]  # This is AUD HIP CA1 CA3 DG

########################
###### hemi data #######
########################

# read data
# df_roi <- read_excel("Spatial_transcriptomics_batch1_AUD_HIP.xlsx", sheet = roi) # read the file in
df_roi <- read_excel(paste0("Spatial_transcriptomics_batch1_batch2_Xenium_extractions_", roi, ".xlsx"), sheet = "gene_density") # read the file in
colnames(df_roi)[1] <- "Gene"

# delete the first 3 and 53th rows
df_roi <- df_roi[c(-3:-1, -53),]

# transpose the data
data2analysis <- t(df_roi) #transpose the data
newnames <- data2analysis[1,] # name the genes
data2analysis <- data2analysis[-1,] # delete the first row

# convert into numerical data
original_dims <- dim(data2analysis)
data2analysis <- as.numeric(data2analysis)
dim(data2analysis) <- original_dims
colnames(data2analysis) <- newnames

# normalization: each gene density divided by the total gene density within each mouse
data2analysis <- data2analysis/rowSums(data2analysis)
data2analysis <- as.data.frame(data2analysis)

# save the data
data2save <- data2analysis

if (roi=="AUD") {
  data2save$hemi <- c(rep("left", 21), rep("right", 21))
  data2save$sex <- c(rep("male", 12), rep("female", 9), rep("male", 12), rep("female", 9))
  data2save$id <- rep(c("M670", "M671", "M672", "M673", "M234", "M253", "M071", "M083", "M650", "M638", "M076", "M236",
                        "F679", "F680", "F682", "F683", "F687", "F688", "F078", "F087", "F090"), 2) 
  data2save <- data2save[, c('id','sex','hemi',setdiff(names(data2save), c('id','sex','hemi')))]
  
} else {
  data2save$hemi <- c(rep("left", 28), rep("right", 28))
  data2save$sex <- c(rep("male", 16), rep("female", 12), rep("male", 16), rep("female", 12))
  data2save$id <- rep(c("M669", "M670", "M671", "M672", "M674", "M676", "M677", "M678", "M234", "M253", "M071", "M083", "M650", "M638", "M076", "M236",
                        "F679", "F680", "F682", "F683", "F685", "F686", "F687", "F688", "F073", "F078", "F087", "F090"), 2) # pay attention to the order compared to sex data
  data2save <- data2save[, c('id','sex','hemi',setdiff(names(data2save), c('id','sex','hemi')))]
}

if (roi=="AUD") {
  write.xlsx(data2save, file = "output/1_Gene_expression/All_genes/Gene_by_gene/region_seperated/batch1_2/Hemi_data2analysis.xlsx", sheetName = roi)
} else {
  wb <- loadWorkbook("output/1_Gene_expression/All_genes/Gene_by_gene/region_seperated/batch1_2/Hemi_data2analysis.xlsx")
  if (roi %in% names(wb)) {
    removeWorksheet(wb, roi)  # Remove the existing sheet
  }
  addWorksheet(wb, roi)
  writeData(wb, roi, data2save)
  saveWorkbook(wb, "output/1_Gene_expression/All_genes/Gene_by_gene/region_seperated/batch1_2/Hemi_data2analysis.xlsx", overwrite=TRUE)
}
