# This code is to do the data curation

# Created 24-Sep-2024; updated 30-09-2024; updated 24-04-2025
# Created by M.-Y. WANG 

# Remove all objects created before to prevent clashing
rm(list = ls())

# Set the working directory to the path where your files are located
# setwd("/Users/joeywang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Mouse_brain") # change it to the file directory
# setwd("/Users/wang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Mouse_brain/")
setwd("/data/workspaces/lag/workspaces/lg-func-asym/working_data/Mengyun")

library(readxl)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(purrr)
library(parallel)
library(future)
plan("multicore", workers = 40) #

args <- commandArgs(trailingOnly=TRUE)
roi <- args[1]  # brain reigon: AUC, CA1, CA3, and DG

#################################Step1: Preprocessing

# define the mice id
mice_id <- list()
if (roi=="AUD") {
  mice_id <- c("M670", "M671", "M672", "M673", "M234", "M253", "M071", "M083",
               "M650", "M638", "M076", "M236", 
               "F679", "F680", "F682", "F683", "F687", "F688", "F078", "F087", 
               "F090")
} else {
  mice_id <- c("M669", "M670", "M671", "M672", "M674", "M676", "M677", "M678", 
               "M234", "M253", "M071", "M083", "M650", "M638", "M076", "M236", 
               "F679", "F680", "F682", "F683", "F685", "F686", "F687", "F688", 
               "F073", "F078", "F087", "F090")
}

# define the functions
read_data <- function (roi, hemi) {
  
  data_list <- lapply(mice_id, function(mice_id) {
    if (roi=="AUD") {
      data_loaded <- read_excel(paste0("output/2_Cell_type/acr_sample/batch1_2/AUD/", mice_id, "_AC_cell_type_harmony.xlsx"), sheet = hemi)
    } else {
      data_loaded <- read_excel(paste0("output/2_Cell_type/acr_sample/batch1_2/", roi, "/", mice_id, "_", roi, "_cell_type_harmony.xlsx"), sheet = hemi)
    }
    data_loaded <- t(data_loaded)
    colume_names <- data_loaded[1, ]
    data_loaded <- data_loaded[-1, ]
    data_loaded <- matrix(as.numeric(data_loaded), nrow = nrow(data_loaded), ncol = ncol(data_loaded))
    colnames(data_loaded) <- colume_names
    data_list <- colSums(data_loaded)
    return(data_list)
  })
  stacked_data <- do.call(rbind, data_list)
  return(stacked_data)
}

save_data <- function(roi, df, cell) {
  
  #save the data
  data2save <- df
  
  if (roi == "AUD") {
    data2save$hemi <- c(rep("left", 21), rep("right", 21))
    data2save$sex <-
      c(rep("male", 12),
        rep("female", 9),
        rep("male", 12),
        rep("female", 9))
    data2save$id <- rep(mice_id, 2)
    data2save <-
      data2save[, c('id', 'sex', 'hemi', setdiff(names(data2save), c('id', 'sex', 'hemi')))]
    
    write.xlsx(data2save, file = paste0("output/3_Cell_density/Cell_type/batch1_2/", cell, "_data2analysis.xlsx"), sheetName = roi)
  } else {
    data2save$hemi <- c(rep("left", 28), rep("right", 28))
    data2save$sex <-
      c(rep("male", 16),
        rep("female", 12),
        rep("male", 16),
        rep("female", 12))
    data2save$id <- rep(mice_id, 2)
    data2save <-
      data2save[, c('id', 'sex', 'hemi', setdiff(names(data2save), c('id', 'sex', 'hemi')))]
    
    wb <-loadWorkbook(paste0("output/3_Cell_density/Cell_type/batch1_2/", cell,"_data2analysis.xlsx"))
    addWorksheet(wb, roi)
    writeData(wb, roi, data2save)
    saveWorkbook(wb,paste0("output/3_Cell_density/Cell_type/batch1_2/", cell,"_data2analysis.xlsx"),overwrite = TRUE)
  }

}

read_save_data <- function (roi) {
  #load the cell counts
  df_cell_counts <- rbind(read_data(roi, "left"), read_data(roi, "right"))
  df_cell_counts <- as.data.frame(df_cell_counts)
  
  # exclude cell type which counts less than 30
  df_cell_count <- df_cell_counts[, colMeans(df_cell_counts, na.rm = TRUE) >= 30]

  
  # exclude cell type CA3 in CA1 and DG
  if (roi=="CA1" || roi=="DG") {
    df_cell_count <- subset(df_cell_count, select = -CA3)
  }
  
  df_cell_total <- rowSums(df_cell_count)
  # # compute the cell density
  # data_loaded <- read_excel(paste0("Mengyun_updated.xlsx"), sheet = "Cell_density")
  # df_data <- data_loaded[, (which(data_loaded[2, ] == "Area"))]
  # df_data <- df_data[c(-1:-2,-21:-26),]
  # colnames(df_data) <- c("AUD", "AUD", "DG", "DG", "CA1", "CA1", "CA3", "CA3")
  # area_size <- rbind(df_data[,c(1,3,5,7)], df_data[,c(2,4,6,8)])
  # area_size <- as.data.frame(area_size)
  # area <- as.numeric(na.omit(area_size[[roi]]))
  # # df_cell_density <- df_cell_count/area
  
  df_cell_density <- df_cell_count/df_cell_total
  
  #save the data
  save_data(roi, df_cell_counts, "counts_all")
  save_data(roi, df_cell_count, "counts_exclude_CA3_30")
  save_data(roi, df_cell_density, "Hemi")

}

 read_save_data(roi)




