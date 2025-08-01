# This code is to do the data curation

# Created 24-Sep-2024; updated 30-09-2024; updated 25-04-2025
# Created by M.-Y. WANG 


# Remove all objects created before to prevent clashing
rm(list = ls())

# Set the working directory to the path where your files are located
# setwd("/Users/joeywang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Mouse_brain") # change it to the file directory
# setwd("/Users/wang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Mouse_brain/")
setwd("/data/workspaces/lag/workspaces/lg-func-asym/working_data/Mengyun/output/3_Cell_density/Cell_type/batch1_2")

library(readxl)
library(openxlsx)

# compute the counts of the hip region
CA1_data <- read_excel("counts_all_data2analysis.xlsx",sheet="CA1")
CA3_data <- read_excel("counts_all_data2analysis.xlsx",sheet="CA3")
DG_data <- read_excel("counts_all_data2analysis.xlsx",sheet="DG")

################################# add them all up
HIP_cell_counts <- CA1_data[, 4:ncol(CA1_data)] + CA3_data[, 4:ncol(CA3_data)] + DG_data[, 4:ncol(DG_data)]

# save the data
data2save_counts <- HIP_cell_counts
data2save_counts$hemi <- CA1_data$hemi
data2save_counts$sex <- CA1_data$sex
data2save_counts$id <- CA1_data$id
data2save_counts <- data2save_counts[, c('id','sex','hemi',setdiff(names(data2save_counts), c('id','sex','hemi')))]

wb <- loadWorkbook("counts_all_data2analysis.xlsx")
if ("HIP_combined" %in% names(wb)) {
  removeWorksheet(wb, "HIP_combined")  # Remove the existing sheet
}
addWorksheet(wb, "HIP_combined")
writeData(wb, "HIP_combined", data2save_counts)
saveWorkbook(wb, "counts_all_data2analysis.xlsx", overwrite=TRUE)


######################## excluding CA3 in CA1 and DG, and cell cunt less than 30
CA1_data$CA3 <- 0
DG_data$CA3 <- 0
# exclude cell type which counts less than 30
HIP_cell_count <- HIP_cell_counts[, colMeans(HIP_cell_counts, na.rm = TRUE) >= 30]

# save the data
data2save_counts <- HIP_cell_count
data2save_counts$hemi <- CA1_data$hemi
data2save_counts$sex <- CA1_data$sex
data2save_counts$id <- CA1_data$id
data2save_counts <- data2save_counts[, c('id','sex','hemi',setdiff(names(data2save_counts), c('id','sex','hemi')))]

wb <- loadWorkbook("counts_exclude_CA3_30_data2analysis.xlsx")
if ("HIP_combined" %in% names(wb)) {
  removeWorksheet(wb, "HIP_combined")  # Remove the existing sheet
}
addWorksheet(wb, "HIP_combined")
writeData(wb, "HIP_combined", data2save_counts)
saveWorkbook(wb, "counts_exclude_CA3_30_data2analysis.xlsx", overwrite=TRUE)


################################################ compute the cell density of HIP
HIP_cell_total <- rowSums(HIP_cell_count)

HIP_cell_density <- HIP_cell_count/HIP_cell_total

data2analysis <- as.data.frame(HIP_cell_density)
data2analysis$hemi <- CA1_data$hemi
data2analysis$sex <- CA1_data$sex
data2analysis$id <- CA1_data$id
data2analysis <- data2analysis[, c('id','sex','hemi',setdiff(names(data2analysis), c('id','sex','hemi')))]

wb <- loadWorkbook("Hemi_data2analysis.xlsx")
if ("HIP_combined" %in% names(wb)) {
  removeWorksheet(wb, "HIP_combined")  # Remove the existing sheet
}
addWorksheet(wb, "HIP_combined")
writeData(wb, "HIP_combined", data2analysis)
saveWorkbook(wb, "Hemi_data2analysis.xlsx", overwrite=TRUE)






