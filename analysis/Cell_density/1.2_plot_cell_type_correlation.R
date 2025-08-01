# This code is plot the cell type correlation matrix

# Created 26-Sep-2024; updated 30-09-2024
# Created by M.-Y. WANG 

# Remove all objects created before to prevent clashing
rm(list = ls())

# Set the working directory to the path where your files are located
# setwd("/Users/joeywang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Mouse_brain") # change it to the file directory
# setwd("/Users/wang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Mouse_brain/")
setwd("/data/workspaces/lag/workspaces/lg-func-asym/working_data/Mengyun/output/3_Cell_density/Cell_type/batch1_2/")

library(readxl)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(purrr)
library(parallel)
args <- commandArgs(trailingOnly=TRUE)
roi <- args[1]  # This is AC CA1 CA3 DG HIP_EX HIP_IN

######--------------------------------------------Step1: load the data
# read data
data2test <- read_excel("Hemi_data2analysis.xlsx", sheet = roi) # read the file in

data2test <- na.omit(data2test[,-3:-1]) # delete the first three columns

# # convert into numerical data:
# original_dims <- dim(data2test)
# data2test <- as.numeric(data2test)
# dim(data2test) <- original_dims
# colnames(data2test) <- newnames
# data2test <- as.data.frame(data2test)

######--------------------------------------------Step 4. Plot the figures
library(reshape2)
library(RColorBrewer)  # For color palettes

### plot the matrix
# data2test <- na.omit(data2test)
correlation_matrix <- cor(scale(data2test))

# Perform hierarchical clustering
hc <- hclust(dist(correlation_matrix))
hc_order <- hc$order

# Reorder the correlation matrix based on the hierarchical clustering
correlation_matrix_ordered <- correlation_matrix[hc_order, hc_order]

# Melt the correlation matrix for visualization
correlation_matrix_melt <- melt(correlation_matrix_ordered)

# Define a reversed RdBu color palette
colors <- rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(100))

# Plotting with ggplot2
plots_qc_matrix <-
  ggplot(data = correlation_matrix_melt, aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = value)) +
  # geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 4.5) + # Add correlation coefficients
  scale_fill_gradientn(
    colors = colors,
    limits = c(-1, 1),
    breaks = c(-0.5, 0, 0.5),
    labels = c("-0.5", "0", "0.5"),
    name = ""
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "") +
  theme(
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 10,
      l = 10,
      unit = "pt"
    ),
    legend.text = element_text(
      size = 12,
      family = "sans",
      face = "bold",
      colour = "black"
    ),
    # plot.title = element_text(hjust = 0.5, size = 18, family = "Arial", face = "bold", colour = "black"),
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1,
      size = 10,
      family = "sans",
      face = "bold",
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 10,
      family = "sans",
      face = "bold",
      colour = "black"
    )
  ) +
  coord_fixed()

ggsave(
  paste0("Plot_cell_density_matrix_", roi, ".png"),
  plot = plots_qc_matrix,
  width = 8,
  height = 8,
  units = 'in',
  dpi = 300
)

