# This code is going to plot the gene by gene correlations

# Created 30-Aug-2024; updated 14-OCT-2024; updated 14-Mar-2025
# Created by M.-Y. WANG 

# Remove all objects created before to prevent clustering
rm(list = ls())

# Set the working directory to the path where your files are located
# setwd("/Users/joeywang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Mouse_brain") # change it to the file directory
# setwd("/Users/wang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Mouse_brain/")
setwd("/data/workspaces/lag/workspaces/lg-func-asym/working_data/Mengyun")

library(readxl)
library(openxlsx)
# library(MVTests)
# library(factoextra)
library(ggplot2)
library(dplyr)
#
args <- commandArgs(trailingOnly=TRUE)
roi <- args[1]  # This is AUD HIP

################################# Preprocessing
# read data
df_roi <- read_excel(paste0("Spatial_transcriptomics_batch1_batch2_Xenium_extractions_", roi, ".xlsx"), sheet = "gene_density") # read the file in
colnames(df_roi)[1] <- "Gene"

# delete the first 3 and 53th rows
df_roi <- df_roi[c(-3:-1, -53),]

# transpose the data
data2test <- t(df_roi) #transpose the data
newnames <- data2test[1,] # name the genes
data2test <- data2test[-1,] # delete the first row

# convert into numerical data:
original_dims <- dim(data2test)
data2test <- as.numeric(data2test)
dim(data2test) <- original_dims
colnames(data2test) <- newnames
# data2test_diff <- data2test[1:16,]-data2test[17:32,]
data2test <- data2test/rowSums(data2test)
data2test <- as.data.frame(data2test)

# do batch effect correction with ComBat
library(sva)
if (roi == "AUD") {
  batch <- factor(c(
    rep(1, 4), rep(2, 8), rep(1, 6), rep(2, 3),  # First 21 samples
    rep(1, 4), rep(2, 8), rep(1, 6), rep(2, 3)   # Next 21 samples
  ))
} else if (roi == "HIP") {
  batch <- factor(c(
    rep(1, 8), rep(2, 8), rep(1, 8), rep(2, 4),  # First 28 samples
    rep(1, 8), rep(2, 8), rep(1, 8), rep(2, 4)   # Next 28 samples
  ))
}

data2test_corrected <- ComBat(dat = t(data2test), batch = batch, par.prior = TRUE, prior.plots = FALSE) %>%
  t(.)%>%
  as.data.frame(.)


### plot the matrix
library(reshape2)
library(RColorBrewer)  # For color palettes
correlation_matrix <- cor(scale(data2test_corrected))


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
plots_qc_matrix <- ggplot(data = correlation_matrix_melt, aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = value)) +
  # geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 4.5) + # Add correlation coefficients
  scale_fill_gradientn(colors = colors, limits = c(-1, 1), breaks = c(-0.8, 0.8), labels = c("-0.8", "0.8"), name = "") +
  theme_minimal() +
  labs(x = "", y = "", title = "") +
  theme(
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
    legend.text = element_text(size = 38, family = "sans", face = "bold", colour = "black"),
    # plot.title = element_text(hjust = 0.5, size = 18, family = "Arial", face = "bold", colour = "black"),
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1,
                               size = 10, family = "sans", face = "bold", colour = "black"),
    axis.text.y = element_text(size = 10, family = "sans", face = "bold", colour = "black")
  ) +
  coord_fixed()


# Print the plot
# plots_qc_matrix
# ggsave(paste0("output/1_Gene_expression/All_genes/Gene_by_gene/region_seperated/batch1_2/Plot_gene_matrix_", roi, ".pdf"),
#        plot = plots_qc_matrix, width = 45, height = 45, units = 'in', dpi = 300, compression = "jpeg")

ggsave(paste0("output/1_Gene_expression/All_genes/Gene_by_gene/region_seperated/batch1_2/Plot_gene_matrix_", roi, ".png"),
       plot = plots_qc_matrix, width = 35, height = 35, units = 'in', dpi = 300)

# ggsave(paste0("output/1_Gene_expression/All_genes/Gene_by_gene/region_seperated/batch1_2/Plot_gene_matrix_", roi, ".pdf"),
#        plot = plots_qc_matrix, width = 35, height = 35, units = 'in')