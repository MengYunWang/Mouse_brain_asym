
# Remove all objects created before to prevent clashing
rm(list = ls())

# Set the working directory to the path where your files are located
# setwd("/Users/joeywang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Mouse_brain") # change it to the file directory
# setwd("/Users/wang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Mouse_brain/")
setwd("/data/workspaces/lag/workspaces/lg-func-asym/working_data/Mengyun/output/1_Gene_expression/All_genes/Gene_by_gene/region_seperated/batch1_2")

# install.packages("ggrain")  # if you haven’t already

library(ggrain)
library(ggplot2)
library(dplyr)
library(readxl)
library(openxlsx)

args <- commandArgs(trailingOnly=TRUE)
roi <- args[1]  # This is AUD HIP CA1 CA3 DG

# read data:
df <- read_excel("Hemi_data2analysis.xlsx", sheet = roi) # read the file in

df2 <- df %>%
  select(id, hemi, Crhr1) %>%
  mutate(
    Crhr1 = as.numeric(gsub(",", ".", Crhr1)),
    hemi   = factor(hemi, levels = c("left", "right"))
  )

# — plot data—
hemi_effect_crhr1 <- ggplot(df2, aes(x = hemi, y = Crhr1, fill = hemi)) +
  geom_rain(
    alpha        = .5,
    rain.side = 'f1x1',
    id.long.var  = "id",    # connect points by id
    boxplot.args = list(outlier.shape = NA)  # drop boxplot outliers
  ) +
  theme_classic() +
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = c("skyblue", "orange")) +
  guides(fill = "none", color = "none")+  
  theme(
    # plot.title = element_text(hjust = 0.5, size = 18, family = "Arial", face = "bold", colour = "black"),
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1,
                               size = 16, family = "Arial", face = "bold", colour = "black"),
    axis.text.y = element_text(size = 16, family = "Arial", face = "bold", colour = "black"),
    axis.title.y = element_text(size = 16, family = "Arial", face = "bold", colour = "black")
  )

ggsave(paste0("crhr1_hemi_effect", "_", roi, ".pdf"),
       plot = hemi_effect_crhr1, device = cairo_pdf, width = 6, height = 7, units = 'in')