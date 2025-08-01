# This code is going to do PCA analysis on gene expression

# Created 28-July-2025
# Created by M.-Y. WANG 

# Remove all objects created before to prevent clashing
rm(list = ls())

# Set the working directory to the path where your files are located
# setwd("/Users/joeywang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Mouse_brain") # change it to the file directory
# setwd("/Users/wang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Mouse_brain/")
setwd("/data/workspaces/lag/workspaces/lg-func-asym/working_data/Mengyun/output/1_Gene_expression/All_genes/PCA/batch1_2")

# install.packages("ggrain")  # if you haven’t already

library(ggrain)
library(ggplot2)
library(dplyr)
library(readxl)
library(openxlsx)

args <- commandArgs(trailingOnly=TRUE)
roi <- args[1]  # This is AUD HIP

# read data:
df <- read_excel("Hemi_data2analysis.xlsx", sheet = roi) # read the file in

###### for PC2
df2 <- df %>%
  select(id, hemi, PC2) %>%
  mutate(
    # PCs = as.numeric(gsub(",", ".", PC2)),
    PC2 = as.numeric(PC2)*-1,
    hemi   = factor(hemi, levels = c("left", "right"))
  )

# — plot data—
hemi_effect_PC2 <- ggplot(df2, aes(x = hemi, y = PC2, fill = hemi)) +
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

ggsave(paste0("PC2_hemi_effect", "_", roi, ".pdf"),
       plot = hemi_effect_PC2, device   = cairo_pdf, width = 6, height = 7)

## for pc9
df9 <- df %>%
  select(id, hemi, PC9) %>%
  mutate(
    # PCs = as.numeric(gsub(",", ".", PC9)),
    hemi   = factor(hemi, levels = c("left", "right"))
  )

# — plot data—
hemi_effect_PC9 <- ggplot(df9, aes(x = hemi, y = PC9, fill = hemi)) +
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

ggsave(paste0("PC9_hemi_effect", "_", roi, ".pdf"),
       plot = hemi_effect_PC9, device   = cairo_pdf, width = 6, height = 7)

## for PC20
df20 <- df %>%
  select(id, hemi, PC20) %>%
  mutate(
    # PCs = as.numeric(gsub(",", ".", PC20)),
    hemi   = factor(hemi, levels = c("left", "right"))
  )

# — plot data—
hemi_effect_PC20 <- ggplot(df20, aes(x = hemi, y = PC20, fill = hemi)) +
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

ggsave(paste0("PC20_hemi_effect", "_", roi, ".pdf"),
       plot = hemi_effect_PC20, device   = cairo_pdf, width = 6, height = 7)
