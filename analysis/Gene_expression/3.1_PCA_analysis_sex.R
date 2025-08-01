# This code is going to do PCA analysis sex effect on gene expression

# Created 30-Aug-2024; updated 09-09-2024; updated 01-April-2025
# Created by M.-Y. WANG 

# Remove all objects created before to prevent clustering
rm(list = ls())

# Set the working directory to the path where your files are located
# setwd("/Users/joeywang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Mouse_brain") # change it to the file directory
# setwd("/Users/wang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Mouse_brain/")
setwd("/data/workspaces/lag/workspaces/lg-func-asym/working_data/Mengyun")


library(readxl)
library(openxlsx)
library(factoextra)
library(ggplot2)
library(dplyr)
#
args <- commandArgs(trailingOnly=TRUE)
roi <- args[1]  # values are AUD HIP

################
# Preprocessing
################

df_roi <- read_excel(paste0("Spatial_transcriptomics_batch1_batch2_Xenium_extractions_", roi, ".xlsx"), sheet = "L+R_gene_densities") # read the file in
colnames(df_roi)[1] <- "Gene"

# delete the first 4 and 54th rows
df_roi <- df_roi[c(-4:-1, -54),]

# transpose the data
data2test <- t(df_roi) #transpose the data
newnames <- data2test[1,] # name the genes
data2test <- data2test[-1,] # delete the first row

# convert into numerical data
original_dims <- dim(data2test)
data2test <- as.numeric(data2test)
dim(data2test) <- original_dims
colnames(data2test) <- newnames

# normalization: each gene density divided by the total gene density within each mouse
data2test <- data2test/rowSums(data2test)
data2test <- as.data.frame(data2test)

# do batch effect correction with ComBat
library(sva)
if (roi == "AUD") {
  batch <- factor(c(
    rep(1, 4), rep(2, 8), rep(1, 6), rep(2, 3)   # Next 21 samples
  ))
} else if (roi == "HIP") {
  batch <- factor(c(
    rep(1, 8), rep(2, 8), rep(1, 8), rep(2, 4)   # Next 28 samples
  ))
}

data2test_corrected <- ComBat(dat = t(data2test), batch = batch, par.prior = TRUE, prior.plots = FALSE) %>%
  t(.)%>%
  as.data.frame(.)


##############
# compute PCA
##############
PCA_roi <- prcomp(data2test_corrected, scale = TRUE)
PCA_roi_batch <- prcomp(data2test, scale = TRUE)

# select the PC 
roi_eg <- get_eigenvalue(PCA_roi) # same thing to PCA_roi$sdev^2
cumu_var_90 <- roi_eg$cumulative.variance.percent[roi_eg$cumulative.variance.percent < 90] # select only the PCs accumulated till 90%
eigenvalues <- roi_eg$eigenvalue[roi_eg$cumulative.variance.percent < 90]
num_cumu_var_90 <- length(cumu_var_90)
# selected_pc <- roi_eg[roi_eg$variance.percent > 5, ]
# eigenvalues <- selected_pc$eigenvalue
# num_selected <- nrow(selected_pc)
# variance_explained <- selected_pc$variance.percent

# project data into PC space
data2analysis <- PCA_roi[["x"]] %>% #same as scale(data2test_corrected)%*%PCA_roi$rotation
  as.data.frame(.)

##############
#save the data
##############
data2save <- data2analysis

if (roi=="AUD") {
  data2save$sex <- c(rep("male", 12), rep("female", 9))
  data2save$id <- c("M670", "M671", "M672", "M673", "M234", "M253", "M071", "M083", "M650", "M638", "M076", "M236",
                        "F679", "F680", "F682", "F683", "F687", "F688", "F078", "F087", "F090")
  data2save <- data2save[, c('id','sex', setdiff(names(data2save), c('id','sex')))]
  
} else {
  data2save$sex <- c(rep("male", 16), rep("female", 12))
  data2save$id <- c("M669", "M670", "M671", "M672", "M674", "M676", "M677", "M678", "M234", "M253", "M071", "M083", "M650", "M638", "M076", "M236",
                        "F679", "F680", "F682", "F683", "F685", "F686", "F687", "F688", "F073", "F078", "F087", "F090")
  data2save <- data2save[, c('id','sex', setdiff(names(data2save), c('id','sex')))]
}

if (roi=="AUD") {
  write.xlsx(data2save, file = "output/1_Gene_expression/All_genes/PCA/batch1_2/Sex_data2analysis.xlsx", sheetName = roi)
} else {
  wb <- loadWorkbook("output/1_Gene_expression/All_genes/PCA/batch1_2/Sex_data2analysis.xlsx")
  if (roi %in% names(wb)) {
    removeWorksheet(wb, roi)  # Remove the existing sheet
  }
  addWorksheet(wb, roi)
  writeData(wb, roi, data2save)
  saveWorkbook(wb, "output/1_Gene_expression/All_genes/PCA/batch1_2/Sex_data2analysis.xlsx", overwrite=TRUE)
}

############################ 
# Analysis
############################

# function to tyde up
test_update <- function (data2updated) {
  p_values <- sapply(data2updated, function(x)
    x$p.value)
  p_adjusted <- p.adjust(p_values[1:num_cumu_var_90], method = "BH")
  data2updated <- mapply(function(x, y) {
    x$p.adjusted <- y
    return(x)
  }, data2updated, p_adjusted, SIMPLIFY = FALSE)
  
  # Create a data frame with the results
  PCA_t_test <- data.frame(
    t = sapply(data2updated, function(x)
      x$statistic),
    df = sapply(data2updated, function(x)
      x$parameter),
    P.Value = sapply(data2updated, function(x)
      x$p.value),
    ajust.P = sapply(data2updated, function(x)
      x$p.adjusted),
    CI.Low = sapply(data2updated, function(x)
      x$conf.int[1]),
    CI.High = sapply(data2updated, function(x)
      x$conf.int[2]),
    eigenvalue = eigenvalues,
    accumulated_variance = cumu_var_90
    # variance_explain = variance_explained
  )
  
  # Add the PC column
  PCA_t_test$PC <- paste0("PC", 1:num_cumu_var_90) # correspond to the number of PCs
  PCA_t_test <- PCA_t_test[, c(
    "PC",
    "t",
    "df",
    "P.Value",
    "ajust.P",
    "CI.Low",
    "CI.High",
    "eigenvalue",
    "accumulated_variance"
    # "variance_explain"
  )]
  return(PCA_t_test)
}

### Perform two sample t-tests on sex
if (roi=="AUD") {
  Two_sample_roi <- lapply(1:num_cumu_var_90, function(i) {
    result <- t.test(data2analysis[1:12,i], data2analysis[13:21,i],
                     paired=FALSE)
    return(result)
  })
} else {
  Two_sample_roi <- lapply(1:num_cumu_var_90, function(i) {
    result <- t.test(data2analysis[1:16,i], data2analysis[17:28,i],  #female
                     paired=FALSE)
    return(result)
  })
}


# Write the results to an Excel file
if (roi=="AUD") {
  write.xlsx(test_update(Two_sample_roi), "output/1_Gene_expression/All_genes/PCA/batch1_2/Sex_PCA_2sample_test.xlsx", sheetName = roi)
} else {
  wb <- loadWorkbook("output/1_Gene_expression/All_genes/PCA/batch1_2/Sex_PCA_2sample_test.xlsx")
  if (roi %in% names(wb)) {
    removeWorksheet(wb, roi)  # Remove the existing sheet
  }
  addWorksheet(wb, roi)
  writeData(wb, roi, test_update(Two_sample_roi))
  saveWorkbook(wb, "output/1_Gene_expression/All_genes/PCA/batch1_2/Sex_PCA_2sample_test.xlsx", overwrite=TRUE)
}


###########################
# some plots about the PCs
###########################

Plot_PC_variance <- fviz_eig(PCA_roi, addlabels = TRUE, main = "", xlab = "Principle Component")

# pca plot for individuals
if (roi=="AUD") {
  labels <- as.character(rep(1:21, times = 2))
  sex <- c(rep("male", 12), rep("female", 9))
} else {
  labels <- as.character(rep(1:28, times = 2))
  sex <- c(rep("male", 16), rep("female", 12))}

Plot_PC_sample <- fviz_pca_ind(PCA_roi,
                               habillage = sex, # Color by the sex
                               col.ind = "cos2", # Color by the quality of representation
                               palette = c("#00AFBB", "#E7B800"),
                               repel = TRUE,     # Avoid text overlapping
                               xlab = "PC1",
                               ylab = "PC2",
                               title = "",
                               legend.title=""
)

Plot_PC_sample_batch <- fviz_pca_ind(PCA_roi_batch,
                               habillage = sex, # Color by the sex
                               col.ind = "cos2", # Color by the quality of representation
                               palette = c("#00AFBB", "#E7B800"),
                               repel = TRUE,     # Avoid text overlapping
                               xlab = "PC1",
                               ylab = "PC2",
                               title = "",
                               legend.title=""
)
# fviz_pca_ind(PCA_roi,
#              col.ind = "cos2", # Color by the quality of representation
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE,     # Avoid text overlapping
#              xlab = "PC1",
#              ylab = "PC2",
#              title = "",
#              legend.title=""
# )
# 
# # pca plot for genes
# fviz_pca_var(PCA_roi,
#              col.var = "contrib", # Color by contributions to the PC
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE     # Avoid text overlapping
# )
# # combined plot
# fviz_pca_biplot(PCA_roi, repel = TRUE,
#                 col.var = "#2E9FDF", # Variables color
#                 col.ind = "#696969",  # Individuals color
#                 xlab = "PC1",
#                 ylab = "PC2"
# )
# Contributions of variables to PC1 to 6
Plot_PC1_contribution <- fviz_contrib(PCA_roi, choice = "var", axes = 1, top = 30, title="PC1")
Plot_PC2_contribution <- fviz_contrib(PCA_roi, choice = "var", axes = 2, top = 30, title="PC2")
Plot_PC3_contribution <- fviz_contrib(PCA_roi, choice = "var", axes = 3, top = 30, title="PC3")
Plot_PC4_contribution <- fviz_contrib(PCA_roi, choice = "var", axes = 4, top = 30, title="PC4")
Plot_PC5_contribution <- fviz_contrib(PCA_roi, choice = "var", axes = 5, top = 30, title="PC5")
Plot_PC6_contribution <- fviz_contrib(PCA_roi, choice = "var", axes = 6, top = 30, title="PC6")

library(patchwork)

Plot_roi_batch_effect<- Plot_PC_sample_batch | Plot_PC_sample

ggsave(paste0("output/1_Gene_expression/All_genes/PCA/batch1_2/Sex_Plot_PCA_", roi, "_batch_effect.png"),
       plot = Plot_roi_batch_effect, width = 10, height = 5, units = 'in', dpi = 300)

Plot_roi<- (Plot_PC_variance  / Plot_PC_sample) |
  (Plot_PC1_contribution / Plot_PC2_contribution /
     Plot_PC3_contribution / Plot_PC4_contribution / 
     Plot_PC5_contribution / Plot_PC6_contribution) 
ggsave(paste0("output/1_Gene_expression/All_genes/PCA/batch1_2/Sex_Plot_PCA_", roi, ".png"),
       plot = Plot_roi, width = 10, height = 9, units = 'in', dpi = 300)


### save the PC 2 4 6 and 8 plots 
# Define a function to generate and save plots
save_contribution_plot <- function(pca, axis_num, top_contrib, roi, pc_name) {
  plot <- fviz_contrib(pca, choice = "var", axes = axis_num, top = top_contrib, title = pc_name) +
    theme(
      axis.title.y = element_text(size = 14, face = "bold", family = "Arial"),
      axis.text.x = element_text(size = 10, face = "bold", family = "Arial"),
      axis.text.y = element_text(size = 10, face = "bold", family = "Arial")
    )
  
  ggsave(paste0("output/1_Gene_expression/All_genes/PCA/batch1_2/Sex_Plot_", pc_name, "_", roi, ".png"),
         plot = plot, width = 14, height = 7, units = 'in', dpi = 300)
}

# Call the function for each PC
save_contribution_plot(PCA_roi, 1, 90, roi, "PC1")
save_contribution_plot(PCA_roi, 2, 90, roi, "PC2")
save_contribution_plot(PCA_roi, 3, 90, roi, "PC3")
save_contribution_plot(PCA_roi, 4, 90, roi, "PC4")


roi.gene <- get_pca_var(PCA_roi)
roi.gene$coord          # Coordinates
roi.gene$contrib        # Contributions to the PCs
roi.gene$cos2.          # Quality of representation 

# Results for individuals
roi.ind <- get_pca_ind(PCA_roi)
roi.ind$coord          # Coordinates
roi.ind$contrib        # Contributions to the PCs
roi.ind$cos2           # Quality of representation 

