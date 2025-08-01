# This code is going to do PCA analysis on gene expression

# Created 30-Aug-2024; updated 09-09-2024
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
library(factoextra)
library(ggplot2)
library(dplyr)
#
args <- commandArgs(trailingOnly=TRUE)
roi <- args[1]  # values are AUD HIP

################
# Preprocessing
################
# read data
df_roi <- read_excel(paste0("Spatial_transcriptomics_batch1_batch2_Xenium_extractions_", roi, ".xlsx"), sheet = "gene_density") # read the file in
colnames(df_roi)[1] <- "Gene"

# delete the first 3 and 53th rows
df_roi <- df_roi[c(-3:-1, -53),]

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
  data2save$hemi <- c(rep("left", 21), rep("right", 21))
  data2save$sex <- c(rep("male", 12), rep("female", 9), rep("male", 12), rep("female", 9))
  data2save$id <- rep(c("M670", "M671", "M672", "M673", "M234", "M253", "M071", "M083", "M650", "M638", "M076", "M236",
                        "F679", "F680", "F682", "F683", "F687", "F688", "F078", "F087", "F090"), 2)
  data2save <- data2save[, c('id','sex','hemi',setdiff(names(data2save), c('id','sex','hemi')))]
  
} else {
  data2save$hemi <- c(rep("left", 28), rep("right", 28))
  data2save$sex <- c(rep("male", 16), rep("female", 12), rep("male", 16), rep("female", 12))
  data2save$id <- rep(c("M669", "M670", "M671", "M672", "M674", "M676", "M677", "M678", "M234", "M253", "M071", "M083", "M650", "M638", "M076", "M236",
                        "F679", "F680", "F682", "F683", "F685", "F686", "F687", "F688", "F073", "F078", "F087", "F090"), 2)
  data2save <- data2save[, c('id','sex','hemi',setdiff(names(data2save), c('id','sex','hemi')))]
}

# if (roi=="AUD") {
#   write.xlsx(data2save, file = "output/1_Gene_expression/All_genes/PCA/batch1_2/Hemi_data2analysis.xlsx", sheetName = roi)
# } else {
#   wb <- loadWorkbook("output/1_Gene_expression/All_genes/PCA/batch1_2/Hemi_data2analysis.xlsx")
#   addWorksheet(wb, roi)
#   writeData(wb, roi, data2save)
#   saveWorkbook(wb, "output/1_Gene_expression/All_genes/PCA/batch1_2/Hemi_data2analysis.xlsx", overwrite=TRUE)
# }


############################ 
# Analysis
############################
# do the multi paired test between left and right hemisphere
# mpair_roi <- Mpaired(T1=data2analysis[1:(nrow(data2analysis)/2),1:num_cumu_var_90], 
#                      T2=data2analysis[(nrow(data2analysis)/2+1):nrow(data2analysis),1:num_cumu_var_90])
# summary(mpair_roi)


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


#### Perform paired t-tests on left and right hemi
paired_roi <- lapply(1:num_cumu_var_90, function(i) {
  result <- t.test(data2analysis[1:(nrow(data2analysis)/2),i], # first half left
                   data2analysis[(nrow(data2analysis)/2+1):nrow(data2analysis),i], #second half right
                   paired=TRUE)
  return(result)
})


# # Write the results to an Excel file
# 
# if (roi=="AUD") {
#   write.xlsx(test_update(paired_roi), "output/1_Gene_expression/All_genes/PCA/batch1_2/Hemi_PCA_paired_test.xlsx", sheetName = roi)
# } else {
#   wb <- loadWorkbook("output/1_Gene_expression/All_genes/PCA/batch1_2/Hemi_PCA_paired_test.xlsx")
#   addWorksheet(wb, roi)
#   writeData(wb, roi, test_update(paired_roi))
#   saveWorkbook(wb, "output/1_Gene_expression/All_genes/PCA/batch1_2/Hemi_PCA_paired_test.xlsx", overwrite=TRUE)
# }

# ### Perform two sample t-tests on sex
# if (roi=="AUD") {
#   Two_sample_roi <- lapply(1:num_cumu_var_90, function(i) {
#     result <- t.test(((data2analysis[1:5,i] + data2analysis[14:18,i])/2),
#                      ((data2analysis[6:13,i] + data2analysis[19:26,i])/2),
#                      paired=FALSE)
#     return(result)
#   })
# } else {
#   Two_sample_roi <- lapply(1:num_cumu_var_90, function(i) {
#     result <- t.test(((data2analysis[1:8,i] + data2analysis[17:24,i])/2), #male
#                      ((data2analysis[9:16,i] + data2analysis[25:32,i])/2),  #female
#                      paired=FALSE)
#     return(result)
#   })
# }
# 
# 
# # Write the results to an Excel file
# if (roi=="CA1") {
#   write.xlsx(test_update(Two_sample_roi), "output/1_Gene_expression/All_genes/PCA/batch1_2/PCA_2sample_test_sex.xlsx", sheetName = roi)
# } else {
#   wb <- loadWorkbook("output/1_Gene_expression/All_genes/PCA/batch1_2/PCA_2sample_test_sex.xlsx")
#   addWorksheet(wb, roi)
#   writeData(wb, roi, test_update(Two_sample_roi))
#   saveWorkbook(wb, "output/1_Gene_expression/All_genes/PCA/batch1_2/PCA_2sample_test_sex.xlsx", overwrite=TRUE)
# }


###########################
# some plots about the PCs
###########################

Plot_PC_variance <- fviz_eig(PCA_roi, addlabels = TRUE, main = "", xlab = "Principle Component") + 
  theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
  )

# pca plot for individuals
if (roi=="AUD") {
  labels <- as.character(rep(1:21, times = 2))
  hemispheres <- rep(c("Left", "Right"), each = 21)
} else {
  labels <- as.character(rep(1:28, times = 2))
  hemispheres <- rep(c("Left", "Right"), each = 28)}

Plot_PC_sample <- fviz_pca_ind(PCA_roi,
                               habillage = hemispheres, # Color by the left/right hemisphere
                               col.ind = "cos2", # Color by the quality of representation
                               palette = c("#00AFBB", "#E7B800"),
                               repel = TRUE,     # Avoid text overlapping
                               xlab = "PC1",
                               ylab = "PC2",
                               title = "",
                               legend.title=""
                               ) +  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

Plot_PC_sample_batch <- fviz_pca_ind(PCA_roi_batch,
                               habillage = hemispheres, # Color by the left/right hemisphere
                               col.ind = "cos2", # Color by the quality of representation
                               palette = c("#00AFBB", "#E7B800"),
                               repel = TRUE,     # Avoid text overlapping
                               xlab = "PC1",
                               ylab = "PC2",
                               title = "",
                               legend.title="") + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
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
Plot_PC1_contribution <- fviz_contrib(PCA_roi, choice = "var", axes = 1, top = 15, title="PC1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Plot_PC2_contribution <- fviz_contrib(PCA_roi, choice = "var", axes = 2, top = 15, title="PC2") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Plot_PC3_contribution <- fviz_contrib(PCA_roi, choice = "var", axes = 3, top = 15, title="PC3") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Plot_PC4_contribution <- fviz_contrib(PCA_roi, choice = "var", axes = 4, top = 15, title="PC4") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Plot_PC5_contribution <- fviz_contrib(PCA_roi, choice = "var", axes = 9, top = 15, title="PC9") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Plot_PC6_contribution <- fviz_contrib(PCA_roi, choice = "var", axes = 20, top = 15, title="PC20") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

library(patchwork)

Plot_roi_batch_effect<- Plot_PC_sample_batch | Plot_PC_sample

ggsave(paste0("output/1_Gene_expression/All_genes/PCA/batch1_2/Hemi_Plot_PCA_", roi, "_batch_effect.png"),
       plot = Plot_roi_batch_effect, width = 10, height = 5, units = 'in', dpi = 300)

Plot_roi<- (Plot_PC_variance  / Plot_PC_sample) |
  (Plot_PC1_contribution / Plot_PC2_contribution /
     Plot_PC3_contribution / Plot_PC4_contribution / 
     Plot_PC5_contribution / Plot_PC6_contribution) 
ggsave(paste0("output/1_Gene_expression/All_genes/PCA/batch1_2/Hemi_Plot_PCA_", roi, ".png"),
       plot = Plot_roi, width = 10, height = 9, units = 'in', dpi = 300)


### save the PC 2 4 6 and 8 plots 
# Define a function to generate and save plots
save_contribution_plot <- function(pca, axis_num, top_contrib, roi, pc_name) {
  plot <- fviz_contrib(pca, choice = "var", axes = axis_num, top = top_contrib, title = pc_name) +
    theme(
      axis.title.y = element_text(size = 14, face = "bold", family = "Arial"),
      axis.text.x = element_text(size = 10, face = "bold", family = "Arial"),
      axis.text.y = element_text(size = 10, face = "bold", family = "Arial"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ggsave(paste0("output/1_Gene_expression/All_genes/PCA/batch1_2/Hemi_Plot_", pc_name, "_", roi, ".pdf"),
         plot = plot, device   = cairo_pdf, width = 14, height = 7)
}

# Call the function for each PC
save_contribution_plot(PCA_roi, 1, 15, roi, "PC1")
save_contribution_plot(PCA_roi, 2, 15, roi, "PC2")
save_contribution_plot(PCA_roi, 9, 15, roi, "PC9")
save_contribution_plot(PCA_roi, 20, 15, roi, "PC20")


roi.gene <- get_pca_var(PCA_roi)
roi.gene$coord          # Coordinates
roi.gene$contrib        # Contributions to the PCs
roi.gene$cos2.          # Quality of representation 

# Results for individuals
roi.ind <- get_pca_ind(PCA_roi)
roi.ind$coord          # Coordinates
roi.ind$contrib        # Contributions to the PCs
roi.ind$cos2           # Quality of representation 

