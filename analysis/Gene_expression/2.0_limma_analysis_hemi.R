# This code is going to do limma on hemi effect

# Created 07-OCT-2024; updated 14-OCT-2024; updated 14-Mar-2025
# Created by M.-Y. WANG 

# Remove all objects created before to prevent clashing
rm(list = ls())

# Set the working directory to the path where your files are located
# setwd("/Users/joeywang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Mouse_brain") # change it to the file directory
# setwd("/Users/wang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Mouse_brain/")
setwd("/data/workspaces/lag/workspaces/lg-func-asym/working_data/Mengyun")

# loading libraries
library(readxl)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(purrr)
library(parallel)
library(limma)
library(future)
plan("multicore", workers = 40) 

args <- commandArgs(trailingOnly=TRUE)
roi <- args[1]  # This is AUD HIP CA1 CA3 DG
######------------------------------ load the data

# define a function to load the data
load_data <- function(roi){
  loaded_data <- read_excel("output/1_Gene_expression/All_genes/Gene_by_gene/region_seperated/batch1_2/Hemi_data2analysis.xlsx", sheet = roi) # read the file in
  data_loaded <- t(loaded_data[,-2:-3]) # get rid of the sex and hemi column
  colsnames <- data_loaded[1,] # extract the first row of the mouse ids, such as M670 ...
  data2analysis <- data_loaded[-1,] # get rid of the first row

  if (roi=="AUD") {
    new_colsnames <- c(paste0("left_", colsnames[1:21]), paste0("right_", colsnames[22:42]))
  } else {
    new_colsnames <- c(paste0("left_", colsnames[1:28]), paste0("right_", colsnames[29:56]))
  }

  colnames(data2analysis) <- new_colsnames
  data2analysis <- as.data.frame(data2analysis)
  names_row <- rownames(data2analysis)

  data2analysis <- as.matrix(apply(data2analysis, 2, as.numeric))
  rownames(data2analysis) <- names_row
  return(data2analysis)
}

#######----------------------------- Permutation of the tests
# define a function to create tibble version of dataframe
creat_shuffle_data <- function(data, roi){
  if (roi=="AUD") {
    data2analysis_test <- tibble(
      id_person = rep(1:21, 2),
      hemisphere = factor(rep(c("Left", "Right"), each = 21)),
      gene = t(data) #
    )
  } else {
    data2analysis_test <- tibble(
      id_person = rep(1:28, 2),
      hemisphere = factor(rep(c("Left", "Right"), each = 28)),
      gene = t(data) #
    )
  }
  return(data2analysis_test)
} 

# define a function to do reshuffle the data
reshuffle <- function(data) {
    df <- data %>%
      group_by(id_person) %>%
      mutate(hemisphere = sample(hemisphere, size = n())) %>% 
      arrange(hemisphere, id_person)%>%
      ungroup()
    left <- filter(df, hemisphere == "Left")$gene
    right <- filter(df, hemisphere == "Right")$gene
    data_shuffled <- rbind(left, right)
    
    return(t(data_shuffled))
}
  
#define a function to do data analysis with limma
limma_test <- function (df, roi) {  
  
  df_scaled <- t(scale(t(df))) 
  # Condition vector (left vs right)
  if (roi=="AUD") {
    hemi <- factor(c(rep("Left", 21), rep("Right", 21)))
    id <- factor(rep(1:21, 2))
  } else {
    hemi <- factor(c(rep("Left", 28), rep("Right", 28)))
    id <- factor(rep(1:28, 2))
  }
  
  # Design matrix 
  design <- model.matrix(~id+hemi)
  
  # Fit the linear model
  fit <- lmFit(df_scaled, design)
  
  # Apply eBayes for statistics
  fit <- eBayes(fit)

  # T-statistics
  return (fit)
}

# define function to extract statistics
permu_limma <- function(data, roi){
  fit <- limma_test(data, roi)
  t.stats <- fit$t[, "hemiRight"]
  largest_t <- max(t.stats, na.rm=TRUE)
  smallest_t <- min(t.stats, na.rm=TRUE)
  median_t <- median(t.stats, na.rm=TRUE)
  positive_count <- sum(t.stats > 0, na.rm=TRUE)
  negative_count <- sum(t.stats < 0, na.rm=TRUE)
  return(list(largest_t, smallest_t, median_t, positive_count, negative_count))
}

##------main analysis
data2analysis <- load_data(roi)
data2analysis_test <- creat_shuffle_data(data2analysis, roi)

results <- mclapply(1:10000,
                    function(x)  permu_limma(reshuffle(data2analysis_test), roi),
                    mc.cores = getOption("mc.cores", 2L)
)

# results <- lapply(1:100, 
#                     function(x)  permu_limma(reshuffle(data2analysis_test), roi)
# )

results_df <- map_df(results, ~{
  tibble(largest_t = .x[[1]], 
         smallest_t = .x[[2]], 
         median_t = .x[[3]], 
         positive_count = .x[[4]], 
         negative_count = .x[[5]])
})

# # save the results into xlsx
# if (roi=="AUD") {
#   write.xlsx(results_df, file = "output/1_Gene_expression/All_genes/Gene_by_gene/region_seperated/batch1_2/Hemi_limma_permutation_distribution.xlsx", sheetName = roi)
# } else {
#   wb <- loadWorkbook("output/1_Gene_expression/All_genes/Gene_by_gene/region_seperated/batch1_2/Hemi_limma_permutation_distribution.xlsx")
#   if (roi %in% names(wb)) {
#     removeWorksheet(wb, roi)  # Remove the existing sheet
#   }
#   addWorksheet(wb, roi)
#   writeData(wb, roi, results_df)
#   saveWorkbook(wb, "output/1_Gene_expression/All_genes/Gene_by_gene/region_seperated/batch1_2/Hemi_limma_permutation_distribution.xlsx", overwrite=TRUE)
# }


#######----------------------------- Get the real data 
t_large_p95 <- quantile(results_df[,1], 0.975, na.rm = TRUE)
t_small_p5 <- quantile(results_df[,2], 0.025, na.rm = TRUE)

# define the function to do limma on real data and extract the statistics
real_limma <- function(data, roi) {
  fit <- limma_test(data, roi)
  logFC <- fit$coefficients[, "hemiRight"]  # Log fold changes
  p.values <- fit$p.value[, "hemiRight"]  # Raw p-values
  prop.TrueNull <- propTrueNull(p.values, method="hist")
  adj.p.values <- p.adjust(p.values, method = "fdr")  # Adjusted p-values (FDR)
  t.stats <- fit$t[, "hemiRight"]
  p_permu <- (t.stats > t_large_p95) | (t.stats < t_small_p5)
  B.stats <- fit$lods[, "hemiRight"]  # B-statistics (log-odds)
  
  fit$adj.p.values <- adj.p.values
  # Create a data frame combining the results
  full_results <- data.frame(
    logFC = logFC,
    P.Value = p.values,
    prop.TrueNull = prop.TrueNull,
    adj.P.Val = adj.p.values,
    t = t.stats,
    t_large_p95 = t_large_p95,
    t_small_p5 = t_small_p5,
    p_permu = p_permu,
    B = B.stats,
    row.names = rownames(fit)  # Retain gene names as row names
  )
  
  colnames(full_results) <- c("logFc","P","prop.TrueNull", "adj.P","t","t_large_p95","t_small_p5", "p_permu", "B")
  
  # 
  # if (roi=="AUD") {
  #   write.xlsx(full_results, file = "output/1_Gene_expression/All_genes/Gene_by_gene/region_seperated/batch1_2/Hemi_limma_gene_by_gene.xlsx", sheetName = roi, rowNames=TRUE, colNames=TRUE)
  # } else {
  #   wb <- loadWorkbook("output/1_Gene_expression/All_genes/Gene_by_gene/region_seperated/batch1_2/Hemi_limma_gene_by_gene.xlsx")
  #   if (roi %in% names(wb)) {
  #     removeWorksheet(wb, roi)  # Remove the existing sheet
  #   }
  #   addWorksheet(wb, roi)
  #   full_results_with_rownames <- cbind(Gene = rownames(full_results), full_results)
  #   writeData(wb, roi, full_results_with_rownames)
  #   saveWorkbook(wb, "output/1_Gene_expression/All_genes/Gene_by_gene/region_seperated/batch1_2/Hemi_limma_gene_by_gene.xlsx", overwrite=TRUE)
  # }
  
  png(filename = paste0("output/1_Gene_expression/All_genes/Gene_by_gene/region_seperated/batch1_2/Hemi_Fig1_residual_", roi, ".png"), 
      width = 6, height = 6, units = 'in', res = 300)
  plotSA(fit, xlab = "Average log-expression", ylab = "sqrt(sigma)", zero.weights = FALSE,
         pch = 16, cex = 0.65, col = c("black","red"))
  dev.off()
  
  
  if (roi == "AUD") {
    xlims  <- c(-0.85, 0.85)
    xticks <- seq(-0.6, 0.6, by = 0.3)
  } else {
    xlims  <- c(-0.65, 0.65)
    xticks <- pretty(xlims)
  }
  p_value_threshold <- 0.05
  cairo_pdf(filename = paste0("output/1_Gene_expression/All_genes/Gene_by_gene/region_seperated/batch1_2/Hemi_Fig2_volcano_", roi, ".pdf"), 
      width = 6.5, height = 9, family = "Arial")
  
  par(family = "Arial", tck = 0.02, cex.axis = 1.2, font.axis = 2)
  
  # volcanoplot(fit, coef = "hemiRight", style = "p-value", xlim = xlims,
  #             highlight = sum(fit$adj.p.values < p_value_threshold), names = rownames(data) , hl.col = "blue",
  #             xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.65, axes=FALSE)
  volcanoplot(fit, coef = "hemiRight", style = "p-value", xlim = xlims,
              xlab = "", ylab = "", pch = 16, cex = 1, axes = FALSE)
  
  axis(1, at = xticks, labels = xticks, cex.axis = 1.2, font.axis = 2)
  axis(2, cex.axis = 1.2, font.axis = 2)
  box(lwd=2)
  
  significant_genes <- fit$adj.p.values < p_value_threshold
  points(fit$coef[significant_genes, "hemiRight"], -log10(fit$p.value[significant_genes, "hemiRight"]), col = "red", pch = 16, cex=1.2)
  
  # # Add custom axis labels in Arial, bold, size 14
  title(xlab = "Log2 Fold Change", font.lab = 2, family = "Arial", cex.lab = 1.4)
  title(ylab = expression(-log[10](italic(p))), font.lab = 2, family = "Arial", cex.lab = 1.4)
  
  dev.off()
  
}


data2analysis <- load_data(roi)
real_limma(data2analysis, roi)
