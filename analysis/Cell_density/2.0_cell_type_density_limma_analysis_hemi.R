# This code is going to do limma on hemi effect

# Created 07-OCT-2024; updated 10-OCT-2024; updated 25-April-2025
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
library(limma)
library(future)
plan("multicore", workers = 40) 

args <- commandArgs(trailingOnly=TRUE)
roi <- args[1]  # This is AUD HIP
######------------------------------ load the data

# load the data
load_data <- function(roi){
  loaded_data <- read_excel("Hemi_data2analysis.xlsx", sheet = roi) # read the file in
  data_loaded <- t(na.omit(loaded_data[,-2:-3]))
  
  colsnames <- data_loaded[1,]
  data2analysis <- data_loaded[-1,]
  
  
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
# create the reshuffle dataset
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

# define function to do reshuffle
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
  
#data analysis with limma
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

results_df <- map_df(results, ~{
  tibble(largest_t = .x[[1]], 
         smallest_t = .x[[2]], 
         median_t = .x[[3]], 
         positive_count = .x[[4]], 
         negative_count = .x[[5]])
})


if (roi=="AUD") {
  write.xlsx(results_df, file = "Hemi_limma_permutation_distribution.xlsx", sheetName = roi)
} else {
  wb <- loadWorkbook("Hemi_limma_permutation_distribution.xlsx")
  if (roi %in% names(wb)) {
    removeWorksheet(wb, roi)  # Remove the existing sheet
  }
  addWorksheet(wb, roi)
  writeData(wb, roi, results_df)
  saveWorkbook(wb, "Hemi_limma_permutation_distribution.xlsx", overwrite=TRUE)
}


#######----------------------------- Get the real data 
t_large_p95 <- quantile(results_df[,1], 0.975, na.rm = TRUE)
t_small_p5 <- quantile(results_df[,2], 0.025, na.rm = TRUE)

# define the function
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
  
  
  if (roi=="AUD") {
    write.xlsx(full_results, file = "Hemi_limma_results_cell_type_density.xlsx", sheetName = roi, rowNames=TRUE, colNames=TRUE)
  } else {
    wb <- loadWorkbook("Hemi_limma_results_cell_type_density.xlsx")
    if (roi %in% names(wb)) {
      removeWorksheet(wb, roi)  # Remove the existing sheet
    }
    addWorksheet(wb, roi)
    full_results_with_rownames <- cbind(Gene = rownames(full_results), full_results)
    writeData(wb, roi, full_results_with_rownames)
    saveWorkbook(wb, "Hemi_limma_results_cell_type_density.xlsx", overwrite=TRUE)
  }
  
  png(filename = paste0("Hemi_Fig1_residual_", roi, ".png"), 
      width = 6, height = 6, units = 'in', res = 300)
  plotSA(fit, xlab = "Average log-expression", ylab = "sqrt(sigma)", zero.weights = FALSE,
         pch = 16, cex = 0.65, col = c("black","red"))
  dev.off()
  
  p_value_threshold <- 0.05
  png(filename = paste0("Hemi_Fig2_volcano_", roi, ".png"), 
      width = 7, height = 8, units = 'in', res = 300)
  volcanoplot(fit, coef = "hemiRight", style = "p-value",
              highlight = sum(fit$adj.p.values < p_value_threshold), names = rownames(data) , hl.col = "blue",
              xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.65)
  significant_genes <- fit$adj.p.values < p_value_threshold
  points(fit$coef[significant_genes, "hemiRight"], -log10(fit$p.value[significant_genes, "hemiRight"]), col = "red", pch = 16)
  dev.off()
  
}

data2analysis <- load_data(roi)
real_limma(data2analysis, roi)





