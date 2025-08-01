# This code is going to do limma on sex effect

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
roi <- args[1]  # This is AUD HIP

######------------------------------ load the data
load_data <- function(roi){
  loaded_data <- read_excel("output/1_Gene_expression/All_genes/Gene_by_gene/region_seperated/batch1_2/Sex_data2analysis.xlsx", sheet = roi) # read the file in
  data_loaded <- t(loaded_data[,-2]) # get rid of the sex column
  colsnames <- data_loaded[1,] # extract the first row of the mouse ids, such as M670 ...
  data2analysis <- data_loaded[-1,] # get rid of the first row
  
  if (roi=="AUD") {
    new_colsnames <- c(paste0("male_", colsnames[1:12]), paste0("female_", colsnames[13:21]))
  } else {
    new_colsnames <- c(paste0("male_", colsnames[1:16]), paste0("female_", colsnames[17:28]))
  }
  
  colnames(data2analysis) <- new_colsnames
  data2analysis <- as.data.frame(data2analysis)
  names_row <- rownames(data2analysis)
  
  data2analysis <- as.matrix(apply(data2analysis, 2, as.numeric))
  rownames(data2analysis) <- names_row
  
  # 
  # data2analysis_tranversed <- t(data2analysis) %>%
  #   as.data.frame(.)
  # 
  # if (roi=="AUD") {
  #   data2analysis_tranversed$batch = rep(c(1, 2, 1, 2), times = c(4,8,6,3))
  # } else {
  #   data2analysis_tranversed$batch = rep(c(1, 2, 1, 2), times = c(8,8,8,4))
  # }
  # 
  # data2analysis <- t(data2analysis_tranversed)
  
  
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
  
  data2test <- ComBat(dat = data2analysis, batch = batch, par.prior = TRUE, prior.plots = FALSE) %>%
    as.data.frame(.)
  
  return(data2test)
}

data2analysis <- load_data(roi)

#######----------------------------- Permutation of the tests

## create tibble format of data
if (roi=="AUD") {
  data2analysis_test <- tibble(
    sex = factor(rep(c("male", "female"), times = c(12,9))),
    gene = t(data2analysis) 
  ) 
} else {
  data2analysis_test <- tibble(
    sex = factor(rep(c("male", "female"), times = c(16,12))),
    gene = t(data2analysis)
  )
}

# define function to do reshuffle the data
reshuffle <- function(data) {
    if (roi == "AUD") {
      data$sex <- sample(c(rep("male", 12), rep("female", nrow(data) - 12)))
    } else {
      data$sex <- sample(c(rep("male", 16), rep("female", nrow(data) - 16)))
    }
    male <- filter(data, sex == "male")$gene
    female <- filter(data, sex == "female")$gene
    data_shuffled <- rbind(male, female)
    return(data_shuffled)
}
  
#defien function to do data analysis with limma
limma_test <- function (df, roi) {  
  
  # df <- as.data.frame(df)
  # batch <- factor(df$batch)
  # 
  # df$batch <- NULL
  
  df_scaled <- scale(t(df))
  
  # Condition vector (male vs female)
  if (roi=="AUD") {
    sex <- factor(rep(c("male", "female"), times = c(12,9)))

  } else {
    sex <- factor(rep(c("male", "female"), times = c(16,12)))

  }
  
  # Design matrix 
  # design <- model.matrix(~batch+sex)
  design <- model.matrix(~sex)
  
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
  t.stats <- fit$t[, "sexmale"]
  largest_t <- max(t.stats, na.rm=TRUE)
  smallest_t <- min(t.stats, na.rm=TRUE)
  median_t <- median(t.stats, na.rm=TRUE)
  positive_count <- sum(t.stats > 0, na.rm=TRUE)
  negative_count <- sum(t.stats < 0, na.rm=TRUE)
  return(list(largest_t, smallest_t, median_t, positive_count, negative_count))
}

##------main analysis

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

# save the results into xlsx
if (roi=="AUD") {
  write.xlsx(results_df, file = "output/1_Gene_expression/All_genes/Gene_by_gene/region_seperated/batch1_2/Sex_limma_permutation_distribution.xlsx", sheetName = roi)
} else {
  wb <- loadWorkbook("output/1_Gene_expression/All_genes/Gene_by_gene/region_seperated/batch1_2/Sex_limma_permutation_distribution.xlsx")
  if (roi %in% names(wb)) {
    removeWorksheet(wb, roi)  # Remove the existing sheet
  }
  addWorksheet(wb, roi)
  writeData(wb, roi, results_df)
  saveWorkbook(wb, "output/1_Gene_expression/All_genes/Gene_by_gene/region_seperated/batch1_2/Sex_limma_permutation_distribution.xlsx", overwrite=TRUE)
}


#######----------------------------- Get the real data 
t_large_p95 <- quantile(results_df[,1], 0.975, na.rm = TRUE)
t_small_p5 <- quantile(results_df[,2], 0.025, na.rm = TRUE)

# define the function to do limma on real data and extract the statistics
real_limma <- function(data, roi) {
  fit <- limma_test(data, roi)
  logFC <- fit$coefficients[, "sexmale"]  # Log fold changes
  p.values <- fit$p.value[, "sexmale"]  # Raw p-values
  prop.TrueNull <- propTrueNull(p.values, method="hist")
  adj.p.values <- p.adjust(p.values, method = "fdr")  # Adjusted p-values (FDR)
  t.stats <- fit$t[, "sexmale"]
  p_permu <- (t.stats > t_large_p95) | (t.stats < t_small_p5)
  B.stats <- fit$lods[, "sexmale"]  # B-statistics (log-odds)
  
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
    write.xlsx(full_results, file = "output/1_Gene_expression/All_genes/Gene_by_gene/region_seperated/batch1_2/Sex_limma_gene_by_gene.xlsx", sheetName = roi, rowNames=TRUE, colNames=TRUE)
  } else {
    wb <- loadWorkbook("output/1_Gene_expression/All_genes/Gene_by_gene/region_seperated/batch1_2/Sex_limma_gene_by_gene.xlsx")
    if (roi %in% names(wb)) {
      removeWorksheet(wb, roi)  # Remove the existing sheet
    }
    addWorksheet(wb, roi)
    full_results_with_rownames <- cbind(Gene = rownames(full_results), full_results)
    writeData(wb, roi, full_results_with_rownames)
    saveWorkbook(wb, "output/1_Gene_expression/All_genes/Gene_by_gene/region_seperated/batch1_2/Sex_limma_gene_by_gene.xlsx", overwrite=TRUE)
  }
  
  png(filename = paste0("output/1_Gene_expression/All_genes/Gene_by_gene/region_seperated/batch1_2/Sex_Fig1_residual_", roi, ".png"), 
      width = 6, height = 6, units = 'in', res = 300)
  plotSA(fit, xlab = "Average log-expression", ylab = "sqrt(sigma)", zero.weights = FALSE,
         pch = 16, cex = 0.65, col = c("black","red"))
  dev.off()
  
  p_value_threshold <- 0.05
  png(filename = paste0("output/1_Gene_expression/All_genes/Gene_by_gene/region_seperated/batch1_2/Sex_Fig2_volcano_", roi, ".png"), 
      width = 7, height = 8, units = 'in', res = 300)
  volcanoplot(fit, coef = "sexmale", style = "p-value",
              highlight = sum(fit$adj.p.values < p_value_threshold), names = rownames(data) , hl.col = "blue",
              xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.65)
  significant_genes <- fit$adj.p.values < p_value_threshold
  points(fit$coef[significant_genes, "sexmale"], -log10(fit$p.value[significant_genes, "sexmale"]), col = "red", pch = 16)
  dev.off()
  
}

data2analysis <- load_data(roi)
real_limma(t(data2analysis), roi)
