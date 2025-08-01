# This code is going to do Permutation analysis on the sex effect on cell type density

# Created 30-Aug-2024
# Created by M.-Y. WANG 

# Remove all objects created before to prevent clusing
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
args <- commandArgs(trailingOnly=TRUE)
roi <- args[1]  # This is CA1 CA3 DG AUD

################################# load the data
# read data
loaded_data <- read_excel("counts_exclude_CA3_30_data2analysis.xlsx", sheet = roi) # read the file in
num_cols <- ncol(loaded_data)

#################### compute the cell density
if (roi=="AUD") {
  
    data_male <- (loaded_data[1:12,4:num_cols] + loaded_data[22:33,4:num_cols])
    data_female <- (loaded_data[13:21,4:num_cols] + loaded_data[34:42,4:num_cols])
    df_cell_count <- rbind(data_male,data_female)
    
    df_cell_total <- rowSums(df_cell_count)
    data2analysis <- t(df_cell_count/df_cell_total)
    
} else {
  data_male <- (loaded_data[1:16, 4:num_cols] + loaded_data[29:44, 4:num_cols]) #male
  data_female <-(loaded_data[17:28 ,4:num_cols] + loaded_data[45:56, 4:num_cols])  #female
  df_cell_count <- rbind(data_male,data_female)
  
  df_cell_total <- rowSums(df_cell_count)
  data2analysis <- t(df_cell_count/df_cell_total)
}

#################### save the data

data2save <- as.data.frame(t(data2analysis))
if (roi == "AUD") {
  data2save$sex <-
    c(rep("male", 12),
      rep("female", 9))
  data2save$id <- loaded_data$id[1:21]
  data2save <-
    data2save[, c('id', 'sex', setdiff(names(data2save), c('id', 'sex')))]
  
  write.xlsx(data2save, file = "Sex_data2analysis.xlsx", sheetName = roi)
} else {
  data2save$sex <-
    c(rep("male", 16),
      rep("female", 12))
  data2save$id <- loaded_data$id[1:28]
  data2save <-
    data2save[, c('id', 'sex', setdiff(names(data2save), c('id', 'sex')))]
  
  wb <-loadWorkbook("Sex_data2analysis.xlsx")
  addWorksheet(wb, roi)
  writeData(wb, roi, data2save)
  saveWorkbook(wb,"Sex_data2analysis.xlsx",overwrite = TRUE)
}


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
  return(t(data_shuffled))
}

#defien function to do data analysis with limma
limma_test <- function (df, roi) {  
  
  # df <- as.data.frame(df)
  # batch <- factor(df$batch)
  # 
  # df$batch <- NULL
  
  df_scaled <- t(scale(t(df)))
  
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
  write.xlsx(results_df, file = "Sex_limma_permutation_distribution.xlsx", sheetName = roi)
} else {
  wb <- loadWorkbook("Sex_limma_permutation_distribution.xlsx")
  if (roi %in% names(wb)) {
    removeWorksheet(wb, roi)  # Remove the existing sheet
  }
  addWorksheet(wb, roi)
  writeData(wb, roi, results_df)
  saveWorkbook(wb, "Sex_limma_permutation_distribution.xlsx", overwrite=TRUE)
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
    write.xlsx(full_results, file = "Sex_limma_results_cell_type_density.xlsx", sheetName = roi, rowNames=TRUE, colNames=TRUE)
  } else {
    wb <- loadWorkbook("Sex_limma_results_cell_type_density.xlsx")
    if (roi %in% names(wb)) {
      removeWorksheet(wb, roi)  # Remove the existing sheet
    }
    addWorksheet(wb, roi)
    full_results_with_rownames <- cbind(Gene = rownames(full_results), full_results)
    writeData(wb, roi, full_results_with_rownames)
    saveWorkbook(wb, "Sex_limma_results_cell_type_density.xlsx", overwrite=TRUE)
  }
  
  png(filename = paste0("Sex_Fig1_residual_", roi, ".png"), 
      width = 6, height = 6, units = 'in', res = 300)
  plotSA(fit, xlab = "Average log-expression", ylab = "sqrt(sigma)", zero.weights = FALSE,
         pch = 16, cex = 0.65, col = c("black","red"))
  dev.off()
  
  p_value_threshold <- 0.05
  png(filename = paste0("Sex_Fig2_volcano_", roi, ".png"), 
      width = 7, height = 8, units = 'in', res = 300)
  volcanoplot(fit, coef = "sexmale", style = "p-value",
              highlight = sum(fit$adj.p.values < p_value_threshold), names = rownames(data) , hl.col = "blue",
              xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.65)
  significant_genes <- fit$adj.p.values < p_value_threshold
  points(fit$coef[significant_genes, "sexmale"], -log10(fit$p.value[significant_genes, "sexmale"]), col = "red", pch = 16)
  dev.off()
  
}

real_limma(data2analysis, roi)

