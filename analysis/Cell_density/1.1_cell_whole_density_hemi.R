# This code is going to do Permutation analysis on the hemi effect on cell density

# Created 23-Sep-2024; updated 15-March-2025
# Created by M.-Y. WANG 

# Remove all objects created before to prevent clushing
rm(list = ls())

# Set the working directory to the path where your files are located
# setwd("/Users/joeywang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Mouse_brain") # change it to the file directory
# setwd("/Users/wang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Mouse_brain/")
setwd("/data/workspaces/lag/workspaces/lg-func-asym/working_data/Mengyun")

library(readxl)
library(openxlsx)
# library(MVTests)
library(ggplot2)
library(dplyr)
library(purrr)
library(parallel)

################################################### Preprocessing
# read data
df_roi_AUD <- read_excel("Spatial_transcriptomics_batch1_batch2_Xenium_extractions_AUD.xlsx", sheet = "cell_density")
df_roi_HIP <- read_excel("Spatial_transcriptomics_batch1_batch2_Xenium_extractions_HIP.xlsx", sheet = "cell_density")
df_roi_CA1 <- read_excel("Spatial_transcriptomics_batch1_batch2_Xenium_extractions_CA1.xlsx", sheet = "cell_density")
df_roi_CA3 <- read_excel("Spatial_transcriptomics_batch1_batch2_Xenium_extractions_CA3.xlsx", sheet = "cell_density")
df_roi_DG <- read_excel("Spatial_transcriptomics_batch1_batch2_Xenium_extractions_DG.xlsx", sheet = "cell_density")

df_data_AUD <- df_roi_AUD[, (which(df_roi_AUD[2, ] == "Density"))]
df_data_HIP <- df_roi_HIP[, (which(df_roi_HIP[2, ] == "Density"))]
df_data_CA1 <- df_roi_CA1[, (which(df_roi_CA1[2, ] == "Density"))]
df_data_CA3 <- df_roi_CA3[, (which(df_roi_CA3[2, ] == "Density"))]
df_data_DG <- df_roi_DG[, (which(df_roi_DG[2, ] == "Density"))]

df_data_AUD <- df_data_AUD[c(-1:-2),]
colnames(df_data_AUD) <- c("AUD", "AUD")
df_data_AUD <- rbind(df_data_AUD[,1], df_data_AUD[,2])
df_data_HIP <- df_data_HIP[c(-1:-2),]
colnames(df_data_HIP) <- c("HIP","HIP")
df_data_HIP <- rbind(df_data_HIP[,1], df_data_HIP[,2])

df_data_CA1 <- df_data_CA1[c(-1:-2),]
colnames(df_data_CA1) <- c("CA1","CA1")
df_data_CA1 <- rbind(df_data_CA1[,1], df_data_CA1[,2])
df_data_CA3 <- df_data_CA3[c(-1:-2),]
colnames(df_data_CA3) <- c("CA3","CA3")
df_data_CA3 <- rbind(df_data_CA3[,1], df_data_CA3[,2])
df_data_DG <- df_data_DG[c(-1:-2),]
colnames(df_data_DG) <- c("DG","DG")
df_data_DG <- rbind(df_data_DG[,1], df_data_DG[,2])

data2analysis <- cbind(df_data_AUD, df_data_HIP, df_data_CA1, df_data_CA3, df_data_DG)
data2analysis <- as.data.frame(data2analysis)


#save the data
data2save <- data2analysis
data2save$hemi <- c(rep("left", 31), rep("right", 31))
data2save$sex <- c(rep("male", 18), rep("female", 13), rep("male", 18), rep("female", 13))
data2save$sample_id <- rep(c("M669", "M670", "M671", "M672", "M673", "M674", "M675", "M676", "M677", 
                      "M678", "M234", "M253", "M071", "M083", "M650", "M638", "M076", "M236", 
                      "F679", "F680", "F681", "F682", "F683", "F685", "F686", "F687", "F688", 
                      "F073", "F078", "F087", "F090"), 2)
data2save$id <- rep(1:31, 2)
data2save <- data2save[, c('id', 'sample_id','sex','hemi',setdiff(names(data2save), c('id','sample_id','sex','hemi')))]

write.xlsx(data2save, file = "output/3_Cell_density/Overall_cell/batch1_2/Hemi_data2analysis_density.xlsx", 
           sheetName="overall_cell_density")

############################################################## Analysis

fit_AUD <- lm(AUD ~ hemi + id, data = data2save, na.action = na.omit)
fit_HIP <- lm(HIP ~ hemi + id, data = data2save, na.action = na.omit)
fit_CA1 <- lm(CA1 ~ hemi + id, data = data2save, na.action = na.omit)
fit_CA3 <- lm(CA3 ~ hemi + id, data = data2save, na.action = na.omit)
fit_DG <- lm(DG ~ hemi + id, data = data2save, na.action = na.omit)

# Define a function to extract stats from an lm() fit
extract_lm_stats <- function(fit) {
  
  s <- summary(fit)
  ctab <- s$coefficients
  idx <- which(rownames(ctab) == "hemiright")
  
  # Extract the relevant stats (assuming there's exactly one match)
  t_value <- ctab[idx, "t value"]
  p_value <- ctab[idx, "Pr(>|t|)"]
  
  # Degrees of freedom for the residuals:
  df_val <- s$df[2]
  
  # Confidence interval
  ci <- confint(fit)
  ci_term <- ci[rownames(ci) %in% rownames(ctab)[idx], , drop = FALSE]
  ci_low <- ci_term[1, 1]
  ci_high <- ci_term[1, 2]
  
  return(list(
    t = t_value,
    df = df_val,
    p = p_value,
    ci_low = ci_low,
    ci_high = ci_high
  ))
}

# Extract stats for each fit
stats_AUD <- extract_lm_stats(fit_AUD)
stats_HIP <- extract_lm_stats(fit_HIP)
stats_CA1 <- extract_lm_stats(fit_CA1)
stats_CA3 <- extract_lm_stats(fit_CA3)
stats_DG <- extract_lm_stats(fit_DG)

# Adjust p-values for multiple testing
raw_pvals <- c(stats_AUD$p, stats_HIP$p, stats_CA1$p, stats_CA3$p, stats_DG$p)
adj_pvals <- p.adjust(raw_pvals, method = "fdr")

# Combine into a single results data frame
results_df <- data.frame(
  region   = c("AUD", "HIP", "CA1", "CA3", "DG"),
  t        = c(stats_AUD$t,    stats_HIP$t,  stats_CA1$t,  stats_CA3$t,  stats_DG$t),
  df       = c(stats_AUD$df,   stats_HIP$df, stats_CA1$df, stats_CA3$df, stats_DG$df),
  P.Value  = c(stats_AUD$p,    stats_HIP$p, stats_CA1$p, stats_CA3$p, stats_DG$p),
  adjust.P = adj_pvals,
  CI.Low   = c(stats_AUD$ci_low,  stats_HIP$ci_low, stats_CA1$ci_low, stats_CA3$ci_low, stats_DG$ci_low),
  CI.High  = c(stats_AUD$ci_high, stats_HIP$ci_high, stats_CA1$ci_high, stats_CA3$ci_high, stats_DG$ci_high)
)

# Write the table to an Excel file
write.xlsx(results_df, "output/3_Cell_density/Overall_cell/batch1_2/Hemi_results_cell_density.xlsx", sheetName = "overall_cell_density")





# 
# # create a tibble version of dataframe
# data2analysis_test <- tibble(
#     id_person = rep(1:31, 2),
#     hemisphere = factor(rep(c("Left", "Right"), each = 31)),
#     cell = data2analysis_corrected)
# 
# # define function to reshuffle the data
# reshuffle <- function(data) {
#   df <- data %>%
#     group_by(id_person) %>%
#     mutate(hemisphere = sample(hemisphere, size = n())) %>% 
#     arrange(hemisphere, id_person)%>%
#     ungroup()
#   
#   left <- filter(df, hemisphere == "Left")$cell
#   right <- filter(df, hemisphere == "Right")$cell
#   
#   data_shuffled <- rbind(left, right)
#   data_shuffled <- as.data.frame(data_shuffled)
#   return(data_shuffled)
# }
# 
# #define function to do the paired t test
# get_stats <- function(data) {
#   
#   data_reshuffled <- reshuffle(data)
#   t_values <- lapply(1:length(data_reshuffled), function(i) {
#     lhemi <- na.omit(data_reshuffled[1:(nrow(data_reshuffled)/2),i])
#     rhemi <- na.omit(data_reshuffled[(nrow(data_reshuffled)/2+1):nrow(data_reshuffled),i])                      
#     result <- t.test(as.numeric(lhemi), 
#                      as.numeric(rhemi),
#                      paired=TRUE)
#     return(result$statistic)
#   })
#   
#   largest_t <- max(unlist(t_values), na.rm = TRUE)
#   smallest_t <- min(unlist(t_values), na.rm = TRUE)
#   
#   return(list(largest_t, smallest_t)) # result as a list
# }
# 
# # Apply get_stats function 10000 times
# results <- mclapply(1:10000, 
#                     function(x) get_stats(data2analysis_test),
#                     mc.cores = getOption("mc.cores", 1L)
# )
# 
# # Convert results to a data frame 
# results_df <- map_df(results, ~{
#   tibble(largest_t = .x[[1]], 
#          smallest_t = .x[[2]])
# })
# 
# # save the results
# write.xlsx(results_df, file = "output/3_Cell_density/Overall_cell/batch1_2/cell_density_permutation_distribution_hemi.xlsx", sheetName = "cell_density")
# 
# 
# # set the alpha threthold
# t_large_p95 <- quantile(results_df[,1], 0.975, na.rm = TRUE)
# t_small_p5 <- quantile(results_df[,2], 0.025, na.rm = TRUE)
# 
# 
# # function to tyde up
# test_update <- function (df) {
#   # get the orig p values
#   p_values <- sapply(df, function(x)
#     x$p.value)
#   
#   # fdr the p values
#   p_adjusted <- p.adjust(p_values, method = "BH")
#   
#   # get the orig t values
#   t_values = sapply(df, function(x)
#     x$statistic)
#   # calculate the p values based on the permutation
#   p_permu <- (t_values > t_large_p95) | (t_values < t_small_p5)
#   
#   df <- mapply(function(x, y, z) {
#     x$p.adjusted <- y
#     x$p.permu <- z
#     return(x)
#   }, df, p_adjusted, p_permu, SIMPLIFY = FALSE)
#   
#   # Create a data frame with the results
#   density_data_t_test <- data.frame(
#     t = sapply(df, function(x)
#       x$statistic),
#     df = sapply(df, function(x)
#       x$parameter),
#     P.Value = sapply(df, function(x)
#       x$p.value),
#     ajust.P = sapply(df, function(x)
#       x$p.adjusted),
#     Permutation = sapply(df, function(x)
#       x$p.permu),
#     CI.Low = sapply(df, function(x)
#       x$conf.int[1]),
#     CI.High = sapply(df, function(x)
#       x$conf.int[2])
#   )
#   
#   # Add the region name column
#   density_data_t_test$region <- c("AUD", "HIP")
#   density_data_t_test <- density_data_t_test[, c(
#     "region",
#     "t",
#     "df",
#     "P.Value",
#     "ajust.P",
#     "Permutation",
#     "CI.Low",
#     "CI.High"
#   )]
#   return(density_data_t_test)
# }
# 
# 
# # Perform paired t-tests on real data left and right hemi
# paired_hemi <- lapply(1:ncol(data2analysis_corrected), function(i) {
#   
#   lhemi <- na.omit(data2analysis_corrected[1:(nrow(data2analysis_corrected)/2),i])
#   rhemi <- na.omit(data2analysis_corrected[(nrow(data2analysis_corrected)/2+1):nrow(data2analysis_corrected),i])                      
#   result <- t.test(as.numeric(lhemi), 
#                    as.numeric(rhemi),
#                    paired=TRUE)
#   return(result)
# })
# 
# 
# write.xlsx(test_update(paired_hemi), "output/3_Cell_density/Overall_cell/batch1_2/cell_density_paired_test_hemi_permu.xlsx", sheetName = "cell_density")



