# Linear-model analysis of hemisphere effects on overall cell density in AUD and HIP.

# Created 23-Sep-2024; updated 15-March-2025; updated 25-Sep-2026
# Created by M.-Y. WANG 

# Remove all objects created before to prevent clushing
rm(list = ls())

# Set the working directory to the path where your files are located
# setwd("/Users/joeywang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Mouse_brain") # change it to the file directory
# setwd("/Users/wang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Mouse_brain/")
setwd("C:/Users/menwan2/Documents/Codex/2026-09-24/creat-a-folder-for-this-project-2/outputs/manuscript-review/data")

library(readxl)
library(openxlsx)
# library(MVTests)
library(ggplot2)
library(dplyr)
library(purrr)
library(parallel)

input_file <- "Hemi_data2analysis.xlsx"
overall_file <- "../../Overall_cell/batch1_2/Hemi_data2analysis_density.xlsx"
output_file <- "Hemi_paired_t_results_cell_type_density.xlsx"

# With no arguments run both regions; optionally specify AUD or HIP_combined.
args <- commandArgs(trailingOnly = TRUE)
rois <- if (length(args) == 0L) c("AUD", "HIP_combined") else unique(args)
if (!all(rois %in% c("AUD", "HIP_combined"))) {
  stop("Only AUD and HIP_combined are supported.")
}

######------------------------------ Load the data

load_data <- function(roi) {
  # Keep numeric columns at full precision, without transposing through text.
  loaded_data <- as.data.frame(read_excel(input_file, sheet = roi,
                                         .name_repair = "check_unique"))
  metadata <- c("id", "sex", "hemi")
  if (!all(metadata %in% names(loaded_data))) {
    stop(roi, ": missing id, sex or hemi metadata.")
  }
  if (anyNA(loaded_data[metadata]) ||
      any(vapply(loaded_data[metadata], function(x) any(trimws(x) == ""), logical(1)))) {
    stop(roi, ": missing sample keys.")
  }
  if (!all(loaded_data$hemi %in% c("left", "right"))) {
    stop(roi, ": hemi must contain only left and right.")
  }
  cell_types <- setdiff(names(loaded_data), c(metadata, "overall"))
  if (length(cell_types) != 15L ||
      !all(vapply(loaded_data[cell_types], is.numeric, logical(1)))) {
    stop(roi, ": expected 15 numeric cell-type proportion columns.")
  }

  overall_data <- as.data.frame(read_excel(overall_file,
                                          sheet = "overall_cell_density",
                                          .name_repair = "check_unique"))
  overall_roi <- if (roi == "HIP_combined") "HIP" else "AUD"
  if (!all(c("sample_id", "sex", "hemi", overall_roi) %in% names(overall_data))) {
    stop(roi, ": overall-density source is missing required columns.")
  }
  source_density <- overall_data[[overall_roi]]
  overall_data$overall <- suppressWarnings(as.numeric(source_density))
  if (any(!is.na(source_density) & is.na(overall_data$overall)) ||
      any(!is.na(overall_data$overall) & !is.finite(overall_data$overall))) {
    stop(roi, ": overall density contains nonnumeric or nonfinite values.")
  }
  # The overall source includes mice absent from this ROI; omit missing densities.
  overall_data <- overall_data[!is.na(overall_data$overall), , drop = FALSE]
  source_metadata <- c("sample_id", "sex", "hemi")
  if (anyNA(overall_data[source_metadata]) ||
      any(vapply(overall_data[source_metadata], function(x) any(trimws(x) == ""), logical(1)))) {
    stop(roi, ": missing sample keys in the overall-density source.")
  }

  # Cell-type id is the M/F mouse label: match it to sample_id, NOT numeric id.
  cell_keys <- paste(loaded_data$id, loaded_data$sex, loaded_data$hemi, sep = "\r")
  source_keys <- paste(overall_data$sample_id, overall_data$sex,
                       overall_data$hemi, sep = "\r")
  if (anyDuplicated(cell_keys) || anyDuplicated(source_keys)) {
    stop(roi, ": duplicate mouse + sex + hemisphere key.")
  }
  if (!setequal(cell_keys, source_keys)) {
    stop(roi, ": cell-type and nonmissing overall-density keys do not fully match.")
  }
  matched_overall <- overall_data$overall[match(cell_keys, source_keys)]
  needs_update <- !"overall" %in% names(loaded_data) ||
    !identical(loaded_data$overall, matched_overall)
  loaded_data$overall <- matched_overall

  # Append/replace only this column (column 19 in the current workbook).
  # An identical numeric column is left untouched; other sheets are preserved.
  if (needs_update) {
    wb <- loadWorkbook(input_file)
    writeData(wb, roi, data.frame(overall = matched_overall),
              startCol = match("overall", names(loaded_data)), startRow = 1,
              colNames = TRUE, rowNames = FALSE)
    saveWorkbook(wb, input_file, overwrite = TRUE)
  }
  return(loaded_data)
}

######------------------------------ Test hemisphere effects

paired_t_test <- function(df, roi) {
  features <- setdiff(names(df), c("id", "sex", "hemi"))
  if (length(features) != 16L || !"overall" %in% features) {
    stop(roi, ": expected 15 cell types plus overall density.")
  }
  left <- df[df$hemi == "left", , drop = FALSE]
  right <- df[df$hemi == "right", , drop = FALSE]
  left_keys <- paste(left$id, left$sex, sep = "\r")
  right_keys <- paste(right$id, right$sex, sep = "\r")
  if (anyDuplicated(left_keys) || anyDuplicated(right_keys) ||
      !setequal(left_keys, right_keys)) {
    stop(roi, ": each mouse + sex must have exactly one left and one right row.")
  }
  # Pair by mouse and sex explicitly; do not assume the spreadsheet row order.
  right <- right[match(left_keys, right_keys), , drop = FALSE]

  results <- lapply(features, function(feature) {
    if (!is.numeric(left[[feature]]) || !is.numeric(right[[feature]])) {
      stop(roi, ": nonnumeric feature ", feature, ".")
    }
    complete_pairs <- is.finite(left[[feature]]) & is.finite(right[[feature]])
    left_values <- left[[feature]][complete_pairs]
    right_values <- right[[feature]][complete_pairs]
    n_pairs <- length(left_values)
    if (n_pairs < 2L) stop(roi, ": fewer than two complete pairs for ", feature, ".")

    # Raw cell-type values remain proportions; overall remains cells/um2.
    # Right is the first argument, so the estimate, CI, t and d_z are right - left.
    fit <- t.test(right_values, left_values, paired = TRUE,
                  alternative = "two.sided", conf.level = 0.95)
    differences <- right_values - left_values
    data.frame(
      Cell_type = feature,
      Measure = if (feature == "overall") "cells/um2" else "proportion",
      n_pairs = n_pairs,
      mean_left = mean(left_values),
      mean_right = mean(right_values),
      mean_right_minus_left = mean(differences),
      CI_low = unname(fit$conf.int[1]),
      CI_high = unname(fit$conf.int[2]),
      t = unname(fit$statistic),
      df = unname(fit$parameter),
      P = fit$p.value,
      d_z = mean(differences) / sd(differences),
      check.names = FALSE
    )
  })
  return(do.call(rbind, results))
}

######------------------------------ Get the real-data results

real_paired_t <- function(data, roi) {
  full_results <- paired_t_test(data, roi)
  # One BH family per ROI: 15 cell-type proportions plus overall density (16 tests).
  full_results$adj.P <- p.adjust(full_results$P, method = "BH")
  full_results <- full_results[c("Cell_type", "Measure", "n_pairs", "mean_left",
                                 "mean_right", "mean_right_minus_left", "CI_low",
                                 "CI_high", "t", "df", "P", "adj.P", "d_z")]
  return(full_results)
}

results <- setNames(lapply(rois, function(roi) {
  data2analysis <- load_data(roi)
  real_paired_t(data2analysis, roi)
}), rois)

# Write the combined workbook once; a single-ROI run preserves the other sheet.
wb <- if (length(rois) == 1L && file.exists(output_file)) {
  loadWorkbook(output_file)
} else {
  createWorkbook()
}
for (roi in rois) {
  if (roi %in% names(wb)) removeWorksheet(wb, roi)
  addWorksheet(wb, roi)
  writeData(wb, roi, results[[roi]], rowNames = FALSE, colNames = TRUE)
  write.table(results[[roi]], file = paste0("Hemi_paired_t_results_", roi, ".tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
saveWorkbook(wb, output_file, overwrite = TRUE)

################################################### Preprocessing
## Read the two main regions only; subregions are outside this analysis.
#df_roi_AUD <- read_excel("Spatial_transcriptomics_batch1_batch2_Xenium_extractions_AUD.xlsx", sheet = "cell_density")
#df_roi_HIP <- read_excel("Spatial_transcriptomics_batch1_batch2_Xenium_extractions_HIP.xlsx", sheet = "cell_density")

#df_data_AUD <- df_roi_AUD[, (which(df_roi_AUD[2, ] == "Density"))]
#df_data_HIP <- df_roi_HIP[, (which(df_roi_HIP[2, ] == "Density"))]
#
#df_data_AUD <- df_data_AUD[c(-1:-2),]
#colnames(df_data_AUD) <- c("AUD", "AUD")
#df_data_AUD <- rbind(df_data_AUD[,1], df_data_AUD[,2])
#df_data_HIP <- df_data_HIP[c(-1:-2),]
#colnames(df_data_HIP) <- c("HIP","HIP")
#df_data_HIP <- rbind(df_data_HIP[,1], df_data_HIP[,2])
#
#data2analysis <- cbind(df_data_AUD, df_data_HIP)
#data2analysis <- as.data.frame(data2analysis)
#
#
##save the data
#data2save <- data2analysis
#data2save$hemi <- c(rep("left", 31), rep("right", 31))
#data2save$sex <- c(rep("male", 18), rep("female", 13), rep("male", 18), rep("female", 13))
#data2save$sample_id <- rep(c("M669", "M670", "M671", "M672", "M673", "M674", "M675", "M676", "M677", 
#                      "M678", "M234", "M253", "M071", "M083", "M650", "M638", "M076", "M236", 
#                      "F679", "F680", "F681", "F682", "F683", "F685", "F686", "F687", "F688", 
#                      "F073", "F078", "F087", "F090"), 2)
## Mouse ID must be categorical: this gives each mouse its own intercept in lm().
## The same ID labels the left and right measurements from that mouse.
#data2save$id <- factor(rep(1:31, 2))
#data2save <- data2save[, c('id', 'sample_id','sex','hemi',setdiff(names(data2save), c('id','sample_id','sex','hemi')))]
#
## Use a distinct filename for the two-region analysis input.
#write.xlsx(data2save, file = "output/3_Cell_density/Overall_cell/batch1_2/Hemi_data2analysis_density.xlsx", 
#           sheetName="overall_cell_density")
#
############################################################## Analysis
#
#fit_AUD <- lm(AUD ~ hemi + id, data = data2save, na.action = na.omit)
#fit_HIP <- lm(HIP ~ hemi + id, data = data2save, na.action = na.omit)
#
## Define a function to extract stats from an lm() fit
#extract_lm_stats <- function(fit) {
#  
#  s <- summary(fit)
#  ctab <- s$coefficients
#  idx <- which(rownames(ctab) == "hemiright")
#  
#  # Extract the relevant stats (assuming there's exactly one match)
#  t_value <- ctab[idx, "t value"]
#  p_value <- ctab[idx, "Pr(>|t|)"]
#  
#  # Degrees of freedom for the residuals:
#  df_val <- s$df[2]
#  
#  # Confidence interval
#  ci <- confint(fit)
#  ci_term <- ci[rownames(ci) %in% rownames(ctab)[idx], , drop = FALSE]
#  ci_low <- ci_term[1, 1]
#  ci_high <- ci_term[1, 2]
#  
#  return(list(
#    t = t_value,
#    df = df_val,
#    p = p_value,
#    ci_low = ci_low,
#    ci_high = ci_high
#  ))
#}
#
## Extract stats for each fit
#stats_AUD <- extract_lm_stats(fit_AUD)
#stats_HIP <- extract_lm_stats(fit_HIP)
#
## BH/FDR correction across AUD and HIP only.
#raw_pvals <- c(stats_AUD$p, stats_HIP$p)
#adj_pvals <- p.adjust(raw_pvals, method = "fdr")
#               
#
## Combine into a single results data frame
#results_df <- data.frame(
#  region   = c("AUD", "HIP"),
#  t        = c(stats_AUD$t, stats_HIP$t),
#  df       = c(stats_AUD$df, stats_HIP$df),
#  P.Value  = c(stats_AUD$p, stats_HIP$p),
#  adjust.P = adj_pvals,
#  CI.Low   = c(stats_AUD$ci_low, stats_HIP$ci_low),
#  CI.High  = c(stats_AUD$ci_high, stats_HIP$ci_high)
#)
#
## Write the table to an Excel file
#write.xlsx(results_df, "output/3_Cell_density/Overall_cell/batch1_2/Hemi_results_cell_density.xlsx", sheetName = "overall_cell_density")





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

