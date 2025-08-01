## This script is to change the niches order and plot the cell type annotation

## created by M.-Y. Wang 19-Sep-2024; updated 07-OCT-2024

# Remove all objects created before to prevent clashing
rm(list = ls())

setwd("/data/workspaces/lag/workspaces/lg-func-asym/working_data/Mengyun")

args <- commandArgs(trailingOnly=TRUE)
roi <- args[1]  # This is the filenames "AC", "hip"


print(roi)


library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(openxlsx)
library(purrr)
library(parallel)
library(future)
plan("multicore", workers = 40) #

if (roi=="AC") {
  mouse_id <- c("M670", "M671", "M672", "M673", "M234", "M253", "M071", "M083", "M650", "M638", "M076", "M236",
                "F679", "F680", "F682", "F683", "F687", "F688", "F078", "F087", "F090")
  # niche_swap_left <- list(
  #   M669 = c(1, 2, 4, 5, 3),
  #   M670 = c(5, 2, 1, 3, 4),
  #   M671 = c(5, 2, 1, 3, 4),
  #   M672 = c(5, 1, 2, 3, 4),
  #   M673 = c(1, 3, 4, 5, 2),
  #   F679 = c(5, 2, 4, 3, 1),
  #   F680 = c(4, 1, 5, 3, 2),
  #   F681 = c(1, 3, 4, 2, 5),
  #   F682 = c(2, 5, 3, 4, 1),
  #   F683 = c(5, 1, 4, 2, 3),
  #   F686 = c(2, 1, 4, 5, 3),
  #   F687 = c(1, 4, 2, 3, 5),
  #   F688 = c(1, 4, 2, 3, 5)
  # )
  # 
  # niche_swap_right <- list(
  #   M669 = c(3,5,4,1,2),
  #   M670 = c(1,3,4,5,1),
  #   M671 = c(3,4,1,2,5),
  #   M672 = c(5,4,2,1,3),
  #   M673 = c(4,3,2,5,1),
  #   F679 = c(2,3,4,1,5),
  #   F680 = c(3,1,4,5,2),
  #   F681 = c(5,4,2,3,1),
  #   F682 = c(5,2,4,1,3),
  #   F683 = c(1,2,3,4,5),
  #   F686 = c(2,3,5,1,4),
  #   F687 = c(3,4,5,2,1),
  #   F688 = c(3,5,2,1,4)
  # )
  
} else {
  mouse_id <- c("M669", "M670", "M671", "M672", "M674", "M676", "M677", "M678", "M234", "M253", "M071", "M083", "M650", "M638", "M076", "M236",
                "F679", "F680", "F682", "F683", "F685", "F686", "F687", "F688", "F073", "F078", "F087", "F090")
  # mouse_id <- "M669"
}

Plot_cell_type <- function(mouse_id){
  mouse_asym_left <- readRDS(file = paste0("intermedia_data/cell_type/",
                                           mouse_id,"_l",roi,"_cell_type_harmony.rds")) 
  
  mouse_asym_right <- readRDS(file = paste0("intermedia_data/cell_type/",
                                            mouse_id, "_r", roi, "_cell_type_harmony.rds"))
# change the niches numbers
  
  # if (roi=="AC") {
  #   mouse_asym_left@meta.data <- mouse_asym_left@meta.data %>%
  #     mutate(
  #       niches_new = case_when(
  #         niches == niche_swap_left[[mouse_id]][1] ~ 1,
  #         niches == niche_swap_left[[mouse_id]][2] ~ 2,
  #         niches == niche_swap_left[[mouse_id]][3] ~ 3,
  #         niches == niche_swap_left[[mouse_id]][4] ~ 4,
  #         niches == niche_swap_left[[mouse_id]][5] ~ 5
  #       ),
  #       niches_new = factor(niches_new, levels = c(1, 2, 3, 4, 5))
  #     )
  #   cell_per_layer <-table(mouse_asym_left$predicted.celltype, mouse_asym_left$niches_new) # cell types in each layer
  #   cell_per_layer_df <- as.data.frame.matrix(cell_per_layer)
  #   cell_per_layer_df$total <- rowSums(cell_per_layer_df)
  # 
  #   cell_type2save <- paste0("output/2_Cell_type/acr_sample/AUD/", mouse_id,"_AC_cell_type_harmony.xlsx")
  #   write.xlsx(cell_per_layer_df, cell_type2save, sheetName = "left", rowNames = TRUE, colNames = TRUE)
  # 
  # 
  #   mouse_asym_right@meta.data <- mouse_asym_right@meta.data %>%
  #     mutate(
  #       niches_new = case_when(
  #         niches == niche_swap_right[[mouse_id]][1] ~ 1,
  #         niches == niche_swap_right[[mouse_id]][2] ~ 2,
  #         niches == niche_swap_right[[mouse_id]][3] ~ 3,
  #         niches == niche_swap_right[[mouse_id]][4] ~ 4,
  #         niches == niche_swap_right[[mouse_id]][5] ~ 5
  #       ),
  #       niches_new = factor(niches_new, levels = c(1, 2, 3, 4, 5))
  #     )
  #   cell_per_layer <- table(mouse_asym_right$predicted.celltype, mouse_asym_right$niches_new) # cell types in each layer
  #   cell_per_layer_df <- as.data.frame.matrix(cell_per_layer)
  #   cell_per_layer_df$total <- rowSums(cell_per_layer_df)
  #   cell_per_layer_df <- cbind(cell = rownames(cell_per_layer_df), cell_per_layer_df)
  # 
  #   cell_type2save <-
  #     paste0("output/2_Cell_type/acr_sample/AUD/",
  #            mouse_id,
  #            "_AC_cell_type_harmony.xlsx")
  #   wb <- loadWorkbook(cell_type2save)
  #   addWorksheet(wb, "right")
  #   writeData(wb, "right", cell_per_layer_df)
  #   saveWorkbook(wb, cell_type2save, overwrite = TRUE)
  # 
  # }
  # Define a unified color palette for shared cell types
  
  # common_colors <- c(
  #   "Astro" = "#DCDCDC",       # Gainsboro (light grey)
  #   "Endo" = "#E2E2E2",        # Bright gold (unchanged)
  #   "Car3" = "#B5B5B5",        # Bright cyan (unchanged)
  #   "Lamp5" = "#D3D3D3",       # Light grey
  #   "Oligo" = "#C0C0C0",       # Silver (lighter grey)
  #   "Pvalb" = "#B0B0B0",       # Medium grey
  #   "Sncg" = "#A9A9A9",        # Dark grey
  #   "Sst" = "#E0E0E0",         # Very light grey
  #   "Vip" = "#D8D8D8",         # Lighter grey
  #   "CR" = "#F0F0F0",          # Very light grey
  #   "Micro-PVM" = "#F5F5F5",   # White smoke (very light grey)
  #   "NP SUB" = "#E8E8E8",      # Very light grey
  #   "SMC-Peri" = "#C8C8C8",    # Light grey with slight green tone
  #   "SUB-ProS" = "#BEBEBE",    # Distinct grey (updated)
  #   "Sst Chodl" = "#E5E5E5"    # Very light grey
  # )
  
  common_colors <- c(
    "Astro"       = "#E41A1C",  # red
    "Endo"        = "#377EB8",  # blue
    "Car3"        = "#4DAF4A",  # green
    "Lamp5"       = "#984EA3",  # purple
    "Oligo"       = "#FF7F00",  # orange
    "Pvalb"       = "#FFFF33",  # yellow
    "Sncg"        = "#A65628",  # brown
    "Sst"         = "#F781BF",  # pink
    "Vip"         = "#999999",  # grey-neutral (just one)
    "CR"          = "#66C2A5",  # teal
    "Micro-PVM"   = "#FC8D62",  # salmon
    "NP SUB"      = "#8DA0CB",  # light blue
    "SMC-Peri"    = "#E78AC3",  # magenta
    "SUB-ProS"    = "#A6D854",  # light green
    "Sst Chodl"   = "#FFD92F"   # gold
  )
  
  # Colors for layers (distinct colors + gradients within layers)
  layer_colors <- c(
    # L2-3 (Blue shades)
    "L2-3 IT CTX" = "#1E90FF",    # Dodger blue for L2-3
    
    # L4-5 (Light orange for L4, Red gradient for L5)
    "L4-5 IT CTX" = "#FFA07A",    # Light salmon (orange) for L4-5 IT
    "L5 PT CTX" = "#FF6347",      # L5 PT: Tomato red
    "L5 IT CTX" = "#FF4500",      # L5 IT: Orange-red
    "L5-6 IT TPE-ENT" = "#FF0000",# L5-6 IT: Pure red
    "L5-6 NP CTX" = "#E9967A",    # L5-6 NP: Dark salmon
    
    # L6 (Purple gradient)
    "L6 CT CTX" = "#9400D3",      # Dark violet for L6 CT
    "L6 IT CTX" = "#9370DB",      # Medium purple for L6 IT (updated)
    "L6b CTX" = "#BA55D3"         # Orchid for L6b (slightly lighter purple)
  )
  
  if (roi == "AC") {
    cell_color_mapping <- c(
      common_colors,  # Use common colors for greyed-out cell types
      layer_colors    # Apply layer colors and gradients
    )
  } else {
    cell_color_mapping <- c(
      common_colors,  # Use common colors for greyed-out cell types
      "CA1-ProS" = "#008000",     # Green for CA1
      "CA2-IG-FC" = "#8A2BE2",    # Blue-violet for CA2 (unchanged)
      "CA3" = "#7FFFD4",          # Aquamarine for CA3
      "DG" = "#0000FF"            # Blue for DG
    )
  }
  
  

  # if (roi=="AC") {
  # cell_color_mapping <- c("Astro" = "#FF0000", "Endo" = "#32CD32", "L2-3 IT CTX" = "#0000FF", "L4-5 IT CTX" = "#EE82EE", "L5 PT CTX" = "#FFFF00",
  #                         "L5 IT CTX" = "#FFA500", "L6 CT CTX" = "#006400", "L6 IT CTX" = "#ADFF2F", "L6b CTX" = "#FF00FF", "Lamp5" = "#A52A2A",
  #                         "Car3" = "#00FFFF", "L5-6 IT TPE-ENT" = "#8A2BE2", "L5-6 NP CTX" = "#7FFFD4", "Oligo" = "#D2691E", "Pvalb" = "#00008B",
  #                         "Sncg" = "#9ACD32", "Sst" = "#4B0082", "Vip" = "#F08080", "Sst Chodl" = "#20B2AA")
  # 
  # } else {
  # cell_color_mapping <- c("Astro" = "#EE82EE", "Endo" = "#FFD700", "Lamp5" = "#A52A2A",
  #                         "Oligo" = "#D2691E", "Pvalb" = "#00008B", "Sncg" = "#9ACD32",
  #                         "Sst" = "#4B0082", "Vip" = "#F08080", "CA1-ProS" =  "#008000",
  #                         "CA2-IG-FC" = "#8A2BE2", "CA3" = "#7FFFD4", "CR" = "#FF0000",
  #                         "DG" = "#0000FF", "Micro-PVM" = "#FFFF00", "NP SUB" = "#FFA500",
  #                         "SMC-Peri" = "#2E8B57", "SUB-ProS" = "#ADFF2F")
  # }
  
  celltype_plot_left <-
    ImageDimPlot(
      mouse_asym_left,
      group.by = "predicted.celltype",
      size = 1.2,
      cols = cell_color_mapping,
      dark.background = F
    ) + 
    ggtitle("")
  left_side <- celltype_plot_left
  
  
  
  celltype_plot_right <-
    ImageDimPlot(
      mouse_asym_right,
      group.by = "predicted.celltype",
      size = 1.2,
      cols = cell_color_mapping,
      dark.background = F
    ) + 
    ggtitle("")
    #  coord_flip()
  
  right_side <- celltype_plot_right
  
  # if (roi=="AC"){
  #   
  #   niche.plot_left <-
  #     ImageDimPlot(
  #       mouse_asym_left,
  #       group.by = "niches",
  #       size = 1.2,
  #       dark.background = F
  #     ) +
  #     scale_fill_manual(values = c("1"="#442288", "2"= "#6CA2EA","3"= "#B5D33D", "4" = "#FED23F","5"= "#EB7D5B")) +
  #     ggtitle("")
  #   
  #   niche.plot_right <-
  #     ImageDimPlot(
  #       mouse_asym_right,
  #       group.by = "niches",
  #       size = 1.2,
  #       dark.background = F
  #     ) +
  #     scale_fill_manual(values = c("1"="#442288", "2"= "#6CA2EA","3"= "#B5D33D", "4" = "#FED23F","5"= "#EB7D5B")) +
  #     ggtitle("")
  # } else {
  #   niche.plot_left <-
  #     ImageDimPlot(
  #       mouse_asym_left,
  #       group.by = "niches",
  #       size = 1.2,
  #       dark.background = F
  #     ) +
  #     scale_fill_manual(values = c("1"="#442288", "2"= "#6CA2EA","3"= "#B5D33D", "4" = "#FED23F","5"= "#EB7D5B")) +
  #     ggtitle("")
  #   
  #   niche.plot_right <-
  #     ImageDimPlot(
  #       mouse_asym_right,
  #       group.by = "niches",
  #       size = 1.2,
  #       dark.background = F
  #     ) +
  #     scale_fill_manual(values = c("1"="#442288", "2"= "#6CA2EA","3"= "#B5D33D", "4" = "#FED23F","5"= "#EB7D5B")) +
  #     ggtitle("")
  # }
  

  
  #plot_annotation <- celltype.plot | niche.plot
  plot_annotation <- right_side | left_side
  
  
  # Define file paths based on roi and save & rotate
  if (roi == "AC") {
    plot2save_left <- paste0("output/2_Cell_type/acr_sample/batch1_2/AUD/Plot_", mouse_id, "_l", roi, "_cell_type_harmony.pdf")
    ggsave(plot2save_left, plot = left_side, device   = cairo_pdf, width = 5, height = 5)
    
    plot2save_right <- paste0("output/2_Cell_type/acr_sample/batch1_2/AUD/Plot_", mouse_id, "_r", roi, "_cell_type_harmony.pdf")
    ggsave(plot2save_right, plot = right_side, device   = cairo_pdf, width = 5, height = 5)
    
  } else {
    plot2save_left <- paste0("output/2_Cell_type/acr_sample/batch1_2/hip/Plot_", mouse_id, "_l", roi, "_cell_type_harmony.pdf")
    ggsave(plot2save_left, plot = left_side, device   = cairo_pdf, width = 5, height = 5)
    
    plot2save_right <- paste0("output/2_Cell_type/acr_sample/batch1_2/hip/Plot_", mouse_id, "_r", roi, "_cell_type_harmony.pdf")
    ggsave(plot2save_right, plot = right_side, device   = cairo_pdf, width = 5, height = 5)
  }
  
  # Define file paths based on roi and save & rotate
#   if (roi == "AC") {
#     plot2save_left <- paste0("output/2_Cell_type/acr_sample/AUD/Plot_", mouse_id, "_l", roi, "_layers_harmony.png")
#     ggsave(plot2save_left, plot = niche.plot_left, width = 5, height = 5, units = 'in', dpi = 300)
#     
#     plot2save_right <- paste0("output/2_Cell_type/acr_sample/AUD/Plot_", mouse_id, "_r", roi, "_layers_harmony.png")
#     ggsave(plot2save_right, plot = niche.plot_right, width = 5, height = 5, units = 'in', dpi = 300)
#     
#   } else {
#     plot2save_left <- paste0("output/2_Cell_type/acr_sample/HIP/Plot_", mouse_id, "_l", roi, "_layers_harmony.png")
#     ggsave(plot2save_left, plot = niche.plot_left, width = 5, height = 5, units = 'in', dpi = 300)
#     
#     plot2save_right <- paste0("output/2_Cell_type/acr_sample/HIP/Plot_", mouse_id, "_r", roi, "_layers_harmony.png")
#     ggsave(plot2save_right, plot = niche.plot_right, width = 5, height = 5, units = 'in', dpi = 300)
#   }
}



mclapply(mouse_id, 
       function (x) Plot_cell_type(x),
       mc.cores = 13
       )




