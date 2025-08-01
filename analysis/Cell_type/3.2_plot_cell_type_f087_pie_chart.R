## This script is to plot the pie chart for the cell type in AUD

## created by M.-Y. Wang 30-July-2025

# Remove all objects created before to prevent clashing
rm(list = ls())

setwd("/data/workspaces/lag/workspaces/lg-func-asym/working_data/Mengyun")


#################### Step 0.  ROI & colour map 
roi <- "AC"          # or "AUD", etc.

common_colors <- c(
  "Astro"       = "#E41A1C",
  "Endo"        = "#377EB8",
  "Car3"        = "#4DAF4A",
  "Lamp5"       = "#984EA3",
  "Oligo"       = "#FF7F00",
  "Pvalb"       = "#FFFF33",
  "Sncg"        = "#A65628",
  "Sst"         = "#F781BF",
  "Vip"         = "#999999",
  "CR"          = "#66C2A5",
  "Micro-PVM"   = "#FC8D62",
  "NP SUB"      = "#8DA0CB",
  "SMC-Peri"    = "#E78AC3",
  "SUB-ProS"    = "#A6D854",
  "Sst Chodl"   = "#FFD92F"
)

layer_colors <- c(
  "L2-3 IT CTX"      = "#1E90FF",
  "L4-5 IT CTX"      = "#FFA07A",
  "L5 PT CTX"        = "#FF6347",
  "L5 IT CTX"        = "#FF4500",
  "L5-6 IT TPE-ENT"  = "#FF0000",
  "L5-6 NP CTX"      = "#E9967A",
  "L6 CT CTX"        = "#9400D3",
  "L6 IT CTX"        = "#9370DB",
  "L6b CTX"          = "#BA55D3"
)

cell_color_mapping <- if (roi == "AC") {
  c(common_colors, layer_colors)
} else {
  c(common_colors,
    "CA1-ProS"  = "#008000",
    "CA2-IG-FC" = "#8A2BE2",
    "CA3"       = "#7FFFD4",
    "DG"        = "#0000FF")
}

#################### Step 1.  Data  
cell_counts <- data.frame(
  cell_type = c("Astro","Car3","Endo","L2-3 IT CTX","L4-5 IT CTX",
                "L5 IT CTX","L5 PT CTX","L5-6 IT TPE-ENT","L5-6 NP CTX",
                "L6 CT CTX","L6 IT CTX","L6b CTX","Lamp5","Oligo",
                "Pvalb","Sncg","Sst","Sst Chodl","Vip"),
  count     = c(354,83,305,664,644,163,125,12,58,
  395,139,51,34,351,157,17,85,4,26), # left one

  # count     = c(340,86,361,651,676,148,109,19,66,
  #               396,165,57,27,383,171,12,89,5,43), # right one
  # 
  stringsAsFactors = FALSE
)

display_name_map <- c(
  "Astro","Car3","Endo","L2-3 IT CTX","L4-5 IT CTX",
  "L5 IT CTX","L5 PT CTX","L5-6 IT TPE-ENT","L5-6 NP CTX",
  "L6 CT CTX","L6 IT CTX","L6b CTX","Lamp5","Oligo",
  "Pvalb","Sncg","Sst","Sst Chodl","Vip"
)

cell_counts$display <- ifelse(cell_counts$cell_type %in% names(display_name_map),
                              display_name_map[cell_counts$cell_type],
                              cell_counts$cell_type)


#################### Step2.  Percentages & legend labels
cell_counts$prop <- cell_counts$count / sum(cell_counts$count)

# Order slices (and legend) by descending percentage
cell_counts <- cell_counts[order(-cell_counts$prop), ]
cell_counts$cell_type <- factor(cell_counts$cell_type, levels = cell_counts$cell_type)

# Legend label: " 9.94%  Astrocyte"
legend_lab <- sprintf("%6.2f%%  %s", 100 * cell_counts$prop, cell_counts$display)


#################### Step 3. Colors (subset to the ones we actually need)
needed_colors <- cell_color_mapping[as.character(cell_counts$cell_type)]


#################### Step 4.  Plot  
library(ggplot2)

p <- ggplot(cell_counts, aes(x = "", y = count, fill = cell_type)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(
    values = needed_colors,
    breaks = levels(cell_counts$cell_type),  # keep our chosen order
    labels = legend_lab,
    name   = NULL
  ) +
  theme_void() +
  theme(
    legend.position   = "right",
    legend.key.width  = unit(10, "pt"),
    legend.key.height = unit(15, "pt"),
    legend.text       = element_text(family = "Arial", size = 10),
    plot.margin       = margin(5.5, 40, 5.5, 5.5)
  ) +
  ggtitle("Cell-type composition")

#################### Step 5.  Save as vector PDF
outfile <- paste0("output/2_Cell_type/acr_sample/batch1_2/AUD/Plot_F087_l", roi, "_cell_type_pie.pdf") # left hemi
# outfile <- paste0("output/2_Cell_type/acr_sample/batch1_2/AUD/Plot_F087_r", roi, "_cell_type_pie.pdf") #right hemi

ggsave(outfile,
       plot   = p,
       device = cairo_pdf,
       width  = 5, height = 5)






