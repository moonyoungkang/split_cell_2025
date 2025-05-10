#Appending metabolome to seurat object
# Open packages --------------
library(tidyverse)
library(stringr)
library(Seurat)
library(readxl)
library(RColorBrewer)
library(viridis)
library(rcartocolor)
library(patchwork)
library(svglite)
library(ggrepel)
library(cowplot)
library(data.table)
library(Matrix)

# Open dataset ------------
seurat_seq <- readRDS(file = "./step2_output/count_Mar2025.rds") #This is seurat object made from count matrices

assay_peakarea <- readRDS(file="./step1_output/targeted_metabolite_peakarea_Mar2025.rds") #This is cell x targeted metabolome table
assay_concentration <- readRDS(file = "./step1_output/targeted_metabolite_concentration_Mar2025.rds") #This is cell x targeted metabolome table

# Append metabolite assay to seurat object----------------------------
# Append targeted metabolome (peak_area) to seurat object 
seurat_seq[["targeted_metabolite_peakarea"]] <- assay_peakarea

# Append targeted metabolome (concentration) to seurat object 
seurat_seq[["targeted_metabolite_concentration"]] <- assay_concentration

# Append cell type annotation by metabolite to seurat object ------------------
# open cell type annotation by metabolite
cell_type_met <- read.csv("./step1_output/cell_type_manual_by_Anh.csv", row.names = 1)
head(cell_type_met)
rownames(cell_type_met) <- cell_type_met$cell_name.cell_name
cell_type_met$Component.Name <- NULL
cell_type_met$cell_name.cell_name <- NULL
colnames(cell_type_met) <- "cell_type_by_targeted_met"
cell_type_met$cell_name <- rownames(cell_type_met)

# append manual cell type to seurat object (As metadata)
metadata <- seurat_seq@meta.data
metadata$cell_name <- rownames(metadata)

metadata <- metadata %>%
  left_join(cell_type_met, by = "cell_name")

head(metadata)
rownames(metadata) <- metadata$cell_name
metadata

seurat_seq@meta.data <- metadata #append metadata back to seurat object


### save seurat object (metabolome appended) --------
saveRDS(seurat_seq, file = "./step3_output/seuratobject_with_metabolome_Mar2025.rds")
