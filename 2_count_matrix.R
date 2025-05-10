#Open packages 
library(tidyverse)
library(Seurat)
library(readxl)
library(RColorBrewer)
library(viridis)
library(rcartocolor)
library(patchwork)
library(svglite)
library(ggrepel)
library(cowplot)

############ Open raw (count) data & filter them #############
#Open count-matrices
counttable <- read.table("./rawdata/expected_count_combined_matrix.txt", header = TRUE)
rownames(counttable) <- counttable[,1]
counttable <- counttable[,-1]

#rename file names
colnames(counttable) <- gsub("X","",colnames(counttable))

#filter tables based on cellcelector pictures
picture <- read_excel("./database/Cells_tobe_removed_byAnh_Feb2025.xlsx")
picture$Cell #they are empty or duplet wells 

names_tobe_removed <- colnames(counttable) %in% picture$Cell
names_tobe_kept <- !names_tobe_removed
counttable <- counttable[,names_tobe_kept] #279 variables kept

############ Import filtered count matrices into Seurat #############
#create seurat object from count matrices
counttable <- CreateSeuratObject(counts = counttable)

############ save RDS ####################
saveRDS(counttable, file="./step2_output/count_Mar2025.rds")
