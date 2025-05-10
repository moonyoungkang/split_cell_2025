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
library(dplyr)

#open seurat object
seurat <- readRDS(file = "./step4_output/seurat_with_metabolome_2025Mar.rds")

#check the number of genes detected
RNA <- seurat@assays$RNA
RNA <- RNA$counts
RNA <- as.data.frame(RNA)
RNA %>%
  filter(if_any(where(is.numeric), ~ . != 0)) %>%
  nrow() #21258 genes detected across the cells that are analyzed

#########Draw lines to connect corresponding cells in DimPlot (Metabolome-guided cell annotation)------------------------------
p1 <- DimPlot(seurat, reduction = "UMAP_RNA", label = TRUE, repel = TRUE, pt.size = 1, group.by = "cell_type_by_targeted_met")
p2 <- DimPlot(seurat, reduction = "UMAP_targMET", label = TRUE, repel = TRUE, pt.size = 1, group.by = "cell_type_by_targeted_met")

wrap_plots(p1, p2, ncol =2)
ggsave("./step5_output/UMAP_seurat_annotationbymet.svg", height = 3, width = 8, bg = "white")

#extract UMAP coordinates for each cell
UMAPRNA_df <- seurat@reductions$UMAP_RNA@cell.embeddings
UMAPtargMET_df <- seurat@reductions$UMAP_targMET@cell.embeddings

UMAPRNA_df <- as.data.frame(UMAPRNA_df)
UMAPtargMET_df <- as.data.frame(UMAPtargMET_df)

#extract cell type annotation information for each cell
cell_type <- as.data.frame(seurat@meta.data)
cell_type <- cell_type[,c("cell_name", "cell_type_by_targeted_met")]
colnames(cell_type) <- c("cell_name","cell_type")
head(cell_type)
cell_type$cell_type[is.na(cell_type$cell_type)] <- "unassigned"
cell_type

#matching size
UMAPtargMET_df[,1] <- UMAPtargMET_df[,1]*0.18
UMAPtargMET_df[,2] <- UMAPtargMET_df[,2]*(1/8)

#adding some values to x axis
UMAPtargMET_df[,1] <- UMAPtargMET_df[,1]+10
UMAPtargMET_df[,2] <- UMAPtargMET_df[,2]+1

#merging dataframe (changing axis names)
head(UMAPRNA_df)
head(UMAPtargMET_df)

colnames(UMAPRNA_df) <- c("UMAP_1","UMAP_2")
colnames(UMAPtargMET_df) <- c("UMAP_1","UMAP_2")

UMAPRNA_df$cell_name <- rownames(UMAPRNA_df)
UMAPtargMET_df$cell_name <- rownames(UMAPtargMET_df)

UMAP_merged <- rbind(UMAPRNA_df, UMAPtargMET_df)

#attaching cell type information
head(UMAP_merged)
head(cell_type)
UMAP_merged <- left_join(UMAP_merged, cell_type, by = "cell_name")
UMAP_merged

#draw plot
ggplot(UMAP_merged, aes(x = UMAP_1, y = UMAP_2, group = cell_name, color = cell_type)) +
  geom_point(size = 1.5, shape = 16) +
  geom_line(size = 0.2, alpha = 0.3) +  
  theme_void()+
  scale_color_manual(values = c("#CCBB44","#4477AA","#EE6677","#BBBBBB"))


ggsave("./step5_output/UMAP_connected_cells_UMAPRNA_andUMAPtargRNA.svg", height = 3, width = 10, bg = "white")

