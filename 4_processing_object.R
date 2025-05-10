###Processing seurat object

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

# Open dataset --------------
seurat_seq <- readRDS("./step3_output/seuratobject_with_metabolome_Mar2025.rds")

############### Filtering seurat object, based on gene numbers------------
#See number of genes and transcripts captured in datasets
VlnPlot(seurat_seq, features = c("nFeature_RNA")) + geom_hline(yintercept = c(1000, 10000), linetype="dashed", color="red", size =1)
ggsave("./step4_output/filtering_criteria_over1000genes.svg", height = 4, width = 3, bg = "white")

#Subsetting dataset based on RNA features
seurat_seq #279 samples within 4 assays

seurat_subset <- subset(seurat_seq, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000)

seurat_subset #193 samples within 4 assays

p1 <- VlnPlot(seurat_subset, features = c("nFeature_RNA"))
p2 <- VlnPlot(seurat_subset, features = c("nCount_RNA"))
p3 <- VlnPlot(seurat_subset, features = c("nFeature_targeted_metabolite_concentration"), assay="targeted_metabolite_concentration")
p4 <- VlnPlot(seurat_subset, features = c("nCount_targeted_metabolite_concentration"), assay="targeted_metabolite_concentration")

wrap_plots(p1,p2,p3,p4, ncol=4)

ggsave("./step4_output/VlnPlot_afterfiltering.svg", height = 4, width = 10, bg = "white")

############### Run log-normalization and dimension reduction of "RNA"####################

seurat <- NormalizeData(
  seurat_subset,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  margin = 1,
  verbose=TRUE,
  assay = "RNA")

#Select variable genes & scale RNA assay & dimension reduction
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 500, assay = "RNA")
seurat <- ScaleData(seurat, verbose = FALSE, assay = "RNA")
seurat <- RunPCA(seurat, features = seurat@assays$RNA@meta.data$var.features, 
                 verbose = FALSE, npcs = 10, assay = "RNA", reduction.name = "pca_rna", reduction.key = "PCRNA_")
ElbowPlot(seurat, ndims = 10, reduction = "pca_rna")
seurat <- FindNeighbors(seurat, dims = 1:10, assay = "RNA", reduction = "pca_rna")
seurat <- FindClusters(seurat, resolution = 1, pc.use = "pca_rna", cluster.name = "RNA_snn_res.1",graph.name = "RNA_snn")
seurat <- RunUMAP(seurat, dims = 1:10, 
                  min.dist = 0.001, respulsion.strength = 1, n.neighbor = 10, spread = 0.5, 
                  assay = "RNA", reduction = "pca_rna",
                  reduction.key = "UMAPRNA_", reduction.name = "UMAP_RNA")

DimPlot(seurat, reduction = "UMAP_RNA", label = TRUE, repel = TRUE, pt.size = 1.5, group.by = "RNA_snn_res.1")
ggsave("./step4_output/UMAP_RNA_over1000genes_npcs10_nosplit.svg", height = 3, width = 3.5, bg = "white")


############### Dimension reduction of targeted metabolite (concentration) ---------------
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 15, assay = "targeted_metabolite_concentration")
seurat <- ScaleData(seurat, verbose = FALSE, assay = "targeted_metabolite_concentration")
seurat <- RunPCA(seurat, verbose = TRUE, features = seurat@assays$targeted_metabolite_concentration@var.features, 
                 npcs = 7, assay = "targeted_metabolite_concentration", reduction.name = "pca_targeted_metabolite", reduction.key = "PCtargMET_")
ElbowPlot(seurat, ndims = 7, reduction = "pca_targeted_metabolite")
seurat <- FindNeighbors(seurat, dims = 1:7, assay = "targeted_metabolite_concentration", reduction = "pca_targeted_metabolite")
seurat <- FindClusters(seurat, resolution = 1, cluster.name = "MET_targ_snn_res.1", graph.name = "targeted_metabolite_concentration_snn",
                       verbose = T)

seurat <- RunUMAP(seurat, dims = 1:7,
                  assay = "targeted_metabolite_concentration", reduction = "pca_targeted_metabolite",
                  reduction.key = "UMAPtargMET_", reduction.name = "UMAP_targMET")

DimPlot(seurat, reduction = "UMAP_targMET", label = TRUE, repel = TRUE, pt.size = 1.5, group.by = "MET_targ_snn_res.1")

#########save seurat object after appending metabolome------------
saveRDS(seurat, file="./step4_output/seurat_with_metabolome_2025Mar.rds")
