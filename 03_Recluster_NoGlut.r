suppressPackageStartupMessages({
library(tidyverse)
library(ggrepel)
library(emmeans)
library(SingleCellExperiment)
library(scater)
library(BiocParallel)
library(ggpubr)
library(speckle)
library(magrittr)
library(broom)
library(muscat)
library(Seurat)
library(clustree)
library(leiden)
library(data.table)
library(cowplot)
library(scDblFinder)
library(BiocSingular)
})
source("utils/Utils.R")

options(future.globals.maxSize = 6000 * 1024^2)

dir.create("output_sct_reclust")


load("output_initial/Brandon_SeuratObj_SCT_30pcs_05res_Slimmed.RData")

# remove glutamatergic neurons from different region. 
seuObject_integrated_subset <- subset(seuObject_integrated_slim, idents = c(7,19,21,24),invert=TRUE)

seuObject_integrated_reclust <- processing_seurat_sctransform(seuObject_integrated_subset, 
															 vars_to_regress = c("nCount_SCT","pMito"), 
															 npcs = 30, 
															 res = 0.5)

# Select the RNA counts slot to be the default assay
DefaultAssay(seuObject_integrated_reclust) <- "RNA"
seuObject_integrated_reclust <- NormalizeData(object = seuObject_integrated_reclust, 
						normalization.method = "LogNormalize", 
						scale.factor = 10000)

save(seuObject_integrated_reclust, file = "output_reclust/Brandon_SCE_SCT_30pcs_05res_Reclust.RData")


# Slimming the files
seuObject_integrated_reclust_slim <- DietSeurat(seuObject_integrated_reclust, 
										counts = TRUE, 
										data = TRUE, 
										scale.data = FALSE,
										assays="RNA",
										dimreducs = c("pca","umap"))

# 
pdf("output_reclust/Data_Integrated_UMAP.pdf", width = 10, height = 6)
p1 <- DimPlot(object = seuObject_integrated_reclust_slim, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seuObject_integrated_reclust_slim, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

save(seuObject_integrated_reclust_slim, file = "output_reclust/Brandon_SCE_SCT_30pcs_05res_Slimmed_Reclust.RData")
