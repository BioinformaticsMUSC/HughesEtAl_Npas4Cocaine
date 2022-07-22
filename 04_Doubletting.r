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
library(scds)
})
source("utils/Utils.R")


load("output_reclust/Brandon_SCE_SCT_30pcs_05res_Slimmed_Reclust.RData")

sce <- as.SingleCellExperiment(seuObject_integrated_reclust_slim)

# Doubleting by genotype
sce <- scDblFinder(sce,
                    samples="Genotype", 
                    BPPARAM=MulticoreParam(3),
                    nfeatures = 3000,
                    dims = 30,
                    dbr.sd = 1)

# add in the column of seurat

seuObject_integrated_reclust_slim@meta.data$Doublets <- sce$scDblFinder.class

seuObject_slim_nodoub <- subset(seuObject_integrated_reclust_slim, subset = Doublets == "singlet")

seuObject_slim_nodoub <- processing_seurat_sctransform(seuObject_slim_nodoub, 
                                                vars_to_regress = c("nCount_RNA","pMito"), 
                                                npcs = 50, 
                                                res = 0.3)

DefaultAssay(seuObject_slim_nodoub) <- "RNA"
seuObject_slim_nodoub <- NormalizeData(object = seuObject_slim_nodoub, 
                        normalization.method = "LogNormalize", 
                        scale.factor = 10000)

save(seuObject_slim_nodoub, file = "output_reclust/Brandon_SCT_30pcs_04res_Reclust_NoDoublet.RData")


# Slimming the files
seuObject_slim_nodoub_slim <- DietSeurat(seuObject_slim_nodoub, 
                                        counts = TRUE, 
                                        data = TRUE, 
                                        scale.data = FALSE,
                                        assays="RNA",
                                        dimreducs = c("pca","umap"))


pdf("output_reclust/Data_SCT_NoDoublet_UMAP.pdf", width = 10, height = 6)
p1 <- DimPlot(object = seuObject_slim_nodoub_slim, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seuObject_slim_nodoub_slim, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

save(seuObject_slim_nodoub_slim, file = "output_reclust/Brandon_SCT_30pcs_04res_Reclust_NoDoublet_Slimmed.RData")



