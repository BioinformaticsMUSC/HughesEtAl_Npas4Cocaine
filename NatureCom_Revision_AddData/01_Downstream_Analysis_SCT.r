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

plan("multiprocess", workers = 10)
options(future.globals.maxSize = 6000 * 1024^2)

dir.create("output_initial")

seuObject_filt <- readRDS("input/jessica_final_filtered.rds")
DefaultAssay(seuObject_filt) <- "RNA"


# Data Integration

seuObject_split <- SplitObject(seuObject_filt, split.by = "sample")

seuObject_split <- seuObject_split[c("control", "drug")]

for (i in 1:length(seuObject_split)) {
    seuObject_split[[i]] <- SCTransform(seuObject_split[[i]], 
				    						vars.to.regress = c("nUMI","pMito"), 
											verbose = FALSE)
    }

integ_features <- SelectIntegrationFeatures(object.list = seuObject_split, 
											nfeatures = 4000) 

seuObject_split <- PrepSCTIntegration(object.list = seuObject_split, 
											anchor.features = integ_features)

integ_anchors <- FindIntegrationAnchors(object.list = seuObject_split, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

seuObject_integrated <- IntegrateData(
								anchorset = integ_anchors,
								new.assay.name = "integrated",
								normalization.method = "SCT",
								dims = 1:30,
								k.weight = 100,
								sd.weight = 1,
								eps = 0.5,
								verbose = TRUE
								)

DefaultAssay(seuObject_integrated) <- "integrated"

seuObject_integrated <- RunPCA(object = seuObject_integrated, 
								features=NULL, 
								weight.by.var = TRUE, 
								ndims.print = 1:5, 
								nfeatures.print = 30, 
								npcs = 30, 
								reduction.name = "pca")

seuObject_integrated <- FindNeighbors(object = seuObject_integrated, 
										reduction = "pca", 
										dims = 1:30, 
										nn.eps = 0.5)

seuObject_integrated <- FindClusters(object = seuObject_integrated, 
										resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2), 
										algorithm = 1,
										n.iter = 1000)

pdf("output_initial/Data_Integrated_Clustree.pdf", width = 12, height = 6)
clustree(seuObject_integrated@meta.data, prefix = "integrated_snn_res.", node_colour = "sc3_stability")
dev.off()

# Select resolution and run UMAP
Idents(object = seuObject_integrated) <- "integrated_snn_res.0.5"

seuObject_integrated <- RunUMAP(object = seuObject_integrated, 
										reduction = "pca", 
										dims = 1:30)

# Select the RNA counts slot to be the default assay
DefaultAssay(seuObject_integrated) <- "RNA"
seuObject_integrated <- NormalizeData(object = seuObject_integrated, 
						normalization.method = "LogNormalize", 
						scale.factor = 10000)


seuObject_integrated@meta.data <- seuObject_integrated@meta.data %>%
                                       rownames_to_column("TMP") %>%
                                       select(TMP,orig.ident,sample,pMito,nCount_RNA,nFeature_RNA,nCount_SCT,nFeature_SCT,seurat_clusters) %>%
                                       column_to_rownames("TMP")

save(seuObject_integrated, file = "output_initial/Jessica_SeuratObj_SCT_30pcs_05res.RData")

# 
pdf("output_initial/Data_Integrated_UMAP.pdf", width = 10, height = 6)
p1 <- DimPlot(object = seuObject_integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seuObject_integrated, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="sample")
plot_grid(p1, p2)
dev.off()

seuObject_integrated_slim <- DietSeurat(seuObject_integrated, 
										counts = TRUE, 
										data = TRUE, 
										scale.data = FALSE,
										assays="RNA",
										dimreducs = c("pca","umap"))

save(seuObject_integrated_slim, file = "output_initial/Jessica_SeuratObj_SCT_30pcs_05res_Slimmed.RData")
