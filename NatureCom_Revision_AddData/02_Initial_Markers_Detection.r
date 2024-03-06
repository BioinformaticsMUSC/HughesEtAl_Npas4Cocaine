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
library(RColorBrewer)
library(MetBrewer)
library(scales)
library(Nebulosa)
})

plan("multiprocess", workers = 10)
options(future.globals.maxSize = 6000 * 1024^2)

load("output_initial/Jessica_SeuratObj_SCT_30pcs_05res.RData")

# Differential Expression
markers<- presto::wilcoxauc(seuObject_integrated, 'seurat_clusters', assay = 'data')

markers$group <- paste("Cluster",markers$group,sep="_")

markers_sign <- markers %>%
                filter(padj < 0.05, pct_in > 20, pct_out < 20) %>% 
                dplyr::rename(Gene = feature, Cluster = group)

openxlsx::write.xlsx(markers_sign, 
                        file = "output_initial/Jessica_Markers_Significant.xlsx", 
                        colNames = TRUE, 
                        borders = "columns",
                        overwrite=T)

write.table(markers_sign,"output_initial/Jessica_Markers_Stats_DGE.txt",sep="\t",quote=F,row.names=F)


df <- markers_sign %>% select(Gene,Cluster)
write.table(df,"output_initial/Jessica_Markers_Stats_DGE.txt",sep="\t",quote=F,row.names=F)


# Run nebulosa for markers genes 
dir.create("output_initial/Nebulosa")

sce_sub <- as.SingleCellExperiment(seuObject_integrated_slim)

gene_list <- c("Slc17a7","Adora2a","Adarb2","Slc17a7","Arc","Calb1","Calcr","Cartpt","Cbln4","Chat","Crym","Cd74","Cxcl14","Dio3","Dlk1",
                "Drd1","Drd2","Drd3","Fos","Gad1","Gad2","Grm8","Junb","Mobp","Npas4","Npy","Pcdh8",
                "Pdgfra","Pdyn","Peg10","Penk","Pnoc","Ppp1r1b","Pvalb","Sst","Tac1","Vip","Cck","Th")

doPlot <- function(sel_name) 
{
PLOT <- plot_density(sce_sub, sel_name, reduction="UMAP")
print(PLOT)
ggsave(sprintf("output_initial/Nebulosa/%s.pdf", sel_name))
}
lapply(unique(gene_list), doPlot)


