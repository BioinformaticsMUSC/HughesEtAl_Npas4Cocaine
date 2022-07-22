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

#plan("multiprocess", workers = 10)
options(future.globals.maxSize = 6000 * 1024^2)

load("output_initial/Brandon_SeuratObj_SCT_30pcs_05res_Slimmed.RData")

# save as Single Cell Experiment for iSEE
sce <- as.SingleCellExperiment(seuObject_integrated_slim)
saveRDS(sce,"output_initial/Brandon_SCE_SCT_30pcs_05res_Slimmed_Reclust.rds")

# Marker using only WT
saline_subset <- subset(seuObject_integrated_slim, subset = Genotype == "SalineCPP")

Integrated.Markers <- FindAllMarkers(object = saline_subset,
            min.cells.gene = 5,
            only.pos = TRUE, 
            logfc.threshold = 0.3)

Integrated.Markers$cluster <- paste("Cluster",Integrated.Markers$cluster,sep="_")

Integrated.Markers.Sign <- Integrated.Markers %>%
                filter(p_val_adj < 0.05)
                
openxlsx::write.xlsx(Integrated.Markers.Sign, file = "output_initial/SalineCPP_Markers_Stats.xlsx", 
                        colNames = TRUE, 
                        borders = "columns",
                        overwrite=T)
write.table(Integrated.Markers,"output_initial/SalineCPP_Markers_Full_Stats.txt",sep="\t",quote=F,row.names=F)

df <- Integrated.Markers.Sign %>% select(gene,cluster) %>% dplyr::rename(Gene = gene)
write.table(df,"output_initial/SalineCPP_Markers_Stats_DGE.txt",sep="\t",quote=F,row.names=F)

# Run nebulosa for markers genes 
dir.create("output_initial/Nebulosa")

sce_sub <- as.SingleCellExperiment(saline_subset)

gene_list <- c("Adora2a","Adarb2","Slc17a7","Arc","Calb1","Calcr","Cartpt","Cbln4","Chat","Crym","Cd74","Cxcl14","Dio3","Dlk1",
                "Drd1","Drd2","Drd3","Egfp","Fos","Gad1","Gad2","Grm8","Junb","Mobp","Npas4","Npy","Pcdh8",
                "Pdgfra","Pdyn","Peg10","Penk","Pnoc","Ppp1r1b","Pvalb","Sst","Tac1","Vip","Cck","Th")

doPlot <- function(sel_name) 
{
PLOT <- plot_density(sce_sub, sel_name, reduction="UMAP")
print(PLOT)
ggsave(sprintf("output_initial/Nebulosa/%s.pdf", sel_name))
}
lapply(unique(gene_list), doPlot)


