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
library(schex)
})

load("output_reclust/Brandon_SCT_30pcs_04res_Reclust_NoDoublet_Slimmed.RData")

saline_subset <- subset(seuObject_slim_nodoub_slim, subset = Genotype == "SalineCPP")


# Differential Expression
Integrated.Markers <- FindAllMarkers(object = saline_subset,
            min.cells.gene = 5,
            only.pos = TRUE, 
            logfc.threshold = 0.3)

Integrated.Markers$cluster <- paste("Cluster",Integrated.Markers$cluster,sep="_")

Integrated.Markers.Sign <- Integrated.Markers %>%
                filter(p_val_adj < 0.05)
openxlsx::write.xlsx(Integrated.Markers.Sign, file = "output_reclust/SalineCPP_Markers_Reclust_Stats.xlsx", 
                        colNames = TRUE, 
                        borders = "columns",
                        overwrite=T)
write.table(Integrated.Markers,"output_reclust/SalineCPP_Markers_Reclust_Full_Stats.txt",sep="\t",quote=F,row.names=F)


df <- Integrated.Markers.Sign %>% select(gene,cluster) %>% dplyr::rename(Gene = gene)
write.table(df,"output_reclust/SalineCPP_Markers_Reclust_Stats_DGE.txt",sep="\t",quote=F,row.names=F)

# Run Enrichment
system("R CMD BATCH --vanilla Fisher_enrich_Markers_NoGlut_CellRep.R")


sce <- as.SingleCellExperiment(seuObject_slim_nodoub_slim)
saveRDS(sce,"output_reclust/Brandon_SCE_SCT_30pcs_05res_Slimmed_Reclust.rds")

# Nebulosa 
dir.create("output_reclust/Nebulosa")

sce_sub <- as.SingleCellExperiment(saline_subset)

gene_list <- c("Adora2a","Arc","Calb1","Calcr","Cartpt","Cbln4","Chat","Crym","Cd74","Cxcl14","Dio3","Dlk1",
                "Drd1","Drd2","Drd3","Egfp","Fos","Gad1","Gad2","Grm8","Junb","Mobp","Npas4","Npy","Pcdh8",
                "Pdgfra","Pdyn","Peg10","Penk","Pnoc","Ppp1r1b","Pvalb","Sst","Tac1","Vip","Cck","Th")

doPlot <- function(sel_name) 
{
PLOT <- plot_density(sce_sub, sel_name, reduction="UMAP")
print(PLOT)
ggsave(sprintf("output_reclust/Nebulosa/%s.pdf", sel_name))
}
lapply(unique(gene_list), doPlot)

# Schex 
dir.create("output_reclust/Schex")

sce_sub <- as.SingleCellExperiment(seuObject_slim_nodoub_slim)

sce_sub <- make_hexbin(sce_sub, nbins = 100, dimension_reduction = "UMAP", use_dims=c(1,2))

gene_list <- c("Adora2a","Arc","Calb1","Calcr","Cartpt","Cbln4","Chat","Crym","Cd74","Cxcl14","Dio3","Dlk1",
                "Drd1","Drd2","Drd3","Egfp","Fos","Gad1","Gad2","Grm8","Junb","Mobp","Npas4","Npy","Pcdh8",
                "Pdgfra","Pdyn","Peg10","Penk","Pnoc","Ppp1r1b","Pvalb","Sst","Tac1","Vip","Cck","Th")

doPlot <- function(sel_name) 
{
PLOT <- plot_hexbin_feature(sce_sub, mod="RNA", type="logcounts", feature=sel_name, action="mean", xlab="UMAP1", ylab="UMAP2",upper_cutoff=0.2,lower_cutoff=0.1)
print(PLOT)
ggsave(sprintf("output_reclust/Schex/%s.pdf", sel_name))
}
lapply(unique(gene_list), doPlot)
