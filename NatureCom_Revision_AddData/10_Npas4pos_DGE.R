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
library(data.table)
library(cowplot)
library(scDblFinder)
library(BiocSingular)
library(Libra)
library(clusterProfiler)
})
source("utils/Utils.R")

#plan("multiprocess", workers = 10)
options(future.globals.maxSize = 6000 * 1024^2)

dir.create("output_dge_Npas4pos")

load("output_reclust/Jessica_SCT_30pcs_05res_Reclust_NoDoublet_Relabel.RData")

Npas4 = GetAssayData(object=seuObject_nodoub_label,assay="RNA",slot="data")["Npas4",]

seuObject_nodoub_label@meta.data$Npas4 <- as.numeric(log2(Npas4+1))

seuObject_nodoub_label@meta.data <- seuObject_nodoub_label@meta.data %>%
                                    mutate(Npas4pos = if_else(Npas4 > 0, "Npas4+","Npas4-"))


# Number of cell with Npas4 > 1 expression
meta <- as.data.frame(seuObject_nodoub_label@meta.data)

meta_filt <- meta %>%
                filter(Npas4 > 1.0) 

x <- table(meta_filt$Cell,meta_filt$Genotype) %>%
        as.data.frame()

write.table(x, "output_dge_Npas4pos/Proportion_Cells_Npas4_More1.txt",sep="\t",quote=F)




drug_subset <- subset(seuObject_nodoub_label, subset = (Genotype == "drug"))
control_subset <- subset(seuObject_nodoub_label, subset = (Genotype == "control"))


# Degs for drug 
tmp <- drug_subset

tmp@meta.data$cell_type <- tmp@meta.data$Cell
tmp@meta.data$replicate <- tmp@meta.data$Npas4pos
tmp@meta.data$label <- tmp@meta.data$Npas4pos

# Run differential expression
DE_drug_Npas4pos <- run_de(tmp, de_family = 'singlecell', de_method = 'wilcox',n_threads = 20)

# Degs for saline
tmp <- control_subset

tmp@meta.data$cell_type <- tmp@meta.data$Cell
tmp@meta.data$replicate <- tmp@meta.data$Npas4pos
tmp@meta.data$label <- tmp@meta.data$Npas4pos

DE_control_Npas4pos <- run_de(tmp, de_family = 'singlecell', de_method = 'wilcox',n_threads = 20)

save(DE_drug_Npas4pos,DE_control_Npas4pos, file = "output_dge_Npas4pos/DGE_Npas4pos_Analysis.RData")


# Filter
DE_drug_Npas4pos_Sign <- DE_drug_Npas4pos %>% filter(p_val_adj < 0.05 & abs(avg_logFC) > 0.2)
DE_control_Npas4pos_Sign <- DE_control_Npas4pos %>% filter(p_val_adj < 0.05 & abs(avg_logFC) > 0.2)


openxlsx::write.xlsx(DE_drug_Npas4pos_Sign, 
                     file = "output_dge_Npas4pos/DE_drug_Npas4pos_Sign.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats",
                     overwrite=T)

openxlsx::write.xlsx(DE_control_Npas4pos_Sign, 
                     file = "output_dge_Npas4pos/DE_control_Npas4pos_Sign.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats",
                     overwrite=T)


write.table(DE_drug_Npas4pos_Sign,"output_dge_Npas4pos/DE_drug_Npas4pos_Sign.txt",sep="\t",quote=F)
write.table(DE_control_Npas4pos_Sign,"output_dge_Npas4pos/DE_control_Npas4pos_Sign.txt",sep="\t",quote=F)

## Gene Onto
dir.create("output_dge_Npas4pos/functional_enrichment_drug/")
DE_drug_Npas4pos_Sign <- DE_drug_Npas4pos %>% filter(p_val_adj < 0.05 & abs(avg_logFC) > 0.2)
l <- split(DE_drug_Npas4pos_Sign, DE_drug_Npas4pos_Sign$cell_type)
l <- l[sapply(l, nrow)>0] #remove objects with less than one gene

GOI <- list()
GeneOnto <- list()

for(i in 1:length(l)){
GOI[[i]] <- bitr(as.character(l[[i]]$gene),  
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Mm.eg.db::org.Mm.eg.db)

GeneOnto[[i]] <- enrichGO(gene = unique(GOI[[i]]$ENTREZID), 
                     keyType = "ENTREZID", 
                     OrgDb = org.Mm.eg.db::org.Mm.eg.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff  = 1, 
                     qvalueCutoff = 1, 
                     readable = TRUE)

openxlsx::write.xlsx(as.data.frame(GeneOnto[[i]]), 
                     file = sprintf("output_dge_Npas4pos/functional_enrichment_drug/%s.GO.xlsx", names(l)[[i]]), 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto", 
                     overwrite = TRUE)


PLOT <- dotplot(GeneOnto[[i]])
print(PLOT)
ggsave(sprintf("output_dge_Npas4pos/functional_enrichment_drug/%s.pdf", names(l)[[i]]))

}

dir.create("output_dge_Npas4pos/functional_enrichment_control/")
DE_control_Npas4pos_Sign <- DE_control_Npas4pos %>% filter(p_val_adj < 0.05 & abs(avg_logFC) > 0.2)
l <- split(DE_control_Npas4pos_Sign, DE_control_Npas4pos_Sign$cell_type)
l <- l[sapply(l, nrow)>0] #remove objects with less than one gene

GOI <- list()
GeneOnto <- list()

for(i in 1:length(l)){
GOI[[i]] <- bitr(as.character(l[[i]]$gene),  
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Mm.eg.db::org.Mm.eg.db)

GeneOnto[[i]] <- enrichGO(gene = unique(GOI[[i]]$ENTREZID), 
                     keyType = "ENTREZID", 
                     OrgDb = org.Mm.eg.db::org.Mm.eg.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff  = 1, 
                     qvalueCutoff = 1, 
                     readable = TRUE)

openxlsx::write.xlsx(as.data.frame(GeneOnto[[i]]), 
                     file = sprintf("output_dge_Npas4pos/functional_enrichment_control/%s.GO.xlsx", names(l)[[i]]), 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto", 
                     overwrite = TRUE)


PLOT <- dotplot(GeneOnto[[i]])
print(PLOT)
ggsave(sprintf("output_dge_Npas4pos/functional_enrichment_control/%s.pdf", names(l)[[i]]))

}


