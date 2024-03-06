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

dir.create("output_dge")

load("output_reclust/Jessica_SCT_30pcs_05res_Reclust_NoDoublet_Relabel.RData")

tmp <- seuObject_nodoub_label

tmp@meta.data$cell_type <- tmp@meta.data$Cell
tmp@meta.data$replicate <- tmp@meta.data$Genotype
tmp@meta.data$label <- tmp@meta.data$Genotype


# Run differential expression
DE_drug <- run_de(tmp, de_family = 'singlecell', de_method = 'wilcox',n_threads = 20)

#DE_drug <- DE_drug %>% mutate(avg_logFC = -1*avg_logFC) # switch sign

save(DE_drug, file = "output_dge/DGE_Mast_Analysis.RData")

# Filter
DE_drug_sign <- DE_drug %>% filter(p_val_adj < 0.05 & abs(avg_logFC) > 0.2)


openxlsx::write.xlsx(DE_drug_sign, 
                     file = "output_dge/DE_drug_sign.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats",
                     overwrite=T)
write.table(DE_drug_sign,"output_dge/DE_drug_sign.txt",sep="\t",quote=F)

#######################
# Convert in Human ID #
#######################
human = biomaRt::useMart(host="https://dec2021.archive.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
mouse = biomaRt::useMart(host="https://dec2021.archive.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")

#######################################
# Load Data Saline CPP vs Cocaine CPP #
#######################################


dge <- DE_drug_sign %>%
        dplyr::select(gene,cell_type) %>%
        dplyr::rename(Gene = gene) %>%
        arrange(cell_type)

MGI = biomaRt::getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = dge$Gene,mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

dge_human <- merge(dge,MGI,by.x="Gene",by.y="MGI.symbol",all=F)

dge_human <- dge_human %>%
                dplyr::select(HGNC.symbol,cell_type) %>%
                dplyr::rename(Gene = HGNC.symbol) %>%
                arrange(cell_type)


write.table(dge,"output_dge/DE_drug_sign_MouseID.txt",sep="\t",quote=F,row.names=F)
write.table(dge_human,"output_dge/DE_drug_sign_HumanID.txt",sep="\t",quote=F,row.names=F)


# Gene Onto
dir.create("output_dge/functional_enrichment/")
dge <- read.table("output_dge/DE_drug_sign_HumanID.txt",header=T)
l <- split(dge, dge$cell_type)
l <- l[sapply(l, nrow)>5] #remove objects with less than one gene

GOI <- list()
GeneOnto <- list()

for(i in 1:length(l)){
GOI[[i]] <- bitr(as.character(l[[i]]$Gene),  
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Hs.eg.db::org.Hs.eg.db)

GeneOnto[[i]] <- enrichGO(gene = unique(GOI[[i]]$ENTREZID), 
                     keyType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff  = 0.2, 
                     qvalueCutoff = 0.2, 
                     readable = TRUE)

openxlsx::write.xlsx(as.data.frame(GeneOnto[[i]]), 
                     file = sprintf("output_dge/functional_enrichment/%s.GO.xlsx", names(l)[[i]]), 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto", 
                     overwrite = TRUE)


PLOT <- dotplot(GeneOnto[[i]])
print(PLOT)
ggsave(sprintf("output_dge/functional_enrichment/%s.pdf", names(l)[[i]]))

}
