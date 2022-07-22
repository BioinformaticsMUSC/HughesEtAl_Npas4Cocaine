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
library(Libra)
})
source("utils/Utils.R")

#plan("multiprocess", workers = 10)
options(future.globals.maxSize = 6000 * 1024^2)

dir.create("output_dge")

load("output_reclust/Brandon_SCT_30pcs_04res_Reclust_NoDoublet_Slimmed_Relabel.RData")

tmp <- seuObject_slim_nodoub_slim

tmp@meta.data$cell_type <- tmp@meta.data$Cell
tmp@meta.data$replicate <- tmp@meta.data$Genotype
tmp@meta.data$label <- tmp@meta.data$Genotype

saline_subset <- subset(tmp, subset = (Genotype == "SalineCPP" | Genotype == "SalineNPAS4"))
cocaine_subset <- subset(tmp, subset = (Genotype == "CocaineCPP" | Genotype == "CocaineNPAS4"))
salcoc_subset <- subset(x = tmp, subset = (Genotype == "SalineCPP" | Genotype == "CocaineCPP"))
salcoc_npas4_subset <- subset(x = tmp, subset = (Genotype == "SalineNPAS4" | Genotype == "CocaineNPAS4"))



# Run differential expression
DE_saline <- run_de(saline_subset, de_family = 'singlecell', de_method = 'MAST',n_threads = 20)
DE_cocaine <- run_de(cocaine_subset, de_family = 'singlecell', de_method = 'MAST',n_threads = 20)
DE_salcoc <- run_de(salcoc_subset, de_family = 'singlecell', de_method = 'MAST',n_threads = 20)
DE_salcoc_npas4 <- run_de(salcoc_npas4_subset, de_family = 'singlecell', de_method = 'MAST',n_threads = 20)

DE_saline <- DE_saline %>% mutate(avg_logFC = -1*avg_logFC) # switch sign
DE_cocaine <- DE_cocaine %>% mutate(avg_logFC = -1*avg_logFC) # switch sign

save(DE_saline,DE_cocaine,DE_salcoc,DE_salcoc_npas4, file = "output_dge/DGE_MultiComp/DGE_Mast_Analysis.RData")

# Filter
DE_saline_sign <- DE_saline %>% filter(p_val_adj < 0.05 & abs(avg_logFC) > 0.2)
DE_cocaine_sign <- DE_cocaine %>% filter(p_val_adj < 0.05 & abs(avg_logFC) > 0.2)
DE_salcoc_sign <- DE_salcoc %>% filter(p_val_adj < 0.05 & abs(avg_logFC) > 0.2)
DE_salcoc_npas4_sign <- DE_salcoc_npas4 %>% filter(p_val_adj < 0.05 & abs(avg_logFC) > 0.2)


openxlsx::write.xlsx(DE_saline_sign, 
                     file = "output_dge/DGE_MultiComp/SalineCPPvsSalineNPAS4.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats",
                     overwrite=T)
write.table(DE_saline_sign,"output_dge/DGE_MultiComp/SalineCPPvsSalineNPAS4.txt",sep="\t",quote=F)

openxlsx::write.xlsx(DE_cocaine_sign, 
                     file = "output_dge/DGE_MultiComp/CocaineCPPvsCocaineNPAS4.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")
write.table(DE_cocaine_sign,"output_dge/DGE_MultiComp/CocaineCPPvsCocaineNPAS4.txt",sep="\t",quote=F)


openxlsx::write.xlsx(DE_salcoc_sign, 
                     file = "output_dge/DGE_MultiComp/SalineCPPvsCocaineCPP.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")
write.table(DE_salcoc_sign,"output_dge/DGE_MultiComp/SalineCPPvsCocaineCPP.txt",sep="\t",quote=F)

openxlsx::write.xlsx(DE_salcoc_npas4_sign, 
                     file = "output_dge/DGE_MultiComp/SalineNPAS4vsCocaineNPAS4.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")
write.table(DE_salcoc_npas4_sign,"output_dge/DGE_MultiComp/SalineNPAS4vsCocaineNPAS4.txt",sep="\t",quote=F)

# Convert in Human ID

human = biomauseMart(host="https://dec2021.archive.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
mouse = useMart(host="https://dec2021.archive.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")

# Load Data Saline vs Cocaine
dge <- read.table("output_dge/DGE_MultiComp/SalineCPPvsCocaineCPP.txt",header=T)

dge <- dge %>%
        dplyr::select(gene,cell_type) %>%
        dplyr::rename(Gene = gene) %>%
        arrange(cell_type)

MGI = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = dge$Gene,mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

dge_human <- merge(dge,MGI,by.x="Gene",by.y="MGI.symbol",all=F)

dge_human <- dge_human %>%
                dplyr::select(HGNC.symbol,cell_type) %>%
                dplyr::rename(Gene = HGNC.symbol) %>%
                arrange(cell_type)


write.table(dge,"output_dge/DGE_MultiComp/SalineCPPvsCocaineCPP_MouseID.txt",sep="\t",quote=F,row.names=F)
write.table(dge_human,"output_dge/DGE_MultiComp/SalineCPPvsCocaineCPP_HumanID.txt",sep="\t",quote=F,row.names=F)



# Load Data SalineNPAS vs CocaineNPAS4
dge <- read.table("output_dge/DGE_MultiComp/SalineNPAS4vsCocaineNPAS4.txt",header=T)

dge <- dge %>%
        dplyr::select(gene,cell_type) %>%
        dplyr::rename(Gene = gene) %>%
        arrange(cell_type)


# Convert to human 
#human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

MGI = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = dge$Gene,mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

dge_human <- merge(dge,MGI,by.x="Gene",by.y="MGI.symbol",all=F)

dge_human <- dge_human %>%
                dplyr::select(HGNC.symbol,cell_type) %>%
                dplyr::rename(Gene = HGNC.symbol) %>%
                arrange(cell_type)


write.table(dge,"output_dge/DGE_MultiComp/SalineNPAS4vsCocaineNPAS4_MouseID.txt",sep="\t",quote=F,row.names=F)
write.table(dge_human,"output_dge/DGE_MultiComp/SalineNPAS4vsCocaineNPAS4_HumanID.txt",sep="\t",quote=F,row.names=F)


# Load Data SalineCCP vs SalineNPAS4
dge <- read.table("output_dge/DGE_MultiComp/SalineCPPvsSalineNPAS4.txt",header=T)

dge <- dge %>%
        dplyr::select(gene,cell_type) %>%
        dplyr::rename(Gene = gene) %>%
        arrange(cell_type)


# Convert to human 
#human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

MGI = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = dge$Gene,mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

dge_human <- merge(dge,MGI,by.x="Gene",by.y="MGI.symbol",all=F)

dge_human <- dge_human %>%
                dplyr::select(HGNC.symbol,cell_type) %>%
                dplyr::rename(Gene = HGNC.symbol) %>%
                arrange(cell_type)


write.table(dge,"output_dge/DGE_MultiComp/SalineCPPvsSalineNPAS4_MouseID.txt",sep="\t",quote=F,row.names=F)
write.table(dge_human,"output_dge/DGE_MultiComp/SalineCPPvsSalineNPAS4_HumanID.txt",sep="\t",quote=F,row.names=F)

# Load Data CocaineCCP vs CocaineNPAS4
dge <- read.table("output_dge/DGE_MultiComp/CocaineCPPvsCocaineNPAS4.txt",header=T)

dge <- dge %>%
        dplyr::select(gene,cell_type) %>%
        dplyr::rename(Gene = gene) %>%
        arrange(cell_type)


# Convert to human 
#human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

MGI = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = dge$Gene,mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

dge_human <- merge(dge,MGI,by.x="Gene",by.y="MGI.symbol",all=F)

dge_human <- dge_human %>%
                dplyr::select(HGNC.symbol,cell_type) %>%
                dplyr::rename(Gene = HGNC.symbol) %>%
                arrange(cell_type)


write.table(dge,"output_dge/DGE_MultiComp/CocaineCPPvsCocaineNPAS4_MouseID.txt",sep="\t",quote=F,row.names=F)
write.table(dge_human,"output_dge/DGE_MultiComp/CocaineCPPvsCocaineNPAS4_HumanID.txt",sep="\t",quote=F,row.names=F)


##############################
# DGE Cocaine NPAS4 Specific #
##############################

load("output_reclust/Brandon_SCT_30pcs_04res_Reclust_NoDoublet_Slimmed_Relabel.RData")

# subset the data
seuObject_slim_nodoub_slim@meta.data <- seuObject_slim_nodoub_slim@meta.data %>% 
                                         select(orig.ident,Genotype,pMito,nCount_RNA,nFeature_RNA,Cell,seurat_clusters)

# subset the for each comparison
a <- subset(seuObject_slim_nodoub_slim, subset = (Genotype == "CocaineNPAS4" | Genotype == "SalineCPP")) %>% as.SingleCellExperiment()

colData(a) <- colData(a) %>% 
                as.data.frame() %>%
                mutate(label = case_when(Genotype == "CocaineNPAS4" ~ "CocaineNPAS4", Genotype != "CocaineNPAS4" ~ "Others")) %>%
                dplyr::rename(cell_type= Cell, replicate = Genotype) %>%
                DataFrame()


b <- subset(seuObject_slim_nodoub_slim, subset = (Genotype == "CocaineNPAS4" | Genotype == "SalineNPAS4")) %>% as.SingleCellExperiment()

colData(b) <- colData(b) %>% 
                as.data.frame() %>%
                mutate(label = case_when(Genotype == "CocaineNPAS4" ~ "CocaineNPAS4", Genotype != "CocaineNPAS4" ~ "Others")) %>%
                dplyr::rename(cell_type= Cell, replicate = Genotype) %>%
                DataFrame()

c <- subset(x = seuObject_slim_nodoub_slim, subset = (Genotype == "CocaineNPAS4" | Genotype == "CocaineCPP")) %>% as.SingleCellExperiment()

colData(c) <- colData(c) %>% 
                as.data.frame() %>%
                mutate(label = case_when(Genotype == "CocaineNPAS4" ~ "CocaineNPAS4", Genotype != "CocaineNPAS4" ~ "Others")) %>%
                dplyr::rename(cell_type= Cell, replicate = Genotype) %>%
                DataFrame()


d <- subset(seuObject_slim_nodoub_slim, subset = (Genotype == "SalineCPP" | Genotype == "CocaineCPP")) %>% as.SingleCellExperiment()

colData(d) <- colData(d) %>% 
                as.data.frame() %>%
                mutate(label = case_when(Genotype == "SalineCPP" ~ "SalineCPP", Genotype != "SalineCPP" ~ "Others")) %>%
                dplyr::rename(cell_type= Cell, replicate = Genotype) %>%
                DataFrame()


e <- subset(seuObject_slim_nodoub_slim, subset = (Genotype == "SalineCPP" | Genotype == "SalineNPAS4")) %>% as.SingleCellExperiment()

colData(e) <- colData(e) %>% 
                as.data.frame() %>%
                mutate(label = case_when(Genotype == "SalineCPP" ~ "SalineCPP", Genotype != "SalineCPP" ~ "Others")) %>%
                dplyr::rename(cell_type= Cell, replicate = Genotype) %>%
                DataFrame()


f <- subset(seuObject_slim_nodoub_slim, subset = (Genotype == "CocaineCPP" | Genotype == "SalineNPAS4")) %>% as.SingleCellExperiment()

colData(f) <- colData(f) %>% 
                as.data.frame() %>%
                mutate(label = case_when(Genotype == "CocaineCPP" ~ "CocaineCPP", Genotype != "CocaineCPP" ~ "Others")) %>%
                dplyr::rename(cell_type= Cell, replicate = Genotype) %>%
                DataFrame()

# Run differetial expression analysis
DE_a <- run_de(a, de_family = 'singlecell', de_method = 'MAST',n_threads = 20)
DE_b <- run_de(b, de_family = 'singlecell', de_method = 'MAST',n_threads = 20)
DE_c <- run_de(c, de_family = 'singlecell', de_method = 'MAST',n_threads = 20)
DE_d <- run_de(d, de_family = 'singlecell', de_method = 'MAST',n_threads = 20)
DE_e <- run_de(e, de_family = 'singlecell', de_method = 'MAST',n_threads = 20)
DE_f <- run_de(f, de_family = 'singlecell', de_method = 'MAST',n_threads = 20)


DE_a <- DE_a %>% 
        select(cell_type, gene, avg_logFC, p_val_adj) %>%
        dplyr::rename(logFC_CocaineNPAS4vsSalineCPP= avg_logFC, FDR_CocaineNPAS4vsSalineCPP = p_val_adj)

DE_b <- DE_b %>% 
        select(cell_type, gene, avg_logFC, p_val_adj) %>%
        dplyr::rename(logFC_CocaineNPAS4vsSalineNPAS4= avg_logFC, FDR_CocaineNPAS4vsSalineNPAS4 = p_val_adj)


DE_c <- DE_c %>% 
        select(cell_type, gene, avg_logFC, p_val_adj) %>%
        dplyr::rename(logFC_CocaineNPAS4vsCocaineCPP= avg_logFC, FDR_CocaineNPAS4vsCocaineCPP = p_val_adj)

DE_d <- DE_d %>% 
        select(cell_type, gene, avg_logFC, p_val_adj) %>%
        dplyr::rename(logFC_SalineCPPvsCocaineCPP= avg_logFC, FDR_SalineCPPvsCocaineCPP = p_val_adj)

DE_e <- DE_e %>% 
        select(cell_type, gene, avg_logFC, p_val_adj) %>%
        dplyr::rename(logFC_SalineCPPvsSalineNPAS4= avg_logFC, FDR_SalineCPPvsSalineNPAS4 = p_val_adj)

DE_f <- DE_f %>% 
        select(cell_type, gene, avg_logFC, p_val_adj) %>%
        dplyr::rename(logFC_CocaineCPPvsSalineNPAS4 = avg_logFC, FDR_CocaineCPPvsSalineNPAS4 = p_val_adj)



# Create a database with all statistics
CocaineNPAS4_DGE_Stats <- list(DE_a,DE_b,DE_c,DE_d,DE_e,DE_f) %>% 
        purrr::reduce(full_join,by=c("cell_type","gene")) %>%
        as.data.frame()


# Filter for CocaineNPAS4 specifically upregulated
CocaineNPAS4_Upreg <- CocaineNPAS4_DGE_Stats %>%
        filter(logFC_CocaineNPAS4vsCocaineCPP > 0.2, logFC_CocaineNPAS4vsSalineNPAS4 > 0, logFC_CocaineNPAS4vsSalineCPP > 0, 
                FDR_CocaineNPAS4vsCocaineCPP < 0.05) %>% 
                #FDR_CocaineNPAS4vsSalineNPAS4 > 0.2, FDR_CocaineNPAS4vsSalineCPP > 0.2, FDR_CocaineCPPvsSalineNPAS4 > 0.2, FDR_SalineCPPvsSalineNPAS4 > 0.2) %>%
        as.data.frame()  %>%
        arrange(cell_type,gene,logFC_CocaineNPAS4vsCocaineCPP)

# Filter for CocaineNPAS4 specifically downregulated
CocaineNPAS4_Downreg <- CocaineNPAS4_DGE_Stats %>%
        filter(logFC_CocaineNPAS4vsCocaineCPP < -0.2, logFC_CocaineNPAS4vsSalineNPAS4 < 0, logFC_CocaineNPAS4vsSalineCPP < 0, 
                FDR_CocaineNPAS4vsCocaineCPP < 0.05) %>% 
                #FDR_CocaineNPAS4vsSalineNPAS4 > 0.2, FDR_CocaineNPAS4vsSalineCPP > 0.2, FDR_CocaineCPPvsSalineNPAS4 > 0.2, FDR_SalineCPPvsSalineNPAS4 > 0.2) %>%
        as.data.frame()  %>%
        arrange(cell_type,gene,logFC_CocaineNPAS4vsCocaineCPP)

# Concatenate and save
CocaineNPAS4_Specific <- rbind(CocaineNPAS4_Upreg,CocaineNPAS4_Downreg)

# Now clean and subset the data
CocaineNPAS4_Specific_Cleaned <- CocaineNPAS4_Specific %>% 
                                   mutate(Direction = case_when(logFC_CocaineNPAS4vsCocaineCPP > 0 ~ "UpReg", logFC_CocaineNPAS4vsCocaineCPP < 0 ~ "DownReg")) %>%
                                   select(gene, cell_type, Direction)

# Convert in Human ID

human = biomaRt::useMart(host="https://dec2021.archive.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
mouse = biomaRt::useMart(host="https://dec2021.archive.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")

MGI = biomaRt::getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = CocaineNPAS4_Specific_Cleaned$gene,mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

CocaineNPAS4_Specific_Cleaned_HumanID <- merge(CocaineNPAS4_Specific_Cleaned,MGI,by.x="gene",by.y="MGI.symbol",all=F)

CocaineNPAS4_Specific_Cleaned_HumanID <- CocaineNPAS4_Specific_Cleaned_HumanID %>%
                dplyr::select(HGNC.symbol,cell_type,Direction) %>%
                dplyr::rename(Gene = HGNC.symbol) %>%
                arrange(cell_type,Direction)


write.table(CocaineNPAS4_Specific_Cleaned,"output_dge/DGE_CocaineNPAS4_Specific/CocaineNPAS4_Specific_MouseID.txt",sep="\t",quote=F,row.names=F)
write.table(CocaineNPAS4_Specific_Cleaned_HumanID,"output_dge/DGE_CocaineNPAS4_Specific/CocaineNPAS4_Specific_HumanID.txt",sep="\t",quote=F,row.names=F)

openxlsx::write.xlsx(CocaineNPAS4_Specific_Cleaned, 
                     file = "output_dge/DGE_CocaineNPAS4_Specific/CocaineNPAS4_Specific.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Dge")

save(CocaineNPAS4_DGE_Stats,CocaineNPAS4_Specific,CocaineNPAS4_Specific_Cleaned,CocaineNPAS4_Specific_Cleaned_HumanID, file = "output_dge/DGE_CocaineNPAS4_Specific/CocaineNPAS4_DEG.RData")


# Gene Onto
dir.create("output_dge/DGE_CocaineNPAS4_Specific/functional_enrichment/")
l <- split(CocaineNPAS4_Specific_Cleaned, CocaineNPAS4_Specific_Cleaned$cell_type)
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
                     pvalueCutoff  = 0.05, 
                     qvalueCutoff = 0.2, 
                     readable = TRUE)

openxlsx::write.xlsx(as.data.frame(GeneOnto[[i]]), 
                     file = sprintf("output_dge/DGE_CocaineNPAS4_Specific/functional_enrichment/%s.GO.xlsx", names(l)[[i]]), 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto", 
                     overwrite = TRUE)


PLOT <- dotplot(GeneOnto[[i]])
print(PLOT)
ggsave(sprintf("output_dge/DGE_CocaineNPAS4_Specific/functional_enrichment/%s.pdf", names(l)[[i]]))

}


# Visualize the gene expression changes of genes. Shisha9 
colore <- c("#e69b00","#17154f","#355828","#b0799a")


Shisa9<- GetAssayData(object=seuObject_slim_nodoub_slim,assay="RNA",slot="data")["Shisa9",] # For Drd1
Cartpt<- GetAssayData(object=seuObject_slim_nodoub_slim,assay="RNA",slot="data")["Cartpt",] # For Drd1

# Add to meta
meta <- seuObject_slim_nodoub_slim@meta.data %>%
            as.data.frame() %>%
            mutate(Shisa9 = log2(Shisa9+1), Cartpt = log2(Cartpt+1))

#my_comparisons <- list( c("SalineCPP", "CocaineCPP"), c("SalineNPAS4", "CocaineNPAS4"), c("CocaineCPP", "CocaineNPAS4"),c("SalineCPP", "SalineNPAS4") )

pdf("output_dge/DGE_CocaineNPAS4_Specific/Shisa9_Expression_Drd2.pdf", width = 2, height = 4)
meta %>%
filter(Cell == "MSN_Drd2+_1") %>%
#filter(Shisa9 > 0) %>%
ggerrorplot(x = "Genotype", y = "Shisa9",
            color = "Genotype", palette = colore, add = "mean", error.plot = "errorbar",
            outlier.shape = NA)+ 
  #stat_compare_means(comparisons = my_comparisons,label.y = c(2.1, 2.3, 2.5, 2.7)) +
        theme(legend.position="none") +
        rotate_x_text(angle = 45) +
        xlab("") + 
        ylab("Expression")
        #ylim(1.6,2.2)
dev.off()

pdf("output_dge/DGE_CocaineNPAS4_Specific/Cartpt_Expression_Drd2.pdf", width = 2, height = 4)
meta %>%
filter(Cell == "MSN_Drd2+_1") %>%
#filter(Shisa9 > 0) %>%
ggerrorplot(x = "Genotype", y = "Cartpt",
            color = "Genotype", palette = colore, add = "mean", error.plot = "errorbar",
            outlier.shape = NA)+ 
  #stat_compare_means(comparisons = my_comparisons,label.y = c(2.1, 2.3, 2.5, 2.7)) +
        theme(legend.position="none") +
        rotate_x_text(angle = 45) +
        xlab("") + 
        ylab("Expression")
        #ylim(1.6,2.2)
dev.off()

