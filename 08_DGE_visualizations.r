suppressPackageStartupMessages({
library(tidyverse)
library(ggrepel)
library(emmeans)
library(SingleCellExperiment)
library(scater)
library(BiocParallel)
library(ggpubr)
library(speckle)
library(muscat)
library(Seurat)
library(clustree)
library(data.table)
library(cowplot)
library(scDblFinder)
library(BiocSingular)
library(RColorBrewer)
library(MetBrewer)
library(scales)
library(scbp)
library(scPred)
library(ComplexHeatmap)
library(schex)
library(ggradar)
library(scGate)
library(UCell)
library(UpSetR)
})

source("utils/Utils.R")

dir.create("output_figures/output_Figure2")

#dge <- read.table("output_sct_reclust/output_mast_manual/SalineCPPvsSalineNPAS4/Final_DEGS_CocaineCPPvsCocaineNPAS4.txt",header=T,sep="\t")
#full <- read.table("output_sct_reclust/output_mast_manual/SalineCPPvsSalineNPAS4/Final_FullStats_CocaineCPPvsCocaineNPAS4.txt",header=T,sep="\t")

load("output_reclust/Brandon_SCT_30pcs_04res_Reclust_NoDoublet_Slimmed_Relabel.RData")
load("output_dge/DGE_MultiComp/DGE_Mast_Analysis.RData")

# Order the cells
order <- c("MSN_Drd1+_1","MSN_Drd1+_2","MSN_Drd1+_3","MSN_Drd2+_1","MSN_Drd2+_2","MSN_Drd3+_1","MSN_Grm8+_1","MSN_Grm8+_2",
  "Interneurons_Chat+","Interneurons_Th+","Interneurons_Sst+","Interneurons_Pvalb+","Interneurons_Pnoc+","Interneurons_Ndnf+",
  "Astrocytes","Oligodendrocytes","OPC","Microglia","Endothelial","Neuroblast")


# Get scale.data
#seuObject_slim_nodoub_slim <- seuObject_slim_nodoub_slim %>%
#                                ScaleData(vars.to.regress = c("nCount_SCT","pMito"), model.use = "linear")


## Visualize number of Npas4 Expressing cells 
Npas4 = GetAssayData(object=seuObject_slim_nodoub_slim,assay="RNA",slot="data")["Npas4",]

seuObject_slim_nodoub_slim@meta.data$Npas4 <- as.numeric(log2(Npas4+1))

umap <- as.data.frame(Embeddings(seuObject_slim_nodoub_slim, reduction = "umap"))
meta <- as.data.frame(seuObject_slim_nodoub_slim@meta.data)

# Npas4 expression violin plot
pdf("output_figures/output_Figure1/Npas4_ViolinPlot.pdf", width = 6, height = 6)
meta %>%
filter(Npas4 > 0) %>%
mutate(Cell=fct_relevel(Cell,order)) %>%
ggviolin(x = "Genotype", y = "Npas4",
            color = "Genotype", 
            palette=c("#17154f","#b0799a","#e69b00","#355828"),
            outlier.shape = NA,
            add = "mean_sd")+ 
        stat_compare_means(comparisons = list(c("SalineCPP", "SalineNPAS4")),label.y = c(3),method = "wilcox.test") +
        theme(legend.position="none") +
        rotate_x_text(angle = 45) +
        xlab("") + 
        ylab("Expression") +
        facet_wrap(.~Cell)
dev.off()

# Add EGFP and filter for cells co-expression both
#EXPR2 = GetAssayData(object=seuObject_slim_nodoub_slim,assay="RNA",slot="data")["Egfp",]
#seuObject_slim_nodoub_slim@meta.data$Egfp <- as.numeric(EXPR2)
#umap <- as.data.frame(Embeddings(seuObject_slim_nodoub_slim, reduction = "umap"))
#meta <- as.data.frame(seuObject_slim_nodoub_slim@meta.data)

# Npas4 expression violin plot in Egfp
#pdf("output_figures/output_Figure1/Npas4_ViolinPlot_OnlyInEGFP.pdf", width = 6, height = 6)
#meta %>%
#filter(Egfp > 0) %>% 
#filter(Npas4 > 0) %>% 
#mutate(Cell=fct_relevel(Cell,order)) %>%
#ggviolin(x = "Genotype", y = "Npas4",
#            color = "Genotype", 
#            palette=c("#17154f","#b0799a","#e69b00","#355828"),
#            outlier.shape = NA,
#            add = "mean_sd")+ 
#        stat_compare_means(comparisons = list(c("SalineCPP", "SalineNPAS4")),label.y = c(3),method = "wilcox.test") +
#        theme(legend.position="none") +
#        rotate_x_text(angle = 45) +
#        xlab("") + 
#        ylab("Expression") +
#        facet_wrap(.~Cell)
#dev.off()

# Egfp expression violin plot
#pdf("output_Figure1/Egfp_ViolinPlot.pdf", width = 6, height = 6)
#meta %>%
#filter(Egfp > 0) %>% 
#mutate(Cell=fct_relevel(Cell,order)) %>%
#ggviolin(x = "Genotype", y = "Egfp",
#            color = "Genotype", 
#            palette=c("#17154f","#b0799a","#e69b00","#355828"),
#            outlier.shape = NA,
#            add = "mean_sd")+ 
#        #stat_compare_means(comparisons = list(c("SalineCPP", "SalineNPAS4")),label.y = c(3),method = "wilcox.test") +
#        theme(legend.position="none") +
#        rotate_x_text(angle = 45) +
#        xlab("") + 
#        ylab("Expression") +
#        facet_wrap(.~Cell)
#dev.off()


# Cell number total
colore <- c("#17154f","#b0799a","#e69b00","#355828")
grouped_data <- table(meta$Cell, meta$Genotype) %>%
                as.matrix() %>%
                as.data.frame() %>%
                arrange(Var1,Var2) %>%
                dplyr::rename(Cell = Var1, Genotype = Var2, Cells = Freq)  %>%
                mutate(Cell=fct_rev(fct_relevel(Cell,order))) %>%
                as.data.frame() %>%
                droplevels()

pdf("output_figures/output_Figure1/Barplot_NumberOfCells.pdf", width = 4, height = 5)
ggbarplot(grouped_data, "Cell", "Cells",
  fill = "Genotype", color = "Genotype", palette = colore,
  label = FALSE, lab.col = "white", lab.pos = "in",
  rotate = TRUE) + 
#scale_y_continuous(labels = dollar_format(suffix = "%", prefix = "")) +
rotate_x_text(angle = 45) +
theme(legend.position="top") +
xlab("")
dev.off()

# NPAS4 only
df <- cbind(umap,meta)%>% 
  group_by(Cell) %>% 
  mutate(Total = n()) %>%
    ungroup() %>% 
    arrange(Genotype) %>%
    filter(Npas4 > 0) %>%
    group_by(Genotype,Cell) %>% 
    mutate(N_npas4 = n()) %>%
    as.data.frame() %>%
    mutate(Perc = (N_npas4*100)/Total)

grouped_data <- table(df$Cell, df$Genotype) %>%
                as.matrix() %>%
                as.data.frame() %>%
                arrange(Var1,Var2) %>%
                dplyr::rename(Cell = Var1, Genotype = Var2, Npas4 = Freq)

total_data <- table(meta$Cell, meta$Genotype) %>%
                as.matrix() %>%
                as.data.frame() %>%
                arrange(Var1,Var2) %>%
                dplyr::rename(Cell = Var1, Genotype = Var2, Total = Freq)

grouped_data$Total <- total_data$Total

grouped_data <- grouped_data %>%
                mutate(Perc = round((Npas4*100)/Total, 1)) %>%
                mutate(Cell=fct_rev(fct_relevel(Cell,order))) %>%
                as.data.frame() 

grouped_data[is.na(grouped_data)] <- 0



pdf("output_figures/output_Figure2/Barplot_Npas4_PercentageByCell.pdf", width = 4, height = 5)
ggbarplot(grouped_data, "Cell", "Perc",
  fill = "Genotype", color = "Genotype", palette = colore,
  label = FALSE, lab.col = "white", lab.pos = "in",
  rotate = TRUE) + 
scale_y_continuous(labels = dollar_format(suffix = "%", prefix = "")) +
rotate_x_text(angle = 45) +
theme(legend.position="top") +
xlab("")
dev.off()

pdf("output_figures/output_Figure2/Barplot_Npas4_TotalByCell.pdf", width = 4, height = 5)
ggbarplot(grouped_data, "Cell", "Npas4",
  fill = "Genotype", color = "Genotype", palette = colore,
  label = FALSE, lab.col = "white", lab.pos = "in",
  rotate = TRUE) + 
#scale_y_continuous(labels = dollar_format(suffix = "%", prefix = "")) +
rotate_x_text(angle = 45) +
theme(legend.position="top") +
xlab("")
dev.off()

# Fold Change Correlations
tmp1 <- DE_salcoc %>% 
        group_by(cell_type) %>%
        mutate(LOG = -log10(p_val_adj), ABS = abs(avg_logFC)) %>% 
        mutate(Threshold = if_else(p_val_adj < 0.05 & ABS > 0.2, "TRUE","FALSE")) %>%
        mutate(Direction = case_when(avg_logFC > 0.2 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.2 & p_val_adj < 0.05 ~ "DownReg")) %>%
        filter(Threshold == "TRUE") %>% 
        dplyr::select(cell_type,gene,avg_logFC) %>%
        dplyr::rename(logFC_SalVsCoc = avg_logFC)

tmp2 <- DE_salcoc_npas4 %>% 
        group_by(cell_type) %>%
        mutate(LOG = -log10(p_val_adj), ABS = abs(avg_logFC)) %>% 
        mutate(Threshold = if_else(p_val_adj < 0.05 & ABS > 0.2, "TRUE","FALSE")) %>%
        mutate(Direction = case_when(avg_logFC > 0.2 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.2 & p_val_adj < 0.05 ~ "DownReg")) %>%
        filter(Threshold == "TRUE") %>% 
        dplyr::select(cell_type,gene,avg_logFC) %>%
        dplyr::rename(logFC_Npas4_SalVsCoc = avg_logFC)

tmp3 <- DE_saline %>% 
        group_by(cell_type) %>%
        mutate(LOG = -log10(p_val_adj), ABS = abs(avg_logFC)) %>% 
        mutate(Threshold = if_else(p_val_adj < 0.05 & ABS > 0.2, "TRUE","FALSE")) %>%
        mutate(Direction = case_when(avg_logFC > 0.2 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.2 & p_val_adj < 0.05 ~ "DownReg")) %>%
        filter(Threshold == "TRUE") %>% 
        dplyr::select(cell_type,gene,avg_logFC) %>%
        dplyr::rename(logFC_SalVsSalNpas4 = avg_logFC)


tmp4 <- DE_cocaine %>% 
        group_by(cell_type) %>%
        mutate(LOG = -log10(p_val_adj), ABS = abs(avg_logFC)) %>% 
        mutate(Threshold = if_else(p_val_adj < 0.05 & ABS > 0.2, "TRUE","FALSE")) %>%
        mutate(Direction = case_when(avg_logFC > 0.2 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.2 & p_val_adj < 0.05 ~ "DownReg")) %>%
        filter(Threshold == "TRUE") %>% 
        dplyr::select(cell_type,gene,avg_logFC) %>%
        dplyr::rename(logFC_CocVsCocNpas4 = avg_logFC)


# Radar plot
df1 <- table(tmp1$cell_type) %>%
    as.data.frame() %>%
    mutate(Class = "DE_salcoc")

df2 <- table(tmp2$cell_type) %>%
    as.data.frame() %>%
    mutate(Class = "DE_salcoc_npas4")

df3 <- table(tmp3$cell_type) %>%
    as.data.frame() %>%
    mutate(Class = "DE_saline")

df4 <- table(tmp4$cell_type) %>%
    as.data.frame() %>%
    mutate(Class = "DE_cocaine")

df <- rbind(df1,df2,df3,df4) %>% 
  pivot_wider(names_from = Var1, values_from = Freq) %>%
  as.data.frame()

labels <- c("Class",order)

df <- df[,match(labels,colnames(df))]
df[is.na(df)] <- 0

pdf("output_figures/output_Figure2/RadarPlot_Degs.pdf", width = 8, height = 8)
ggradar(
  df, 
  values.radar = c("0", "500", "1000"),
  grid.min = 0, grid.mid = 500, grid.max = 1000,
  # Polygons
  group.line.width = 1, 
  group.point.size = 3,
  group.colours = c("#04a3bd","#f0be3d","#da7901","#247d3f"),
  # Background and grid lines
  background.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = "bottom"
  )
dev.off()

# Radar only Saline and Cocaine vs CvsCnpas4
df1 <- table(tmp1$cell_type) %>%
    as.data.frame() %>%
    mutate(Class = "DE_salcoc")

df4 <- table(tmp4$cell_type) %>%
    as.data.frame() %>%
    mutate(Class = "DE_cocaine")

df <- rbind(df1,df4) %>% 
  pivot_wider(names_from = Var1, values_from = Freq) %>%
  as.data.frame()

labels <- c("Class",order)

df <- df[,match(labels,colnames(df))]
df[is.na(df)] <- 0

pdf("output_figures/output_Figure2/RadarPlot_Degs_SvsC_CvsCnpas4.pdf", width = 8, height = 8)
ggradar(
  df, 
  values.radar = c("0", "500", "1000"),
  grid.min = 0, grid.mid = 500, grid.max = 1000,
  # Polygons
  group.line.width = 1, 
  group.point.size = 3,
  group.colours = c("#eacc62","#924099"),
  # Background and grid lines
  background.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = "bottom"
  )
dev.off()


# MSN only 
pdf("output_figures/output_Figure2/RadarPlot_Degs_SvsC_CvsCnpas4_MSN.pdf", width = 5, height = 5)
ggradar(
  df[c(1:9)], 
  values.radar = c("0", "250", "500"),
  grid.min = 0, grid.mid = 250, grid.max = 500,
  # Polygons
  group.line.width = 1, 
  group.point.size = 3,
  group.colours = c("#eacc62","#924099"),
  # Background and grid lines
  background.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = "bottom"
  )
dev.off()


###########################################
# Vulcano plots Saline CPP Vs Cocaine CPP #
###########################################

df1 <- DE_salcoc %>% 
    group_by(cell_type) %>%
        mutate(LOG = -log10(p_val_adj), ABS = abs(avg_logFC)) %>% 
        mutate(Threshold = if_else(p_val_adj < 0.05 & ABS > 0.2, "TRUE","FALSE")) %>%
        mutate(Direction = case_when(avg_logFC > 0.2 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.2 & p_val_adj < 0.05 ~ "DownReg")) 

# cell colors original
col <- data.table::fread("output_figures/output_Figure1/Colors_Used_Initial_Data.txt") %>% as.data.frame()
col <- col[match(order,col$Cells),]

colors <- as.character(col$Colors)


pdf("output_figures/output_Figure2/Vulcano_SalCPPVsCocCPP_AllCells.pdf", width = 4, height = 4)
df1 %>%
ggplot(aes(x=avg_logFC, y=LOG)) +
ggrastr::geom_point_rast(aes(colour = cell_type),size=0.5) +
      xlab("log2(Fold Change)")+ 
      ylab("-log10(FDR)")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_vline(xintercept = 0.2, colour = "grey",linetype="dashed",size=1,alpha=0.5) +
      geom_vline(xintercept = -0.2, colour = "grey",linetype="dashed",size=1,alpha=0.5) + 
      geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      theme(legend.position="none")+
       scale_colour_manual(values = colors) +   # customized color palette
      ylim(0,130) + xlim(-3,3)
dev.off()

pdf("output_figures/output_Figure2/Vulcano_SalCPPVsCocCPP_Drd2_1.pdf", width = 4, height = 4)
top_labelled <- df1 %>% 
                  filter(cell_type == "MSN_Drd2+_1", gene %in% c("Calm2", "Slc8a1", "Camk2a", "Penk", "Scg2", "Ptn","Ctnna3")) 

df1 %>%
  filter(cell_type == "MSN_Drd2+_1") %>%
  ggscatter( 
            x = "avg_logFC", 
            y = "LOG",
            color = "Threshold",
            palette=c("grey","red"),
            size = 1,
            alpha=0.3,
            shape=19)+
      xlab("log2(Fold Change)")+ 
      ylab("-log10(FDR)")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_vline(xintercept = 0.2, colour = "grey",linetype="dashed",size=1,alpha=0.5) +
      geom_vline(xintercept = -0.2, colour = "grey",linetype="dashed",size=1,alpha=0.5) + 
      geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_text_repel(data = top_labelled, 
                        mapping = aes(label = gene), 
                        size = 5,
                        nudge_x = .15,
                        box.padding = 0.5,
                        nudge_y = 1,
                        segment.curvature = -0.1,
                        segment.ncp = 3,
                        segment.angle = 20)+
      theme(legend.position="none")+
      ylim(0,130) + xlim(-1,1)
dev.off()


pdf("output_figures/output_Figure2/Vulcano_SalCPPVsCocCPP_Drd1_1.pdf", width = 4, height = 4)
top_labelled <- df1 %>% 
                  filter(cell_type == "MSN_Drd1+_1", gene %in% c("Pcdh15", "Tac1", "Reln", "Htr4", "Pcbp2", "Atp2b1","Ntm")) 

df1 %>%
  filter(cell_type == "MSN_Drd1+_1") %>%
  ggscatter( 
            x = "avg_logFC", 
            y = "LOG",
            color = "Threshold",
            palette=c("grey","red"),
            size = 1,
            alpha=0.3,
            shape=19)+
      xlab("log2(Fold Change)")+ 
      ylab("-log10(FDR)")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_vline(xintercept = 0.2, colour = "grey",linetype="dashed",size=1,alpha=0.5) +
      geom_vline(xintercept = -0.2, colour = "grey",linetype="dashed",size=1,alpha=0.5) + 
      geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_text_repel(data = top_labelled, 
                        mapping = aes(label = gene), 
                        size = 5,
                        nudge_x = .15,
                        box.padding = 0.5,
                        nudge_y = 1,
                        segment.curvature = -0.1,
                        segment.ncp = 3,
                        segment.angle = 20)+
      theme(legend.position="none")+
      ylim(0,130) + xlim(-1,1)
dev.off()




##############################################
# Vulcano plots Cocaine CPP Vs Cocaine NPAS4 #
##############################################
df2 <- DE_cocaine %>% 
    group_by(cell_type) %>%
        mutate(LOG = -log10(p_val_adj), ABS = abs(avg_logFC)) %>% 
        mutate(Threshold = if_else(p_val_adj < 0.05 & ABS > 0.2, "TRUE","FALSE")) %>%
        mutate(Direction = case_when(avg_logFC > 0.2 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.2 & p_val_adj < 0.05 ~ "DownReg")) 


pdf("output_figures/output_Figure2/Vulcano_CocCPPVsCocNPAS4_AllCells.pdf", width = 4, height = 4)
df2 %>%
ggplot(aes(x=avg_logFC, y=LOG)) +
ggrastr::geom_point_rast(aes(colour = cell_type),size=0.5) +
      xlab("log2(Fold Change)")+ 
      ylab("-log10(FDR)")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_vline(xintercept = 0.2, colour = "grey",linetype="dashed",size=1,alpha=0.5) +
      geom_vline(xintercept = -0.2, colour = "grey",linetype="dashed",size=1,alpha=0.5) + 
      geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      theme(legend.position="none")+
       scale_colour_manual(values = colors) +   # customized color palette
      ylim(0,130) + xlim(-3,3)
dev.off()


pdf("output_figures/output_Figure2/Vulcano_CocCPPVsCocNpas4_Drd2_1.pdf", width = 4, height = 4)
top_labelled <- df2 %>% 
                  filter(cell_type == "MSN_Drd2+_1", gene %in% c("Penk", "Oprm1", "Cartpt", "Tac1", "Robo2", "Reln","Calm1","Celf3")) 

#Ppp1r1b, Calb1, Calm1, Penk, Slc1a2, Kcna1, Cartpt, Tac1, Camkv, Camk2n1


df2 %>%
  filter(cell_type == "MSN_Drd2+_1") %>%
  ggscatter( 
            x = "avg_logFC", 
            y = "LOG",
            color = "Threshold",
            palette=c("grey","red"),
            size = 1,
            alpha=0.3,
            shape=19)+
      xlab("log2(Fold Change)")+ 
      ylab("-log10(FDR)")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_vline(xintercept = 0.2, colour = "grey",linetype="dashed",size=1,alpha=0.5) +
      geom_vline(xintercept = -0.2, colour = "grey",linetype="dashed",size=1,alpha=0.5) + 
      geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_text_repel(data = top_labelled, 
                        mapping = aes(label = gene), 
                        size = 5,
                        nudge_x = .15,
                        box.padding = 0.5,
                        nudge_y = 1,
                        segment.curvature = -0.1,
                        segment.ncp = 3,
                        segment.angle = 20)+
      theme(legend.position="none")+
      ylim(0,130) + xlim(-1,1)
dev.off()



pdf("output_figures/output_Figure2/Vulcano_CocCPPVsCocNpas4_Drd1_1.pdf", width = 4, height = 4)
top_labelled <- df2 %>% 
                  filter(cell_type == "MSN_Drd1+_1", gene %in% c("Syne1", "Tac1", "Clk1", "Grm1", "Tcf4", "Kcnk1","Pdzd2","Actb")) 

#Ppp1r1b, Calb1, Calm1, Penk, Slc1a2, Kcna1, Cartpt, Tac1, Camkv, Camk2n1


df2 %>%
  filter(cell_type == "MSN_Drd1+_1") %>%
  ggscatter( 
            x = "avg_logFC", 
            y = "LOG",
            color = "Threshold",
            palette=c("grey","red"),
            size = 1,
            alpha=0.3,
            shape=19)+
      xlab("log2(Fold Change)")+ 
      ylab("-log10(FDR)")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_vline(xintercept = 0.2, colour = "grey",linetype="dashed",size=1,alpha=0.5) +
      geom_vline(xintercept = -0.2, colour = "grey",linetype="dashed",size=1,alpha=0.5) + 
      geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_text_repel(data = top_labelled, 
                        mapping = aes(label = gene), 
                        size = 5,
                        nudge_x = .15,
                        box.padding = 0.5,
                        nudge_y = 1,
                        segment.curvature = -0.1,
                        segment.ncp = 3,
                        segment.angle = 20)+
      theme(legend.position="none")+
      ylim(0,130) + xlim(-1,1)
dev.off()


############################################
# Vulcano plots Saline CPP Vs Saline NPAS4 #
############################################

df3 <- DE_saline %>% 
    group_by(cell_type) %>%
        mutate(LOG = -log10(p_val_adj), ABS = abs(avg_logFC)) %>% 
        mutate(Threshold = if_else(p_val_adj < 0.05 & ABS > 0.2, "TRUE","FALSE")) %>%
        mutate(Direction = case_when(avg_logFC > 0.2 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.2 & p_val_adj < 0.05 ~ "DownReg")) 


pdf("output_figures/output_Figure2/Vulcano_SalCPPVsSalNPAS4_AllCells.pdf", width = 4, height = 4)
df3 %>%
ggplot(aes(x=avg_logFC, y=LOG)) +
ggrastr::geom_point_rast(aes(colour = cell_type),size=0.5) +
      xlab("log2(Fold Change)")+ 
      ylab("-log10(FDR)")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_vline(xintercept = 0.2, colour = "grey",linetype="dashed",size=1,alpha=0.5) +
      geom_vline(xintercept = -0.2, colour = "grey",linetype="dashed",size=1,alpha=0.5) + 
      geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      theme(legend.position="none")+
       scale_colour_manual(values = colors) +   # customized color palette
      ylim(0,130) + xlim(-3,3)
dev.off()


pdf("output_figures/output_Figure2/Vulcano_SalCPPVsSalNPAS4_Drd2_1.pdf", width = 4, height = 4)
top_labelled <- df3 %>% 
                  filter(cell_type == "MSN_Drd2+_1", gene %in% c("Penk", "Tac1", "Cacna1h", "Calm1", "Drd1", "Grm1", "Slc32a1", "Pde10a")) 

#Ppp1r1b, Calb1, Calm1, Penk, Slc1a2, Kcna1, Cartpt, Tac1, Camkv, Camk2n1

df3 %>%
  filter(cell_type == "MSN_Drd2+_1") %>%
  ggscatter( 
            x = "avg_logFC", 
            y = "LOG",
            color = "Threshold",
            palette=c("grey","red"),
            size = 1,
            alpha=0.3,
            shape=19)+
      xlab("log2(Fold Change)")+ 
      ylab("-log10(FDR)")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_vline(xintercept = 0.2, colour = "grey",linetype="dashed",size=1,alpha=0.5) +
      geom_vline(xintercept = -0.2, colour = "grey",linetype="dashed",size=1,alpha=0.5) + 
      geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_text_repel(data = top_labelled, 
                        mapping = aes(label = gene), 
                        size = 5,
                        nudge_x = .15,
                        box.padding = 0.5,
                        nudge_y = 1,
                        segment.curvature = -0.1,
                        segment.ncp = 3,
                        segment.angle = 20)+
      theme(legend.position="none")+
      ylim(0,130) + xlim(-1,1)
dev.off()


pdf("output_figures/output_Figure2/Vulcano_SalCPPVsSalNPAS4_Drd1_1.pdf", width = 4, height = 4)
top_labelled <- df3 %>% 
                  filter(cell_type == "MSN_Drd1+_1", gene %in% c("Actn2", "Ptn", "Syt11", "Syp1", "Pcdh7", "Acvr1c", "Cntn5", "Nsg1")) 

#Ppp1r1b, Calb1, Calm1, Penk, Slc1a2, Kcna1, Cartpt, Tac1, Camkv, Camk2n1


df3 %>%
  filter(cell_type == "MSN_Drd1+_1") %>%
  ggscatter( 
            x = "avg_logFC", 
            y = "LOG",
            color = "Threshold",
            palette=c("grey","red"),
            size = 1,
            alpha=0.3,
            shape=19)+
      xlab("log2(Fold Change)")+ 
      ylab("-log10(FDR)")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_vline(xintercept = 0.2, colour = "grey",linetype="dashed",size=1,alpha=0.5) +
      geom_vline(xintercept = -0.2, colour = "grey",linetype="dashed",size=1,alpha=0.5) + 
      geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_text_repel(data = top_labelled, 
                        mapping = aes(label = gene), 
                        size = 5,
                        nudge_x = .15,
                        box.padding = 0.5,
                        nudge_y = 1,
                        segment.curvature = -0.1,
                        segment.ncp = 3,
                        segment.angle = 20)+
      theme(legend.position="none")+
      ylim(0,130) + xlim(-1,1)
dev.off()



###############################################
# Vulcano plots Saline NPAS4 Vs Cocaine NPAS4 #
###############################################
df4 <- DE_salcoc_npas4 %>% 
    group_by(cell_type) %>%
        mutate(LOG = -log10(p_val_adj), ABS = abs(avg_logFC)) %>% 
        mutate(Threshold = if_else(p_val_adj < 0.05 & ABS > 0.2, "TRUE","FALSE")) %>%
        mutate(Direction = case_when(avg_logFC > 0.2 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.2 & p_val_adj < 0.05 ~ "DownReg")) 


pdf("output_figures/output_Figure2/Vulcano_SalNPAS4VsCocNPAS4_AllCells.pdf", width = 4, height = 4)
df4 %>%
ggplot(aes(x=avg_logFC, y=LOG)) +
ggrastr::geom_point_rast(aes(colour = cell_type),size=0.5) +
      xlab("log2(Fold Change)")+ 
      ylab("-log10(FDR)")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_vline(xintercept = 0.2, colour = "grey",linetype="dashed",size=1,alpha=0.5) +
      geom_vline(xintercept = -0.2, colour = "grey",linetype="dashed",size=1,alpha=0.5) + 
      geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      theme(legend.position="none")+
       scale_colour_manual(values = colors) +   # customized color palette
      ylim(0,130) + xlim(-3,3)
dev.off()

pdf("output_figures/output_Figure2/Vulcano_SalNPAS4VsCocNPAS4_Drd2_1.pdf", width = 4, height = 4)
top_labelled <- df4 %>% 
                  filter(cell_type == "MSN_Drd2+_1", gene %in% c("Ppp1r1b", "Calb1", "Calm1", "Penk", "Slc1a2", "Kcna1", "Cartpt", "Camkv")) 

#Ppp1r1b, Calb1, Calm1, Penk, Slc1a2, Kcna1, Cartpt, Tac1, Camkv, Camk2n1


df4 %>%
  filter(cell_type == "MSN_Drd2+_1") %>%
  ggscatter( 
            x = "avg_logFC", 
            y = "LOG",
            color = "Threshold",
            palette=c("grey","red"),
            size = 1,
            alpha=0.3,
            shape=19)+
      xlab("log2(Fold Change)")+ 
      ylab("-log10(FDR)")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_vline(xintercept = 0.2, colour = "grey",linetype="dashed",size=1,alpha=0.5) +
      geom_vline(xintercept = -0.2, colour = "grey",linetype="dashed",size=1,alpha=0.5) + 
      geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_text_repel(data = top_labelled, 
                        mapping = aes(label = gene), 
                        size = 5,
                        nudge_x = .15,
                        box.padding = 0.5,
                        nudge_y = 1,
                        segment.curvature = -0.1,
                        segment.ncp = 3,
                        segment.angle = 20)+
      theme(legend.position="none")+
      ylim(0,130) + xlim(-1,1)
dev.off()

pdf("output_figures/output_Figure2/Vulcano_SalNPAS4VsCocNPAS4_Drd1_1.pdf", width = 4, height = 4)
top_labelled <- df4 %>% 
                  filter(cell_type == "MSN_Drd1+_1", gene %in% c("Eef1a2", "Pcbp2", "Atp1a1", "Syt4", "Calm2", "Kcna1", "Kcna1", "Usp50")) 

#Ppp1r1b, Calb1, Calm1, Penk, Slc1a2, Kcna1, Cartpt, Tac1, Camkv, Camk2n1


df4 %>%
  filter(cell_type == "MSN_Drd1+_1") %>%
  ggscatter( 
            x = "avg_logFC", 
            y = "LOG",
            color = "Threshold",
            palette=c("grey","red"),
            size = 1,
            alpha=0.3,
            shape=19)+
      xlab("log2(Fold Change)")+ 
      ylab("-log10(FDR)")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_vline(xintercept = 0.2, colour = "grey",linetype="dashed",size=1,alpha=0.5) +
      geom_vline(xintercept = -0.2, colour = "grey",linetype="dashed",size=1,alpha=0.5) + 
      geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_text_repel(data = top_labelled, 
                        mapping = aes(label = gene), 
                        size = 5,
                        nudge_x = .15,
                        box.padding = 0.5,
                        nudge_y = 1,
                        segment.curvature = -0.1,
                        segment.ncp = 3,
                        segment.angle = 20)+
      theme(legend.position="none")+
      ylim(0,130) + xlim(-1,1)
dev.off()


##$$################
# Intersection Drd2#
####################
df1 <- tmp1 %>%
    as.data.frame() %>%
    filter(cell_type == "MSN_Drd2+_1") %>%
    mutate(Class = "SalVsCoc_CPP") %>%
    select(gene, Class)

df2 <- tmp2 %>%
    as.data.frame() %>%
    filter(cell_type == "MSN_Drd2+_1") %>%
    mutate(Class = "SalVsCoc_NPAS4") %>%
    select(gene, Class)

df3 <- tmp3 %>%
    as.data.frame() %>%
    filter(cell_type == "MSN_Drd2+_1") %>%
    mutate(Class = "SalVsSal") %>%
    select(gene, Class)

df4 <- tmp4 %>%
    as.data.frame() %>%
    filter(cell_type == "MSN_Drd2+_1") %>%
    mutate(Class = "CocVsCoc") %>%
    select(gene, Class)

pdf("output_figures/output_Figure2/UpsetPlot_Drd2.pdf",width=5,height=3,useDingbats=FALSE)
tmpUpset <- rbind(df1,df2,df3,df4)
l <- split(as.character(tmpUpset$gene),tmpUpset$Class)
Class <- names(l)
ToTGene <- as.numeric(sapply(l, length))
metadata <- as.data.frame(cbind(Class, ToTGene))
names(metadata) <- c("Class", "ToTGene")
metadata$ToTGene <- as.numeric(as.character(metadata$ToTGene))
upset(fromList(l),,nsets = 6, set.metadata = list(data = metadata, plots = list(list(type = "hist", 
    column = "ToTGene", assign = 20), 
    list(type = "matrix_rows", column = "sets", colors = c(CocVsCoc = "#04a3bd", SalVsCoc_CPP = "#f0be3d", SalvsCoc_Npas4 = "#da7901",SalVsSal = "#247d3f"), 
    alpha = 0.5))))
dev.off()


#####################
# Intersection Drd1 #
#####################
df1 <- tmp1 %>%
    as.data.frame() %>%
    filter(cell_type == "MSN_Drd1+_1") %>%
    mutate(Class = "SalVsCoc_CPP") %>%
    select(gene, Class)

df2 <- tmp2 %>%
    as.data.frame() %>%
    filter(cell_type == "MSN_Drd1+_1") %>%
    mutate(Class = "SalVsCoc_NPAS4") %>%
    select(gene, Class)

df3 <- tmp3 %>%
    as.data.frame() %>%
    filter(cell_type == "MSN_Drd1+_1") %>%
    mutate(Class = "SalVsSal") %>%
    select(gene, Class)

df4 <- tmp4 %>%
    as.data.frame() %>%
    filter(cell_type == "MSN_Drd1+_1") %>%
    mutate(Class = "CocVsCoc") %>%
    select(gene, Class)

pdf("output_figures/output_Figure2/UpsetPlot_Drd1.pdf",width=5,height=3,useDingbats=FALSE)
tmpUpset <- rbind(df1,df2,df3,df4)
l <- split(as.character(tmpUpset$gene),tmpUpset$Class)
Class <- names(l)
ToTGene <- as.numeric(sapply(l, length))
metadata <- as.data.frame(cbind(Class, ToTGene))
names(metadata) <- c("Class", "ToTGene")
metadata$ToTGene <- as.numeric(as.character(metadata$ToTGene))
upset(fromList(l),,nsets = 6, set.metadata = list(data = metadata, plots = list(list(type = "hist", 
    column = "ToTGene", assign = 20), 
    list(type = "matrix_rows", column = "sets", colors = c(CocVsCoc = "#04a3bd", SalVsCoc_CPP = "#f0be3d", SalvsCoc_Npas4 = "#da7901",SalVsSal = "#247d3f"), 
    alpha = 0.5))))
dev.off()


################
# Heatmap Drd2 #
################
genes <- unique(unlist(l))


df1 <- DE_salcoc %>%
    as.data.frame() %>%
    filter(cell_type == "MSN_Drd2+_1") %>%
    select(gene, avg_logFC)

df2 <- DE_salcoc_npas4 %>%
    as.data.frame() %>%
    filter(cell_type == "MSN_Drd2+_1") %>%
    select(gene, avg_logFC)

df3 <- DE_saline %>%
    as.data.frame() %>%
    filter(cell_type == "MSN_Drd2+_1") %>%
    select(gene, avg_logFC)

df4 <- DE_cocaine %>%
    as.data.frame() %>%
    filter(cell_type == "MSN_Drd2+_1") %>%
    select(gene, avg_logFC)

listdata <- list(df1,df2,df3,df4)

new_names <- c("SalVsCoc_CPP","SalVsCoc_NPAS4","SalVsSal","CocVsCoc")
names(listdata) <- new_names


pdf("output_figures/output_Figure2/FoldChange_Heatmap_Drd2.pdf", width=2, height=5)
map(listdata, ~ (.x %>% select(gene, avg_logFC))) %>% 
      map2(new_names, ~setnames(.x, 'avg_logFC', .y)) %>%
      purrr::reduce(left_join, by = "gene") %>% 
      filter(gene %in% genes) %>% 
      column_to_rownames("gene") %>% 
      pheatmap::pheatmap(
          color         = viridis::viridis(100), #colorRampPalette(c("navy", "white", "firebrick3"))(50),
          border_color  = NA,
          show_colnames = TRUE,
          show_rownames = FALSE,
          scale="row"
        )
dev.off()



################
# Heatmap Drd1 #
################
genes <- unique(unlist(l))


df1 <- DE_salcoc %>%
    as.data.frame() %>%
    filter(cell_type == "MSN_Drd1+_1") %>%
    select(gene, avg_logFC)

df2 <- DE_salcoc_npas4 %>%
    as.data.frame() %>%
    filter(cell_type == "MSN_Drd1+_1") %>%
    select(gene, avg_logFC)

df3 <- DE_saline %>%
    as.data.frame() %>%
    filter(cell_type == "MSN_Drd1+_1") %>%
    select(gene, avg_logFC)

df4 <- DE_cocaine %>%
    as.data.frame() %>%
    filter(cell_type == "MSN_Drd1+_1") %>%
    select(gene, avg_logFC)

listdata <- list(df1,df2,df3,df4)

new_names <- c("SalVsCoc_CPP","SalVsCoc_NPAS4","SalVsSal","CocVsCoc")
names(listdata) <- new_names


pdf("output_figures/output_Figure2/FoldChange_Heatmap_Drd1.pdf", width=2, height=5)
map(listdata, ~ (.x %>% select(gene, avg_logFC))) %>% 
      map2(new_names, ~setnames(.x, 'avg_logFC', .y)) %>%
      purrr::reduce(left_join, by = "gene") %>% 
      filter(gene %in% genes) %>% 
      column_to_rownames("gene") %>% 
      pheatmap::pheatmap(
          color         = viridis::viridis(100), #colorRampPalette(c("navy", "white", "firebrick3"))(50),
          border_color  = NA,
          show_colnames = TRUE,
          show_rownames = FALSE,
          scale="row"
        )
dev.off()



#############################
# Intersection Drd1 vs Drd2 #
#############################

# SalVsCoc_CPP

a <- tmp1 %>%
    as.data.frame() %>%
    filter(cell_type == "MSN_Drd1+_1") %>%
    mutate(Class = "Drd1") %>%
    select(gene, Class)

b <- tmp1 %>%
    as.data.frame() %>%
    filter(cell_type == "MSN_Drd2+_1") %>%
    mutate(Class = "Drd2") %>%
    select(gene, Class)


pdf("output_figures/output_Figure2/UpsetPlot_Drd1vsDrd2_SalCPPVsCocCPP.pdf",width=5,height=3,useDingbats=FALSE)
tmpUpset <- rbind(a,b)
l <- split(as.character(tmpUpset$gene),tmpUpset$Class)
Class <- names(l)
ToTGene <- as.numeric(sapply(l, length))
metadata <- as.data.frame(cbind(Class, ToTGene))
names(metadata) <- c("Class", "ToTGene")
metadata$ToTGene <- as.numeric(as.character(metadata$ToTGene))
upset(fromList(l),nsets = 6, set.metadata = list(data = metadata, plots = list(list(type = "hist", 
    column = "ToTGene", assign = 20), 
    list(type = "matrix_rows", column = "sets", colors = c(Drd1 = "#67E5C8", Drd2 = "#D85BA6"), 
    alpha = 0.5))))
dev.off()

# DE_salcoc_npas4

a <- tmp2 %>%
    as.data.frame() %>%
    filter(cell_type == "MSN_Drd1+_1") %>%
    mutate(Class = "Drd1") %>%
    select(gene, Class)

b <- tmp2 %>%
    as.data.frame() %>%
    filter(cell_type == "MSN_Drd2+_1") %>%
    mutate(Class = "Drd2") %>%
    select(gene, Class)


pdf("output_figures/output_Figure2/UpsetPlot_Drd1vsDrd2_SalNPAS4VsCocNPAS4.pdf",width=5,height=3,useDingbats=FALSE)
tmpUpset <- rbind(a,b)
l <- split(as.character(tmpUpset$gene),tmpUpset$Class)
Class <- names(l)
ToTGene <- as.numeric(sapply(l, length))
metadata <- as.data.frame(cbind(Class, ToTGene))
names(metadata) <- c("Class", "ToTGene")
metadata$ToTGene <- as.numeric(as.character(metadata$ToTGene))
upset(fromList(l),nsets = 6, set.metadata = list(data = metadata, plots = list(list(type = "hist", 
    column = "ToTGene", assign = 20), 
    list(type = "matrix_rows", column = "sets", colors = c(Drd1 = "#67E5C8", Drd2 = "#D85BA6"), 
    alpha = 0.5))))
dev.off()

# DE_saline

a <- tmp3 %>%
    as.data.frame() %>%
    filter(cell_type == "MSN_Drd1+_1") %>%
    mutate(Class = "Drd1") %>%
    select(gene, Class)

b <- tmp3 %>%
    as.data.frame() %>%
    filter(cell_type == "MSN_Drd2+_1") %>%
    mutate(Class = "Drd2") %>%
    select(gene, Class)


pdf("output_figures/output_Figure2/UpsetPlot_Drd1vsDrd2_SalCPPVsSalNPAS4.pdf",width=5,height=3,useDingbats=FALSE)
tmpUpset <- rbind(a,b)
l <- split(as.character(tmpUpset$gene),tmpUpset$Class)
Class <- names(l)
ToTGene <- as.numeric(sapply(l, length))
metadata <- as.data.frame(cbind(Class, ToTGene))
names(metadata) <- c("Class", "ToTGene")
metadata$ToTGene <- as.numeric(as.character(metadata$ToTGene))
upset(fromList(l),nsets = 6, set.metadata = list(data = metadata, plots = list(list(type = "hist", 
    column = "ToTGene", assign = 20), 
    list(type = "matrix_rows", column = "sets", colors = c(Drd1 = "#67E5C8", Drd2 = "#D85BA6"), 
    alpha = 0.5))))
dev.off()

# DE_cocaine

a <- tmp4 %>%
    as.data.frame() %>%
    filter(cell_type == "MSN_Drd1+_1") %>%
    mutate(Class = "Drd1") %>%
    select(gene, Class)

b <- tmp4 %>%
    as.data.frame() %>%
    filter(cell_type == "MSN_Drd2+_1") %>%
    mutate(Class = "Drd2") %>%
    select(gene, Class)


pdf("output_figures/output_Figure2/UpsetPlot_Drd1vsDrd2_CocCPPVsCocNPAS4.pdf",width=5,height=3,useDingbats=FALSE)
tmpUpset <- rbind(a,b)
l <- split(as.character(tmpUpset$gene),tmpUpset$Class)
Class <- names(l)
ToTGene <- as.numeric(sapply(l, length))
metadata <- as.data.frame(cbind(Class, ToTGene))
names(metadata) <- c("Class", "ToTGene")
metadata$ToTGene <- as.numeric(as.character(metadata$ToTGene))
upset(fromList(l),nsets = 6, set.metadata = list(data = metadata, plots = list(list(type = "hist", 
    column = "ToTGene", assign = 20), 
    list(type = "matrix_rows", column = "sets", colors = c(Drd1 = "#67E5C8", Drd2 = "#D85BA6"), 
    alpha = 0.5))))
dev.off()

# Scatter plot fold changes Drd1 vs Drd2
# Fold Change Correlations
tmp1 <- DE_salcoc %>% 
        dplyr::select(cell_type,gene,avg_logFC) %>%
        dplyr::rename(logFC_SalCPPVsCocCPP = avg_logFC) %>%
        filter(cell_type %in% c("MSN_Drd2+_1","MSN_Drd1+_1")) %>%
        pivot_wider(names_from = cell_type, values_from = logFC_SalCPPVsCocCPP) %>%
        mutate(Class = "SalCPPVsCocCPP")

tmp2 <- DE_salcoc_npas4 %>% 
        dplyr::select(cell_type,gene,avg_logFC) %>%
        dplyr::rename(logFC_SalNpas4VsCocNpas4 = avg_logFC) %>%
        filter(cell_type %in% c("MSN_Drd2+_1","MSN_Drd1+_1")) %>%
        pivot_wider(names_from = cell_type, values_from = logFC_SalNpas4VsCocNpas4) %>%
        mutate(Class = "SalNpas4VsCocNpas4")

tmp3 <- DE_saline %>% 
        dplyr::select(cell_type,gene,avg_logFC) %>%
        dplyr::rename(logFC_SalCppVsSalNpas4 = avg_logFC) %>%
        filter(cell_type %in% c("MSN_Drd2+_1","MSN_Drd1+_1")) %>%
        pivot_wider(names_from = cell_type, values_from = logFC_SalCppVsSalNpas4) %>%
        mutate(Class = "SalCppVsSalNpas4")


tmp4 <- DE_cocaine %>% 
        dplyr::select(cell_type,gene,avg_logFC) %>%
        dplyr::rename(logFC_CocCppVsCocNpas4 = avg_logFC) %>%
        filter(cell_type %in% c("MSN_Drd2+_1","MSN_Drd1+_1")) %>%
        pivot_wider(names_from = cell_type, values_from = logFC_CocCppVsCocNpas4) %>%
        mutate(Class = "CocCppVsCocNpas4")


combined <- rbind(tmp4,tmp1,tmp3,tmp2)
colnames(combined) <- c("gene","Drd2","Drd1","Class")


pdf("output_figures/output_Figure2/Scatter_FoldChange.pdf",width=5,height=5,useDingbats=FALSE)
combined %>%
mutate(Class=fct_relevel(Class,c("SalCPPVsCocCPP","CocCppVsCocNpas4","SalCppVsSalNpas4","SalNpas4VsCocNpas4"))) %>%
ggplot(aes(x=Drd1, y=Drd2)) +
ggrastr::geom_point_rast(aes(colour = Class),size=0.5) +
      xlab("Drd1 log2(Fold Change)")+ 
      ylab("Drd2 log2(Fold Change)")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_hline(yintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      theme(legend.position="none") +
    facet_wrap(.~Class) + 
    xlim(-1.5,1.5) + 
    ylim(-1.5,1.5) + 
  stat_cor(aes(color = Class), method = "spearman", label.y = 1.5) +
  scale_colour_manual(values = c("#F0BE3D","#04A3BD","#247D3F","#DA7901"))   # customized color palette
dev.off()
