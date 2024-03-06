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

load("output_reclust/Jessica_SCT_30pcs_05res_Reclust_NoDoublet_Relabel.RData")
load("output_dge/DGE_Mast_Analysis.RData")

# Order the cells
order <- c("MSN_Drd1+_1","MSN_Drd1+_2","MSN_Drd1+_3","MSN_Drd2+_1","MSN_Drd2+_2","MSN_Drd3+_1","MSN_Grm8+_1","MSN_Grm8+_2",
  "Interneurons_Chat+","Interneurons_Th+","Interneurons_Sst+","Interneurons_Pvalb+","Interneurons_Pnoc+","Interneurons_Ndnf+",
  "Astrocytes","Oligodendrocytes","OPC","Microglia","Endothelial","Neuroblast")



## Visualize number of Npas4 Expressing cells 
Npas4 = GetAssayData(object=seuObject_nodoub_label,assay="RNA",slot="data")["Npas4",]

seuObject_nodoub_label@meta.data$Npas4 <- as.numeric(log2(Npas4+1))

umap <- as.data.frame(Embeddings(seuObject_nodoub_label, reduction = "umap"))
meta <- as.data.frame(seuObject_nodoub_label@meta.data)

# Npas4 expression violin plot
pdf("output_figures/output_Figure1/Npas4_ViolinPlot.pdf", width = 8, height = 8)
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
        stat_compare_means(aes(label = ..p.signif..),label.x = 1.5) +
        facet_wrap(.~Cell)
dev.off()



# Cell number total
colore <- c("#17154f","#b0799a")
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
                dplyr::rename(Cell = Var1, Genotype = Var2, Total = Freq) %>%
                filter(Cell %in% grouped_data$Cell)

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
tmp1 <- DE_drug %>% 
        group_by(cell_type) %>%
        mutate(LOG = -log10(p_val_adj), ABS = abs(avg_logFC)) %>% 
        mutate(Threshold = if_else(p_val_adj < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
        mutate(Direction = case_when(avg_logFC > 0.2 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.2 & p_val_adj < 0.05 ~ "DownReg")) %>%
        filter(Threshold == "TRUE") %>% 
        dplyr::select(cell_type,gene,avg_logFC)

# Radar plot
df1 <- table(tmp1$cell_type) %>%
    as.data.frame() %>%
    mutate(Class = "DE_drug")

df <- df1 %>% 
  pivot_wider(names_from = Var1, values_from = Freq) %>%
  as.data.frame()

labels <- c("Class",order)

df <- df[,match(labels,colnames(df))]
df[is.na(df)] <- 0

pdf("output_figures/output_Figure2/RadarPlot_Degs.pdf", width = 8, height = 8)
ggradar(
  df, 
  values.radar = c("0", "750", "1500"),
  grid.min = 0, grid.mid = 750, grid.max = 1500,
  # Polygons
  group.line.width = 1, 
  group.point.size = 3,
  group.colours = c("#04a3bd"),
  # Background and grid lines
  background.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = "bottom"
  )
dev.off()

###########################################
# Vulcano plots Saline CPP Vs Cocaine CPP #
###########################################

df1 <- DE_drug %>% 
    group_by(cell_type) %>%
        mutate(LOG = -log10(p_val_adj), ABS = abs(avg_logFC)) %>% 
        mutate(Threshold = if_else(p_val_adj < 0.05 & ABS > 0.2, "TRUE","FALSE")) %>%
        mutate(Direction = case_when(avg_logFC > 0.2 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.2 & p_val_adj < 0.05 ~ "DownReg")) 

# cell colors original
col <- read.csv("utils/Colors_Used_Initial_Data_Brandon.txt", sep = "\t", skip = 0, header = TRUE, comment.char = "", check.names = FALSE)
col <- col[match(order, col$Cells),]
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
      ylim(0,130) + xlim(-8,8)
dev.off()

pdf("output_figures/output_Figure2/Vulcano_Drd2_1.pdf", width = 4, height = 4)
top_labelled <- df1 %>% 
                  filter(cell_type == "MSN_Drd2+_1") %>%
                  group_by(Direction) %>% 
                  na.omit() %>%
                  arrange(p_val_adj) %>%
                  top_n(n = 5, wt = LOG) 

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
      ylim(0,100) + xlim(-3,3)
dev.off()


pdf("output_figures/output_Figure2/Vulcano_Drd1_1.pdf", width = 4, height = 4)
top_labelled <- df1 %>% 
                  filter(cell_type == "MSN_Drd1+_1") %>%
                  group_by(Direction) %>% 
                  na.omit() %>%
                  arrange(p_val_adj) %>%
                  top_n(n = 5, wt = LOG) 


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
      ylim(0,100) + xlim(-3,3)
dev.off()




