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
library(scbp)
library(scPred)
library(ComplexHeatmap)
})


# create dir for initial figure 1
dir.create("output_figures")
dir.create("output_figures/output_Figure1")

#plan("multiprocess", workers = 10)
options(future.globals.maxSize = 6000 * 1024^2)

# Transfer label
load("output_reclust/Jessica_SCT_30pcs_05res_Reclust_NoDoublet.RData")
jessica <- seuObject_nodoub

load("utils/Brandon_SCT_30pcs_04res_Reclust_NoDoublet_Slimmed_Relabel.RData")
brandon <- seuObject_slim_nodoub_slim

DefaultAssay(jessica) <- "RNA"
DefaultAssay(brandon) <- "RNA"


brandon <- NormalizeData(brandon)
jessica <- NormalizeData(jessica)

brandon <- FindVariableFeatures(brandon)
jessica <- FindVariableFeatures(jessica)

anchors <- FindTransferAnchors(reference = brandon, query = jessica)

# transfer labels
predictions <- TransferData(anchorset = anchors, refdata = brandon$Cell)
jessica <- AddMetaData(object = jessica, metadata = predictions)



seuObject_nodoub_label <- jessica

save(seuObject_nodoub_label,file = "output_reclust/Jessica_SCT_30pcs_05res_Reclust_NoDoublet_Relabel.RData")


umap <- as.data.frame(Embeddings(seuObject_nodoub_label, reduction = "umap"))
meta <- as.data.frame(seuObject_nodoub_label@meta.data)

df <- cbind(umap,meta)%>% 
  dplyr::rename(Cell = predicted.id,Genotype = sample) %>%
  group_by(Cell) %>% 
  mutate(N = n()) %>%
    ungroup() %>% 
    arrange(Genotype) %>%
    as.data.frame()  


label_2 <- df %>% 
  group_by(Cell) %>% 
  summarize(umap_1 = median(umap_1), umap_2 = median(umap_2),N = n()) %>% 
  as.data.frame() 


# Colors
col <- read.csv("utils/Colors_Used_Initial_Data_Brandon.txt", sep = "\t", skip = 0, header = TRUE, comment.char = "", check.names = FALSE)
col <- col[match(label_2$Cell, col$Cells),]
colors <- as.character(col$Colors)
df$Cell <- factor(df$Cell, levels=col$Cells)

pdf("output_figures/output_Figure1/Data_Integrated_UMAP_Labels_Alternative.pdf", width = 8, height = 8)
#n <- 60
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',] 
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ggplot(df, aes(x=umap_1, y=umap_2)) +
ggrastr::geom_point_rast(aes(colour = Cell),size=0.5) +
ggrepel::geom_text_repel(data = label_2, aes(label = Cell),
              color = "black",
              #fontface = 'bold',
              segment.colour = "grey60",
                box.padding = unit(0.25, "lines"),
                point.padding = unit(0.5, "lines"),
                nudge_x = .15,
                nudge_y = 1,
                size = 6) + 
    #scale_color_viridis(discrete=TRUE,option="inferno")
    #scale_colour_manual(values = colorRampPalette(brewer.pal(13, "Paired"))(31))
    #scale_colour_manual(values = met.brewer(name="Klimt",n=13))+    
    scale_colour_manual(values = colors)+
    theme_void() +
    theme(legend.position="none")

dev.off()


# Proportion by Genotype

order <- c("MSN_Drd1+_1","MSN_Drd1+_2","MSN_Drd1+_3","MSN_Drd2+_1","MSN_Drd2+_2","MSN_Drd3+_1","MSN_Grm8+_1","MSN_Grm8+_2",
  "Interneurons_Chat+","Interneurons_Th+","Interneurons_Sst+","Interneurons_Pvalb+","Interneurons_Pnoc+","Interneurons_Ndnf+",
  "Astrocytes","Oligodendrocytes","OPC","Microglia","Endothelial","Neuroblast")

colore <- c("#17154f","#b0799a")

prop_cell <-  df %>%
  as.data.frame() %>%
  mutate(Cell=fct_relevel(Cell,order)) %>%
  arrange(Cell) %>%
    group_by(Cell,Genotype) %>% 
    summarise(cnt = n()) %>% 
    group_by(Cell) %>%  
     mutate(percent = 100*round(cnt/sum(cnt), 2)) %>% 
     group_by(Cell) %>%
      ggbarplot("Cell", "percent",
    fill = "Genotype", color = "Genotype", palette = colore,
    label = TRUE,lab.pos = "in", lab.col = "white") +
  scale_y_continuous(labels = dollar_format(suffix = "%", prefix = "")) +
rotate_x_text(angle = 45) +
theme(legend.position="right") +
xlab("")
ggsave("output_figures/output_Figure1/CellProportion_By_Genotype.pdf", plot = prop_cell, width = 9, height = 5, units = "in", dpi = 150)

# UMAP by Genotype
pdf("output_figures/output_Figure1/Data_Integrated_ByGenotype.pdf", width = 6, height = 6)
#n <- 60
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',] 
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ggplot(df, aes(x=umap_1, y=umap_2)) +
ggrastr::geom_point_rast(aes(colour = Genotype),size=0.5) + 
    #scale_color_viridis(discrete=TRUE,option="inferno")
    #scale_colour_manual(values = colorRampPalette(brewer.pal(13, "Paired"))(31))
    #scale_colour_manual(values = met.brewer(name="Klimt",n=13))+    
    scale_colour_manual(values = colore)+
    theme_void() +
    theme(legend.position="none") +
    facet_wrap(.~Genotype)
dev.off()


#Fetch Gene expression for markers
Drd1<- FetchData(seuObject_nodoub_label, vars = "Drd1")
Drd2<- FetchData(seuObject_nodoub_label, vars = "Drd2") 
Drd3<- FetchData(seuObject_nodoub_label, vars = "Drd3") 
Snap25<- FetchData(seuObject_nodoub_label, vars = "Snap25") 
Mobp<- FetchData(seuObject_nodoub_label, vars = "Mobp") 
C1qa<- FetchData(seuObject_nodoub_label, vars = "C1qa") 
Pdgfra<- FetchData(seuObject_nodoub_label, vars = "Pdgfra") 
Gja1<- FetchData(seuObject_nodoub_label, vars = "Gja1")
Vtn<- FetchData(seuObject_nodoub_label, vars = "Vtn") 
Flt1<-FetchData(seuObject_nodoub_label, vars = "Flt1")
Col4a6<-FetchData(seuObject_nodoub_label, vars = "Col4a6")
Igfbpl1<-FetchData(seuObject_nodoub_label, vars = "Igfbpl1")


GeneExp_Data <- df %>% 
      mutate(Drd1 = Drd1$Drd1, Drd2 = Drd2$Drd2, Drd3 = Drd3$Drd3,
             Snap25 = Snap25$Snap25, Mobp = Mobp$Mobp, C1qa = C1qa$C1qa,
             Gja1 = Gja1$Gja1, Pdgfra = Pdgfra$Pdgfra, Vtn = Vtn$Vtn, 
             Flt1 = Flt1$Flt1, Col4a6 = Col4a6$Col4a6, Igfbpl1 = Igfbpl1$Igfbpl1)


a <- ggplot(GeneExp_Data, aes(x=umap_1, y=umap_2)) +
ggrastr::geom_point_rast(aes(colour = Snap25),size=0.1, alpha=0.1) + 
    scale_colour_gradient(low = "gray90", high = "red") +
    #scale_color_viridis(option="plasma") +
    theme_void() +
    theme(legend.position="none") + 
    ggtitle("Neuron: Snap25")

b <- ggplot(GeneExp_Data, aes(x=umap_1, y=umap_2)) +
ggrastr::geom_point_rast(aes(colour = Drd1),size=0.1, alpha=0.1) + 
    scale_colour_gradient(low = "gray90", high = "red") +
    #scale_color_viridis(option="plasma") +
    theme_void() +
    theme(legend.position="none") + 
    ggtitle("Drd1")

c <- ggplot(GeneExp_Data, aes(x=umap_1, y=umap_2)) +
ggrastr::geom_point_rast(aes(colour = Drd2),size=0.1, alpha=0.1) + 
    scale_colour_gradient(low = "gray90", high = "red") +
    #scale_color_viridis(option="plasma") +
    theme_void() +
    theme(legend.position="none") + 
    ggtitle("Drd2")

d <- ggplot(GeneExp_Data, aes(x=umap_1, y=umap_2)) +
ggrastr::geom_point_rast(aes(colour = Drd3),size=0.1, alpha=0.1) + 
    scale_colour_gradient(low = "gray90", high = "red") +
    #scale_color_viridis(option="plasma") +
    theme_void() +
    theme(legend.position="none") + 
    ggtitle("Drd3")

e <- ggplot(GeneExp_Data, aes(x=umap_1, y=umap_2)) +
ggrastr::geom_point_rast(aes(colour = Mobp),size=0.1, alpha=0.1) + 
    scale_colour_gradient(low = "gray90", high = "red") +
    #scale_color_viridis(option="plasma") +
    theme_void() +
    theme(legend.position="none") + 
    ggtitle("Oligo: Mobp")

f <- ggplot(GeneExp_Data, aes(x=umap_1, y=umap_2)) +
ggrastr::geom_point_rast(aes(colour = C1qa),size=0.1, alpha=0.1) + 
    scale_colour_gradient(low = "gray90", high = "red") +
    #scale_color_viridis(option="plasma") +
    theme_void() + 
    theme(legend.position="none") + 
    ggtitle("Microglia: C1qa")

g <- ggplot(GeneExp_Data, aes(x=umap_1, y=umap_2)) +
ggrastr::geom_point_rast(aes(colour = Pdgfra),size=0.1, alpha=0.1) + 
    scale_colour_gradient(low = "gray90", high = "red") +
    #scale_color_viridis(option="plasma") +
    theme_void() +
    theme(legend.position="none") + 
    ggtitle("OPC: Pdgfra")

h <- ggplot(GeneExp_Data, aes(x=umap_1, y=umap_2)) +
ggrastr::geom_point_rast(aes(colour = Vtn),size=0.1, alpha=0.1) + 
    scale_colour_gradient(low = "gray90", high = "red") +
    #scale_color_viridis(option="plasma") +
    theme_void() +
    theme(legend.position="none") + 
    ggtitle("Mural: Vtn")

i <- ggplot(GeneExp_Data, aes(x=umap_1, y=umap_2)) +
ggrastr::geom_point_rast(aes(colour = Flt1),size=0.1, alpha=0.1) + 
    scale_colour_gradient(low = "gray90", high = "red") +
    #scale_color_viridis(option="plasma") +
    theme_void() +
    theme(legend.position="none") + 
    ggtitle("Endothelial: Flt1")

l <- ggplot(GeneExp_Data, aes(x=umap_1, y=umap_2)) +
ggrastr::geom_point_rast(aes(colour = Gja1),size=0.1, alpha=0.1) + 
    scale_colour_gradient(low = "gray90", high = "red") +
    #scale_color_viridis(option="plasma") +
    theme_void() +
    theme(legend.position="none") + 
    ggtitle("Astrocytes: Gja1")


m <- ggplot(GeneExp_Data, aes(x=umap_1, y=umap_2)) +
ggrastr::geom_point_rast(aes(colour = Igfbpl1),size=0.1, alpha=0.1) + 
    scale_colour_gradient(low = "gray90", high = "red") +
    #scale_color_viridis(option="plasma") +
    theme_void() +
    theme(legend.position="none") + 
    ggtitle("Neuroblast: Igfbpl1")


plotting <- cowplot::plot_grid(a,b,c,e,f,g,h,i,l,m, ncol = 2, nrow=5)
cowplot::save_plot("output_figures/output_Figure1/FeaturePlot_Basic_Markers.pdf", plotting, 
                     base_height = 6, base_width = 3)
 

# Dot plot
Idents(seuObject_nodoub_label) <- factor(Idents(seuObject_nodoub_label), levels = order)

pdf("output_figures/output_Figure1/Dot_Plot_Markers.pdf", width = 7, height = 6)
cd_genes <- c("Drd1","Drd2","Drd3","Grm8","Chat","Th","Sst","Pvalb","Elavl2","Pde3a","Gja1","Mobp",
              "Pdgfra","C1qa","Flt1","Igfbpl1",
              "Snap25","Ppp1r1b","Foxp2","Syt1")
DotPlot(object = seuObject_nodoub_label, features = cd_genes) + 
      theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
          theme(legend.position="top") 
dev.off()


# Do heatmap
markers<- presto::wilcoxauc(seuObject_nodoub_label, 'predicted.id', assay = 'data')
markers<- presto::top_markers(markers, n = 20, auc_min = 0.5, pct_in_min = 20, pct_out_max = 20)

all_markers<- markers %>%
  select(-rank) %>% 
  unclass() %>% 
  stack() %>%
  pull(values) %>%
  unique() %>%
   .[!is.na(.)]

cols <- viridis::viridis(100)[c(1, 50, 100)]

downsamples <- subset(seuObject_nodoub_label, downsample = 300)


mat<- downsamples[["RNA"]]@data[all_markers, ] %>% as.matrix()
mat<- t(scale(t(mat)))
cluster_anno<- downsamples@meta.data$predicted.id

cluster_anno <- factor(cluster_anno, levels = order)

# Colors
col_fun = circlize::colorRamp2(c(-1, 0, 3), c("#440154FF", "#2D708EFF", "#FDE725FF"))

pdf("output_figures/output_Figure1/heatmap_markers.pdf", width = 8, height = 3)
ComplexHeatmap::Heatmap(mat, name = "Expression",  
        column_split = cluster_anno,
        cluster_columns = FALSE,
        show_column_dend = FALSE,
        cluster_column_slices = TRUE,
        column_title_gp = grid::gpar(fontsize = 8),
        column_gap = unit(0.5, "mm"),
        cluster_rows = FALSE,
        show_row_dend = FALSE,
        col = col_fun,
        row_names_gp = grid::gpar(fontsize = 4),
        column_title_rot = 90,
        top_annotation = ComplexHeatmap::HeatmapAnnotation(foo = anno_block(gp = grid::gpar(fill = colors))),
        show_column_names = FALSE,
        use_raster = TRUE,
        show_row_names = FALSE,
        raster_quality = 4)
dev.off()


# Relabel and slim the data

seuObject_nodoub_label@meta.data <- seuObject_nodoub_label@meta.data %>%
                                      dplyr::rename(Cell = predicted.id,Genotype = sample) %>%
                                      select(orig.ident,Genotype,pMito,nCount_RNA,nFeature_RNA, nCount_SCT,nFeature_SCT,seurat_clusters,SCT_snn_res.0.5,Cell)


save(seuObject_nodoub_label,file = "output_reclust/Jessica_SCT_30pcs_05res_Reclust_NoDoublet_Relabel.RData")




