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

load("output_reclust/Brandon_SCT_30pcs_04res_Reclust_NoDoublet_Slimmed.RData")

labels <- read.table("output_reclust/Labels_Clusters.txt",header=T,sep="\t")

new <- labels %>% 
                separate(variable, c("Cell", "Class","Definition"),"_",remove = FALSE) %>%
                separate(Rows, c("ID", "Cluster"),"_",remove = FALSE) %>%
                unite(Ident, c("Cell","Definition"),sep = "_",remove = FALSE,na.rm = TRUE) %>%
                unite(Ident2, c("Cell","Cluster","Definition"),sep = "_",remove = FALSE,na.rm = TRUE)


new <- new[order(as.numeric(as.character(new$Ident2))), ]

current.cluster.ids <- new$Cluster
new.cluster.ids <- as.character(new$variable)

seuObject_slim_nodoub_slim@active.ident <- plyr::mapvalues(x = seuObject_slim_nodoub_slim@active.ident, 
                                                     from = current.cluster.ids, 
                                                     to = new.cluster.ids)

# UMAP plot with new labels
seuObject_slim_nodoub_slim@meta.data$Cell <- seuObject_slim_nodoub_slim@active.ident

save(seuObject_slim_nodoub_slim,file = "output_reclust/Brandon_SCT_30pcs_04res_Reclust_NoDoublet_Slimmed_Relabel.RData")


umap <- as.data.frame(Embeddings(seuObject_slim_nodoub_slim, reduction = "umap"))
meta <- as.data.frame(seuObject_slim_nodoub_slim@meta.data)

df <- cbind(umap,meta)%>% 
  group_by(Cell) %>% 
  mutate(N = n()) %>%
    ungroup() %>% 
    arrange(Genotype) %>%
    as.data.frame()  


label <- data.frame(Cell=levels(df$Cell),label=levels(df$Cell))

label_2 <- df %>% 
  group_by(Cell) %>% 
  summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2),N = n()) %>% 
  left_join(label) %>%
  as.data.frame() 

col <- randomcoloR::distinctColorPalette(20)
write.table(col, "output_figures/output_Figure1/Colors_Used_Initial_Data.txt", sep="\t",quote=F)
#col <- data.table::fread("output_sct_reclust/Figure1/Colors_Used_Initial_Data.txt")
#col <- as.character(col$Colors)

pdf("output_figures/output_Figure1/Data_Integrated_UMAP_Labels_Alternative.pdf", width = 8, height = 8)
#n <- 60
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',] 
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ggplot(df, aes(x=UMAP_1, y=UMAP_2)) +
ggrastr::geom_point_rast(aes(colour = Cell),size=0.5) +
ggrepel::geom_text_repel(data = label_2, aes(label = label),
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
    scale_colour_manual(values = col)+
    theme_void() +
    theme(legend.position="none")

dev.off()

# Color table
colors_used <- data.frame(Cells = unique(df$Cell), Colors = col)
write.table(colors_used, "output_figures/output_Figure1/Colors_Used_Initial_Data.txt",sep="\t",quote=F)

# Proportion by Genotype

order <- c("MSN_Drd1+_1","MSN_Drd1+_2","MSN_Drd1+_3","MSN_Drd2+_1","MSN_Drd2+_2","MSN_Drd3+_1","MSN_Grm8+_1","MSN_Grm8+_2",
  "Interneurons_Chat+","Interneurons_Th+","Interneurons_Sst+","Interneurons_Pvalb+","Interneurons_Pnoc+","Interneurons_Ndnf+",
  "Astrocytes","Oligodendrocytes","OPC","Microglia","Endothelial","Neuroblast")

colore <- c("#17154f","#b0799a","#e69b00","#355828")

prop_cell <-  seuObject_slim_nodoub_slim@meta.data %>%
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
ggplot(df, aes(x=UMAP_1, y=UMAP_2)) +
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
Drd1<- FetchData(seuObject_slim_nodoub_slim, vars = "Drd1")
Drd2<- FetchData(seuObject_slim_nodoub_slim, vars = "Drd2") 
Drd3<- FetchData(seuObject_slim_nodoub_slim, vars = "Drd3") 
Snap25<- FetchData(seuObject_slim_nodoub_slim, vars = "Snap25") 
Mobp<- FetchData(seuObject_slim_nodoub_slim, vars = "Mobp") 
C1qa<- FetchData(seuObject_slim_nodoub_slim, vars = "C1qa") 
Pdgfra<- FetchData(seuObject_slim_nodoub_slim, vars = "Pdgfra") 
Gja1<- FetchData(seuObject_slim_nodoub_slim, vars = "Gja1")
Vtn<- FetchData(seuObject_slim_nodoub_slim, vars = "Vtn") 
Flt1<-FetchData(seuObject_slim_nodoub_slim, vars = "Flt1")
Col4a6<-FetchData(seuObject_slim_nodoub_slim, vars = "Col4a6")
Igfbpl1<-FetchData(seuObject_slim_nodoub_slim, vars = "Igfbpl1")


GeneExp_Data <- cbind(umap,meta)%>% 
      select(-integrated_snn_res.0.6,-SCT_snn_res.0.5,-SCT_snn_res.0.3,-Doublets,-nCount_RNA,-nFeature_RNA) %>%
      mutate(Drd1 = Drd1$Drd1, Drd2 = Drd2$Drd2, Drd3 = Drd3$Drd3,
             Snap25 = Snap25$Snap25, Mobp = Mobp$Mobp, C1qa = C1qa$C1qa,
             Gja1 = Gja1$Gja1, Pdgfra = Pdgfra$Pdgfra, Vtn = Vtn$Vtn, 
             Flt1 = Flt1$Flt1, Col4a6 = Col4a6$Col4a6, Igfbpl1 = Igfbpl1$Igfbpl1)


a <- ggplot(GeneExp_Data, aes(x=UMAP_1, y=UMAP_2)) +
ggrastr::geom_point_rast(aes(colour = Snap25),size=0.1, alpha=0.1) + 
    scale_colour_gradient(low = "gray90", high = "red") +
    #scale_color_viridis(option="plasma") +
    theme_void() +
    theme(legend.position="none") + 
    ggtitle("Neuron: Snap25")

b <- ggplot(GeneExp_Data, aes(x=UMAP_1, y=UMAP_2)) +
ggrastr::geom_point_rast(aes(colour = Drd1),size=0.1, alpha=0.1) + 
    scale_colour_gradient(low = "gray90", high = "red") +
    #scale_color_viridis(option="plasma") +
    theme_void() +
    theme(legend.position="none") + 
    ggtitle("Drd1")

c <- ggplot(GeneExp_Data, aes(x=UMAP_1, y=UMAP_2)) +
ggrastr::geom_point_rast(aes(colour = Drd2),size=0.1, alpha=0.1) + 
    scale_colour_gradient(low = "gray90", high = "red") +
    #scale_color_viridis(option="plasma") +
    theme_void() +
    theme(legend.position="none") + 
    ggtitle("Drd2")

d <- ggplot(GeneExp_Data, aes(x=UMAP_1, y=UMAP_2)) +
ggrastr::geom_point_rast(aes(colour = Drd3),size=0.1, alpha=0.1) + 
    scale_colour_gradient(low = "gray90", high = "red") +
    #scale_color_viridis(option="plasma") +
    theme_void() +
    theme(legend.position="none") + 
    ggtitle("Drd3")

e <- ggplot(GeneExp_Data, aes(x=UMAP_1, y=UMAP_2)) +
ggrastr::geom_point_rast(aes(colour = Mobp),size=0.1, alpha=0.1) + 
    scale_colour_gradient(low = "gray90", high = "red") +
    #scale_color_viridis(option="plasma") +
    theme_void() +
    theme(legend.position="none") + 
    ggtitle("Oligo: Mobp")

f <- ggplot(GeneExp_Data, aes(x=UMAP_1, y=UMAP_2)) +
ggrastr::geom_point_rast(aes(colour = C1qa),size=0.1, alpha=0.1) + 
    scale_colour_gradient(low = "gray90", high = "red") +
    #scale_color_viridis(option="plasma") +
    theme_void() + 
    theme(legend.position="none") + 
    ggtitle("Microglia: C1qa")

g <- ggplot(GeneExp_Data, aes(x=UMAP_1, y=UMAP_2)) +
ggrastr::geom_point_rast(aes(colour = Pdgfra),size=0.1, alpha=0.1) + 
    scale_colour_gradient(low = "gray90", high = "red") +
    #scale_color_viridis(option="plasma") +
    theme_void() +
    theme(legend.position="none") + 
    ggtitle("OPC: Pdgfra")

h <- ggplot(GeneExp_Data, aes(x=UMAP_1, y=UMAP_2)) +
ggrastr::geom_point_rast(aes(colour = Vtn),size=0.1, alpha=0.1) + 
    scale_colour_gradient(low = "gray90", high = "red") +
    #scale_color_viridis(option="plasma") +
    theme_void() +
    theme(legend.position="none") + 
    ggtitle("Mural: Vtn")

i <- ggplot(GeneExp_Data, aes(x=UMAP_1, y=UMAP_2)) +
ggrastr::geom_point_rast(aes(colour = Flt1),size=0.1, alpha=0.1) + 
    scale_colour_gradient(low = "gray90", high = "red") +
    #scale_color_viridis(option="plasma") +
    theme_void() +
    theme(legend.position="none") + 
    ggtitle("Endothelial: Flt1")

l <- ggplot(GeneExp_Data, aes(x=UMAP_1, y=UMAP_2)) +
ggrastr::geom_point_rast(aes(colour = Gja1),size=0.1, alpha=0.1) + 
    scale_colour_gradient(low = "gray90", high = "red") +
    #scale_color_viridis(option="plasma") +
    theme_void() +
    theme(legend.position="none") + 
    ggtitle("Astrocytes: Gja1")


m <- ggplot(GeneExp_Data, aes(x=UMAP_1, y=UMAP_2)) +
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
Idents(seuObject_slim_nodoub_slim) <- factor(Idents(seuObject_slim_nodoub_slim), levels = order)

pdf("output_figures/output_Figure1/Dot_Plot_Markers.pdf", width = 7, height = 6)
cd_genes <- c("Drd1","Drd2","Drd3","Grm8","Chat","Th","Sst","Pvalb","Elavl2","Pde3a","Gja1","Mobp",
              "Pdgfra","C1qa","Flt1","Igfbpl1",
              "Snap25","Ppp1r1b","Foxp2","Syt1")
DotPlot(object = seuObject_slim_nodoub_slim, features = cd_genes) + 
      theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
          theme(legend.position="top") 
dev.off()


# Do heatmap
markers<- presto::wilcoxauc(seuObject_slim_nodoub_slim, 'Cell', assay = 'data')
markers<- presto::top_markers(markers, n = 20, auc_min = 0.5, pct_in_min = 20, pct_out_max = 20)

all_markers<- markers %>%
  select(-rank) %>% 
  unclass() %>% 
  stack() %>%
  pull(values) %>%
  unique() %>%
   .[!is.na(.)]

cols <- viridis::viridis(100)[c(1, 50, 100)]

downsamples <- subset(seuObject_slim_nodoub_slim, downsample = 300)


mat<- downsamples[["RNA"]]@data[all_markers, ] %>% as.matrix()
mat<- t(scale(t(mat)))
cluster_anno<- downsamples@meta.data$Cell

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
        top_annotation = ComplexHeatmap::HeatmapAnnotation(foo = anno_block(gp = grid::gpar(fill = col))),
        show_column_names = FALSE,
        use_raster = TRUE,
        show_row_names = FALSE,
        raster_quality = 4)
dev.off()


# Predict with Chen's data
load("utils/NAc_Mouse_Chen.RData")

downsamples@meta.data$Data <- "Brandon"
chen_seuobj@meta.data$Data <- "Chen"
chen_seuobj@meta.data$Cell <- chen_seuobj@meta.data$orig.ident

# scPred using Chen Data
chen_seuobj <- chen_seuobj %>% 
                FindVariableFeatures(selection.method = "vst", nfeatures = 3000, verbose = FALSE) %>% 
                ScaleData(vars.to.regress = c("percent.mito","nUMI"), model.use = "linear") %>%
                RunPCA(features=NULL, weight.by.var = TRUE, npcs = 30, reduction.name = "pca")  


reference <- getFeatureSpace(chen_seuobj, "CellType")
reference <- trainModel(reference,model = "mda")
#downsamples <- subset(seuObject_slim_nodoub_slim, downsample = 300)
query <- scPredict(seuObject_slim_nodoub_slim, reference, max.iter.harmony = 30)

mat <- crossTab(query, "Cell", "scpred_prediction", output = "prop") %>% 
                as.data.frame() %>%
                rownames_to_column("Cell") %>%
                 pivot_longer(!Cell, names_to = "Cluster", values_to = "Pred") %>%
                 group_by(Cluster) %>%
                 filter(!(Cell %in% "unassigned")) %>%
            filter(Pred == max(Pred)) %>% 
            arrange(Cluster) %>% 
            as.data.frame()

# Colors
col <- data.table::fread("output_figures/output_Figure1/Colors_Used_Initial_Data.txt") %>%
        as.data.frame()

mat <- merge(mat,col,by.x="Cluster",by.y="Cells")

mat$Cluster <- factor(mat$Cluster, levels = order)


pdf("output_figures/output_Figure1/Chen_Pred.pdf", width = 4, height = 5)
ggbarplot(mat, x = "Cluster", y = "Pred",
          fill = "Cluster",           # change fill color by mpg_level
          color = "black",            # Set bar border colors to white
          palette = mat$Colors,            # jco journal color palett. see ?ggpar
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 45,          # Rotate vertically x axis texts
          ylab = "Chen et al, Predicted %",
          rotate = TRUE,
          ggtheme = theme_classic()
          ) +
    theme(legend.position="none")
dev.off()



mat2 <- crossTab(query, "Cell", "scpred_prediction", output = "prop") %>% 
                as.matrix()
mat2 <- mat2[-10,]

labels2 <- c("D1","D2","IN","Astro","Oligo","OPC","Micro","Endo","NB")


mat2 <- mat2[match(labels2,rownames(mat2)),match(order,colnames(mat2))]


pdf("output_figures/output_Figure1/Chen_Pred_Heatmap.pdf", width = 8, height = 5)
pheatmap::pheatmap(
  mat               = mat2,
  color             = viridis::viridis(n = 256, alpha = 1, begin = 0, end = 1, option = "viridis"),
  border_color      = NA,
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  drop_levels       = TRUE,
  cluster_cols      = FALSE,
  cluster_rows      = FALSE,
  fontsize          = 12
  )
dev.off()

openxlsx::write.xlsx(mat, file = "output_figures/output_Figure1/ScPred_ChenData.xlsx", 
                        colNames = TRUE, 
                        borders = "columns",
                        overwrite=T)
write.table(mat,"output_figures/output_Figure1/ScPred_ChenData.txt",sep="\t",quote=F,row.names=F)


