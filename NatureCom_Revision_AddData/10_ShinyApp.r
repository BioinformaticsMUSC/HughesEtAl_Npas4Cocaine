library(rsconnect)
library(ShinyCell)
library(Seurat)
library(tidyverse)


load("output_reclust/Jessica_SCT_30pcs_05res_Reclust_NoDoublet_Relabel.RData")


seu_filt <- seuObject_nodoub_label

seu_filt@meta.data <- seu_filt@meta.data %>% 
                            select(Genotype, pMito, nCount_RNA, nFeature_RNA, seurat_clusters, Cell)


scConf = createConfig(seu_filt)
makeShinyApp(seu_filt, scConf, gene.mapping = TRUE,
             shiny.title = "Jessica - NAc, Cocaine, and Npas4") 

rsconnect::setAccountInfo(name='bioinformatics-musc',
              token='838A2925138D0715F2D093E909823204',
              secret='suQLpnyt0hEXRa+LXT9IQxaYNCSIITnH2vBymGPJ')

options(repos = BiocManager::repositories())

rsconnect::deployApp("shinyApp/",
    appName = 'Jessica_NAc_Cocaine_Npas4', 
    account = 'bioinformatics-musc', 
    server = 'shinyapps.io')