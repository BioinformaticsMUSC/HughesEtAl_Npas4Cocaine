NPAS4 controls cell type-specific circuit adaptations underlying drug seeking behavior 
==========================

![](snRNA-seq Experimental Design_HughesBW.jpf)

This repository contains analysis code for the single nuclei RNA-seq project carried out by researchers at the [Cowan lab, MUSC](https://medicine.musc.edu/departments/neuroscience/research/cowan) and [Berto Lab, MUSC](https://bertolab.org/)

## Cite this

If you use anything in this repository please cite the following publication:

Pre-print URL: 

## Files

| directory | contents | code |
| --------- | -------- | -------- |
| [`input`](input/) | Input/Output data of the initial processing and quality check. | 01_Downstream_Analysis_SCT.r|
| [`output_initial`](output_initial/) | Output data of the initial clustering and integration. | 01_Downstream_Analysis_SCT.r \ 02_Initial_Markers_Detection.r|
| [`output_reclust`](output_reclust/) | Output data of the reclustering and subset analyses. | 03_Recluster_NoGlut.r \ 04_Doubletting.r \ 05_Markers_Detection_NoGlut.r \ 06_Relabel_And_InitialViz.r|
| [`output_dge`](output_dge/) | Output data of the DGE analyses. | 07_Differential_Expression.r |
| [`output_figures`](output_figures/) | Output for figures and additional visualizations. | 08_DGE_visualizations.r |
| [`shinyApp`](shinyApp/) | Output of the ShinyApp. | 09_ShinyApp.r|
