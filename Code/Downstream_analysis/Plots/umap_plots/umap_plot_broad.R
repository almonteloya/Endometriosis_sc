#### UMAP plots of broad celltypes
### April 2025

#### All cells


library(ggplot2)
library(Seurat)
library(dittoSeq)
library(dplyr)

dir_in='/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/seurat_objects/'
dir_out<- '/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/Fig2/'

dir.create(dir_out)

s_obj = readRDS(file=paste0(dir_in,"merged_data_celltype_all_2025.RDS"))

pdf(paste0(dir_out,"umap_celltype.pdf"),height =8, width = 12)
Idents(s_obj)<- s_obj$labels_singleR_reduced
DimPlot(s_obj,cols=dittoColors()[c(9,15,13,10,12,7)])+ theme(legend.text=element_text(size=20))
DimPlot(s_obj,cols=dittoColors()[c(9,15,13,10,12,7)])+ (NoLegend())
dev.off()

pdf(paste0(dir_out,"umap_phase.pdf"),height =8, width = 12)
Idents(s_obj)<- s_obj$Phase_menstrual_reduced
DimPlot(s_obj,cols=dittoColors()[c(1,2)])+ theme(legend.text=element_text(size=20))
DimPlot(s_obj,cols=dittoColors()[c(1,2)])+ NoLegend()
dev.off()


pdf(paste0(dir_out,"umap_condition.pdf"),height =8, width = 12)
Idents(s_obj)<- s_obj$condition_control
DimPlot(s_obj,cols=dittoColors()[c(3,7,4)]) + theme(legend.text=element_text(size=20))
DimPlot(s_obj,cols=dittoColors()[c(3,7,4)]) + NoLegend()
dev.off()



