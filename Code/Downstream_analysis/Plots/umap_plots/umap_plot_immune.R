### UMAP for immune cells
dir_out <- '/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/Fig5/'
dir2='/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/'

s_obj = readRDS(file=paste0(dir2,"Immune_obj_2025.RDS"))

pdf(paste0(dir_out,"umap_immune_celltype_1.pdf"),height =8, width = 12)
Idents(s_obj)<- s_obj$immune_agregate1_c
DimPlot(s_obj,cols=dittoColors()[c(5,6,8,9,10,11,12,13,14,15,16,17)]) + theme(legend.text=element_text(size=20))
dev.off()

dir_out <- '/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/SFig15/'
pdf(paste0(dir_out,"umap_immune_celltype_2.pdf"),height =8, width = 12)
Idents(s_obj)<- s_obj$NK_annotation
DimPlot(s_obj,cols=dittoColors()) + theme(legend.text=element_text(size=20))
dev.off()

