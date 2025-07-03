#### UMAP PLOTS 
### April 2024

#### STROMAL

dir_in <- '/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/'
dir_out<-"/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/Fig3/"

s_obj <- readRDS(sprintf("%s/Stromal_obj_2025.RDS",dir_in))


pdf(paste0(dir_out,"umap_celltype2.pdf"),height =8, width = 7)
Idents(s_obj)<- s_obj$labels_singleR_reduced
DimPlot(s_obj, cols = dittoColors()[c(15, 12)]) + 
  theme( legend.position = "bottom",      # Move the legend to the bottom
         legend.text = element_text(size = 20))  # Increase legend text size)
dev.off()


pdf(paste0(dir_out,"umap_celltype.pdf"),height =8, width = 12)
Idents(s_obj)<- s_obj$annotation_updated
DimPlot(s_obj,cols=dittoColors()) + theme(legend.text=element_text(size=20))
dev.off()

pdf(paste0(dir_out,"umap_celltype_color.pdf"),height =8, width = 12)
Idents(s_obj)<- s_obj$annotation_updated
DimPlot(s_obj,cols=dittoColors()[c(5,6,8,9,10,11,12,13)]) + theme(legend.text=element_text(size=20))
dev.off()



