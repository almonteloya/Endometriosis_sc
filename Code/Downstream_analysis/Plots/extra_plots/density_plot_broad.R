#### UMAP density plot

library(ggplot2)
library(Seurat)
library(dittoSeq)
library(dplyr)

dir_in='/seurat_objects/'
dir_out<- '/outdir/'

dir.create(dir_out)

s_obj = readRDS(file=paste0(dir_in,"merged_data_celltype_all_2025.RDS"))

meta_w_umap = cbind(s_obj@meta.data,s_obj@reductions$umap@cell.embeddings)

pdf(paste0(dir_out,"umap_density_condition.pdf"),height =6, width = 12)
ggplot(meta_w_umap, aes(x=UMAP_1, y=UMAP_2))+
  geom_density_2d_filled(contour_var = "ndensity") +
  #scale_fill_brewer(palette="Oranges")+
  scale_fill_manual(values=c("white", '#fff5eb','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801','#a63603','#7f2704'))+
  theme_bw()+
  theme(panel.grid=element_blank())+
  labs(fill="percentile\n(per plot)")+
  facet_wrap(vars(condition_control), nrow=1)
dev.off()


pdf(paste0(dir_out,"umap_density_phase.pdf"),height =6, width = 10)
ggplot(meta_w_umap, aes(x=UMAP_1, y=UMAP_2))+
  geom_density_2d_filled(contour_var = "ndensity") +
  #scale_fill_brewer(palette="Oranges")+
  scale_fill_manual(values=c("white", '#fff5eb','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801','#a63603','#7f2704'))+
  theme_bw()+
  theme(panel.grid=element_blank())+
  labs(fill="percentile\n(per plot)")+
  facet_wrap(vars(Phase_menstrual_reduced), nrow=1)
dev.off()
