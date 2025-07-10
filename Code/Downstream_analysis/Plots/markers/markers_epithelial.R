## Markers for epithelial


library(ggplot2)
library(Seurat)
library(dittoSeq)
library(dplyr)

dir_in_marker='epithelial/object/dir'
s_obj=readRDS('/Epithelial_obj_2025.RDS')
##Markers from find markers
markers<- readRDS( paste0(dir_in_marker,"epithelial_filter_markers.RDS"))

dir_out<-"/Fix2025/Fig4/"

markers %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC) -> top3

genes_plot <- unique(top3$gene)
genes_plot<- c(genes_plot,"MUC5B","FOXJ1","PIFO") ## ciliated markers
pdf(paste0(dir_out,"epithelial_filter_markers_find.pdf"),height =4, width=10)
Idents(s_obj)<-s_obj$annotation
print(DotPlot(s_obj, features=genes_plot, dot.scale=6, cols=c("yellow", "red")) +
        theme(axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5, hjust=1))) 
dev.off()

