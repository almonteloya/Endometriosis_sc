

library(ggplot2)
library(Seurat)
library(dittoSeq)
library(dplyr)

dir='/krummellab/data1/immunox/XREP1a/10x/merged_XREP/subset_celltype/Mesenchymal/'

dir_in <- '/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/'
dir_out<-"/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/Fig3/"
s_obj <- readRDS(sprintf("%s/Stromal_obj_2025.RDS",dir_in))


markers <- readRDS(paste0(dir,"filtered_markers.RDS"))
markers %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC) -> top5

genes_plot <- unique(top5$gene)
pdf(paste0(dir_out,"stromal_filter_markers_find.pdf"),height =4, width=10)
Idents(s_obj)<-s_obj$annotation
print(DotPlot(s_obj, features=genes_plot, dot.scale=6, cols=c("yellow", "red")) +
        theme(axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5, hjust=1))) 
dev.off()

