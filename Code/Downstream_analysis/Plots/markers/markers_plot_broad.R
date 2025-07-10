## Markers for broad celltypes

library(ggplot2)
library(Seurat)
library(dittoSeq)
library(dplyr)

dir_in='/seurat/objects'
dir_out<- '/Sfig1/'

dir.create(dir_out)
s_obj = readRDS(file=paste0(dir_in,"merged_data_celltype_all_2025.RDS"))




## Markers based on literature, see methods 

Stromal <- c("LUM")
Fibroblast <- c("COL1A1", "COL3A1", "COL1A2")
Endothelia <- c("PECAM1","ENG")
Smooth_muscle <- "ACTA2"
Immune <- c("CD14","PTPRC")
Epithelial <- c("EPCAM","KRT19", "KRT18")
Ciliated <-  c("SNTN", "FOXJ1", "TP73", "PIFO")

All_markers<-c(Stromal,Fibroblast,Smooth_muscle,Epithelial,Ciliated,Endothelia,Immune)

markers_order <- c("Immune", "Endothelia", "Ciliated Epithelia", "Unciliated_Epithelia", "Smooth_muscle", "Stromal_fibroblast")

s_obj$labels_singleR_reduced <- factor(s_obj$labels_singleR_reduced, levels = markers_order)
Idents(s_obj) <- s_obj$labels_singleR_reduced


pdf(paste0(dir_out,"Markers_broad.pdf"),height = 5)
DotPlot( s_obj, features=All_markers,cols = c("yellow", "red")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),  
    axis.text.y = element_text(size = 16),  
    axis.title = element_text(size = 16),  
    plot.title = element_text(size = 16)    
  )
dev.off()
