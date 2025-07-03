library(Seurat)
library(tidyverse)
library(dplyr)
library(dittoSeq)
library(viridis)
library(tibble)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)

#dir1 <- ("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/pristine_disease/broad/Stromal_fibroblast_proliferative_DEG.RDS")
#dir2 <- ("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/pristine_disease/broad/Stromal_fibroblast_secretory_DEG.RDS")
#dir_out <- "/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/SFig8/"

dir1 <- ("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/pristine_disease/broad/Unciliated_Epithelia_proliferative_DEG.RDS")
dir2 <- ("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/pristine_disease/broad/Unciliated_Epithelia_secretory_DEG.RDS")
dir_out <- "/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/Sfig11/"

files<-list(dir1,dir2)

DEG_genes <- lapply(files, readRDS)

r_c <- function(x,dataname){
  x <- x %>% filter (abs(avg_log2FC) > .25 & p_val_adj < 0.05)
  # x <- x %>% filter ( p_val_adj < 0.05)
  x<- tibble::rownames_to_column(x, "Gene")
  x<- dplyr::select(x, c(Gene,avg_log2FC,p_val_adj))
  x
}
DEG_genes <- lapply(DEG_genes, r_c)

merged<- DEG_genes %>% purrr::reduce(full_join, by='Gene') %>% na.exclude()

new_colnames <- c( "Proliferative","Secretory")

new_colnames_lfc <- paste0("LFC_",new_colnames)
new_colnames_p <- paste0("pvalue_",new_colnames)
new_colnames<-c(rbind(new_colnames_lfc, new_colnames_p))
new_colnames<- c("Gene",new_colnames)
colnames(merged) <-new_colnames

mat <- as.matrix(merged %>% dplyr::select(new_colnames_lfc))
rownames(mat)<-merged$Gene

heatmap <- ComplexHeatmap::Heatmap(mat, name = "logFC", 
                                   cluster_rows = T, cluster_columns = T, show_row_names = T,
                                    row_title = "Gene Id", 
                                   row_names_gp = grid::gpar(fontsize = 6), row_title_side = "left", column_names_gp = grid::gpar(fontsize = 8),
                                   column_title_side = "top", column_title_gp = gpar(fontsize = 3), row_dend_reorder = T)


pdf(paste0(dir_out, "heatmap_overlap.pdf"), height = 14, width = 5)
heatmap
dev.off()

