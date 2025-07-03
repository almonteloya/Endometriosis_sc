
## Hestmap diases vs control

library(Seurat)
library(tidyverse)
library(dplyr)
library(ComplexHeatmap)



celltype_heatmap <- function(celltype){
  files <- all_files[grepl(celltype, all_files)]
  all_markers <- lapply(files, readRDS)
  pattern <- ".*\\/([^/]+)\\.RDS$"
  names_file <- sub(pattern, "\\1", files)
  names(all_markers)<- names_file
  
  
  
  filter_mark <- function(f_markers){
    {{f_markers}}  %>%  tibble::rownames_to_column(.,"Gene")  %>% filter (abs(avg_log2FC) > .25) %>% dplyr::filter (p_val_adj < 0.05) %>% dplyr::select("Gene",avg_log2FC)
      }
  
  
  filtered_markers <- lapply(all_markers, filter_mark)
  
  lapply(filtered_markers, nrow)
  
  mat<-(filtered_markers %>% purrr::reduce(full_join, by='Gene')) 
  
  mat <- mat[rowSums(is.na(mat[, colnames(mat)])) < 2,]
  
  words <- c("disease", "Fibroid", "Pristine")
  names_file_short <- words[str_detect(names_file, words)]
  
  colnames(mat)<- c("Gene",names_file_short)
  
  df_long <-  as.data.frame(mat) %>%
    pivot_longer(cols = names_file_short, names_to = "variable", values_to = "value")
  
  
  mat <-as.matrix (mat %>% remove_rownames %>% 
                     column_to_rownames(var="Gene")) 
  
  
  
  yb<-colorRampPalette(c("goldenrod","white","purple"))(50)
  yb<-colorRampPalette(c("#56B4E9","white","darkorange"))(50)
  
  #yb<-colorRampPalette(c("purple","white","darkgreen"))(50)
  
  
  
  heatmap <- Heatmap(mat, name = "logFC", col=yb,
                     cluster_rows = F, cluster_columns = T, show_row_names = TRUE,
                     column_title = "celltype", row_title = "Gene", 
                     row_names_gp = gpar(fontsize = 8), row_title_side = "left",
                     column_title_side = "top", column_title_gp = gpar(fontsize = 8), row_dend_reorder = F)      
  heatmap
  
  
  heatmap <- Heatmap(mat, name = "logFC", col=yb,
                     cluster_rows = T, cluster_columns = T, show_row_names = F,
                     column_title = celltype, row_title = "Gene", 
                     row_names_gp = gpar(fontsize = 8), row_title_side = "left",
                     column_title_side = "top", column_title_gp = gpar(fontsize = 8), row_dend_reorder = T)      
  heatmap
  
  col_dend = dendsort::dendsort(hclust(dist(t(mat))))
  
  heatmap <- Heatmap(mat, name = "logFC", col=yb, cluster_columns = col_dend,
                     cluster_rows = T, show_row_names = F,
                     column_title = celltype, row_title = "Gene", 
                     row_names_gp = gpar(fontsize = 8), row_title_side = "left",
                     column_title_side = "top", column_title_gp = gpar(fontsize = 8), row_dend_reorder = T)      
  heatmap
}


dir_in1= "/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/menstrual_phase/"
patt <- paste0(".*", ".*RDS")

all_files <- list.files( path= dir_in1, pattern = patt, full.names = TRUE)

pdf("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/SPhase/heatmap_unciliated.pdf",width = 2)
print(celltype_heatmap("Unciliated"))
dev.off()

pdf("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/SPhase/heatmap_fibro.pdf",width = 2)
print(celltype_heatmap("fibroblast"))
dev.off()
