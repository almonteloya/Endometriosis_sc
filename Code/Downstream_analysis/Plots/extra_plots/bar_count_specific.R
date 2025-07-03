### CLEANup pelase
## Barplot of number ofr DEG upregulated and downregulated by cell type

library(dplyr)
library(ggplot2)
library(ggrepel)



data_generate <- function(dir, phase) {
  ## Generate the tables with the number of genes by each phase
  
  ## Get the RDS files of the DEG 
  files <- list.files( path= dir,pattern = "RDS", full.names = TRUE)
  filtered_markers <- lapply(files, readRDS)
  names_file <- list.files( path= dir,pattern = "RDS", full.names = F)
  names(filtered_markers)<-names_file
  
  ## Count genes  
  dim_cat_fil_pos<- function(f_markers){
    nrow({{f_markers}} %>% dplyr::filter (p_val_adj < 0.05 & abs(avg_log2FC) > .25  & avg_log2FC > 0))
  }
  dim_cat_fil_neg<- function(f_markers){
    nrow({{f_markers}} %>% dplyr::filter (p_val_adj < 0.05 & abs(avg_log2FC) > .25  & avg_log2FC < 0))
  }
  
  ## Apply and divide by upregulated and downregulated
  DEG_num_log_pos <- unlist(sapply(filtered_markers,dim_cat_fil_pos))
  DEG_num_log_neg <- unlist(sapply(filtered_markers,dim_cat_fil_neg))
  #Concat 
  DEG_count<-data.frame(celltype=names_file, Upregulated = DEG_num_log_pos,Downregulated=DEG_num_log_neg)
  DEG_count
  
  ## Plot
  ggplot(data=DEG_count, aes(x=celltype, y=Negative, fill=Negative)) +
    geom_bar(stat="identity", position=position_dodge())
  DEG_count$group <- phase
  DEG_count
}



### Stromal cells

dir <- ("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/pristine_disease/stromal_specific/secretory/")
secretory_genes<- data_generate(dir,"secretory")
dir <- ("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/pristine_disease/stromal_specific/proliferative/")
proliferative_genes<- data_generate(dir,"proliferative")
melted_df_pro <- reshape2::melt(proliferative_genes,  id.vars = "celltype",measure.vars = c("Upregulated", "Downregulated"))

melted_df_sec <- reshape2::melt(secretory_genes,  id.vars = "celltype",measure.vars = c("Upregulated", "Downregulated"))



#### Select genes

#cell_selected<-c( "CD8+-T-cells-_DEG.RDS" ,"Macrophages_DEG.RDS", "NK-cells_DEG.RDS" ,       "Monocytes_DEG.RDS")
#cell_selected<-c( "Glandular_1.RDS" ,"EMT_epithelia_2.RDS", "EMT_epithelia_1.RDS" ,"Ciliated_epithelia.RDS")
#cell_selected <- c("HSP_stromal_fibroblast.RDS", "Proliferative_fibroblast.RDS", "SFRP4.RDS", "Stromal_fibroblast_1.RDS", "Stromal_fibroblast_2.RDS", "Stromal_fibroblast_3.RDS", "vascular_smooth_muscle.RDS")

#melted_df_sec<-melted_df_sec[melted_df_sec$celltype %in% cell_selected,]
pdf("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/Fig3/barplot_secretory_stromal.pdf", height = 4,width = 8)
#pdf("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/DEG2/viz_2024/barplot_secretory_epithelial.pdf", height = 4,width = 8)

ggplot(melted_df_sec, aes(x = celltype, y = value,fill=variable)) +
  geom_bar(stat = "identity") + scale_fill_manual(values=c('red','blue'))  + ylim(0,150)+
  theme_minimal() +  coord_flip() + scale_x_discrete(position = "top")  + theme(text = element_text(size = 20),legend.position="bottom")
dev.off()

#melted_df_pro<-melted_df_pro[melted_df_pro$celltype %in% cell_selected,]

pdf("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/Fig3/barplot_proliferative_stromal.pdf",  height = 4,width = 8)
#pdf("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/DEG2/viz_2024/barplot_proliferative_epithelial.pdf", height = 4,width = 8)
 
ggplot(melted_df_pro, aes(x = celltype, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c('red', 'blue')) +
  theme_minimal() +
  coord_flip() +
  scale_x_discrete(position = "top") +
  theme(text = element_text(size = 20), legend.position = "bottom") +
  scale_y_reverse(limits = c(150, 0))

dev.off()




### Epithelial 

dir <- ("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/pristine_disease/epithelial_specific/secretory/")
secretory_genes<- data_generate(dir,"secretory")
dir <- ("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/pristine_disease/epithelial_specific/proliferative/")
proliferative_genes<- data_generate(dir,"proliferative")
melted_df_pro <- reshape2::melt(proliferative_genes,  id.vars = "celltype",measure.vars = c("Upregulated", "Downregulated"))

melted_df_sec <- reshape2::melt(secretory_genes,  id.vars = "celltype",measure.vars = c("Upregulated", "Downregulated"))



#### Select genes

#cell_selected<-c( "CD8+-T-cells-_DEG.RDS" ,"Macrophages_DEG.RDS", "NK-cells_DEG.RDS" ,       "Monocytes_DEG.RDS")
cell_selected<-c( "Glandular-1_secretory_DEG.RDS", "EMT-epithelia-1_secretory_DEG.RDS" ,"Ciliated-epithelia_secretory_DEG.RDS")

melted_df_sec<-melted_df_sec[melted_df_sec$celltype %in% cell_selected,]
pdf("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/Fig4/barplot_secretory_epi.pdf", height = 4,width = 8)

ggplot(melted_df_sec, aes(x = celltype, y = value,fill=variable)) +
  geom_bar(stat = "identity") + scale_fill_manual(values=c('red','blue'))  + ylim(0,900)+
  theme_minimal() +  coord_flip() + scale_x_discrete(position = "top")  + theme(text = element_text(size = 20),legend.position="bottom")
dev.off()


cell_selected<-c( "Glandular-1_proliferative_DEG.RDS", "EMT-epithelia-1_proliferative_DEG.RDS" ,
                  "Ciliated-epithelia_proliferative_DEG.RDS")
melted_df_pro<-melted_df_pro[melted_df_pro$celltype %in% cell_selected,]

pdf("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/Fig4/barplot_proliferative_epi.pdf",  height = 4,width = 8)

ggplot(melted_df_pro, aes(x = celltype, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c('red', 'blue')) +
  theme_minimal() +
  coord_flip() +
  scale_x_discrete(position = "top") +
  theme(text = element_text(size = 20), legend.position = "bottom") +
  scale_y_reverse(limits = c(900, 0))

dev.off()




#### Immune


dir <- ("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/pristine_disease/immune_broad/secretory/")
secretory_genes<- data_generate(dir,"secretory")
dir <- ("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/pristine_disease/immune_broad/proliferative/")
proliferative_genes<- data_generate(dir,"proliferative")


## Adding a row to make it simetrical 
extrainfo<- data.frame(celltype="MAST_secretory_DEG.RDS",Upregulated=0,Downregulated=0,group="secretory")
secretory_genes<-rbind(secretory_genes,extrainfo)


melted_df_pro <- reshape2::melt(proliferative_genes,  id.vars = "celltype",measure.vars = c("Upregulated", "Downregulated"))
melted_df_sec <- reshape2::melt(secretory_genes,  id.vars = "celltype",measure.vars = c("Upregulated", "Downregulated"))


#### Select genes


cell_selected<-c( "CD8+_T_cells_secretory_DEG.RDS" , "Monocytes_secretory_DEG.RDS", 
                  "NK_cells_secretory_DEG.RDS",  "MAST_secretory_DEG.RDS","Macrophages_secretory_DEG.RDS")

melted_df_sec<-melted_df_sec[melted_df_sec$celltype %in% cell_selected,]
pdf("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/Fig5/barplot_secretory_immune.pdf", height = 4,width = 8)

ggplot(melted_df_sec, aes(x = celltype, y = value,fill=variable)) +
  geom_bar(stat = "identity") + scale_fill_manual(values=c('red','blue'))  + ylim(0,120)+
  theme_minimal() +  coord_flip() + scale_x_discrete(position = "top")  + theme(text = element_text(size = 20),legend.position="bottom")
dev.off()



cell_selected<-c( "CD8+_T_cells_proliferative_DEG.RDS" ,"Monocytes_proliferative_DEG.RDS", 
                  "NK_cells_proliferative_DEG.RDS",  "MAST_proliferative_DEG.RDS","Macrophages_proliferative_DEG.RDS")

melted_df_pro<-melted_df_pro[melted_df_pro$celltype %in% cell_selected,]

pdf("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/Fig5/barplot_proliferative_immune.pdf",  height = 4,width = 8)

ggplot(melted_df_pro, aes(x = celltype, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c('red', 'blue')) +
  theme_minimal() +
  coord_flip() +
  scale_x_discrete(position = "top") +
  theme(text = element_text(size = 20), legend.position = "bottom") +
  scale_y_reverse(limits = c(120, 0))

dev.off()

