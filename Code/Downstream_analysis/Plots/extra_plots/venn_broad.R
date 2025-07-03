## Comapring overlap broad
## Secretory vs proliferative
library(ggplot2)
library(Seurat)
library(dittoSeq)
library(ggrepel)
library(dplyr)

dir_endo <-"/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/pristine_disease/broad/"

r_filter <- function(df){
  x<- df %>% dplyr::filter(avg_log2FC > .25 & p_val_adj < 0.05 )
  genesup <- paste0(rownames(x),"_UP")
  x<- df %>% dplyr::filter(avg_log2FC < -0.25 & p_val_adj < 0.05 )
  genesdown <- paste0(rownames(x),"_DOWN")
  return(c(genesup, genesdown))
}


files_endo_secretory <- list.files( path=dir_endo,pattern = "(_|/)secretory(_DEG\\.RDS)", full.names = TRUE)
files_endo_proliferative <- list.files( path=dir_endo,pattern = "(_|/)proliferative(_DEG\\.RDS)", full.names = TRUE)

names(files_endo_secretory)<-basename(files_endo_secretory)
names(files_endo_proliferative)<-basename(files_endo_proliferative)


files_list <- list(files_endo_secretory,files_endo_proliferative)

rds_objects <- lapply(files_list, function(group) lapply(group, readRDS))
significant <- lapply(rds_objects, function(group) lapply(group, r_filter))


pdf("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/SFBroadVen/venn_broad.pdf",width = 4, height = 4)
for (i in seq_len(length(significant[[1]]))) {
  x<-c(Proliferative=significant[[2]][i],Secretory=significant[[1]][i])
  gplot<-ggVennDiagram(x) + 
    scale_fill_gradient(low = "#F4FAFE", high = "#F4FAFE") +  
    scale_color_manual(values = c("1" = "lightblue", "2" = "orange"))+
  ggtitle(names_celltypes_sec[i])
  print(gplot)
}


dev.off()
