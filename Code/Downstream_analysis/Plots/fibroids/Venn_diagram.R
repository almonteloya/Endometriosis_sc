

## Comparing overlap
## Pristine vs Fibroid
## Pristine vs Disease
library(ggplot2)
library(Seurat)
library(dittoSeq)
library(ggrepel)
library(dplyr)

dir_fibroids <-"DEG/fibroid_pristine/"
dir_endo <-"DEG/pristine_disease/broad/"

r_filter <- function(df){
  x<- df %>% dplyr::filter(avg_log2FC > .25 & p_val_adj < 0.05 )
  genesup <- paste0(rownames(x),"_UP")
  x<- df %>% dplyr::filter(avg_log2FC < -0.25 & p_val_adj < 0.05 )
  genesdown <- paste0(rownames(x),"_DOWN")
  return(c(genesup, genesdown))
}


files_fibroid_secretory <- list.files( path=dir_fibroids,pattern = "(_|/)secretory(_DEG\\.RDS)", full.names = TRUE)
files_fibroid_proliferative <- list.files( path=dir_fibroids,pattern = "(_|/)proliferative(_DEG\\.RDS)", full.names = TRUE)
files_endo_secretory <- list.files( path=dir_endo,pattern = "(_|/)secretory(_DEG\\.RDS)", full.names = TRUE)
files_endo_proliferative <- list.files( path=dir_endo,pattern = "(_|/)proliferative(_DEG\\.RDS)", full.names = TRUE)

names(files_fibroid_secretory)<-basename(files_fibroid_secretory)
names(files_fibroid_proliferative)<-basename(files_fibroid_proliferative )
names(files_endo_secretory)<-basename(files_endo_secretory)
names(files_endo_proliferative)<-basename(files_endo_proliferative)


## Names 
cleaned <- sub("_DEG\\.RDS$", "", basename(files_fibroid_secretory)) ## remove DEG
# Replace underscores with spaces
names_celltypes_sec <- gsub("_", " ", cleaned)

## Names 
cleaned <- sub("_DEG\\.RDS$", "", basename(files_fibroid_proliferative)) ## remove DEG
# Step 2: Replace underscores with spaces
names_celltypes_pro <- gsub("_", " ", cleaned)


files_list <- list(files_fibroid_secretory,files_fibroid_proliferative,
            files_endo_secretory,files_endo_proliferative)
            
rds_objects <- lapply(files_list, function(group) lapply(group, readRDS))
significant <- lapply(rds_objects, function(group) lapply(group, r_filter))


pdf("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/SFibroids/Venn.pdf",
    height = 3, width = 3)
# Secreotry
for (i in seq_len(length(significant[[1]]))) {
  x<-list(Fibroids=significant[[1]][i][[1]], Endometriosis=significant[[3]][i][[1]])
  gplot<-ggVennDiagram(x)+ 
    scale_fill_gradient(low = "#F4FAFE", high = "#F4FAFE") + ggtitle(names_celltypes_sec[i])
  print(gplot)
  }

# Proliferative

for (i in seq_len(length(significant[[1]]))) {
  x<-list(Fibroids=significant[[2]][i][[1]], Endometriosis=significant[[4]][i][[1]])
  gplot<-ggVennDiagram(x)+ 
    scale_fill_gradient(low = "#F4FAFE", high = "#F4FAFE") + ggtitle(names_celltypes_pro[i])
  print(gplot)
}
dev.off()
