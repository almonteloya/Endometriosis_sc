## DEG all stromal, epi and endothelia

library(ggplot2)
library(Seurat)
library(dittoSeq)
library(ggrepel)
library(dplyr)


DEG_cells<-function(cell_type,m_condition){
  Idents(S_object) <- S_object$labels_singleR_reduced
  sobjs_clust <- subset(S_object,idents= cell_type)
  ## Dividing it by phase cycle 
  Idents(sobjs_clust)<-sobjs_clust$condition_control
  sobjs_clust <- subset(sobjs_clust,idents= m_condition)
  ## Checking we are running the correct subset
  print(unique(sobjs_clust$labels_singleR_reduced))
  print(unique(sobjs_clust$condition_control))
  ## DEG by condition
  Idents(sobjs_clust)<-sobjs_clust$Phase_menstrual_reduced
  DEG<-FindMarkers(sobjs_clust ,
                   test.use="MAST",
                   logfc.threshold =0,
                   latent.vars= c("LIBRARY","Patient_id"),
                   ident.1 = "secretory", ident.2 = "proliferative")
  
  saveRDS(DEG,paste0(dir_out,cell_type,"_",m_condition,"_DEG.RDS"))
}

main_deg<-function(S_object){
  # Filter Mitocondrial
  S_object <- S_object[!grepl("^MT-", rownames(S_object)), ]
  Idents(S_object) <- S_object$labels_singleR_reduced
  celltypes <- unique(Idents(S_object))
  condition<-(unique(S_object$condition_control))
  for (i in seq_along(condition)) {
    for (j in seq_along(celltypes)) {
      DEG_cells(celltypes[[j]],condition[[i]])  
    }
  }
}


dir_out<-'path/'
dir_in <- 'path/'


## Epihtelial
S_object <- readRDS(sprintf("%s/Epithelial_obj_2025.RDS",dir_in))
main_deg(S_object)



## Stromal 
S_object <- readRDS(sprintf("%s/Stromal_obj_2025.RDS",dir_in))
main_deg(S_object)


# Endothelia
S_object <- readRDS(sprintf("%s/merged_data_celltype_all_2025.RDS",dir_in))
Idents(S_object)<-S_object$labels_singleR_reduced
S_endothelia <- subset(S_object, idents="Endothelia")
main_deg(S_endothelia)
