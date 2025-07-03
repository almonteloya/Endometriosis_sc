## DEG all stromal, epi and endothelia

library(ggplot2)
library(Seurat)
library(dittoSeq)
library(ggrepel)
library(dplyr)


DEG_cells<-function(cell_type,mphase,mseverity){
  Idents(S_object) <- S_object$labels_singleR_reduced
  sobjs_clust <- subset(S_object,idents= cell_type)
  ## Dividing it by phase cycle 
  Idents(sobjs_clust)<-sobjs_clust$Phase_menstrual_reduced
  sobjs_clust <- subset(sobjs_clust,idents= mphase)
  ## Checking we are running the correct subset
  print("Calculating:")
  print(unique(sobjs_clust$labels_singleR))
  print(unique(sobjs_clust$Phase_menstrual_reduced))
  print(mseverity)
  ## DEG by condition
  Idents(sobjs_clust)<-sobjs_clust$condition_full
  DEG<-FindMarkers(sobjs_clust ,
                   test.use="MAST",
                   logfc.threshold =0,
                   latent.vars= c("LIBRARY","Patient_id"),
                   ident.1 = mseverity, ident.2 = "Pristine_control")
  
  saveRDS(DEG,paste0(dir_out,cell_type,"_",mphase,"_",mseverity,"_DEG.RDS"))
}

main_deg<-function(S_object,dir_out){
  # Filter Mitocondrial
  S_object <- S_object[!grepl("^MT-", rownames(S_object)), ]
  Idents(S_object) <- S_object$labels_singleR_reduced
  celltypes <- unique(Idents(S_object))
  phase<-(unique(S_object$Phase_menstrual_reduced))
  severity<-c("I-II","III-IV")
  for (i in seq_along(phase)) {
    for (j in seq_along(celltypes)) {
      for (h in seq_along(severity)) {
      DEG_cells(celltypes[[j]],phase[[i]],severity[[h]])  
    }
  }
}
}

dir_in <- 'path/'
dir_out<-'path/'
## Epithelial
S_object <- readRDS(sprintf("%s/Epithelial_obj_2025.RDS",dir_in))
main_deg(S_object,dir_out)



## Stromal 
S_object <- readRDS(sprintf("%s/Stromal_obj_2025.RDS",dir_in))
main_deg(S_object,dir_out)


# Endothelia
S_object <- readRDS(sprintf("%s/merged_data_celltype_all_2025.RDS",dir_in))
Idents(S_object)<-S_object$labels_singleR_reduced
S_endothelia <- subset(S_object, idents="Endothelia")
main_deg(S_endothelia,dir_out)
