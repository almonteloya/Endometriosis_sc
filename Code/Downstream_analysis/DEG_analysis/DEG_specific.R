## DEG all stromal and muscle

library(ggplot2)
library(Seurat)
library(dittoSeq)
library(ggrepel)
library(dplyr)


DEG_cells<-function(cell_type,mphase){
  Idents(S_object) <- S_object$annotation
  sobjs_clust <- subset(S_object,idents= cell_type)
  ## Dividing it by phase cycle 
  Idents(sobjs_clust)<-sobjs_clust$Phase_menstrual_reduced
  sobjs_clust <- subset(sobjs_clust,idents= mphase)
  ## Checking we are running the correct subset
  print(unique(sobjs_clust$labels_singleR))
  print(unique(sobjs_clust$Phase_menstrual_reduced))
  ## DEG by condition
  Idents(sobjs_clust)<-sobjs_clust$condition_control
  DEG<-FindMarkers(sobjs_clust ,
                   test.use="MAST",
                   logfc.threshold =0,
                   latent.vars= c("LIBRARY","Patient_id"),
                   ident.1 = "disease", ident.2 = "Pristine_control")
  
  saveRDS(DEG,paste0(dir_out,cell_type,"_",mphase,"_DEG.RDS"))
}

main_deg<-function(S_object,dir_out){
  
  s_deg<-function(cellt){
    print(cellt)
    sobjs_clust <- subset(S_object,idents= cellt)
    x <- table(sobjs_clust$condition_control)
    if (length(x)==3 & x[3]>3){
      DEG_plot(sobjs_clust,cellt)
    }
  }
  
  # Filter Mitocondrial
  S_object <- S_object[!grepl("^MT-", rownames(S_object)), ]
  Idents(S_object)<- S_object$condition_full
  Idents(S_object) <- S_object$annotation
  celltypes <- unique(Idents(S_object))
  phase<-(unique(S_object$Phase_menstrual_reduced))
  
  for (i in seq_along(phase)) {
    for (j in seq_along(celltypes)) {
      ## Check for enough cells
      x<-table(S_object$condition_control,S_object$annotation,S_object$Phase_menstrual_reduced)
      x<-as.data.frame(x)
      Fre<-filter(x,Var3==phase[i] & Var2==celltypes[j] & Var1 %in% c("disease","Pristine_control"))$Freq
      print(Fre)
      if (all(Fre  > 4)){
        DEG_cells(celltypes[[j]],phase[[i]]) 
      }
    }
  }
}

## Epihtelial
dir_in <- 'path/'
dir_out<-'path/'
S_object <- readRDS(sprintf("%s/Epithelial_obj_2025.RDS",dir_in))


main_deg(S_object,dir_out)

## Stromal
dir_in <- 'path/'
dir_out<-'path/'
S_object <- readRDS(sprintf("%s/Stromal_obj_2025.RDS",dir_in))

main_deg(S_object,dir_out)


