### DEG immune cells proliferative

library(dplyr)
library(ggplot2)
library(Seurat)


DEG_cells<-function(cell_type,mphase,S_object_sub){
  Idents(S_object_sub) <-S_object_sub[["NK_annotation"]]
  S_object_sub <- subset(S_object_sub,idents= cell_type)
  print(paste("Processing",unique(Idents(S_object_sub))))
  Idents(S_object_sub)<-S_object_sub[['Phase_menstrual_reduced']]
  S_object_sub <- subset(S_object_sub,idents= mphase)
  print(paste("Processing",unique(Idents(S_object_sub))))
  Idents(S_object_sub)<-S_object_sub[['condition_full_reduced']]
  ## this time woth the standard threshold
  DEG<-FindMarkers(S_object_sub ,
                   test.use="MAST",
                   logfc.threshold = 0,
                   latent.vars= c("LIBRARY","Patient_id"),
                   ident.1 = "Disease", ident.2 = "Pristine_control")
  saveRDS(DEG,paste0(dir_out,cell_type,"_",mphase,"_DEG.RDS"))
  print("DONE***")
  DEG
}


main_deg<-function(S_object,dir_out){
  # Filter Mitocondrial
  S_object <- S_object[!grepl("^MT-", rownames(S_object)), ]
  Idents(S_object) <- S_object[["NK_annotation"]]
  celltypes <- unique(Idents(S_object))
  Idents(S_object) <- S_object[["Phase_menstrual_reduced"]]
  phase<- unique(Idents(S_object))
  for (i in seq_along(phase)) {
    for (j in seq_along(celltypes)) {
      ## Check for enough cells
      x<-table(S_object$condition_full_reduced,S_object$NK_annotation,
               S_object$Phase_menstrual_reduced)
      x<-as.data.frame(x)
      Fre<-filter(x,Var3==phase[i] & Var2==celltypes[j] & Var1 %in% c("Disease","Pristine_control"))$Freq
      print(Fre)
      if (all(Fre  > 4)){
        DEG_cells(celltypes[[j]],phase[[i]],S_object) 
      }
    }
  }
}




dir2='path/'
dir_out <- 'path/'
dir.create(dir_out)
S_object = readRDS(file=paste0(dir2,"Immune_obj_2025.RDS"))
main_deg(S_object,dir_out)

