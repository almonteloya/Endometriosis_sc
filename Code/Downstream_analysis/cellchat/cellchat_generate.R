## Cellchat 
## Process and generate cellchat objects to compare 

library(CellChat)
library(Seurat)

merged_data_sub<-readRDS("merged_data_cellchat.RDS") ##Object with non-immune and immnune celltypes
dir_out<-("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/cellchat/")

## some of the Celltypes are not present in all of the data sets, it's easier to remove
Idents(merged_data_sub)<-merged_data_sub$immune_annotation
merged_data_sub <- subset(merged_data_sub, idents= c("cDC1","MAST"), invert=TRUE)


CellChatDB <- CellChatDB.human # loading the data
showDatabaseCategory(CellChatDB)
## we can explore what the database has
CellChatDB.use <- CellChatDB # simply use the default CellChatDB


cell_chat_obj <- function(seurat_obj,name_obj){
  ### For each seurat obj create a cellchat object
  data.input <- GetAssayData(seurat_obj, assay = "RNA", slot = "data") # normalized data matrix
  Idents(seurat_obj)<-seurat_obj$immune_annotation
  labels <- Idents(seurat_obj)
  meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
  
  ## Adding metadata information
  meta$condition<-seurat_obj@meta.data$Condition
  meta$patient_id<-seurat_obj@meta.data$Patient_id
  meta$phase<-seurat_obj@meta.data$Phase_menstrual_reduced
  cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")
  #cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
  levels(cellchat@idents) # show factor levels of the cell labels
  
  
  # set the used database in the object
  cellchat@DB <- CellChatDB.use
  
  
  # Preprocessing the expression data for cell-cell communication analysis
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  ## Compute the communication probability and infer cellular communication network
  cellchat <- computeCommunProb(cellchat)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  saveRDS(cellchat,paste0(dir_out, name_obj,"_cellchat_obj.RDS"))
}



cellchat_object <- function(){
  ## Stratified by menstrual phase and disease condition 
  disease_obj<-subset(merged_data_sub, subset = Condition =="disease")
  disease_sec<-subset(disease_obj, subset = Phase_menstrual_reduced =="secretory")
  cell_chat_obj(disease_sec,  "disease_secretory")
  
  disease_pro <-subset(disease_obj, subset = Phase_menstrual_reduced=="proliferative")
  cell_chat_obj(disease_pro,  "disease_proliferative")
  rm(disease_obj)
  
  contr_obj<-subset(merged_data_sub, subset = condition_control=="Pristine_control")
  rm(merged_data_sub)
  control_sec <- subset(contr_obj, subset = Phase_menstrual_reduced=="secretory")
  cell_chat_obj(control_sec,  "control_sec")
  control_pro <-subset(contr_obj, subset = Phase_menstrual_reduced=="proliferative")
  cell_chat_obj(control_pro,  "control_pro")
  rm(contr_obj)
}

cellchat_object()

