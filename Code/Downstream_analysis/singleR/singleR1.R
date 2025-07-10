
### Single R annotation
## Using Wanxin's et al as reference 
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111976

library(SingleR)
library(scRNAseq)
library(ggplot2)
library(Seurat)
library(dittoSeq)
library(dplyr)
library(scuttle)

dir1='/merged_XREP/'
dir2="matrix_singleR/" ## Matrix for downloaded data
dir_out='single_R/'


wanxin_data<-function(){
  #### Get Wanxin data and convert to appropiate format
wan_data<-readRDS("/path/GSE111976_ct_endo_10x")
## This is only the matrix, we need to add extra information and convert it
## This are the labels
## they are in the same order
summary_wan_data<- read.csv('/GSE111976_summary_10x_day_donor_ctype.csv')
summary_wan_data<- summary_wan_data$cell_type
## we convert it 
wan_data_sce <- SingleCellExperiment(list(counts=wan_data),
                                     colData=DataFrame(label=summary_wan_data))
wan_data_sce <- wan_data_sce[,!is.na(wan_data_sce$label)]
# SingleR() expects reference datasets to be normalized and log-transformed.
wan_data_sce <- logNormCounts(wan_data_sce)
wan_data_sce
}


conver_pred <- function(seurat_data,ref_data,data_name){
  ### convert the data to a single cell experiment 
  ## Assing labels using SingleR
  
  ##-- Seurat Data: our data
  ##-- Ref data: Wanxin's data
  #-- Data name: name for our output
  ## output:Our data with the new lables
  Seurat_Object_Diet <- DietSeurat(seurat_data, graphs = "umap")
  sce_ct <- as.SingleCellExperiment(Seurat_Object_Diet, assay="RNA")
  colData(sce_ct)
  sce_ct <- sce_ct[,colSums(counts(sce_ct)) > 0] # Remove libraries with no counts.
  sce_ct <- logNormCounts(sce_ct) 
  sce_ct  
  pred.grun <- SingleR(test=sce_ct, ref=ref_data, 
                       labels=ref_data$label, de.method="wilcox",
                       clusters=sce_ct$seurat_clusters)
  print(table(pred.grun$labels))
  print(plotScoreHeatmap(pred.grun))
  Idents(seurat_data) <- seurat_data$seurat_clusters
  names(pred.grun$first.labels) <- levels(seurat_data)
  seurat_data <- RenameIdents(seurat_data, pred.grun$first.labels)
  print(DimPlot(seurat_data, reduction = "umap", label = TRUE, pt.size = 0.5) )
  saveRDS(seurat_data,paste0(dir_out,data_name,".RDS"))
} 

## This is our data 
merged_data <-  readRDS(file=paste0(dir1,"merged_data_celltype_res8_filtered_3_7.RDS"))
##Creating a refrence data
reference_data <- wanxin_data()
conver_pred(merged_data,reference_data,"Wanxin")
#### we can also see by sub types
cell_type_single <- function(cell_type){
  print(paste("processing",cell_type))
  sobjs_clust <- readRDS(paste0(dir1,"subset_celltype/", cell_type, "_obj.RDS"))
  conver_pred(sobjs_clust,reference_data,paste0(cell_type , "Wanxin"))
}
sapply(c("stromal","immune","endothelium","epithelium"), cell_type_single)
dev.off() 
