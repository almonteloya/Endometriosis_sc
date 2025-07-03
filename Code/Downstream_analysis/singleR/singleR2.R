## Adding the labels from Wanxin data to ours
## First need to create object (SingleR1)
library(Seurat)
library(dplyr)

dir='/dir/in'
dir_single = paste0(dir,"single_R/") ## Input for seurat object
dir_out="dir/out/"


merged_data <-  readRDS(file=paste0(dir,"merged_data_celltype_filter.RDS"))
wanxin_annoation <- readRDS(paste0(dir_single,"Wanxin.RDS")) ## Produced in previous step
DimPlot(wanxin_annoation)
wanxin_annoation<- AddMetaData(wanxin_annoation,Idents(wanxin_annoation),"labels_singleR")


## We want broad annotation
wanxin_annoation$labels_singleR <- recode(wanxin_annoation$labels_singleR,
                                          "Lymphocytes" = "Immune", "Macrophages" = "Immune", 
                                          "Endothelia" = "Endothelia", "Stromal fibroblasts" ="Stromal_fibroblast",
                                          "Smooth muscle cells" = "VSmooth_muscle",
                                          "Unciliated epithelia 1" = "Glandular_Epithelia", "Unciliated epithelia 2" = "Glandular_Epithelia")



Idents(wanxin_annoation)<-wanxin_annoation$labels_singleR
DimPlot(wanxin_annoation,label=TRUE)

saveRDS(wanxin_annoation,paste0(dir,"merged_data_celltype_annotation.RDS"))

