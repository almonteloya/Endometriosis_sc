## Subsetting Stromal population
## Cleaning possible doublets

library(ggplot2)
library(Seurat)
library(dittoSeq)
library(dplyr)
library(harmony)


dir='path/here'
dir_out='/path/here'
dir.create(dir_out) 


## Reading the object with all of the cells
merged_data = readRDS(file=paste0(dir,"merged_data_celltype_all.RDS"))

## Labels from general annotation
Idents(merged_data)<-merged_data$labels_singleR

## We are interested in this two celltypes

cell_type<-c("Stromal_fibroblast","Smooth_muscle")

subset_plot_markers<- function(cell_type){
  ## Subset and recluster
  ## Plot UMAP
  ## Save object
  ## Run FindMarkers
  print(paste("Processing", cell_type))
  subset_cell <- subset(merged_data, idents = cell_type )
  diet_sub <- DietSeurat(subset_cell)
  diet_sub = NormalizeData(diet_sub)
  diet_sub <- FindVariableFeatures(diet_sub, selection.method = "vst", nfeatures = 3000)
  diet_sub <- ScaleData(diet_sub)
  diet_sub <- RunPCA(diet_sub)
  pdf(paste0(dir_out, 'harmony_convergence.pdf'))
  
  ElbowPlot(diet_sub)
  diet_sub <- RunHarmony(diet_sub,
                         "LIBRARY",
                         assay.use='RNA',
                         plot_convergence = TRUE,
                         max.iter.harmony=20,
                         max.iter.cluster=30)
  dev.off()
  
  diet_sub <- RunUMAP(diet_sub,
                      reduction = "harmony",
                      dims = 1:30,  # Num PCs to use
                      n.neighbors = 30,  # Default. Controls how UMAP balances local (low) versus global (large) structure in the data
                      min.dist = 0.3,   # Default. Controls the size of the clusters. Should be smaller than spread
                      spread = 1,  # Default. Controls the inter-cluster distances to some extent. Should be larger than min_dist
                      a = NULL,  # Default. Can be used with b instead of using min.dist/spread
                      b = NULL,  # Default. Can be used with a instead of using min.dist/spread
                      verbose = FALSE)
  
  diet_sub <- FindNeighbors(diet_sub,
                               dims = 1:30,  # Num PCs to use
                               reduction='harmony',
                               k.param = 20,  # k for the knn algorithm
                               verbose = FALSE)
  diet_sub <- FindClusters(diet_sub, resolution = 0.3)
  pdf(paste0(dir_out,"Umap.pdf"))
  print(DimPlot(diet_sub,label = TRUE) + ggtitle("Mesenchymal"))
  dev.off()
  saveRDS(diet_sub,paste0(dir_out,"seurat_obj.RDS"))
  
  markers <- FindAllMarkers(diet_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
  saveRDS(markers, paste0(dir_out,"markers.RDS"))
  markers %>%
    group_by(cluster) %>%
    top_n(n = 8, wt = avg_log2FC) -> top8
  pdf(paste0(dir_out,"markers.pdf"),width=18, height = 10)
  
  print(DotPlot(diet_sub, features=unique(top8$gene), dot.scale=6, cols="RdBu") +
          theme(axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5, hjust=1)) +
          ggtitle(.3))
  
  dev.off()
}
subset_plot_markers(cell_type)

#### Exploring stromal cells output


sub_obj <- readRDS(paste0(dir_out,"seurat_obj.RDS"))

###---- exploring some populations ----------
explore_immune <- function(){
  ### I notice that there's one cluster which has PTRC so most likely immune doublets. 
  ## we are going to filter them
  VlnPlot(sub_obj, features = c("PTPRC"), ncol = 3,pt.size=0)
  sub_obj<-subset(sub_obj,seurat_clusters!=7)
  saveRDS(sub_obj,paste0(dir_out,"seurat_obj_no_immune.RDS"))
}


## Output: Stromal and VSM with no imunne subclusters
