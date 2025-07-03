# re-cluster after filtering immune cells
## Annotate based in literature markers

library(ggplot2)
library(Seurat)
library(dittoSeq)
library(dplyr)
library(harmony)
dir='path/here'
dir_out='/path/here'
dir.create(dir_out) 

## S.F and VSM no immune
sub_obj <- readRDS(paste0(dir_out,"seurat_obj_no_immune.RDS")) 

subset_plot_markers<- function(s_object){
  ##Cluster again
  ## Save object 
  ## FindMarkers and annotate
  diet_sub <- DietSeurat(s_object)
  diet_sub <- NormalizeData(diet_sub)
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
  diet_sub <- FindClusters(diet_sub, resolution = 0.2)
  pdf(paste0(dir_out,"Umap_filtered.pdf"))
  print(DimPlot(diet_sub,label = TRUE) + ggtitle("Stromal"))
  dev.off()
  saveRDS(diet_sub,paste0(dir_out,"seurat_filtered_obj.RDS"))
  
  markers <- FindAllMarkers(diet_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
  saveRDS(markers, paste0(dir_out,"filtered_markers.RDS"))
  markers %>%
    group_by(cluster) %>%
    top_n(n = 8, wt = avg_log2FC) -> top8
  pdf(paste0(dir_out,"filtered_markers.pdf"),width=18, height = 10)
  
  print(DotPlot(diet_sub, features=unique(top8$gene), dot.scale=6, cols="RdBu") +
          theme(axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5, hjust=1)) +
          ggtitle(.3))
  
  dev.off()
}
subset_plot_markers(sub_obj)



#### Adding the cell names for each cluster
dir_out='/krummellab/data1/immunox/XREP1a/10x/merged_XREP/subset_celltype/Mesenchymal/'

obj <- readRDS(paste0(dir_out,"seurat_filtered_obj.RDS"))
DimPlot(obj)
new.cluster.ids <- c("Activated-fibroblast","IGF1+-fibroblast","HSP-stromal-fibroblast",
                     "SFRP4+-stromal-fibroblast","MMP+-fibroblast","Proliferating-fibroblast",
                     "vascular-smooth-muscle","WFC2-Stromal-fibroblast")
names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)
DimPlot(obj, reduction = "umap", pt.size = 0.5)
obj@meta.data$annotation <- Idents(obj)
saveRDS(obj@meta.data, paste0(dir_out,"metadata_annotation.RDS"))
