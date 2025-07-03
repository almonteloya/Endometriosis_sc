
### Filtering on cluster 7 and 3 after merging

library(Seurat)
library(dplyr)
library(ggplot2)

dir='/insert/path/in/'
dir_out="/insert/path/in/"

merged_data <- readRDS(paste0(dir,"merged_data_celltype.RDS"))


Idents(merged_data)<-merged_data@meta.data$louvain_res0.8
sobj.markers <- FindAllMarkers(merged_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

sobj.markers <- readRDS("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/lovain8_markers.RDS")#

sobj.markers %>%
  group_by(cluster) %>%
  top_n(n = 8, wt = avg_log2FC) -> top5

#### look at the markers
Idents(merged_data)<-merged_data@meta.data$louvain_res0.8



Idents(merged_data)<-merged_data$louvain_res0.8
VlnPlot(merged_data, features="nCount_RNA",pt.size=0.0)
VlnPlot(merged_data, features="nFeature_RNA",pt.size=0.0)
VlnPlot(merged_data, features="percent.mt",pt.size=0.0)

#### looking at cluster 3 must genes seem to be mithocondrial , so I'll filter it
#### looking at cluster 7 looks like low count so I'll subset it 

merged_filter <-  subset(merged_data,idents=c(3,7),invert=TRUE)

DimPlot(merged_filter,label=TRUE)#

## the general rule: if the cluster is either (1) in the center
## or (2) more than 10% of the cells --> re-run PC, NN, cluster, UMAP#

merged_filter <- RunPCA(merged_filter)
merged_filter <- RunUMAP(merged_filter,
                         dims = 1:30,  # Num PCs to use
                         n.neighbors = 30,  # Default. Controls how UMAP balances local (low) versus global (large) structure in the data
                         min.dist = 0.3,   # Default. Controls the size of the clusters. Should be smaller than spread
                         spread = 1,  # Default. Controls the inter-cluster distances to some extent. Should be larger than min_dist
                         a = NULL,  # Default. Can be used with b instead of using min.dist/spread
                         b = NULL,  # Default. Can be used with a instead of using min.dist/spread
                         verbose = FALSE)#

merged_filter <- FindNeighbors(merged_filter, dims = 1:20)
merged_filter <- FindClusters(merged_filter)

saveRDS(merged_filter,paste0(dir,"merged_data_celltype_filter.RDS"))
