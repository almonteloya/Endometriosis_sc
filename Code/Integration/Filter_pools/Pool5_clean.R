### Cleaning pool 5
## We are going to filter everything that overlaps with cluster 4,
## which according to the heatmap is from a different sample
library(Seurat)
library(dplyr)
project="XREP1a"
x=5
print("LOADING SAMPLE")
sample=paste0('XREP1-POOL-SCG',x)
in_dir= "/path/here/"
load(file.path(in_dir, paste0(sample,'.RData')))
sobj <- get(sample)
rm(list=sample)

print("CLEANING FREEMUXLET")
pdf(paste0(in_dir,sample,'clean_freemuxlet.pdf'))
print(DimPlot(sobj,label=TRUE))

Idents(sobj)<-sobj@meta.data$SAMPLE.by.SNPs
print(Seurat::DimPlot(sobj, reduction = "umap",label = TRUE))

table (sobj@meta.data$SAMPLE.by.SNPs)

## Cluster 4 accroding to freemuxlet is the "other sample"
sobj = subset(x = sobj, idents = c(0,1,2,3))
print(Seurat::DimPlot(sobj, reduction = "umap",label = TRUE))
sobj <- RunUMAP(sobj,
                      dims = 1:30,  # Num PCs to use
                      n.neighbors = 30,  # Default. Controls how UMAP balances local (low) versus global (large) structure in the data
                      min.dist = 0.3,   # Default. Controls the size of the clusters. Should be smaller than spread
                      spread = 1,  # Default. Controls the inter-cluster distances to some extent. Should be larger than min_dist
                      a = NULL,  # Default. Can be used with b instead of using min.dist/spread
                      b = NULL,  # Default. Can be used with a instead of using min.dist/spread
                      verbose = FALSE)
print(Seurat::DimPlot(sobj, reduction = "umap",label = TRUE))
# Calculate the neighborhood graph
sobj <- FindNeighbors(sobj,
                            dims = 1:30,  # Num PCs to use
                            k.param = 20,  # k for the knn algorithm
                            verbose = FALSE)
# Use the neighborhood graph to cluster the data
sobj <- FindClusters(sobj, verbose = TRUE,
                           algorithm = 4)  # Use Leiden

Idents(sobj)<-sobj@meta.data$seurat_clusters

print(Seurat::DimPlot(sobj, reduction = "umap",label = TRUE))
dev.off()

saveRDS(sobj,file=paste0(in_dir,sample,'.RDS'))

sobj.markers <- FindAllMarkers(sobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
##. We look for the top  genes
sobj.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
## Here we can do a heatmap of the most important genes
sonj_by_ident = subset(x = sobj, idents = c(1,2,3,4,5,6))
pdf(paste0(in_dir,sample,"_heatmap_violin_markers.pdf"))
print(DoHeatmap(sonj_by_ident, features = top5$gene[top5$cluster %in% c(1,2,3,4,5,6)],raster=FALSE) + NoLegend())
dev.off()


### Looks fine, no further filtering 
