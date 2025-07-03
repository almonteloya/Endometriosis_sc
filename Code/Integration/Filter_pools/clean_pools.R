## Go into each pool, get the markers, print the cluster heatmap 
## and the count filter if its junk 
## Pool 5 and pool 10 were filtered differently

pools<-c(1,2,3,4,6,7,8,9,11,12,13)

for (x in pools) {
  print(paste("Processing pool",x))
  
  library(Seurat)
  library(tidyverse)
  library(plyr)
  library(cowplot)
  
  project="XREP1a"
  sample=paste0('XREP1-POOL-SCG',x) ## Pool prefix
  in_dir="insert/path/here/" ## Pathway for all seurat objects
  load(file.path(in_dir, paste0(sample,'.RData')))
  sobj <- get(sample)
  rm(list=sample)
  ##--- Finding which clusters are junk cells
  
  sobj.markers <- FindAllMarkers(sobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  #saveRDS(sobj.markers,file=paste0(in_dir,sample,"_markers.RDS"),compress = FALSE)
  
  ##. We look for the top 5 genes
  sobj.markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC) -> top5
  
  ## Here we can do a heatmap of the most important genes
  sonj_by_ident = subset(x = sobj, idents = c(1,2,3,4,5,6))
  pdf(paste0(in_dir,sample,"_heatmap_violin_markers.pdf"))
 print(DoHeatmap(sonj_by_ident, features = top5$gene[top5$cluster %in% c(1,2,3,4,5,6)],raster=FALSE) + NoLegend())
  dev.off()
  
  ## check ncount
pdf(paste0(in_dir,sample,"_violin.pdf"))
Idents(object = sobj) <- sobj@meta.data$seurat_clusters
print(VlnPlot(sobj,features='nCount_RNA',pt.size=0))
dev.off()
  
  
  ### Manually looking at the plots there are some pools who have a cluster of likely dead cells
  ## Pool 3 cluster 1
  ## Pool 6 cluster 1
  ## Pool 7 cluster 1
  ## Pool 11 cluster 6
  
  ## But the rest of them we are also converting them from Rdata to Rbject
  no_filter<-c(1,2,4,8,9,12,13)
  
  if (x %in% no_filter){
    saveRDS(sobj, file = paste0(in_dir,sample,".RDS"))
    pdf(paste0(in_dir,sample,"umap_filtered.pdf"))
    print(DimPlot(sobj,reduction = "umap",label = TRUE))
    dev.off()
  }
  
  if (x %in% c(3,6,7)){
    #number of clusters
    nclusters=length(levels(sobj$seurat_clusters))
    sobj<-subset(sobj, idents = c(2:nclusters))
    sobj <- RunUMAP(sobj, dims = 1:25)
    saveRDS(sobj, file = paste0(in_dir,sample,".RDS"))
    pdf(paste0(in_dir,sample,"umap_filtered.pdf"))
    print(DimPlot(sobj,reduction = "umap",label = TRUE))
    dev.off()
  }
  
  if (x == 11){
    #number of clusters
    all_no6<-levels(sobj$seurat_clusters)[levels(sobj$seurat_clusters)!=6]
    sobj<-subset(sobj, idents = all_no6)
    sobj <- RunUMAP(sobj, dims = 1:25)
    saveRDS(sobj, file = paste0(in_dir,sample,".RDS"))
    pdf(paste0(in_dir,sample,"umap_filtered.pdf"))
    print(DimPlot(sobj,reduction = "umap",label = TRUE))
    dev.off()
  }
  ### 

  rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
  gc() #free up memrory and report the memory usage.a
}


