library(dittoSeq)
library(Seurat)
library(harmony)
library(dplyr)

##### Since we saw a lot of overlap in clusters 0,1,and 5
## We are going to subset each cluster and look for low quality cells

subset_plot_markers<- function(cell_type,res,nn){
  ## Subset by clusters
  ## Re-cluster
  print(paste("Processing", cell_type))
  subset_cell <- subset(epi_data, idents = cell_type )
  diet_sub <- DietSeurat(subset_cell)
  diet_sub = NormalizeData(diet_sub)
  diet_sub <- FindVariableFeatures(diet_sub, selection.method = "vst", nfeatures = 3000)
  diet_sub <- ScaleData(diet_sub)
  diet_sub <- RunPCA(diet_sub,features = VariableFeatures(object = diet_sub))
  pdf(paste0(dir_out, 'harmony_convergence.pdf'))

  print(ElbowPlot(diet_sub))
  diet_sub <- RunHarmony(diet_sub,
                         "LIBRARY",
                         assay.use='RNA',
                         dims.use = 1:15,
                         plot_convergence = TRUE,
                         max.iter.harmony=20,
                         max.iter.cluster=30)
  dev.off()


  diet_sub <- RunUMAP(diet_sub,
                      dims = 1:15,  # Num PCs to use
                      reduction='pca',
                      n.neighbors = nn,  # Default. Controls how UMAP balances local (low) versus global (large) structure in the data
                      min.dist = 0.2,   # Default. Controls the size of the clusters. Should be smaller than spread
                      spread = .8,  # Default. Controls the inter-cluster distances to some extent. Should be larger than min_dist
                      a = NULL,  # Default. Can be used with b instead of using min.dist/spread
                      b = NULL,  # Default. Can be used with a instead of using min.dist/spread
                      verbose = FALSE)

  diet_sub <- FindNeighbors(diet_sub, dims = 1:15,   prune.SNN = 1/10) ## k.param = 20
  diet_sub <- FindClusters(diet_sub, resolution = .15)

  pdf(paste0(dir_out,"umap_",res,"_",nn,".pdf"))
  print(DimPlot(diet_sub,label = TRUE) + ggtitle(cell_type))
  print(density_plotter(diet_sub, "seurat_clusters"))

  dev.off()
  saveRDS(diet_sub,paste0(dir_out,cell_type,"_01.obj.RDS"))

  markers <- FindAllMarkers(diet_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
  saveRDS(markers, paste0(dir_out,cell_type,"_markers_0.1.RDS"))
  markers %>%
    group_by(cluster) %>%
    top_n(n = 8, wt = avg_log2FC) -> top8
  pdf(paste0(dir_out,cell_type,"_markers_0.1.pdf"),width=18, height = 10)

  print(DotPlot(diet_sub, features=unique(top8$gene), dot.scale=6, cols="RdBu") +
          theme(axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5, hjust=1)) +
          ggtitle(cell_type))

  dev.off()
}


cell_type = c(0,1,5)

for (i in cell_type) {
  subset_plot_markers(i, 0.05 ,25)
}




########## Now let's explore each subclusters

dir_out='path/here'

files <- list.files( path=dir_out, pattern = "obj.RDS", full.names = TRUE)
obj <- lapply(files, readRDS) ## Each one has a subcluster of cluster 0,1,5


plotUMAP <- function(s_obj){
  print(DimPlot(s_obj))
  print(density_plotter(s_obj, "seurat_clusters"))
}

sapply(obj, plotUMAP)

###-------- To clean things up I'm removing some cells, with high mitochondrial content

## For the cluster 0, subcluster 0 has high MT
mt_cells_0 <- names(obj[[1]]$seurat_clusters[obj[[1]]$seurat_clusters==0])

## For the cluster 1, subcluster sucluster 0 has high MT
mt_cells_1 <- names(obj[[2]]$seurat_clusters[obj[[2]]$seurat_clusters==0])
mt_cells<- c(mt_cells_0,mt_cells_1)

## Cluster 5 : all subclusters look OK. We are not filtering 

### ----- Now let's filter those cells 

epi_data <- readRDS(paste0(dir_in,"seurat_obj_no_immune.RDS"))

epi_data <- subset(epi_data, subset = nCount_RNA < 30000 & percent.mt < 15 )
epi_data <- epi_data[,!colnames(epi_data) %in% mt_cells]
DimPlot(epi_data)





## Now we'll recluster again -----


subset_plot_markers<- function(nn){
  diet_sub <- DietSeurat(epi_data)
  diet_sub = NormalizeData(diet_sub)
  diet_sub <- FindVariableFeatures(diet_sub, selection.method = "vst", nfeatures = 3000)
  diet_sub <- ScaleData(diet_sub)
  diet_sub <- RunPCA(diet_sub,features = VariableFeatures(object = diet_sub))
  pdf(paste0(dir_out, 'harmony_convergence.pdf'))
  
  ElbowPlot(diet_sub)
  diet_sub <- RunHarmony(diet_sub,
                         "LIBRARY",
                         assay.use='RNA',
                         dims.use = 1:15,
                         plot_convergence = TRUE,
                         max.iter.harmony=20,
                         max.iter.cluster=30)
  dev.off()
  
  diet_sub <- FindNeighbors(diet_sub, dims = 1:15,   prune.SNN = 1/10) ## k.param = 20
  diet_sub <- FindClusters(diet_sub, resolution = .09)
  diet_sub <- RunUMAP(diet_sub,
                       dims = 1:15,  # Num PCs to use
                       reduction="harmony",
                       n.neighbors = 17,  # Default 20. Controls how UMAP balances local (low) versus global (large) structure in the data
                       min.dist = 0.2,   # Default. Controls the size of the clusters. Should be smaller than spread
                       spread = .8,  # Default. Controls the inter-cluster distances to some extent. Should be larger than min_dist
                       verbose = FALSE)
  saveRDS(diet_sub, paste0(dir_out,"epithelial_filter_obj.RDS"))
  
  pdf(paste0(dir_out,"umap_",".pdf"))
  print(DimPlot(diet_sub,label = TRUE))
  dev.off()
  markers <- FindAllMarkers(diet_sub, only.pos = TRUE, min.pct = 0.25)
  saveRDS(markers, paste0(dir_out,"epithelial_filter_markers.RDS"))
  markers %>%
    group_by(cluster) %>%
    top_n(n = 8, wt = avg_log2FC) -> top5
  pdf(paste0(dir_out,"epithelial_filter_markers.pdf"),width=18, height = 10)
  
  print(DotPlot(diet_sub, features=unique(top5$gene), dot.scale=6, cols="RdBu") +
          theme(axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5, hjust=1))) 
  
  dev.off()
  
  ##Explore overlap
  
  clus <- "seurat_clusters"
  ggsave(
    paste0(dir_out, clus,"epithelial_filter_hex.png"),
    dittoSeq::dittoDimHex(diet_sub, split.by = clus, color.panel = rep("black", 4000),
                main = clus)  +
      guides(fill=guide_colourbar(barheight=30)),
    units = "in",
    width = 15,
    height = 15
  )
  
}

subset_plot_markers(0.1)





#----- Now we can explore some markers
## Plus based on the literature we can name some clusters

VlnPlot(epi_data, features = c("MUC1","EPCAM"), ncol = 3,pt.size=0)

ciliated <- c("FOXJ1","PIFO","TP73")
FeaturePlot(epi_data,ciliated)

new.cluster.ids <- c("Glandular-1","EMT-epithelia-1","Glandular-2","EMT-epithelia-2","Ciliated-epithelia",
                     "Glandular-3","Glandular-4","Progenitor-MUC5B+TFF3+")
names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)
obj@meta.data$annotation <- Idents(obj)

saveRDS(obj@meta.data, paste0(dir_out,"metadata_annotation.RDS"))


