library(ggplot2)
library(Seurat)
library(dittoSeq)
library(dplyr)
library(harmony)
dir='/krummellab/data1/immunox/XREP1a/10x/merged_XREP/'
dir_out='/krummellab/data1/immunox/XREP1a/10x/merged_XREP/subset_celltype/epithelia/'
merged_data = readRDS(file=paste0(dir,"merged_data_celltype_all.RDS"))


Idents(merged_data)<-merged_data$labels_singleR_reduced

cell_type<-"Epithelia"

subset_plot_markers<- function(cell_type,res,nn){
  print(paste("Processing", cell_type))
  subset_cell <- subset(merged_data, idents = cell_type )
  diet_sub <- DietSeurat(subset_cell)
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

  
  diet_sub <- RunUMAP(diet_sub,
                      dims = 1:15,  # Num PCs to use
                      reduction='harmony',
                      n.neighbors = nn,  # Default. Controls how UMAP balances local (low) versus global (large) structure in the data
                      min.dist = 0.2,   # Default. Controls the size of the clusters. Should be smaller than spread
                      spread = .8,  # Default. Controls the inter-cluster distances to some extent. Should be larger than min_dist
                      a = NULL,  # Default. Can be used with b instead of using min.dist/spread
                      b = NULL,  # Default. Can be used with a instead of using min.dist/spread
                      verbose = FALSE)
  
  diet_sub <- FindNeighbors(diet_sub, dims = 1:15)
  diet_sub <- FindClusters(diet_sub, resolution = res)
  
pdf(paste0(dir_out,"umap_",res,"_",nn,".pdf"))
print(DimPlot(diet_sub,label = TRUE) + ggtitle(cell_type))
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


subset_plot_markers("Epithelia",0.1,25) ## using smaller resolution and neighbors: less cells


dir_out='/krummellab/data1/immunox/XREP1a/10x/merged_XREP/subset_celltype/epithelia/'
 epi_data = readRDS(file=paste0(dir_out,"Epithelia_01.obj.RDS"))
 
 sub_markers <- readRDS(paste0(dir_out,"Epithelia_markers_0.1.RDS"))
 sub_markers %>%
   group_by(cluster) %>%
   top_n(n = 10, wt = avg_log2FC) -> top8
 
 ###---- exploring some populations ----------
 explore_immune <- function(){
   ### I notice that there's one cluster which has PTRC so most likely immune doublets. 
   ## we are going to filter them
   VlnPlot(epi_data, features = c("PTPRC"), ncol = 3,pt.size=0)
   epi_data<-subset(epi_data,seurat_clusters!=6)
   saveRDS(epi_data,paste0(dir_out,"seurat_obj_no_immune.RDS"))
 }
 
 
 ## Output Seurat object with epithelial cells (ciliated and unciliated) after filtering non-immune
 

