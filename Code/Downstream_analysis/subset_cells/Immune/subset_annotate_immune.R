## subsetiing and annotating immune cells

library(Seurat)
library(tidyverse)

dir='path/here/'
dir_out = "path/here"
dir.create(dir_out)

merged_data <- readRDS(paste0(dir,"merged_data.RDS"))
                       
                      
### subset by immune 

# --------- Lymphatic cells cluster close to immune subtypes 

cell_names <- colnames(subset(merged_data,louvain_res0.4 ==18))
Idents(merged_data, cells = cell_names) <- "immune"
merged_data$cell_type_general <- Idents(merged_data)



Idents(merged_data)<-merged_data$cell_type_general

cell_type<-"immune"
subset_plot_markers<- function(cell_type){
  print(paste("Processing", cell_type))
  subset_cell <- subset(merged_data, idents = cell_type )
  diet_sub <- DietSeurat(subset_cell)
  diet_sub = NormalizeData(diet_sub)
  diet_sub <- FindVariableFeatures(diet_sub, selection.method = "vst", nfeatures = 3000)
  diet_sub <- ScaleData(diet_sub)
  diet_sub <- RunPCA(diet_sub)
  
  diet_sub <- RunUMAP(diet_sub,
                      dims = 1:30,  # Num PCs to use
                      n.neighbors = 30,  # Default. Controls how UMAP balances local (low) versus global (large) structure in the data
                      min.dist = 0.3,   # Default. Controls the size of the clusters. Should be smaller than spread
                      spread = 1,  # Default. Controls the inter-cluster distances to some extent. Should be larger than min_dist
                      a = NULL,  # Default. Can be used with b instead of using min.dist/spread
                      b = NULL,  # Default. Can be used with a instead of using min.dist/spread
                      verbose = FALSE)
  
  diet_sub <- FindNeighbors(diet_sub, dims = 1:20)
  diet_sub <- FindClusters(diet_sub, resolution = 0.4)
  pdf(paste0(dir_out,cell_type,"_2_umap.pdf"))
  print(DimPlot(diet_sub,label = TRUE) + ggtitle(cell_type))
  dev.off()
  saveRDS(diet_sub,paste0(dir_out,cell_type,"_2_obj.RDS"))
  
  markers <- FindAllMarkers(diet_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  saveRDS(markers, paste0(dir_out,cell_type,"_2_markers.RDS"))
  markers %>%
    group_by(cluster) %>%
    top_n(n = 8, wt = avg_log2FC) -> top8
  pdf(paste0(dir_out,cell_type,"_2_markers.pdf"),width=18, height = 10)
  
  print(DotPlot(diet_sub, features=unique(top8$gene), dot.scale=6, cols="RdBu") +
          theme(axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5, hjust=1)) +
          ggtitle(cell_type))
  
  dev.off()
}
subset_plot_markers("immune")



immune_object<-  readRDS(paste0(dir_out,cell_type,"_2_obj.RDS"))
DimPlot(immune_object)


new.cluster.ids <- c("Stressed_Tcells","NK_cells","CD8+_Tcells ",
"Macrophages","Stressed_Tcells","NK_cells","NK_cells",
"Monocytes","NK_cells ",
"Monocytes","Stressed_Tcells", "Lymphatic",
"NK_cells","T_reg","Macrophage","Innate_lymphoid","Plasma_cells",
"MAST","cDC1","PD C")
names(new.cluster.ids) <- levels(immune_object)
immune_object <- RenameIdents(immune_object, new.cluster.ids)


saveRDS(immune_object, paste0(dir,"immune_annotation.RDS"))

