#### Nov 17 2022
## When talking with alexis we realized that the clusters on NK cells were separated
## I'll run harmony again after the DNA

#### Dec 8 2022
## run with new object

## Dec 14 2022
## run with new object

## March 2023
## There's a subset of cells that are lymphatic cells according to wanxin


library(Seurat)
library(tidyverse)
library(dittoSeq)

dir='/krummellab/data1/immunox/XREP1a/10x/merged_XREP/'
dir_out = "/krummellab/data1/immunox/XREP1a/10x/merged_XREP/subset_celltype/"

## This is produced by markers_louvan8_2

merged_data <- readRDS(paste0(dir,"merged_data_celltype_res8_filtered_3_7.RDS"))
### Linda have us more info on this patient 
merged_data@meta.data <- merged_data@meta.data %>% 
  mutate(Condition = ifelse(Patient_id=='STF59',"control",Condition),
         condition_full = ifelse(Patient_id=='STF59',"control",condition_full))

##lymphatic cells are immune 
Idents(merged_data)<-merged_data$cell_type_general
cell_names <- colnames(subset(merged_data,louvain_res0.4 ==18))
Idents(merged_data, cells = cell_names) <- "immune"

## smooth muscle is not stromal
## based on wanxins
cell_names <- colnames(subset(merged_data,louvain_res0.4 ==11))
Idents(merged_data, cells = cell_names) <- "Smooth_muscle"

### subset by general type
merged_data$cell_type_general <-Idents(merged_data)
## for immune we do a smaller clustering
#general_cell_type<-c("epithelium","stromal","immune","endothelium")

### we should save the plot
pdf(paste0(dir,"/plots/cell_type_general.pdf"))
dittoDimPlot(merged_data,"cell_type_general",size=.3)
dev.off()

general_cell_type<-c("epithelium","stromal","endothelium")

subset_plot_markers<- function(cell_type){
  print(paste("Processing",cell_type))
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
  diet_sub <- FindClusters(diet_sub, resolution = 0.2)
  pdf(paste0(dir_out,cell_type,"_umap.pdf"))
  print(DimPlot(diet_sub) + ggtitle(cell_type))
  dev.off()
  saveRDS(diet_sub,paste0(dir_out,cell_type,"_obj.RDS"))

# markers<-readRDS(paste0(dir_out,cell_type,"_markers.RDS"))
  markers <- FindAllMarkers(diet_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  markers %>%
  group_by(cluster) %>%
  top_n(n = 8, wt = avg_log2FC) -> top8
  pdf(paste0(dir_out,cell_type,"_markers.pdf"),width=18, height = 10)
  
  print(DotPlot(diet_sub, features=unique(top8$gene), dot.scale=6, cols="RdBu") +
        theme(axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5, hjust=1)) +
        ggtitle(cell_type))
  
  dev.off()
}

sapply(general_cell_type,subset_plot_markers)
