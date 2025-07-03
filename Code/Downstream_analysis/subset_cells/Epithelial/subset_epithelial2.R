
## Re-clustering epithelia after filtering immune subset
## Explore cluster outputs

library(dittoSeq)
library(Seurat)
library(harmony)
library(dplyr)

dir_out='/path/here/'
dir_in ='path/here'
epi_data <- readRDS(paste0(dir_in,"seurat_obj_no_immune.RDS"))


### Re cluster and process
subset_plot_markers<- function(res,nn){
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
  
  
  diet_sub <- RunUMAP(diet_sub,
                      dims = 1:15,  # Num PCs to use
                      reduction='harmony',
                      n.neighbors = nn,  # Controls how UMAP balances local (low) versus global (large) structure in the data
                      min.dist = 0.2,   # Default. Controls the size of the clusters. Should be smaller than spread
                      spread = .8,  # Default. Controls the inter-cluster distances to some extent. Should be larger than min_dist
                      a = NULL,  # Default. Can be used with b instead of using min.dist/spread
                      b = NULL,  # Default. Can be used with a instead of using min.dist/spread
                      verbose = FALSE)
  
  diet_sub <- FindNeighbors(diet_sub, dims = 1:15, prune.SNN = 1/10 )
  diet_sub <- FindClusters(diet_sub, resolution = res)
  
  pdf(paste0(dir_out,"umap_",res,"_",nn,".pdf"))
  print(DimPlot(diet_sub,label = TRUE))
  dev.off()
  saveRDS(diet_sub,paste0(dir_out,res,".obj.RDS"))
  
   markers <- FindAllMarkers(diet_sub, only.pos = TRUE, min.pct = 0.25)
  saveRDS(markers, paste0(dir_out,res,"markers.RDS"))
  markers %>%
   group_by(cluster) %>%
     top_n(n = 5, wt = avg_log2FC) -> top5
  pdf(paste0(dir_out,res,"markers_0.1.pdf"),width=18, height = 10)
   
  print(DotPlot(diet_sub, features=unique(top5$gene), dot.scale=6, cols="RdBu") +
          theme(axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5, hjust=1))) 
  
  dev.off()

}


subset_plot_markers(0.1,25)
  

  
  
## Visualize --- 

Dimplot(diet_sub,"seurat_clusters")  


## Some of the clusters seem to be overlapping between each other
## Might be because of low quality cells?
## Plotting hexplot to see overlap

  legend_breaks <- c(
    0, 0.001, 0.005,
    0.01, 0.05,
    0.1, 0.2, 1)
  lightblue <- dittoSeq::dittoColors()[2]
  yellow <- dittoSeq::dittoColors()[4]
  orangered <- dittoSeq::dittoColors()[6]
  legend_colors <- c(
    "grey95", Lighten(lightblue, 0.6), Lighten(lightblue, 0.3),
    Lighten(yellow, 0.3), Lighten(orangered, 0.3),
    orangered, Darken(orangered, 0.3), Darken(orangered, 0.75))
  
  clus <- "seurat_clusters"
  ggsave(
    paste0(dir_out, clus, "__hex__faceted_umap.png"),
    dittoDimHex(diet_sub, split.by = clus, color.panel = rep("black", 4000),
                main = clus) +
      scale_fill_gradientn(
        name = "Cells",
        values = legend_breaks, colors = legend_colors,
        breaks = function(limits) {c(0, 50, 200, 500, scales::breaks_extended(5)(limits))}
      ) +
      guides(fill=guide_colourbar(barheight=30)),
    units = "in",
    width = 15,
    height = 15
  )
  
## Since there's some overlap we'll filter some cells


