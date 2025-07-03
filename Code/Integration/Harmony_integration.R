### Aug 22 2022, 
## Harmony all the samples filtered


library(cowplot)
library(tidyverse)
library(Seurat)
library(harmony)

setwd("/path/here/")
prefix='merged_XREP'
dir.create(dirname(prefix), showWarnings=T)

IN_DIR ="path/here/input/"
OUT_DIR = "path/here/output"

### LOAD DATA ###

list_samples = c("XREP1-POOL-SCG1","XREP1-POOL-SCG2","XREP1-POOL-SCG3",
                 "XREP1-POOL-SCG7","XREP1-POOL-SCG8","XREP1-POOL-SCG4",
                 "XREP1-POOL-SCG5","XREP1-POOL-SCG6","XREP1-POOL-SCG9","XREP1-POOL-SCG10","XREP1-POOL-SCG11",
                 "XREP1-POOL-SCG12","XREP1-POOL-SCG13")

### normalize each dataset ###
set.seed(1208)
sobjs = lapply(list_samples, FUN=function(SAMPLE_NAME){
  sobj = readRDS(paste0(IN_DIR,SAMPLE_NAME,"/automated_processing_filter/",SAMPLE_NAME,"_scTransformed_processed_filtered.RDS"))
  DefaultAssay(sobj) = "RNA" # need to remove SCT data
  sobj[["SCT"]] = NULL
  sobj = NormalizeData(sobj)
  sobj@meta.data$LIBRARY = sobj@meta.data$orig.ident
  return(sobj)
})
names(sobjs) = list_samples
print("Samples loaded")

gene_intersection <- lapply(sobjs, row.names) %>% Reduce(intersect, .)

first_sobj <- sobjs[[names(sobjs)[1]]]
sobjs[[names(sobjs)[1]]] <- NULL
print("abt to run merge")
merged_data <- merge(x = first_sobj, y = unname(sobjs))
print("ran merge")
merged_data <- FindVariableFeatures(merged_data, selection.method = "vst", nfeatures = 3000)
merged_data <- ScaleData(merged_data, vars.to.regress=c('percent.mt', 'percent.ribo', 'S.Score', 
                                                        'G2M.Score'), verbose = FALSE)
merged_data <- RunPCA(merged_data, npcs = 30, verbose = FALSE)

save(merged_data, file=paste0(prefix, '_merged_temp.RData'))
# plot harmony convergence
pdf(paste0(prefix, '_harmony_convergence.pdf'))
merged_data <- RunHarmony(merged_data,
                          "LIBRARY",
                          assay.use='RNA',
                          plot_convergence = TRUE,
                          max.iter.harmony=20,
                          max.iter.cluster=30)
dev.off()

merged_data <- RunUMAP(merged_data,
                       dims = 1:30,  # Num PCs to use
                       reduction='harmony',
                       n.neighbors = 30,  # Default. Controls how UMAP balances local (low) versus global (large) structure in the data
                       min.dist = 0.3,   # Default. Controls the size of the clusters. Should be smaller than spread
                       spread = 1,  # Default. Controls the inter-cluster distances to some extent. Should be larger than min_dist
                       a = NULL,  # Default. Can be used with b instead of using min.dist/spread
                       b = NULL,  # Default. Can be used with a instead of using min.dist/spread
                       verbose = FALSE)


# Calculate the neighborhood graph
merged_data <- FindNeighbors(merged_data,
                             dims = 1:30,  # Num PCs to use
                             reduction='harmony',
                             k.param = 20,  # k for the knn algorithm
                             verbose = FALSE
)

png(filename=paste0(prefix, '_library_umap.png'), width = 5, height = 5, units = "in", res = 300)
print(DimPlot(merged_data, group.by='LIBRARY') + NoLegend())
dev.off()

png(filename=paste0(prefix, '_split_library_umap.png'), width = 12, height = 6, units = "in", res = 300)
print(DimPlot(merged_data, split.by='LIBRARY', ncol=4) + theme(legend.position="none", axis.title=element_blank(), axis.text=element_blank()))
dev.off()


# run clustering on at of different resolutions

for (res in c(0.2,0.6,0.4, 0.8,1)){
  if (paste0('louvain_res', res) %in% colnames(merged_data@meta.data)){
    next
  }
  merged_data <- FindClusters(merged_data, verbose = TRUE,
                              algorithm = 1,
                              resolution = res)
  merged_data@meta.data[[paste0('louvain_res', res)]] <- merged_data@meta.data$seurat_clusters
  
  png(filename=paste0(prefix, '_louvain_res_', res, '.png'), width = 5, height = 5, units = "in", res = 300)
  print(DimPlot(merged_data, group.by=paste0('louvain_res', res), label=T) + NoLegend())
  dev.off()
}

save(merged_data, file=paste0(prefix, '_merged_processed.RData'))

# merge the metadata
metadata <- merged_data@meta.data
for (redn in c('pca', 'umap')){
  cell_embeddings <- as.data.frame(merged_data@reductions[[redn]]@cell.embeddings[, c(1: min(5, dim(merged_data@reductions[[redn]]@cell.embeddings)[2]))])
  metadata <- merge(metadata,
                    cell_embeddings,
                    by=0)
  rownames(metadata) <- metadata$Row.names
  metadata$Row.names <- NULL
}

write.table(metadata,
            file=paste(prefix, 'metadata.tsv', sep='_'),
            sep='\t',
            row.names=T,
            col.names=T,
            quote=F)



