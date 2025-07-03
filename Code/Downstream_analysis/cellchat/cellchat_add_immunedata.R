## For cellchat we want to connect immune and non-immune
## Add the immune data into the non_immune one


library(Seurat)

## preparing object
dir_obj="path/"

merged_data = readRDS(file=paste0(dir_obj,"merged_data"))

### We also need the immune cell part 
immune_data <- readRDS("immune_annotation.RDS")
ident_number <- levels(Idents(immune_data))

for (i in 0:length(ident_number)) {
  cell_names <- colnames(subset(immune_data,idents=ident_number[i]))
  Idents(merged_data, cells = cell_names) <- ident_number[i]
}

## This is the final annotation object
saveRDS(merged_data,paste0(dir_obj,"merged_data_celltype_all.RDS"))
