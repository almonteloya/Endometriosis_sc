library(dplyr)
library(ggplot2)
library(ggrepel)


DEG_dir<-"/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/pristine_disease/broad/"
dir_out<-"/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/Fig2/"
data1 <- readRDS(paste0(DEG_dir,"Stromal_fibroblast_proliferative_DEG.RDS"))
data2 <- readRDS(paste0(DEG_dir,"Stromal_fibroblast_secretory_DEG.RDS"))

data3 <- readRDS(paste0(DEG_dir,"Unciliated_Epithelia_proliferative_DEG.RDS"))
data4 <- readRDS(paste0(DEG_dir,"Unciliated_Epithelia_secretory_DEG.RDS"))
                                      
                              
data5 <- readRDS(paste0(DEG_dir,"Ciliated_Epithelia_proliferative_DEG.RDS"))
data6 <- readRDS(paste0(DEG_dir,"Ciliated_Epithelia_secretory_DEG.RDS"))
                                  
data7 <- readRDS(paste0(DEG_dir,"Smooth_muscle_proliferative_DEG.RDS"))
data8 <- readRDS(paste0(DEG_dir,"Smooth_muscle_secretory_DEG.RDS"))

data9 <- readRDS(paste0(DEG_dir,"Endothelia_proliferative_DEG.RDS"))
data10 <- readRDS(paste0(DEG_dir,"Endothelia_secretory_DEG.RDS"))


library(ggplot2)
library(ComplexUpset)

###### DEG -- both up and down
## Proliferative
listInput <- list(data1,data3,data5,data7,data9)

r_filter <- function(df){
  x<- df %>% dplyr::filter(avg_log2FC > .25 & p_val_adj < 0.05 )
  genesup <- paste0(rownames(x),"_UP")
  x<- df %>% dplyr::filter(avg_log2FC < -0.25 & p_val_adj < 0.05 )
  genesdown <- paste0(rownames(x),"_DOWN")
  return(c(genesup, genesdown))
}

significant <- lapply(listInput,r_filter)

names(significant) <- c('Fibroblast','Unciliated Epithelia',
                        "Ciliated Epithelia","V.S.M",
                        "Endothelia")


all_genes <- unique(unlist(significant))

# Step 3: Create an empty dataframe with gene IDs as rows
gene_presence_df <- data.frame(GeneID = all_genes)

# Step 4: Fill the dataframe with TRUE/FALSE values
for (vector_name in names(significant)) {
  gene_presence_df[[vector_name]] <- gene_presence_df$GeneID %in% significant[[vector_name]]
}

# View the dataframe
gene_presence_df <- gene_presence_df %>% tibble::column_to_rownames("GeneID")


gene_presence_df$sign <- sub(".*_(.*)", "\\1", rownames(gene_presence_df))

selected_col <- colnames(gene_presence_df)[-length(colnames(gene_presence_df))]

pdf(paste0(dir_out,"proliferative_upset.pdf"),width = 8,height = 6)
ComplexUpset::upset(gene_presence_df, selected_col,set_sizes=(upset_set_size() + theme(axis.text.x=element_text(angle=45))
),width_ratio=0.2,min_size=3,    themes=upset_modify_themes(
  list(
    'intersections_matrix'=theme(text=element_text(size=20)),
    'overall_sizes'=theme(axis.text.x=element_text(angle=90)))),
base_annotations = 
  list(  'Number of DEG'=(intersection_size(bar_number_threshold=.3, counts=FALSE,
                                            mapping=aes(fill=sign)) + 
                            scale_fill_manual(values=c("UP"="red","DOWN"="blue"))+
                            ylim(c(0, 900)) ))
)

dev.off()

#### Secretroy 
listInput <- list(data2,data4,data6,data8,data10)


significant <- lapply(listInput,r_filter)

names(significant) <- c('Fibroblast','Unciliated_Epithelia',
                        "Ciliated_Epithelia","V.S.M",
                        "Endothelia")


all_genes <- unique(unlist(significant))

# Step 3: Create an empty dataframe with gene IDs as rows
gene_presence_df <- data.frame(GeneID = all_genes)

# Step 4: Fill the dataframe with TRUE/FALSE values
for (vector_name in names(significant)) {
  gene_presence_df[[vector_name]] <- gene_presence_df$GeneID %in% significant[[vector_name]]
}

# View the dataframe
gene_presence_df <- gene_presence_df %>% tibble::column_to_rownames("GeneID")


gene_presence_df$sign <- sub(".*_(.*)", "\\1", rownames(gene_presence_df))

selected_col <- colnames(gene_presence_df)[-length(colnames(gene_presence_df))]

pdf(paste0(dir_out,"secretory_upset.pdf"),width = 8,height = 6)
ComplexUpset::upset(gene_presence_df, selected_col,set_sizes=(upset_set_size() + theme(axis.text.x=element_text(angle=45))
),width_ratio=0.2,min_size=3,    themes=upset_modify_themes(
  list(
    'intersections_matrix'=theme(text=element_text(size=20)),
    'overall_sizes'=theme(axis.text.x=element_text(angle=90)))),
base_annotations = 
  list(  'Number of DEG'=(intersection_size(bar_number_threshold=.3, counts=FALSE,
                                            mapping=aes(fill=sign)) + 
                            scale_fill_manual(values=c("UP"="red","DOWN"="blue"))+
                            ylim(c(0, 900)) ))
)
dev.off()
