

## Comapring overlap between phases
## 
library(ggplot2)
library(Seurat)
library(dittoSeq)
library(ggrepel)
library(dplyr)

dir_phase <-"/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/menstrual_phase/"
all_files <- list.files( path= dir_phase, pattern = paste0(".*", ".*RDS"), full.names = TRUE)

r_filter <- function(df){
  x<- df %>% dplyr::filter(avg_log2FC > .25 & p_val_adj < 0.05 )
  genesup <- paste0(rownames(x),"_UP")
  x<- df %>% dplyr::filter(avg_log2FC < -0.25 & p_val_adj < 0.05 )
  genesdown <- paste0(rownames(x),"_DOWN")
  return(c(genesup, genesdown))
}


celltypes = c("Ciliated","Unciliated","Endothelia","fibroblast","Smooth-Muscle")

vennplot <- function(celltype){
  files <- all_files[grepl(celltype, all_files)]
  all_markers <-lapply(files, readRDS)
  type_dis <- basename(files)
  matched_words <- str_extract(type_dis, "(disease|Pristine|Fibroid)")
  names(all_markers)<-matched_words
  genes_filter<-lapply(all_markers,r_filter)
    gplot<-ggVennDiagram(genes_filter)+ 
      scale_fill_gradient(low = "#F4FAFE", high = "#F4FAFE") +
      scale_x_continuous(expand = expansion(mult = .1)) + ggtitle(celltype)
    print(gplot)
}

pdf("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/SPhase/venn.pdf", height = 4,width = 5)
lapply(celltypes, vennplot)
dev.off()
