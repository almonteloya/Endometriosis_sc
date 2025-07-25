### GSEA
## Looking into multiple subtypes and generating GSEA files
## Using DEG results from DEG_specific.R 

library(ggplot2)
library(Seurat)
library(dittoSeq)
library(clusterProfiler)
library(fgsea)

organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

dir = "/dir/DEG/output"

dir_out = paste0(dir,"GSEA/")
dir.create(dir_out)

##-- Reading the gene information

pat<-"*RDS$"
files <- list.files( path= dir, pattern = pat, full.names = TRUE)
file_name <- basename(files)
filtered_markers <- lapply(files, readRDS)
names(filtered_markers)<- basename(files)



GSEA_plot<- function(markers,cell_type){
  markers_pvalue2 <- (markers %>% filter(p_val_adj<.05) %>% 
                        arrange(desc((avg_log2FC))) %>% filter(abs(avg_log2FC) > 0.25))
  print(dim(markers_pvalue2))
  if ((nrow(markers_pvalue2) > 20)) { ## having at least 20 DEG
    genelist<-markers_pvalue2$avg_log2FC
    names(genelist)<- rownames(markers_pvalue2)
    se <- gseGO(geneList= genelist,  ont ="BP", keyType = "SYMBOL",
                OrgDb= org.Hs.eg.db,pvalueCutoff = 0.05,eps=0)
    
    require(DOSE)
    dim(se@result)
    saveRDS(se@result, paste0(dir_out,cell_type,"_GSEA.RDS"))
  }
  
}

for (i in c(1:length(filtered_markers))) {
  GSEA_plot(filtered_markers[[i]], basename(files)[[i]] )
}


