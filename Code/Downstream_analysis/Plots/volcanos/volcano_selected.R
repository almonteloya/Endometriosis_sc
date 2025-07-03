## Volcano plots with hilighthing genes
## In Stromal genes 

library(ggplot2)
library(dittoSeq)
library(clusterProfiler)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

volcanoplot<-function(DEG,gene_list){
  DEG$diffexpressed <- "NO"
  ##UP VS DOWN
  DEG$diffexpressed[DEG$avg_log2FC > .25 & DEG$avg_log2FC > 0 & DEG$p_val_adj < 0.05] <- "UP"
  DEG$diffexpressed[DEG$avg_log2FC < -0.25  & DEG$p_val_adj < 0.05] <- "DOWN"
  DEG$delabel <- NA
  ## If its either up or down get the name of the genes
  DEG$delabel[DEG$diffexpressed != "NO"] <- rownames(DEG)[DEG$diffexpressed != "NO"]
  DEG$fontface <- ifelse(DEG$delabel %in% gene_list, "bold", "plain")
  print(ggplot(data=DEG, aes(x=avg_log2FC, y=-log10(p_val), col=diffexpressed, label=delabel))  +
          geom_point(alpha = .5) + theme_minimal()  +  geom_text_repel(aes(fontface = fontface),force=.2,size=6,max.overlaps =20) +
          scale_color_manual(values=c("blue", "gray", "red")) + 
          labs(x = "Downregulated endometriosis     Avg Log2FC       Upregulated endometriosis"))    +
    theme(axis.title.x = element_text(size = 12))
}



###------------------ STROMAL

##---- Secretory
dir_GSEA = "/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/pristine_disease/stromal_specific/secretory/GSEA/"
dir_out="/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/Fig3/"
##-- Reading the GSEA
pat<-"*RDS$"
files_GSEA<- list.files( path= dir_GSEA, pattern = pat, full.names = TRUE)
filtered_markers <- lapply(files_GSEA, readRDS)
names(filtered_markers)<-  basename(files_GSEA)

####------ VSM
DEG_VSM<-readRDS("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/pristine_disease/stromal_specific/secretory/vascular-smooth-muscle_secretory_DEG.RDS")
gsea<-filtered_markers[[1]]

path1<-gsea[(gsea$Description=="regulation of cell death"),]$core_enrichment
path2<-gsea[(gsea$Description=="response to temperature stimulus"),]$core_enrichment
path3<-gsea[(gsea$Description=="cellular response to stress"),]$core_enrichment

genes_selected<-c(path1,path2,path3)
split_vector <- strsplit(genes_selected, "/")
gene_list<-unique(unlist(split_vector))
## Printing volcano
pdf(paste0(dir_out,"VSM_volcano.pdf"))
volcanoplot(DEG_VSM,gene_list)
dev.off()

### WFC2 -- we are using the same pathways since they are shared
## Printing volcano
DEG_WFC2<-readRDS("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/pristine_disease/stromal_specific/secretory/WFDC2+-fibroblast_secretory_DEG.RDS")
pdf(paste0(dir_out,"WFC2_volcano.pdf"))
volcanoplot(DEG_WFC2,gene_list)
dev.off()



### ----- Proliferative
gsea = readRDS("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/pristine_disease/stromal_specific/proliferative/GSEA/HSP-stromal-fibroblast_proliferative_DEG.RDS_GSEA.RDS")
dir_out="/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/Fig3/"
##-- Reading the GSEA
#### HSP
DEG_HSP<-readRDS("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/pristine_disease/stromal_specific/proliferative/HSP-stromal-fibroblast_proliferative_DEG.RDS")

path1<-gsea[(gsea$Description=="immune system development"),]$core_enrichment
path2<-gsea[(gsea$Description=="cellular response to stress"),]$core_enrichment

genes_selected<-c(path1,path2)
split_vector <- strsplit(genes_selected, "/")
gene_list<-unique(unlist(split_vector))
## Printing volcano
pdf(paste0(dir_out,"HSP_volcano.pdf"))
volcanoplot(DEG_VSM,gene_list)
dev.off()


###------------------ Epithelial

##---- Secretory
dir_GSEA = "/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/pristine_disease/epithelial_specific/GSEA/"
dir_out="/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/Fig4/"
##-- Reading the GSEA
pat<-"*RDS$"
files_GSEA<- list.files( path= dir_GSEA, pattern = pat, full.names = TRUE)
filtered_markers <- lapply(files_GSEA, readRDS)
names(filtered_markers)<-  basename(files_GSEA)

####------ EMT
DEG_EMT<-readRDS("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/pristine_disease/epithelial_specific/secretory/EMT-epithelia-1_secretory_DEG.RDS")
gsea<-filtered_markers[[2]]

path1<-gsea[(gsea$Description=="regulation of locomotion"),]$core_enrichment
path2<-gsea[(gsea$Description=="cell migration"),]$core_enrichment
path3<-gsea[(gsea$Description=="cell adhesion"),]$core_enrichment

genes_selected<-c(path1,path2,path3)
split_vector <- strsplit(genes_selected, "/")
gene_list<-unique(unlist(split_vector))
## Printing volcano
pdf(paste0(dir_out,"EMT_volcano.pdf"))
volcanoplot(DEG_EMT,gene_list)
dev.off()

####------ Ciliated
DEG_Cil<-readRDS("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/pristine_disease/epithelial_specific/secretory/Ciliated-epithelia_secretory_DEG.RDS")
gsea<-filtered_markers[[1]]
genes_selected<-gsea[(gsea$Description=="programmed cell death"),]$core_enrichment
split_vector <- strsplit(genes_selected, "/")
gene_list<-unique(unlist(split_vector))
## Printing volcano
pdf(paste0(dir_out,"Cil_volcano.pdf"))
volcanoplot(DEG_Cil,gene_list)
dev.off()



##---- Proliferative
dir_GSEA = "/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/pristine_disease/epithelial_specific/GSEA/"
dir_out="/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/Fig4/"
##-- Reading the GSEA
pat<-"*RDS$"
files_GSEA<- list.files( path= dir_GSEA, pattern = pat, full.names = TRUE)
filtered_markers <- lapply(files_GSEA, readRDS)
names(filtered_markers)<-  basename(files_GSEA)

####------ G1
DEG_G1<-readRDS("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/pristine_disease/epithelial_specific/proliferative/Glandular-1_proliferative_DEG.RDS")
gsea<-filtered_markers[[3]]

path1<-gsea[(gsea$Description=="female pregnancy"),]$core_enrichment

genes_selected<-c(path1)
split_vector <- strsplit(genes_selected, "/")
gene_list<-unique(unlist(split_vector))
## Printing volcano
pdf(paste0(dir_out,"Glandular1_volcano.pdf"))
volcanoplot(DEG_G1,gene_list)
dev.off()
