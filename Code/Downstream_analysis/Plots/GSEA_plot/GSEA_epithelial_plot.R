### reading and plotting GSEA
## Epithelia -- secretory

library(ggplot2)
library(clusterProfiler)
library(fgsea)

dir_secretory<-("data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/pristine_disease/epithelial_specific//GSEA/")
dir_out<-("data1/immunox/XREP1a/10x/merged_XREP/Fix2025/Fig4/")

pat<-"*RDS$"
files <- list.files( path= dir_secretory, pattern = pat, full.names = TRUE)
#files<-c(files_sec,files_pro)
cellNames<-basename(files)

gc_list<-lapply(files, readRDS)
names(gc_list) <- cellNames

filtered_gc <- gc_list[sapply(gc_list, function(df) nrow(df) > 0)]

result_access <- function(x){
   x<-x %>% filter(enrichmentScore>0,p.adjust<0.05) %>% arrange(desc(enrichmentScore),p.adjust)
  #x<-x %>% filter(p.adjust<0.05) %>% arrange(desc(enrichmentScore),p.adjust)
  x <- x %>% arrange(desc(NES))
  print(x)
  x<-x[1:20,]
  x
}

top_10 <- lapply(filtered_gc, result_access)


for (i in 1:length(top_10)) {
  top_10[[i]]$gene <- names(top_10)[i] 
}


### In this case I only want to plot 2
#common_df<-rbind(top_10[[1]],top_10[[2]],top_10[[3]])
common_df<-rbind(top_10[[1]],top_10[[2]],top_10[[3]])
common_df$gene_f <- factor(common_df$gene, levels = c(
  "Glandular-1_proliferative_DEG.RDS_GSEA.RDS",
  "Ciliated-epithelia_secretory_DEG.RDS_GSEA.RDS",
  "EMT-epithelia-1_secretory_DEG.RDS_GSEA.RDS"
))

S1<- ggplot(common_df, aes(x= NES, y=reorder(Description,enrichmentScore), size=setSize, color=p.adjust,label=NA)) + geom_point(alpha = 0.8)+
  scale_color_gradient(low="red", high="blue") + geom_text(size=15) + labs(y="Description")+ theme_minimal() +
  theme(axis.text = element_text(size = 10),axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=.5)) 
S1 <- S1  + facet_wrap(~ gene_f)

pdf(paste0(dir_out,"GSEA.pdf"),width = 7)
print(S1)
dev.off()

