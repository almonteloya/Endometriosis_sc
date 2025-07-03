### Reading and plotting immune

library(ggplot2) 
library(Seurat)
library(dittoSeq)

dir_out="//krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/Fig5/"
NK_GSEA = readRDS("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/pristine_disease/immune_broad/secretory/GSEA2/NK_cells_secretory_DEG.RDS_GSEA.RDS")

## Top 10 and filter

NK_GSEA<-NK_GSEA %>% filter(p.adjust<0.05) %>% arrange(desc(enrichmentScore),p.adjust)
NK_GSEA<- NK_GSEA %>% arrange(desc(NES))
NK_GSEA



### Plotting 

plot_gsea<- ggplot(NK_GSEA, aes(x= NES, y=reorder(Description,enrichmentScore), size=setSize, color=p.adjust,label=NA)) + geom_point(alpha = 0.8)+
    scale_color_gradient(low="red", high="blue") + geom_text(size=15) + labs(y="Description")+ theme_bw() +
    theme(axis.text = element_text(size = 13),axis.text.x = element_text(size=12,angle = 45, vjust = 0.5, hjust=.5)) +
    ggtitle("NK Secretory")
  
  pdf(paste0(dir_out,"NK_GSEA.pdf"),height =4,width = 6)
  print(plot_gsea)
  dev.off()


