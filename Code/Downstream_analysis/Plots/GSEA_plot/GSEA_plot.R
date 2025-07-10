### Reading and plotting GSEA

library(ggplot2) 
library(Seurat)
library(dittoSeq)

dir_out="broad/GSEA/"
dir = "DEG/pristine_disease/broad/"

pat<-"*RDS$"
files <- list.files( path= dir_out, pattern = pat, full.names = TRUE)
pattern <- ".*\\/([^/]+)_GSEA\\.RDS$"
names_file <- sub(pattern, "\\1", files)
gc_list <- lapply(files, readRDS)
names(gc_list)<- names_file

filtered_gc <- gc_list[sapply(gc_list, function(df) nrow(df) > 0)]
names_file<-names(filtered_gc)

## Top 10 and filter
top10filter<-function(x){
  x<-x %>% filter(p.adjust<0.05) %>% arrange(desc(enrichmentScore),p.adjust)
  x <- x %>% arrange(desc(NES))
  print(x)
  x<-x[1:10,]
  x
}


### Plotting 

for (i in (1:length(filtered_gc))) {
fgsea<- top10filter(filtered_gc[[i]])
plot_gsea<- ggplot(fgsea, aes(x= NES, y=reorder(Description,enrichmentScore), size=setSize, color=p.adjust,label=NA)) + geom_point(alpha = 0.8)+
  scale_color_gradient(low="red", high="blue") + geom_text(size=15) + labs(y="Description")+ theme_bw() +
  theme(axis.text = element_text(size = 13),axis.text.x = element_text(size=12,angle = 45, vjust = 0.5, hjust=.5)) +
  ggtitle(names_file[i])

print(plot_gsea)
pdf(paste0(dir_out,names_file[i],".pdf"),height =4,width = 6)
print(plot_gsea)
dev.off()
}

