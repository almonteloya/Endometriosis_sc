

### Differentil abundance for subpopulations
## Pristine vs disease 
## Stratified by proliferative
library(ggplot2)
library(Seurat)
library(dittoSeq)
library(dplyr)
library(ggpubr)
library(ggforce)

get_info <- function(meta.data, cluster_type, condition1){
  ##Get fraction of cells per each cell type per contion
  frac_group_ind = {{meta.data}} %>% 
    dplyr::select({{cluster_type}}, Patient_id, {{condition1}}) %>%
    group_by(Patient_id) %>%
    dplyr::mutate(tot=n()) %>%
    ungroup() %>%
    group_by({{condition1}}, Patient_id, tot, {{cluster_type}}) %>%
    dplyr::count() %>%
    mutate(frac=n/tot) %>%
    dplyr::arrange(desc(frac))
  frac_group_ind
}


plot_box <- function(s_obj,PHASE){
  # Plot boxplot
    Idents(s_obj)<-s_obj$condition_full_reduced
  s_obj<-subset(s_obj,idents=c("Disease","Pristine_control"))
  Idents(s_obj)<-s_obj@meta.data$Phase_menstrual_reduced
  s_data <- subset(s_obj, idents= PHASE)
  data_ind<-get_info (s_data@meta.data, annotation, condition_full_reduced)
  
  ## need at least two per condition to perfom t-test
  filtered_data <- data_ind %>%
    group_by(annotation, condition_full_reduced) %>%
    filter(n() > 2) %>% 
    ungroup() %>%
    group_by(annotation) %>%
    filter(n_distinct(condition_full_reduced) == 2)
  
  ## Perform t-test 
  stat.test <- filtered_data %>%
    group_by(annotation) %>%
    t_test(frac ~ condition_full_reduced) %>%
    adjust_pvalue() %>%
    mutate(
      y.position = 1,
      significance = case_when(
        p.adj < 0.001 ~ "***",
        p.adj < 0.01  ~ "**",
        p.adj < 0.05  ~ "*",
        TRUE            ~ "ns"
      )
    )
  
  ## Plot and add p-value
  grid <- data_ind %>% na.omit() %>% 
    ggplot(aes(x=condition_full_reduced, y=frac, color=condition_full_reduced))+
    geom_boxplot(width=0.1,outlier.alpha=0)+ geom_jitter(alpha=0.8)+    
    theme_light()+ theme(axis.text.x=element_blank()) +
    ylab("Fraction of cells") +
    scale_color_manual(values=dittoColors()[c(3,1,4)])+
    facet_grid(as.formula(paste0("~ ","annotation"))) +
    scale_y_continuous(expand = expansion(mult = c(.2, 0.2))) + theme(legend.position = "none") 
  
  p1<- grid +stat_pvalue_manual(stat.test,label = "p.adj") + ggtitle(PHASE)
  p2<-grid +stat_pvalue_manual(stat.test,label = "significance") + ggtitle(PHASE)
  
  return(list(p1,p2))
}

dir_in=("dir/seurat/object")
DA_dir=('out/dir')

#### Epithelial:
s_obj = readRDS(file=paste0(dir,"Epithelial_obj_2025.RDS"))
DimPlot(s_obj)
p_pro<-plot_box(s_obj,"proliferative")
p_sec<-plot_box(s_obj,"secretory")

pdf(paste0(DA_dir,"epithelial_padj.pdf"),height = 2)
print(p_pro[[1]])
print(p_pro[[2]])
print(p_sec[[1]])
print(p_sec[[2]])
dev.off()


## STROMAL
s_obj <- readRDS(sprintf("%s/Stromal_obj_2025.RDS",dir_in))
DimPlot(s_obj)
p_pro<-plot_box(s_obj,"proliferative")
p_sec<-plot_box(s_obj,"secretory")

pdf(paste0(dir_out,"stromal_padj.pdf"),height = 2)
print(p_pro[[1]])
print(p_pro[[2]])
print(p_sec[[1]])
print(p_sec[[2]])
dev.off()


## Immune
s_obj <- readRDS(paste0(dir_in,"Immune_obj_2025.RDS"))
s_obj$annotation<- s_obj$immune_agregate1_c
Idents(s_obj)<-s_obj$annotation
DimPlot(s_obj)
p_pro<-plot_box(s_obj,"proliferative")
p_sec<-plot_box(s_obj,"secretory")

pdf(paste0(DA_dir,"immune_padj.pdf"),height = 2)
print(p_pro[[1]])
print(p_pro[[2]])
print(p_sec[[1]])
print(p_sec[[2]])
dev.off()
