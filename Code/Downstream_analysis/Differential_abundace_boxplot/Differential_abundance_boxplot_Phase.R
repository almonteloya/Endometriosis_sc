
library(ggplot2)
library(Seurat)
library(dittoSeq)
library(dplyr)
library(ggpubr)
library(ggforce)

DA_dir<-"/krummellab/data1/immunox/XREP1a/10x/merged_XREP/DA/"

get_info <- function(meta.data, cluster_type, condition1){
  ##Get fraction of cells 
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


plot_box <- function(s_obj,D_S){
  s_obj@meta.data <- s_obj@meta.data %>% mutate(condition_full_reduced=recode(condition_full,
                                                                              "III-IV"="Disease",
                                                                              "I-II"="Disease",
                                                                              "Pathological_control" = "Pathological_control",
                                                                              "Pristine_control"="Pristine_control"))
  
  
  Idents(s_obj)<-s_obj$condition_full_reduced
  s_obj<-subset(s_obj,idents=D_S)
  data_ind<-get_info (s_obj@meta.data, annotation, Phase_menstrual_reduced)
  filtered_data <- data_ind %>%
    group_by(annotation, Phase_menstrual_reduced) %>%
    filter(n() > 2) %>%
    ungroup() %>%
    group_by(annotation) %>%
    filter(n_distinct(Phase_menstrual_reduced) == 2)
  
  stat.test <- filtered_data %>%
    group_by(annotation) %>%
    t_test(frac ~ Phase_menstrual_reduced) %>%
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
  
  
  grid <- data_ind %>% na.omit() %>% 
    ggplot(aes(x=Phase_menstrual_reduced, y=frac, color=Phase_menstrual_reduced))+
    geom_boxplot(width=0.1,outlier.alpha=0)+ geom_jitter(alpha=0.8)+    
    theme_light()+ theme(axis.text.x=element_blank()) +
    ylab("Fraction of cells") +
    scale_color_manual(values=dittoColors()[c(1,2)])+
    facet_grid(as.formula(paste0("~ ","annotation"))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + theme(legend.position = "none") 
  
  p1<- grid +stat_pvalue_manual(stat.test,label = "p.adj") + ggtitle(D_S)
  p2<-grid +stat_pvalue_manual(stat.test,label = "significance") + ggtitle(D_S)
  
  return(list(p1,p2))
}


#### Epithelial:
dir='/krummellab/data1/immunox/XREP1a/10x/merged_XREP/'
dir_out='/krummellab/data1/immunox/XREP1a/10x/merged_XREP/subset_celltype/epithelia/difuse_fix/recluster/'
s_obj = readRDS(file=paste0(dir_out,"epithelial_filter_obj.RDS"))
metadata<-readRDS(file=paste0(dir_out,"metadata_annotation.RDS"))
s_obj@meta.data <- metadata
DimPlot(s_obj)
p_pro<-plot_box(s_obj,"Disease")
p_sec<-plot_box(s_obj,"Pristine_control")

pdf(paste0(DA_dir,"epithelial_padj_DS.pdf"),height = 2)
print(p_pro[[1]])
print(p_pro[[2]])
print(p_sec[[1]])
print(p_sec[[2]])
dev.off()


## STROMAL
dir_in <- '/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/'
dir_out<-"/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/SFig8/"
s_obj <- readRDS(sprintf("%s/Stromal_obj_2025.RDS",dir_in))
DimPlot(s_obj)
p_pro<-plot_box(s_obj,"Disease")
p_sec<-plot_box(s_obj,"Pristine_control")

pdf(paste0(dir_out,"stromal_padj_DS.pdf"),height = 2)
print(p_pro[[1]])
print(p_pro[[2]])
print(p_sec[[1]])
print(p_sec[[2]])
dev.off()


## Immune
dir_out='/krummellab/data1/immunox/XREP1a/10x/merged_XREP/subset_celltype/immune/'
s_obj <- readRDS(paste0(dir_out,"immune_annotation_all.RDS"))
meta_annotation <- readRDS(file=paste0(dir_out,"immune_aggregate.RDS"))
s_obj@meta.data <- meta_annotation
s_obj$annotation<- s_obj$immune_agregate1_c
Idents(s_obj)<-s_obj$annotation
DimPlot(s_obj)

p_pro<-plot_box(s_obj,"Disease")
p_sec<-plot_box(s_obj,"Pristine_control")

pdf(paste0(DA_dir,"immune_padj_phase.pdf"),height = 2,width = 10)
print(p_pro[[1]])
print(p_pro[[2]])
print(p_sec[[1]])
print(p_sec[[2]])
dev.off()
