library(ggplot2)
library(Seurat)
library(dittoSeq)
library(dplyr)

dir_in='/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/seurat_objects/'
s_obj = readRDS(file=paste0(dir_in,"Immune_obj_2025.RDS"))
s_obj@meta.data <- s_obj@meta.data %>% mutate(condition_control=recode(condition_full,
                                                                       "III-IV" ="disease",
                                                                       "I-II" ="disease",
                                                                       "Pathological_control"= "Pathological_control" ,
                                                                       "Pristine_control" = "Pristine_control"))
subset_cell<-subset(x = s_obj, subset = condition_control != "Pathological_control")

dir_out<- '/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/Sfig14/'
dir.create(dir_out)

##### Barplots
make_meta_barplot = function(sobj, feature, cluster, col_vec){
  sobj@meta.data %>% na.omit() %>%
    dplyr::select({{feature}}, {{cluster}}) %>%
    dplyr::rename(cluster={{cluster}}) %>%
    ggplot(aes(x=cluster, fill={{feature}}))+
    geom_bar(position="fill")+
    theme_bw()+
    scale_fill_manual(values=dittoColors(length(unique(sobj@meta.data %>% pull({{feature}})))))+
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2),
          axis.title.x=element_blank())+
    ylab("fraction of cells")
}

pdf(paste0(dir_out,"bar_Condition_full_nf.pdf"),height =8, width = 8)
make_meta_barplot(subset_cell, condition_full, immune_agregate1, dittoColors())+
  scale_fill_manual(values=dittoColors()[c(3,11,4)])+
  ylab("fraction of cluster") + theme(legend.text = element_text(size=15),
                                      axis.text=element_text(size=15),)
dev.off()

pdf(paste0(dir_out,"bar_Condition_control_nf.pdf"),height =8, width = 8)
make_meta_barplot(subset_cell, condition_control, immune_agregate1, dittoColors())+
  scale_fill_manual(values=dittoColors()[c(3,4)])+
  ylab("fraction of cluster") + theme(legend.text = element_text(size=15),
                                      axis.text=element_text(size=15),)
dev.off()

pdf(paste0(dir_out,"bar_Phase_nf.pdf"),height =8, width = 8)
make_meta_barplot(subset_cell, Phase_menstrual_reduced, immune_agregate1, dittoColors())+
  scale_fill_manual(values=dittoColors()[c(2,1)])+
  ylab("fraction of cluster") + theme(legend.text = element_text(size=15),
                                      axis.text=element_text(size=15),)
dev.off()

pdf(paste0(dir_out,"bar_Library_nf.pdf"),height =8, width = 8)
make_meta_barplot(subset_cell, LIBRARY, immune_agregate1, dittoColors())+
  scale_fill_manual(values=dittoColors())+
  ylab("fraction of cluster") + theme(legend.text = element_text(size=15),
                                      axis.text=element_text(size=15),)
dev.off()

pdf(paste0(dir_out,"bar_patient_nf.pdf"),height =8, width = 8)
make_meta_barplot(subset_cell, Patient_id, immune_agregate1, dittoColors())+
  scale_fill_manual(values=dittoColors())+
  ylab("fraction of cluster") + theme(legend.text = element_text(size=15),
                                      axis.text=element_text(size=15),)
dev.off()






### ------------ NK Cells

dir_out<- '/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/SFig15/'
  pdf(paste0(dir_out,"bar_patient.pdf"),height =8, width = 8)
  make_meta_barplot(subset_cell, Patient_id, NK_annotation, dittoColors())+
    scale_fill_manual(values=dittoColors())+
    ylab("fraction of cluster") + theme(legend.text = element_text(size=18),
                                        axis.text=element_text(size=12),)
  dev.off()
  
  pdf(paste0(dir_out,"bar_Condition.pdf"),height =8, width = 8)
  make_meta_barplot(subset_cell, Condition, NK_annotation, dittoColors())+
    scale_fill_manual(values=dittoColors()[c(7,3)])+
    ylab("fraction of cluster") 
  dev.off()
  
  pdf(paste0(dir_out,"bar_Condition_full.pdf"),height =8, width = 8)
  make_meta_barplot(subset_cell, condition_full, NK_annotation, dittoColors())+
    scale_fill_manual(values=dittoColors()[c(3,11,4)])+
    ylab("fraction of cluster") 
  dev.off()
  
  pdf(paste0(dir_out,"bar_Condition_control.pdf"),height =8, width = 8)
  make_meta_barplot(subset_cell, condition_control, NK_annotation, dittoColors())+
    scale_fill_manual(values=dittoColors()[c(3,4)])+
    ylab("fraction of cluster") 
  dev.off()
  
  pdf(paste0(dir_out,"bar_Phase.pdf"),height =8, width = 8)
  make_meta_barplot(subset_cell, Phase_menstrual_reduced, NK_annotation, dittoColors())+
    scale_fill_manual(values=dittoColors()[c(2,1)])+
    ylab("fraction of cluster") 
  dev.off()
  
  pdf(paste0(dir_out,"bar_Library.pdf"),height =8, width = 8)
  make_meta_barplot(subset_cell, LIBRARY, NK_annotation, dittoColors())+
    scale_fill_manual(values=dittoColors())+
    ylab("fraction of cluster") 
  dev.off()





