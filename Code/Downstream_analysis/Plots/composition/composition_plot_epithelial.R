## cluster compostion 
## Broad celltypes

library(ggplot2)
library(Seurat)
library(dittoSeq)
library(dplyr)

dir_in='/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/seurat_objects/'
dir_out<- '/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/Sfig11/'
dir.create(dir_out)

s_obj = readRDS(file=paste0(dir_in,"Epithelial_obj_2025.RDS"))

### Now cluster compostions


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

## Exclude fibroids from analysis 
subset_cell<-subset(x = s_obj, subset = condition_control != "Fibroid_control")

pdf(paste0(dir_out,"bar_Condition_full_nf.pdf"),height =8, width = 8)
make_meta_barplot(subset_cell, condition_full, annotation, dittoColors())+
  scale_fill_manual(values=dittoColors()[c(3,11,4)])+
  ylab("fraction of cluster") + theme(legend.text = element_text(size=15),
                                      axis.text=element_text(size=15),)
dev.off()

pdf(paste0(dir_out,"bar_Condition_control_nf.pdf"),height =8, width = 8)
make_meta_barplot(subset_cell, condition_control, annotation, dittoColors())+
  scale_fill_manual(values=dittoColors()[c(3,4)])+
  ylab("fraction of cluster") + theme(legend.text = element_text(size=15),
                                      axis.text=element_text(size=15),)
dev.off()

pdf(paste0(dir_out,"bar_Phase_nf.pdf"),height =8, width = 8)
make_meta_barplot(subset_cell, Phase_menstrual_reduced, annotation, dittoColors())+
  scale_fill_manual(values=dittoColors()[c(2,1)])+
  ylab("fraction of cluster") + theme(legend.text = element_text(size=15),
                                      axis.text=element_text(size=15),)
dev.off()

pdf(paste0(dir_out,"bar_Library_nf.pdf"),height =8, width = 8)
make_meta_barplot(subset_cell, LIBRARY, annotation, dittoColors())+
  scale_fill_manual(values=dittoColors())+
  ylab("fraction of cluster") + theme(legend.text = element_text(size=15),
                                      axis.text=element_text(size=15),)
dev.off()

pdf(paste0(dir_out,"bar_patient_nf.pdf"),height =8, width = 8)
make_meta_barplot(subset_cell, Patient_id, annotation, dittoColors())+
  scale_fill_manual(values=dittoColors())+
  ylab("fraction of cluster") + theme(legend.text = element_text(size=15),
                                      axis.text=element_text(size=15),)
dev.off()

