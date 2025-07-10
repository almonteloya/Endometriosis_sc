## Markers based on literature

dir_out <- 'Fig5/'
dir2='seurat/object/immune'

s_obj = readRDS(file=paste0(dir2,"Immune_obj_2025.RDS"))


#### markers for this based on literature

tcell<-c("CD3D","CD3G","CD3E")
cd8_tcell <-c( "CD8A")
T_reg<-c("CD2","CD4","FOXP3")
monocyte <-c("CD68", "MS4A4A", "MS4A7")
macrophages<-c("CD14", "CD16", "CCR2")
cdc1<-c("HLA−DRA", "HLA−DQA1","CLEC9A")
innate_lymphoid<-c("ALDOC","IL4I1")
lymphatic<-c("PECAM1", "PDPN")
mast<-c("CD117", "TPSB2", "TPSAB1")
uNK <- c( "KLRC1","NCAM1", "FCGRT","CD3D" , "CD3G")
PDC<-c( "IL3RA", "CLEC4C", "LILRA4","RIMKLB","CLDN4",
        "TM4SF1", "SLPI", "SCGB2A1", "WFDC")
PlasmaB<-c("JCHAIN","IL3RA")
Allimm<-c("PTPRC","CD74")


immune_markers<-c(Allimm,PlasmaB,macrophages,monocyte,tcell,T_reg,cd8_tcell,lymphatic,PDC,
                  innate_lymphoid,uNK,cdc1,mast)
Idents(s_obj)<-s_obj$immune_agregate1_c

pdf(paste0(dir_out,"immune_markers.pdf"),height =4, width=10)
DotPlot( s_obj, features=unique(immune_markers),cluster.idents=T,cols = c("yellow", "red")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
