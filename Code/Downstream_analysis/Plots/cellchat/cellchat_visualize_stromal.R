## Cellchat
### Secretory exploration 
## Selecting pathways 

library(dplyr)
library(CellChat)


CellChatDB.use <- CellChatDB.human


read_process <- function(x){
  print(paste0(dir_pristine, x,"_cellchat_obj.RDS"))
  cellchat <- readRDS(paste0(dir_pristine, x,"_cellchat_obj.RDS"))
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat) 
  ## Extra analysys
  cellchat <- netAnalysis_computeCentrality(cellchat)
  cellchat
}


dir_out="/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/Fig6/"

## Reading the files
dir_pristine="/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/cellchat/"
condi<-c("disease_secretory","control_sec")
cellchat_objs <- lapply(condi,read_process)
## fixing a problem with the columns
common_names <- c("labels","Condition","Patient_id", "Phase")
colnames(cellchat_objs[[2]]@meta) <- common_names
colnames(cellchat_objs[[1]]@meta) <- common_names
names(cellchat_objs) <-condi 


cellchat<- mergeCellChat(cellchat_objs, add.names = condi, cell.prefix = TRUE)


pdf(paste0(dir_out,"netVisual_diffInteraction_Secretory.pdf"))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",top=.10, 
                          comparison = c(2,1), arrow.width=8)
dev.off()


## Getting the data
communication_data_specific <- subsetCommunication(cellchat,thresh=.05,slot.name = "netP")
d_s<-communication_data_specific$disease_secretory
c_s<-communication_data_specific$control_sec


### Plotting stromal
st_c<-c_s %>% filter (source == "Smooth_muscle") 
st_c$conditon<-"Control"

st_d<-d_s %>% filter (source == "Smooth_muscle") 
st_d$conditon<-"Disease"

d_s<- data.frame(rbind(st_c,st_d))

#### 

#Merging both dataframes by Target and Pathway
merged_df <- merge(st_d, st_c, by = c('target', 'pathway_name'), suffixes = c('_df1', '_df2'), all = TRUE)

# Replace NA with 0 for the difference calculation
merged_df[is.na(merged_df)] <- 0

## Addding a pseudocount 
merged_df$prob_df1_pseudo<-merged_df$prob_df1 + .000001
merged_df$prob_df2_pseudo<-merged_df$prob_df2 + .000001


merged_df$Probability_Log2_Ratio <- log2(merged_df$prob_df1_pseudo / merged_df$prob_df2_pseudo)


bplot <- ggplot(merged_df, aes(x = target, y = pathway_name, color =Probability_Log2_Ratio)) +
  geom_point(alpha = 0.7,size=4) +  # Adjust the transparency of the bubbles
  #scale_color_gradient(low = "lightblue", high = "red") +  # Gradient color from blue to red
  labs(x = "Target", y = "Pathway Name", color = "Log2 Strenght") +  # Correct usage of labs for axis and legend titles
  theme_minimal() +  # Use a minimal theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis text
    axis.text.y = element_text(size = 10),  # Adjust y-axis text size
    plot.title = element_text(hjust = 0.5)) +scale_color_gradient2(low="blue", mid="white", high="red")


pdf(paste0(dir_out,"SmoothM_sec_1.pdf"),height = 5,width = 4.5)
bplot
dev.off()


### Adding pathway annotation
Com <- subsetCommunication(cellchat,thresh=.05,slot.name = "net")
path_annotation<- data.frame(pathway_name=Com$disease_secretory$pathway_name,
                             annotation=Com$disease_secretory$annotation)
path_annotation2<- data.frame(pathway_name=Com$control_sec$pathway_name,
                              annotation=Com$control_sec$annotation)

path_anno<-rbind(path_annotation,path_annotation2)
path_anno<-path_anno[!duplicated(path_anno), ]

merged_df2<-left_join(merged_df,
                      path_anno,
                      by = c("pathway_name"))

bplot <- ggplot(merged_df2, aes(x = target, y = pathway_name, color =Probability_Log2_Ratio)) +
  geom_point(alpha = 0.7,size=4) +  # Adjust the transparency of the bubbles
  labs(x = "Target", y = "Pathway Name", size = "Frequency", color = "Log2FC Prob.") +  # Correct usage of labs for axis and legend titles
  theme_minimal() +  # Use a minimal theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis text
    axis.text.y = element_text(size = 10),  # Adjust y-axis text size
    plot.title = element_text(hjust = 0.5)) +scale_color_gradient2(midpoint=0, low="blue", mid="white", high="red") +
  facet_grid(annotation ~ ., scales = "free", space = "free")


pdf(paste0(dir_out,"SmoothM_sec_2.pdf"),height = 6,width = 4.5)
bplot
dev.off()


### ----- Specific L and R pathways -- Not used

VB<-netVisual_bubble(cellchat, sources.use = "Smooth_muscle",comparison = c(2,1), angle.x = 90,
                     thresh=0.05,return.data = T)


VB_df<-VB$communication
VB_df$target_dataset <- paste(VB_df$target,VB_df$dataset)

VB_df<-VB_df[VB_df$pathway_name %in% "COLLAGEN",]

VB_df<-VB_df[VB_df$pathway_name %in% c("LAMININ","COLLAGEN","FN1","THBS"),]


bplot <- ggplot(VB_df, aes(x = target_dataset, y = interaction_name_2, color =prob)) +
  geom_point(alpha = 0.7,size=3) +  # Adjust the transparency of the bubbles
  #scale_color_gradient(low = "lightblue", high = "red") +  # Gradient color from blue to red
  labs(x = "Target", y = "Pathway Name", color = "Prob.") +  # Correct usage of labs for axis and legend titles
  theme_minimal() +  # Use a minimal theme
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis text
    axis.text.y = element_text(size = 10),  # Adjust y-axis text size
    plot.title = element_text(hjust = 0.5)) +
  facet_grid( .~ target, scales = "free_x") 

#pdf(paste0(dir_out,"Innate_Lymphoid_LR2.pdf"),height = 8,width = 8)
bplot
dev.off()



