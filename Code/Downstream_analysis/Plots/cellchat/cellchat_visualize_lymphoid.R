## Cellchat
### pristine controls vs endometriosis 
## Proliferative  
## Net differentila visualization
## Innate lymphoid visualization 

library(dplyr)
library(CellChat)


CellChatDB.use <- CellChatDB.human
read_process <- function(x){
  ## Read objects and finish some extra processing 
  print(paste0(dir_pristine, x,"_cellchat_obj.RDS"))
  cellchat <- readRDS(paste0(dir_pristine, x,"_cellchat_obj.RDS"))
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat) 
  ## Extra analysys
  cellchat <- netAnalysis_computeCentrality(cellchat)
  cellchat
}


dir_out="/out/dir"
dir_pristine="dir/with/cellchat/objects"

condi<-c("disease_proliferative","control_pro")
cellchat_objs <- lapply(condi,read_process)
## fixing a problem with the columns
common_names <- c("labels","Condition","Patient_id", "Phase")
colnames(cellchat_objs[[2]]@meta) <- common_names
colnames(cellchat_objs[[1]]@meta) <- common_names
names(cellchat_objs) <-condi 

## Merge
cellchat<- mergeCellChat(cellchat_objs, add.names = condi, cell.prefix = TRUE)


pdf(paste0(dir_out,"netVisual_diffInteraction_Proliferative.pdf"))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",top=.10, 
                          comparison = c(2,1), arrow.width=8)
dev.off()
## Main trend: Downregulation of Innate lymphoid

## Getting the data for ploting specific pathways
communication_data_specific <- subsetCommunication(cellchat,thresh=0.05,slot.name = "netP")
d_s<-communication_data_specific$disease_proliferative
c_s<-communication_data_specific$control_pro



### Plotting Innate lymphoid as a source and cheking pathwys
## Control
st_c<-c_s %>% filter (source == "Innate lymphoid") 
st_c$conditon<-"Control"
##Disease
st_d<-d_s %>% filter (source == "Innate lymphoid") 
st_d$conditon<-"Disease"

d_s<- data.frame(rbind(st_c,st_d))


#### 

#Merging both dataframes by Target and Pathway
merged_df <- merge(st_d, st_c, by = c('target', 'pathway_name'), suffixes = c('_df1', '_df2'), all = TRUE)

# Replace NA with 0 for the missinng pathwys
merged_df[is.na(merged_df)] <- 0

## Addding a pseudocount 
merged_df$prob_df1_pseudo<-merged_df$prob_df1 + .000001
merged_df$prob_df2_pseudo<-merged_df$prob_df2 + .000001

## Calculate logF change
merged_df$Probability_Log2_Ratio <- log2(merged_df$prob_df1_pseudo / merged_df$prob_df2_pseudo)


bplot <- ggplot(merged_df, aes(x = target, y = pathway_name, color =Probability_Log2_Ratio)) +
  geom_point(alpha = 0.7,size=4) +  # Adjust the transparency of the bubbles
  #scale_color_gradient(low = "lightblue", high = "red") +  # Gradient color from blue to red
  labs(x = "Target", y = "Pathway Name", color = "Log2FC Prob") +  # Correct usage of labs for axis and legend titles
  theme_minimal() +  # Use a minimal theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis text
    axis.text.y = element_text(size = 10),  # Adjust y-axis text size
    plot.title = element_text(hjust = 0.5)) +scale_color_gradient2(low="blue", mid="white", high="red")

bplot
pdf(paste0(dir_out,"Innate_Lymphoid_pro_1.pdf"),height = 5,width = 4.5)
bplot
dev.off()

### Adding the pathway classification

Com <- subsetCommunication(cellchat,thresh=.05,slot.name = "net")
path_annotation<- data.frame(pathway_name=Com$disease_proliferative$pathway_name,
                             annotation=Com$disease_proliferative$annotation)
path_annotation2<- data.frame(pathway_name=Com$control_pro$pathway_name,
                              annotation=Com$control_pro$annotation)

path_anno<-rbind(path_annotation,path_annotation2)
path_anno<-path_anno[!duplicated(path_anno), ]

merged_df2<-left_join(merged_df,
                      path_anno,
                      by = c("pathway_name"))

bplot <- ggplot(merged_df2, aes(x = target, y = pathway_name, color =Probability_Log2_Ratio)) +
  geom_point(alpha = 0.7,size=4) +  # Adjust the transparency of the bubbles
  labs(x = "Target", y = "Pathway Name", color = "Log2FC Prob.") +  # Correct usage of labs for axis and legend titles
  theme_minimal() +  # Use a minimal theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis text
    axis.text.y = element_text(size = 10),  # Adjust y-axis text size
    plot.title = element_text(hjust = 0.5)) +scale_color_gradient2(midpoint=0, low="blue", mid="white", high="red") +
  facet_grid(annotation ~ ., scales = "free", space = "free")


pdf(paste0(dir_out,"Innate_Lymphoid_pro_2.pdf"),height = 5,width = 4.5)
bplot
dev.off()


####------ L-R interactions

VB<-netVisual_bubble(cellchat, sources.use = "Innate lymphoid",comparison = c(2,1), angle.x = 90,
                 thresh=0.05,return.data = T)


pdf(paste0(dir_out,"Innate_Lymphoid_LR1.pdf"),height = 8,width = 8)
VB$gg.obj
dev.off()


VB_df<-VB$communication
VB_df$target_dataset <- paste(VB_df$target,VB_df$dataset)


bplot <- ggplot(VB_df, aes(x = target_dataset, y = interaction_name_2, color =prob)) +
  geom_point(alpha = 0.7,size=3) +  # Adjust the transparency of the bubbles
  #scale_color_gradient(low = "lightblue", high = "red") +  # Gradient color from blue to red
  labs(x = "Target", y = "Pathway Name", color = "Prob.") +  # Correct usage of labs for axis and legend titles
  theme_minimal() +  # Use a minimal theme
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis text
    axis.text.y = element_text(size = 10),  # Adjust y-axis text size
    plot.title = element_text(hjust = 0.5)) +
  facet_grid(. ~ target, scales = "free_x", space = "free_x") 

pdf(paste0(dir_out,"Innate_Lymphoid_LR2.pdf"),height = 8,width = 8)
  bplot
dev.off()



VB_df_filtered<-VB_df[VB_df$pathway_name=="TNF",]
