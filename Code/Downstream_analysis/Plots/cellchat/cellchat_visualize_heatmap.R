### comparing pristine cellchat
## Plotting based on cellchat output

library(Seurat)
library(CellChat)


## Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human # loading the data
showDatabaseCategory(CellChatDB)
## we can explore what the database has
#CellChatDB.use <- CellChatDB # simply use the default CellChatDB

read_process <- function(x){
  print(paste0(dir_obj, x,"_cellchat_obj.RDS"))
  cellchat <- readRDS(paste0(dir_out, x,"_cellchat_obj.RDS"))
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat) 
  ## Extra analysys
  cellchat <- netAnalysis_computeCentrality(cellchat)
  cellchat
}

## preparing object
dir_obj="/dir/cellchat/output"
dir_out = "out/dir"

read_cell_chat <- function(condi){
  cellchat_objs <- lapply(condi,read_process)
  ## fixing a problem with the columns
  common_names <- c("labels","Condition","Patient_id", "Phase")
  colnames(cellchat_objs[[2]]@meta) <- common_names
  colnames(cellchat_objs[[1]]@meta) <- common_names
  names(cellchat_objs) <-condi 
  cellchat_objs
}

phase <- c("disease_secretory","control_sec")
cellchat_objs_sec<-read_cell_chat(phase)
#cellchat_secretory <- mergeCellChat(cellchat_objs, add.names = phase, cell.prefix = TRUE)

phase <- c("disease_proliferative","control_pro")
cellchat_objs_pro<-read_cell_chat(phase)
# cellchat_proliferative <- mergeCellChat(cellchat_objs, add.names = phase, cell.prefix = TRUE)



#### Read all 
objects_name <- c("disease_secretory","control_sec","disease_proliferative","control_pro")
cellchat_objs_all <- read_cell_chat(objects_name)

cellchat_all <- mergeCellChat(cellchat_objs_all, add.names = objects_name , cell.prefix = TRUE)

gg1 <- compareInteractions(cellchat_all, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_all, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2


gg1 <- rankNet(cellchat_all, mode = "comparison", stacked = T, do.stat = TRUE, comparison = c(1,2),measure = 'weight', cutoff.pvalue = 0.001)$data
gg2 <- rankNet(cellchat_all, mode = "comparison", stacked = T, do.stat = TRUE, comparison = c(3,4), measure = 'weight', cutoff.pvalue = 0.001)$data
gg3 <- rankNet(cellchat_all, mode = "comparison", stacked = T, do.stat = TRUE, comparison = c(1,3), measure = 'weight', cutoff.pvalue = 0.001)$data
gg4 <- rankNet(cellchat_all, mode = "comparison", stacked = T, do.stat = TRUE, comparison = c(2,4), measure = 'weight', cutoff.pvalue = 0.001)$data

p1 <- gg1 %>% 
  group_by(name) %>% 
  mutate(fold = contribution[group == 'disease_secretory']/contribution[group == 'control_sec']) %>%
  mutate(logfold = log2(fold)) %>%
  group_by(name) %>% 
  mutate(t = sum(contribution)) %>% 
  mutate(prop1 = contribution/t) %>%
  #filter(abs(logfold) > 1) %>%
  filter(pvalues < 0.001) %>%
  ggplot(aes(x=name, y=contribution, fill = group)) + 
  geom_bar(stat="identity",width = 0.7, position ="fill") +
  coord_flip() +
  CellChat_theme_opts()


p2 <- gg2 %>% 
  group_by(name) %>% 
  mutate(fold = contribution[group == 'disease_proliferative']/contribution[group == 'control_pro']) %>%
  mutate(logfold = log2(fold)) %>%
  group_by(name) %>% 
  mutate(t = sum(contribution)) %>% 
  mutate(prop1 = contribution/t) %>%
  #filter(abs(logfold) > 1) %>%
  filter(pvalues < 0.001) %>%
  ggplot(aes(x=name, y=contribution, fill = group)) + 
  geom_bar(stat="identity",width = 0.7, position ="fill") +
  coord_flip() +
  CellChat_theme_opts()

p3 <- gg3 %>% 
  group_by(name) %>% 
  mutate(fold = contribution[group == 'disease_secretory']/contribution[group == 'disease_proliferative']) %>%
  mutate(logfold = log2(fold)) %>%
  group_by(name) %>% 
  mutate(t = sum(contribution)) %>% 
  mutate(prop1 = contribution/t) %>%
  #filter(abs(logfold) > 1) %>%
  filter(pvalues < 0.001) %>%
  ggplot(aes(x=name, y=contribution, fill = group)) + 
  geom_bar(stat="identity",width = 0.7, position ="fill") +
  coord_flip() +
  CellChat_theme_opts()

p4 <- gg4 %>% 
  group_by(name) %>% 
  mutate(fold = contribution[group == 'control_sec']/contribution[group == 'control_pro']) %>%
  mutate(logfold = log2(fold)) %>%
  group_by(name) %>% 
  mutate(t = sum(contribution)) %>% 
  mutate(prop1 = contribution/t) %>%
  #filter(abs(logfold) > 1) %>%
  filter(pvalues < 0.001) %>%
  ggplot(aes(x=name, y=contribution, fill = group)) + 
  geom_bar(stat="identity",width = 0.7, position ="fill") +
  coord_flip() +
  CellChat_theme_opts()

p1 + p2 + p3 + p4

data1 <- p1$data
data2 <- p2$data
data3 <- p3$data
data4 <- p4$data

data1$Comparison <- 'Disease sec vs Control sec'
data2$Comparison <- 'Disease pro vs Control pro'
data3$Comparison <- 'Disease sec vs Disease pro'
data4$Comparison <- 'Control sec vs control pro'
data <- do.call('rbind', list(data1, data2, data3, data4))
data <- data %>% group_by(name) %>% mutate(count = n()) %>% ungroup() %>% arrange(name) %>% arrange(desc(count)) %>% mutate(name = as.factor(name))
t <- unique(data$name)
t <- factor(t, levels = t)
# data$Comparison <- factor(data$Comparison, levels = c('RA vs control', 'High DA vs control', 'High DA vs low DA', 'Low DA vs control'))
data$name <- factor(data$name, levels = t)

data$Presence <- ifelse(data$logfold == 'Inf', 'Only in group 1', ifelse(data$logfold == '-Inf', 'Only in group 2', 'In both groups'))
data$size <- data$logfold
data$size <- ifelse(data$size == 'Inf', 2, ifelse(data$logfold == '-Inf', 2, ''))
data
toDelete <- seq(0, nrow(data), 2)
toDelete
data <-  data[-toDelete, ]
data
library(scales)



plot3 <- ggplot(data, aes(Comparison, name, fill = logfold, colour = Presence)) +
  geom_point(shape = 21, stroke = 1, size = 3) +
  scale_fill_gradientn(colours = c('blue', 'white', 'red'), values = rescale(c(-8,0,3)), guide = 'colorbar', na.value = 'white') +
  scale_colour_manual(values = c('black', "red", 'blue')) +
  theme_bw() +  
  theme(legend.key.size = unit(1, "cm"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        plot.caption = element_text(hjust = 0, face = 'bold', size = 10),
        axis.text = element_text(size = 10, face = 'bold')) +
  ylab('') +
  labs(title = 'Logfold change per group and pathway')
plot3 


ggplot(data, aes(Comparison, name, fill = logfold, colour = Presence)) +
  geom_point(shape = 21, stroke = 1, size = 3) +
  scale_fill_gradientn(colours = c('blue', 'white', 'red'), values = rescale(c(-8,0,3)), guide = 'colorbar', na.value = 'white') +
  scale_colour_manual(values = c('black', "red", 'blue')) +
  theme_bw() +  
  theme(legend.key.size = unit(1, "cm"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        plot.caption = element_text(hjust = 0, face = 'bold', size = 10),
        axis.text = element_text(size = 9),
        axis.title.y = element_text(size = 10)) +
  ylab('') +
  labs(title = 'Logfold change per group and pathway')

data1 <- data %>% 
  filter(Comparison %in% c("Disease pro vs Control pro","Disease sec vs Control sec"))

ggplot(data1, aes(Comparison, name, fill = logfold, colour = Presence)) +
  geom_point(shape = 21, stroke = 1, size = 3) +
  scale_fill_gradientn(colours = c('blue', 'white', 'red'), values = rescale(c(-8,0,3)), guide = 'colorbar', na.value = 'white') +
  scale_colour_manual(values = c('black', "red", 'blue')) +
  theme_bw() +  
  theme(legend.key.size = unit(1, "cm"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        plot.caption = element_text(hjust = 0, face = 'bold', size = 10),
        axis.text = element_text(size = 9),
        axis.title.y = element_text(size = 10)) +
  ylab('') +
  labs(title = 'Logfold change per group and pathway')

data1$logfold
library(dplyr)
data1 <- data1 %>%
  mutate(logfold2 = case_when(
    is.infinite(logfold) & logfold > 0 ~ 10,  # Check for positive infinity
    is.infinite(logfold) & logfold < 0 ~ -10, # Check for negative infinity
    TRUE ~ logfold  # Else leave it as it is
  ))


library(ggplot2)
library(scales)

ggplot(data1, aes(Comparison, name, fill = logfold2)) +
  geom_point(shape = 21, stroke = 1, size = 3) +
  #scale_fill_gradientn(colours = c('blue', 'white', 'red'), values = rescale(c(-8, 0, 3)), guide = 'colorbar', na.value = 'white') +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, na.value = 'white') +
  theme_bw() +
  theme(
    legend.key.size = unit(1, "cm"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.caption = element_text(hjust = 0, face = 'bold', size = 10),
    axis.text = element_text(size = 9),
    axis.title.y = element_text(size = 10)
  ) +
  ylab('') +
  labs(title = 'Logfold change per group and pathway')