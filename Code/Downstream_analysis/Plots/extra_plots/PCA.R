##PCA of cell compositions
## Because of how the data was processed take a different approach:
## We take the fraction of each cell type across every major celltype. 

library(dplyr)

dir_out<-"/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/Fig2/"

increase_theme <- theme(
  axis.text.x = element_text(size = 14),    
  axis.text.y = element_text(size = 14),    
  axis.title.x = element_text(size = 16),  
  axis.title.y = element_text(size = 16),   
  plot.title = element_text(size = 20),     
  legend.title = element_text(size = 18),   
  legend.text = element_text(size = 18)     
)

percentage_generate<-function(metadata){
  # Step 1: Count the occurrences of each celltype per Patient
  df_counts <- metadata %>%
    group_by(Patient_id, annotation) %>%
    summarise(count = n(), .groups = 'drop')
  
  # Step 2: Calculate total count of cells per Patient
  df_totals <- df_counts %>%
    group_by(Patient_id) %>%
    summarise(total_count = sum(count), .groups = 'drop')
  
  # Step 3: Merge the counts with the total counts and calculate the percentage
  df_result <- df_counts %>%
    left_join(df_totals, by = "Patient_id") %>%
    mutate(percentage = (count / total_count))
  
  result <- df_result %>% select(Patient_id,annotation,percentage) %>%
    pivot_wider(names_from = Patient_id, values_from = percentage)
  
  result[is.na(result)] <- 0
  
  as.data.frame(result)
}


percentage_generate_imm<-function(metadata){
  # Step 1: Count the occurrences of each celltype per Patient
  df_counts <- metadata %>%
    group_by(Patient_id, immune_agregate1_c) %>%
    summarise(count = n(), .groups = 'drop')
  
  # Step 2: Calculate total count of cells per Patient
  df_totals <- df_counts %>%
    group_by(Patient_id) %>%
    summarise(total_count = sum(count), .groups = 'drop')
  
  # Step 3: Merge the counts with the total counts and calculate the percentage
  df_result <- df_counts %>%
    left_join(df_totals, by = "Patient_id") %>%
    mutate(percentage = (count / total_count))
  
  # Format the result
  result<-df_result  %>% select(Patient_id,immune_agregate1_c,percentage) %>%
    pivot_wider(names_from = Patient_id, values_from = percentage)
  
  result[is.na(result)] <- 0
  
  as.data.frame(result)
}

## STROMAL
s_obj <- readRDS('/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/Stromal_obj_2025.RDS')
metadata_stromal<- s_obj@meta.data
stromal_per<-percentage_generate(metadata_stromal)


summary_df <- metadata_stromal %>%
  select(Patient_id, Phase_menstrual_reduced, condition_full_reduced) %>%
  distinct()

summary_df <- summary_df %>%
  mutate(
    condition_full_reduced = case_when(
      condition_full_reduced == "Pristine_control" ~ "Pristine controls",
      condition_full_reduced == "Pathological_control" ~ "Fibroid controls",
      condition_full_reduced == "Disease" ~ "Endometriosis",
      TRUE ~ condition_full_reduced  # Retain original value if it doesn't match any of the conditions
    )
  )


## Epithelial
s_obj <- readRDS('/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/Epithelial_obj_2025.RDS')
  metadata_epi <- s_obj@meta.data
epi_per<-percentage_generate(metadata_epi)

## Immune
s_obj <- readRDS('/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/Immune_obj_2025.RDS')
metadata_immune <- s_obj@meta.data
imm_per<-percentage_generate_imm(metadata_immune)
colnames(imm_per)[1]<-"annotation"

### Rbinding everything 
total_counts<-rbind(epi_per,stromal_per,imm_per)
rownames(total_counts) <- total_counts$annotation
total_counts$annotation<-NULL

pca_result <- prcomp(t(total_counts), scale. = TRUE)
summary_pca <- summary(pca_result)
var1<-summary_pca$importance[2,][1]
var1 <- sprintf("(%.2f%%)", var1 * 100)
var2<-summary_pca$importance[2,][2]
var2 <- sprintf("(%.2f%%)", var2 * 100)

library(ggplot2)

dtp <- data.frame('Patient_id' = colnames(total_counts), pca_result$x[,1:2]) 


dtp<-merge(summary_df,dtp,by="Patient_id")

pdf(paste0(dir_out,"PCA.pdf")) 
ggplot(data = dtp) + 
  geom_point(aes(x = PC1, y = PC2, col = condition_full_reduced),size=3) + 
  theme_minimal() +
  scale_color_manual(values=dittoColors()[c(3,7,4)])  + 
  labs(title = "Condition\n", color = "Condition\n", x=paste("PC1",var1),y=paste("PC2",var2)) + increase_theme

ggplot(data = dtp) + 
  geom_point(aes(x = PC1, y = PC2, col = Phase_menstrual_reduced),size=3) + 
  theme_minimal()  + 
  scale_color_manual(values=dittoColors()[c(2,1)]) +
  labs(title = "Phase\n", color = "Phase\n", x=paste("PC1",var1),y=paste("PC2",var2)) + increase_theme

dev.off()
