###----- Density plots differences between LFC 

library(dplyr)
library(ggplot2)
library(ggrepel)


dir <- "/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/severity/"
dir2<- "/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/pristine_disease/broad/"
sev <- list.files( path=dir,pattern = ".RDS", full.names = TRUE)
all <-list.files( path=dir2,pattern = ".RDS", full.names = TRUE)

markers_sev <- lapply(sev, readRDS)
markers_all <- lapply(all, readRDS)


filter_genes <- function(df){
  x <- df 
  x <- x %>% mutate(Significant = ifelse((avg_log2FC < -0.25 | avg_log2FC > 0.25) & p_val_adj < 0.05, "Significant", "Non-Significant"))
  x<- tibble::rownames_to_column(x, "Gene")
  x<- dplyr::select(x, c(Gene,avg_log2FC,p_val_adj,Significant))
  x
}

significant_sev <- lapply(markers_sev,filter_genes)
significant_all <- lapply(markers_all,filter_genes)

names(significant_sev) <- basename(sev)
names(significant_all) <- basename(all)

celltypes = c("Ciliated","Unciliated","Endothelia","fibroblast","Smooth")
phase<-c("secretory","proliferative")


pdf("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/SFseverity/correlation.pdf", height = 5)
for (i in c(1:2)) {
  for (j in c(1:5)) {
    print(phase[i])
    print(celltypes[j])
phase_sev <- significant_sev[grepl(phase[i], names(significant_sev))]
phase_sig <- significant_all[grepl(phase[i], names(significant_all))]
sev_comp <- phase_sev[grepl(celltypes[j], names(phase_sev))]
all_comp <- phase_sig[grepl(celltypes[j], names(phase_sig))]

  for (z in c(1:2)) {
Stromal_secretory<- list(sev_comp[[z]],all_comp[[1]])
merged_stromal <- Stromal_secretory %>% purrr::reduce(full_join, by='Gene') 
merged_stromal_na <- merged_stromal %>% na.exclude()
  
merged_stromal_na<- merged_stromal_na%>% mutate(Significant = ifelse(Significant.x == "Significant" & Significant.y == "Significant", "Significant Both", 
                                                                       ifelse(Significant.x == "Significant", "Significant group 1", 
                                                                              ifelse(Significant.y == "Significant", "Significant group 2", "Non-Significant"))))
correlation <- cor(merged_stromal_na$avg_log2FC.x,merged_stromal_na$avg_log2FC.y)
  # Create a ggplot with scatter plot and add the correlation value
  p<-ggplot(merged_stromal_na, aes(x = avg_log2FC.x, y = avg_log2FC.y,color=Significant)) +
    geom_point(alpha = 0.3) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    annotate("text", x = min(merged_stromal_na$avg_log2FC.x), y = max(merged_stromal_na$avg_log2FC.y), 
             label = paste("Correlation = ", round(correlation, 2)),
             hjust = 0, vjust = 1, size = 5, color = "red") +labs(y= "Stratified my severity", x = "All") +
    scale_color_manual(values = c("Significant Both" = "black", "Significant group 1" = "darkblue", "Significant group 2" = "darkblue", "Non-Significant" = "gray")) +
    theme_minimal()  + ggtitle(paste(names(sev_comp)[z], "with",  names(all_comp)[1]))+
    theme(
      legend.text = element_text(size = 12),  # Adjust the size of the legend text
      legend.title = element_text(size = 14), # Adjust the size of the legend title
      legend.key.size = unit(1.5, "lines")    # Adjust the size of the legend keys
    )
  

  print(p)
}
  }
}
dev.off()

