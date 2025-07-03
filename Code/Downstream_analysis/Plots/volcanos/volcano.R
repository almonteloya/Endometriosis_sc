
### April 2024
## Making a better visualization with the volcano plot
library(ggplot2)
library(Seurat)
library(dittoSeq)
library(ggrepel)


dir_in= "/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/fibroid_pristine/"
dir_in= "/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/pristine_disease/broad/"
dir_in= "/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/menstrual_phase/"
dir_in= "/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/pristine_disease/immune_specific/"

dir_out= paste0(dir_in,"volcano/")
dir.create(dir_out)


files <- list.files( path=dir_in,pattern = ".RDS", full.names = TRUE)
files
markers <- lapply(files, readRDS)



pattern <- ".*\\/([^/]+)\\.RDS$"
names_file <- sub(pattern, "\\1", files)
names_file <- sub("_DEG$", "", names_file)

names(markers)<- names_file




DEG_plot<-function(cell_type){
  DEG <-markers[[cell_type]]
  
  #### VOLCANO PLOT
  
  pdf(paste0(dir_out,cell_type,"_volcano.pdf"))
  #https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html
  
  DEG$diffexpressed <- "NO"
  DEG$diffexpressed[DEG$avg_log2FC > 0.25 & DEG$p_val_adj < 0.05] <- "UP"
  DEG$diffexpressed[DEG$avg_log2FC < -0.25 & DEG$p_val_adj < 0.05] <- "DOWN"
  
  # Convert diffexpressed to a factor and set levels
  DEG$diffexpressed <- factor(DEG$diffexpressed, levels = c("DOWN", "NO", "UP"))
  
  # Create a column for labels
  DEG$delabel <- NA
  DEG$delabel[DEG$diffexpressed != "NO"] <- rownames(DEG)[DEG$diffexpressed != "NO"]
  
  # Create ggplot
  p <- ggplot(data = DEG, aes(x = avg_log2FC, y = -log10(p_val), col = diffexpressed, label = delabel)) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    geom_text_repel(force = 0.5, size = 6) +
    scale_color_manual(values = c("DOWN" = "blue", "NO" = "gray", "UP" = "red")) +
    ggtitle(cell_type) +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 12)) +
    #labs(x = "Down in Endometriosis     Avg Log2FC    Up in Endometriosis ") 
    labs(x = "Proliferative     Avg Log2FC       Secretory") 
 
  
  # Print the plot
  print(p)
  
  # Close the device
  dev.off()
}

lapply(names_file, DEG_plot)




