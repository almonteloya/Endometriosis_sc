## Barplot Pristine vs Fibroid
library(ggplot2)
library(Seurat)
library(dittoSeq)
library(dplyr)

dir_fibroids <-"/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/DEG/fibroid_pristine/"
dir_out <- "/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/SFibroids/"
r_filter <- function(df){
  x<- df %>% dplyr::filter(avg_log2FC > .25 & p_val_adj < 0.05 )
  genesup <- paste0(rownames(x),"_UP")
  x<- df %>% dplyr::filter(avg_log2FC < -0.25 & p_val_adj < 0.05 )
  genesdown <- paste0(rownames(x),"_DOWN")
  return(c(genesup, genesdown))
}

files_fibroid_secretory <- list.files( path=dir_fibroids,pattern = "(_|/)secretory(_DEG\\.RDS)", full.names = TRUE)
names(files_fibroid_secretory)<-basename(files_fibroid_secretory)

rds_objects <- lapply(files_fibroid_secretory, readRDS)
significant <- lapply(rds_objects, r_filter)


gene_summary <- lapply(significant, function(x) {
  data.frame(
    Status = c("UP", "DOWN"),
    Count = c(sum(grepl("_UP$", x)), sum(grepl("_DOWN$", x)))
  )
})

# Add list names as a column and combine
gene_summary_df <- bind_rows(
  lapply(names(gene_summary), function(name) {
    df <- gene_summary[[name]]
    df$List <- name
    df
  })
)

# Reorder columns
gene_summary_df <- gene_summary_df %>% select(List, Status, Count)

# Use sub to replace the last column values
gene_summary_df$List <- sub("_secretory.*", "", gene_summary_df$List)



sec= ggplot(gene_summary_df, aes(x = List, y = Count, fill = Status)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("UP" = "red", "DOWN" = "blue")) +
  theme_minimal() +
  labs(title = "DEG genes in secretory", y = "Gene Count", x="Cell type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12))


pdf(paste0(dir_out,"Sec_barplot.pdf"))
sec
dev.off()

#### ------Now proliferative

files_fibroid_proliferative <- list.files( path=dir_fibroids,pattern = "(_|/)proliferative(_DEG\\.RDS)", full.names = TRUE)
names(files_fibroid_proliferative)<-basename(files_fibroid_proliferative)

rds_objects <- lapply(files_fibroid_proliferative, readRDS)
significant <- lapply(rds_objects, r_filter)


gene_summary <- lapply(significant, function(x) {
  data.frame(
    Status = c("UP", "DOWN"),
    Count = c(sum(grepl("_UP$", x)), sum(grepl("_DOWN$", x)))
  )
})

# Add list names as a column and combine
gene_summary_df <- bind_rows(
  lapply(names(gene_summary), function(name) {
    df <- gene_summary[[name]]
    df$List <- name
    df
  })
)

# Reorder columns
gene_summary_df <- gene_summary_df %>% select(List, Status, Count)

# Use sub to replace the last column values
gene_summary_df$List <- sub("_proliferative.*", "", gene_summary_df$List)



pro = ggplot(gene_summary_df, aes(x = List, y = Count, fill = Status)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("UP" = "red", "DOWN" = "blue")) +
  theme_minimal() +
  labs(title = "DEG genes in proliferative", y = "Gene Count", x="Cell type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12))



pdf(paste0(dir_out,"Pro_barplot.pdf"))
pro
dev.off()
