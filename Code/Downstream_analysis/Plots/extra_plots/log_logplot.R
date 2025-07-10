### Log log plot 
## Stromal cells

library (dplyr)
library(ggplot2)
library(ggrepel)


data1 <- readRDS("Path") ### Proliferative DEG list
data2 <- readRDS("Path")## Secretory DEG list
celltype<-"Unciliated_epithelia" ## Or stromal fibroblast


merged_df <- merge(data1, data2, by =0,all=T)


merged_df$gene <- merged_df$Row.names

merged_df <- merged_df %>%
  mutate(significance = case_when(
    abs(avg_log2FC.x) > 0.25 & p_val_adj.x < 0.05 & (is.na(p_val_adj.y)) ~ "Unique_proliferative",
    abs(avg_log2FC.y) > 0.25 & p_val_adj.y < 0.05 & (is.na(p_val_adj.x)) ~ "Unique_secretory",
    abs(avg_log2FC.x) > 0.25 & p_val_adj.x < 0.05 & (abs(avg_log2FC.y) < 0.25 | p_val_adj.y > 0.05) ~ "Unique_proliferative",
    abs(avg_log2FC.y) > 0.25 & p_val_adj.y < 0.05 & (abs(avg_log2FC.x) < 0.25 | p_val_adj.x >= 0.05) ~ "Unique_secretory",
    (avg_log2FC.x) > 0.25 & p_val_adj.x < 0.05 & (avg_log2FC.y) > 0.25 & p_val_adj.y < 0.05 ~ "Upregulated_in_endometriosis",
    (avg_log2FC.x) < -0.25 & p_val_adj.x < 0.05 & (avg_log2FC.y) < -0.25 & p_val_adj.y < 0.05 ~ "Downregulted_in_endometriosis",
    (avg_log2FC.x)  > 0.25 & p_val_adj.x < 0.05 & (avg_log2FC.y) < -0.25 & p_val_adj.y < 0.05 ~ "Opposite_direction",
    (avg_log2FC.x)  < -0.25 & p_val_adj.x < 0.05 & (avg_log2FC.y) > 0.25 & p_val_adj.y < 0.05 ~ "Opposite_direction",
    TRUE ~ "None"
  ))


library(ggplot2)
library(ggrepel)

# Ensure `significance` is a factor
merged_df$significance <- factor(merged_df$significance, levels = c("None", "Opposite_direction", "Upregulated_in_endometriosis", "Downregulted_in_endometriosis", "Unique_proliferative", "Unique_secretory"))

# Define colors for each category
color_map <- c(None = "gray", Opposite_direction = "purple", Upregulated_in_endometriosis = "red", Downregulted_in_endometriosis = "blue", 
               Unique_proliferative = "#56B4E9", Unique_secretory = "darkorange")

# Create the legend labels
legend_plot <- as.data.frame(table(merged_df$significance))
legend_labels <- paste0(legend_plot$Var1, " (n=", legend_plot$Freq, ")")

# Ensure that the names of the colors match the levels of the factor
names(color_map) <- levels(merged_df$significance)

## Chnage gene_labels:

merged_df <- merged_df %>%
  mutate(gene_label = ifelse(significance %in% c("Upregulated_in_endometriosis", "Downregulted_in_endometriosis"), gene, NA))

# Sort data by the sum of absolute values of `avg_log2FC.x` and `avg_log2FC.y`
sorted_df <- merged_df %>%
  mutate(abs_sum_log2FC = abs(avg_log2FC.x) + abs(avg_log2FC.y)) %>%
  arrange(desc(abs_sum_log2FC))

# Plot the data
pdf("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/Fig4/log_log.pdf",width =8)
#pdf("/krummellab/data1/immunox/XREP1a/10x/merged_XREP/Fix2025/Fig3/log_log.pdf",width =8)
ggplot(merged_df, aes(x = avg_log2FC.x, y = avg_log2FC.y, color = significance)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(labels = legend_labels, values = color_map) +
  labs(
    title = celltype,
    x = "Log Fold Change Proliferative",
    y = "Log Fold Change Secretory",
    color = "Significance"
  ) +
  theme_minimal() +
  geom_text_repel(data = subset(sorted_df, significance %in% c("Upregulated_in_endometriosis", "Downregulted_in_endometriosis")), 
                  aes(label = gene_label), size = 5, box.padding = 0.2, point.padding = 0.2, segment.color = 'grey50')
dev.off()