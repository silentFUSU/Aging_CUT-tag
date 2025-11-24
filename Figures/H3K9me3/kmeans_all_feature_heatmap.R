rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)
library(rtracklayer)
library(gridExtra)
library(grid)  
library(GenomicRanges)
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggsignif)
options(bitmapType="cairo")  
directory_path <- "/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/H3K9me3/recursion_peaks_diff_table/feature_table_combined_chrY/"
csv_files <- list.files(path = directory_path, pattern = "\\.csv$", full.names = TRUE)
data_list <- lapply(csv_files, function(file){
  data <- read.csv(file)
  names(data)[2] <- "cluster"
  data <- data[,-1]
  return(data)
})
names(data_list) <- tools::file_path_sans_ext(basename(csv_files))
merged_data <- Reduce(function(x, y) merge(x, y, by = "cluster", all = TRUE), data_list)
colnames(merged_data)[-1] <- names(data_list)

to_plot <- merged_data
rownames(to_plot) <- to_plot$cluster
to_plot <- to_plot[,-1]
colnames(to_plot) <- c("H3K27me3 signal","H3K9me3 signal","peak length","compartmentB","DNA methylation","ERV1","ERVK","LTR","phastCons score")
to_plot <- to_plot[c("kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"),
                   c("H3K9me3 signal","H3K27me3 signal","DNA methylation","compartmentB","phastCons score","LTR","ERV1","ERVK")]
to_plot <- as.data.frame(t(to_plot))
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
pheatmap::pheatmap(to_plot,scale="row",cluster_rows = F,cluster_cols = F,color = color_palette,filename = "result/figures/H3K9me3_all_feature.pdf",width = 6,height = 4)

color <- setNames(c("#e64b35","#4dbbd5","#00a087","#3c5488","#f39b7f","grey","black"),c("kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"))
to_plot <- as.data.frame(t(to_plot))
to_plot$cluster <- rownames(to_plot)
summary <- to_plot
p_list <- list()
for(i in c(1:8)){
  to_plot <- summary[,c(i,9)]  
  colnames(to_plot)[1] <- "feature"
  to_plot$feature <- as.numeric(to_plot)
  if(i == 1){
    p_list[[i]] <- ggplot(to_plot, aes(x = cluster, y = feature, fill=cluster)) +
      geom_violin(color = "black") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12,face = "bold", color = "black"),  
        axis.text.y = element_text(size = 12,face = "bold", color = "black"),  
        axis.title.x = element_text(size = 14,face = "bold", color = "black"), 
        axis.title.y = element_text(size = 14,face = "bold", color = "black"), 
        legend.text = element_text(size = 12),  
        panel.background = element_blank(),  
        panel.grid.major = element_line(size = 0.1, linetype = 'solid', color = "grey"),
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid', color = "lightgrey"),
        panel.border = element_rect(color = "black", fill = NA, size = 1) 
      ) + 
      xlim(20,100)+
      labs(x = "DNA Methylation", y= NULL)
  }
}

