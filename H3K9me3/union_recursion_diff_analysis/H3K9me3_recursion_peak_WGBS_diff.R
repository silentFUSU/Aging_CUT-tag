rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
library(ggrepel)
library(grid)
library(ggsignif)
library(data.table)
tissue <- "mammarygland"
tissue_label_change <- function(tissue){
  if(tissue=="brain"){
    tissue_label <- "Cortex"
  }else if(tissue == "Hip"){
    tissue_label <- "Hippocampus"
  }else if(tissue == "CB"){
    tissue_label <- "Cerebellum"
  }else{
    tissue_label <- str_to_title(tissue)
    if(tissue_label == "Bonemarrow"){
      tissue_label <- "Bone Marrow"
    }else if(tissue_label == "Bat"){
      tissue_label <- "BAT"
    }else if(tissue_label=="Mammarygland"){
      tissue_label <- "Mammary Gland"
    }else if(tissue_label=="Iwat"){
      tissue_label <- "IWAT"
    }
  }
  return(tissue_label)
}
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
                  "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT","ileum")) 
tissue_summary <- data.frame()
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  t_search_table <- search_table[which(search_table$tissue==tissue),]
  H3K9me3_peaks <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation.csv")
  split_chr <- strsplit(as.character(H3K9me3_peaks$X), ":")  
  chr_column <- sapply(split_chr, `[[`, 1)  
  split_start_end <- strsplit(sapply(split_chr, `[[`, 2), "-")  
  start_column <- sapply(split_start_end, `[[`, 1)  
  end_column <- sapply(split_start_end, `[[`, 2)  
  H3K9me3_peaks <- data.frame(chr = chr_column,start = start_column, end = end_column, cluster=H3K9me3_peaks$cluster)
  if(tissue %in% c("ovary","mammarygland","uterus")){
    H3K9me3_peaks <- H3K9me3_peaks[-which(H3K9me3_peaks$chr == "chrY"),]
  }
  H3K9me3_peaks$start <- as.numeric(H3K9me3_peaks$start)
  H3K9me3_peaks$end <- as.numeric(H3K9me3_peaks$end)
  H3K9me3_peaks$label <- paste0(H3K9me3_peaks$chr,":",H3K9me3_peaks$start,"-",H3K9me3_peaks$end)
  H3K9me3_peaks <- as.data.table(H3K9me3_peaks)
  setDT(H3K9me3_peaks)
  setkey(H3K9me3_peaks,chr,start,end)
  summary <- data.frame()
  for(sample in t_search_table$sample_name){
    df <- fread(paste0("data/samples/WGBS/",tissue,"/bdg/",sample,"_CpG.bdg"),sep = "\t")
    setDT(df)
    setkey(df,V1,V2,V3)  
    overlaps <- foverlaps(df, H3K9me3_peaks, type = "any", nomatch = 0L)  
    
    result <- overlaps[, .(V4_sum = sum(V4), V5_sum = sum(V5)), by = label]
    result <- as.data.frame(result)
    result$methylation <- result$V4_sum/result$V5_sum
    result <- result[,c("label","methylation")]
    colnames(result)[2] <- sample
    if(nrow(summary)==0){
      summary <- result
    }else{
      summary <- merge(summary,result,by="label")
    }
  }
  young_summary <- summary[,c("label",t_search_table$sample_name[which(t_search_table$age=="3M")])]
  old_summary <- summary[,c("label",t_search_table$sample_name[which(t_search_table$age=="24M")])]
  
  young_summary$young_methylation <- rowMeans(young_summary[,-1])
  old_summary$old_methylation <- rowMeans(old_summary[,-1])
  t_tissue_summary <- merge(young_summary,old_summary,by="label")
  t_tissue_summary$delta <- log2(t_tissue_summary$old_methylation/t_tissue_summary$young_methylation)
  t_tissue_summary <- t_tissue_summary[,c("label","delta")]
  colnames(t_tissue_summary)[2] <- tissue_label_change(tissue)
  if(nrow(tissue_summary)==0){
    tissue_summary <- t_tissue_summary
  }else{
    tissue_summary <- merge(tissue_summary,t_tissue_summary,by="label",all=T)
  }
}
# write.csv(tissue_summary,"data/samples/WGBS/all_tissues_delta_in_H3K9me3_recursion_peaks.csv")
annotation <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation.csv")
colnames(annotation)[1]<-"label"
to_plot <- merge(tissue_summary,annotation,by="label")
rownames(to_plot) <- to_plot$label
to_plot <- to_plot[order(to_plot$cluster),]
to_plot <- to_plot[,-which(colnames(to_plot)%in% c("label","cluster"))]
rownames(annotation) <- annotation$label
annotation <- annotation[,-1,drop =F]
annotation$cluster <- as.character(annotation$cluster)
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- c(seq(-0.5, -0.11, length.out = 40), seq(-0.1, 0.1, length.out = 20), seq(0.11, 0.5, length.out = 40))
tissue_order <- c("Kidney","Muscle","Skin","Bladder","Stomach","Heart","Hippocampus","Uterus","Liver","Aorta","Testis","Cortex","Tongue","Cerebellum","BAT","Lung","Mammary Gland","Pancreas","Bone Marrow","iWAT","Cecum","Colon","Jejunum","Spleen","Thymus","Ileum","Ovary")
to_plot <- to_plot[,tissue_order]
pheatmap::pheatmap(to_plot,cluster_rows = F,cluster_cols = F,annotation_row = annotation,breaks = breaks,color = color_palette,show_rownames = F,main = "CpG methylation")
H3K9me3_tissue_order <- c("Lung","Cerebellum","BAT","Muscle","Heart","Aorta","Skin","Kidney","Hippocampus","Cortex",
                          "Liver","Tongue","Bladder","Pancreas","Cecum","Spleen","Stomach","Colon","Bone Marrow","Jejunum",
                          "iWAT","Thymus","Ileum")
to_plot_H3K9me3_order <- to_plot[,H3K9me3_tissue_order]
pheatmap::pheatmap(to_plot_H3K9me3_order,cluster_rows = F,cluster_cols = F,annotation_row = annotation,breaks = breaks,color = color_palette,show_rownames = F,main = "CpG methylation")

to_plot_H3K9me3_order_long <- to_plot_H3K9me3_order
to_plot_H3K9me3_order_long$label <- rownames(to_plot_H3K9me3_order_long)
to_plot_H3K9me3_order_long <- reshape2::melt(to_plot_H3K9me3_order_long)
to_plot_H3K9me3_order_long$variable <- factor(to_plot_H3K9me3_order_long$variable,levels = H3K9me3_tissue_order)
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(colnames(to_plot)))
ggplot(to_plot_H3K9me3_order_long, aes(x = variable, y = value,fill=variable)) +
  geom_boxplot() +
  scale_fill_manual(values = color, name = "Tissue") +
  labs(title = "CpG methylation log2(Fold change)",
       x = NULL,
       y = "log2(old/young)") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 14),  # 调整图例标题大小
    legend.text = element_text(size = 12),   # 调整图例项的字符大小
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12)  # 调整x轴文本角度和大小
  )

to_plot_H3K9me3_order_long$condition <- "Large changes"
to_plot_H3K9me3_order_long$condition[which(to_plot_H3K9me3_order_long$variable %in% c("Ileum","Bladder","Testis","Tongue","Cecum","Pancreas","Stomach","Bone Marrow","Colon","Spleen","Thymus","Jejunum","iWAT"))] <- "Small changes"

ggplot(to_plot_H3K9me3_order_long, aes(x = condition, y = value,fill=condition)) +
  geom_boxplot(outliers =F) +
  # scale_fill_manual(values = color, name = "Tissue") +
  labs(title = "CpG methylation log2(Fold change) in H3K9me3 peaks",
       x = NULL,
       y = "log2(old/young)") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 14),  # 调整图例标题大小
    legend.text = element_text(size = 12),   # 调整图例项的字符大小
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12)  # 调整x轴文本角度和大小
  )+
  geom_signif(comparisons = list(c("Large changes", "Small changes")),
              textsize = 4,test = "t.test",
              map_signif_level = TRUE,
              y_position = c(0.07),
              tip_length = c(1/100))



