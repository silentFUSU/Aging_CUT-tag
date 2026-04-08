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
options(scipen = 999)
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
      tissue_label <- "iWAT"
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
  # H3K9me3_peaks <- read.table(paste0("data/samples/",tissue,"/H3K9me3/bed/H3K9me3_young_old_merge-W5000-G10000-E100.bed"))
  H3K9me3_peaks <- read.table(paste0("data/samples/all/H3K9me3/bed/H3K9me3_young_old_merge-W5000-G10000-E100_recursion.bed"))
  H3K9me3_peaks <- H3K9me3_peaks[which(H3K9me3_peaks$V3 - H3K9me3_peaks$V2 +1 >= 200000),]
  if(tissue %in% c("ovary","mammarygland","uterus")){
    H3K9me3_peaks <- H3K9me3_peaks[which(H3K9me3_peaks$V1 %in% paste0("chr",c(1:19,"X"))),]
  }
  setDT(H3K9me3_peaks)
  setkey(H3K9me3_peaks,V1,V2,V3)
  summary <- data.frame()
  for(sample in t_search_table$sample_name){
    df <- fread(paste0("data/samples/WGBS/",tissue,"/bdg/",sample,"_CpG.bdg"),sep = "\t")
    df[, label := paste(V1, V2, sep = "-")]
    setDT(df)
    setkey(df,V1,V2,V3)  
    overlaps <- foverlaps(df, H3K9me3_peaks, type = "any", nomatch = 0L)  
    df[, peak_condition := "out"]
    df[label %in% overlaps$label, peak_condition := "peak"]
    
    result <- df[, .(V4_sum = sum(V4), V5_sum = sum(V5)), by = peak_condition]
    result <- as.data.frame(result)
    result$methylation <- result$V4_sum/result$V5_sum
    result <- result[,c("peak_condition","methylation")]
    colnames(result)[2] <- sample
    if(nrow(summary)==0){
      summary <- result
    }else{
      summary <- merge(summary,result,by="peak_condition")
    }
  }
  young_summary <- summary[,c("peak_condition",t_search_table$sample_name[which(t_search_table$age=="3M")])]
  old_summary <- summary[,c("peak_condition",t_search_table$sample_name[which(t_search_table$age=="24M")])]
  young_summary$young_methylation <- rowMeans(young_summary[,-1])
  old_summary$old_methylation <- rowMeans(old_summary[,-1])
  t_tissue_summary <- merge(young_summary,old_summary,by="peak_condition")
  t_tissue_summary$delta <- t_tissue_summary$old_methylation - t_tissue_summary$young_methylation
  t_tissue_summary <- t_tissue_summary[,c("peak_condition","delta")]
  colnames(t_tissue_summary)[2] <- tissue_label_change(tissue)
  if(nrow(tissue_summary)==0){
    tissue_summary <- t_tissue_summary
  }else{
    tissue_summary <- merge(tissue_summary,t_tissue_summary,by="peak_condition",all=T)
  }
}
# write.csv(tissue_summary,"data/samples/WGBS/all/DNA_methylation_change_in_H3K9me3_recursion_peaks.csv")
tissue_summary <- read.csv("data/samples/WGBS/all/DNA_methylation_change_in_H3K9me3_recursion_peaks.csv",row.names = 1)
rownames(tissue_summary) <- tissue_summary$peak_condition
tissue_summary <- tissue_summary[,-1]


tissue_order <- c("Mammary.Gland","Cecum","Thymus","Uterus","iWAT","Stomach","Skin","Spleen","Muscle","Bone.Marrow","Liver",
                  "Ileum","Testis","Cortex","Jejunum","Tongue","Hippocampus","Colon","Bladder","Aorta","Cerebellum","Lung",
                  "Heart","Kidney","BAT","Ovary","Pancreas")
tissue_summary <- as.data.frame(t(tissue_summary))

to_plot <- tissue_summary[tissue_order,]
to_plot <- to_plot * 100
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)

breaks <- c(seq(-5, -0.61, length.out = 40), seq(-0.6, 0.6, length.out = 20), seq(0.61, 5, length.out = 40))
pheatmap::pheatmap(to_plot,cluster_rows = F,cluster_cols = F,breaks = breaks, color = color_palette)

to_plot <- tissue_summary
to_plot$tissue <- rownames(tissue_summary)
to_plot <- reshape2::melt(to_plot)
to_plot$value <- to_plot$value * 100
to_plot <- to_plot[which(to_plot$tissue %in% c("Mammary.Gland","Cecum","Uterus","iWAT","Stomach","Skin","Spleen")),]

p_value_test <- tissue_summary[which(rownames(tissue_summary) %in% c("Mammary.Gland","Cecum","Uterus","iWAT","Stomach","Skin","Spleen")),]
p_value_test <- p_value_test * 100
wilcox.test(p_value_test$peak, p_value_test$out, paired = TRUE, alternative = "less")
t.test(p_value_test$peak, p_value_test$out, paired = TRUE, alternative = "less")
differences <- p_value_test$out - p_value_test$peak
t.test(differences, mu = 0, alternative = "greater")

color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(tissue_order))
p <- ggplot(to_plot, aes(x = variable, y =value, fill = variable)) +  
  geom_boxplot(alpha = 0.7,outliers = F) +  
  scale_fill_brewer(palette = "Pastel1") +
  geom_point(aes(color = tissue), size = 2) + 
  scale_color_manual(values = color) +
  labs(  
    title = paste0("DNA methylation delta in H3K9me3 peaks"),  
    x = "Cluster",  
    y="abs(Delta)"
  ) + 
  theme_minimal() +  
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 14,angle = 45,hjust = 1),  
    axis.text.y = element_text(size = 14),  
    plot.title = element_text(size = 16, face = "bold")
  )
ggsave("result/figures/DNA_methylation_change_in_H3K9me3_peaks_boxplot.pdf",p,width = 4,height = 5)
