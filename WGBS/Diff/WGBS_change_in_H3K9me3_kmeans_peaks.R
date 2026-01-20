rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)  
library(ggsignif)
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
  H3K9me3_peaks <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation_add_stable.csv",row.names = 1)
  split_chr <- strsplit(as.character(H3K9me3_peaks$label), ":")  
  chr_column <- sapply(split_chr, `[[`, 1)  
  split_start_end <- strsplit(sapply(split_chr, `[[`, 2), "-")  
  start_column <- sapply(split_start_end, `[[`, 1)  
  end_column <- sapply(split_start_end, `[[`, 2)  
  H3K9me3_peaks <- data.frame(chr = chr_column,start = start_column, end = end_column, cluster=H3K9me3_peaks$cluster)
  random_regions1 <- read.table("data/samples/all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY.bed")
  random_regions1$cluster <- "random whole genome"
  colnames(random_regions1)[1:3] <- c("chr","start","end")
  H3K9me3_peaks <- rbind(H3K9me3_peaks,random_regions1)
  
  random_regions2 <- read.table("data/samples/all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY_out_recursion_peaks.bed")
  random_regions2$cluster <- "random out of peak"
  colnames(random_regions2)[1:3] <- c("chr","start","end")
  H3K9me3_peaks <- rbind(H3K9me3_peaks,random_regions2)
  
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
    
    result <- overlaps[, .(V4_sum = sum(V4), V5_sum = sum(V5)), by = cluster]
    result <- as.data.frame(result)
    result$methylation <- result$V4_sum/result$V5_sum *100
    result <- result[,c("cluster","methylation")]
    colnames(result)[2] <- sample
    if(nrow(summary)==0){
      summary <- result
    }else{
      summary <- merge(summary,result,by="cluster")
    }
  }
  young_summary <- summary[,c("cluster",t_search_table$sample_name[which(t_search_table$age=="3M")])]
  old_summary <- summary[,c("cluster",t_search_table$sample_name[which(t_search_table$age=="24M")])]
  young_summary$young_methylation <- rowMeans(young_summary[,-1])
  old_summary$old_methylation <- rowMeans(old_summary[,-1])
  
  t_tissue_summary <- merge(young_summary,old_summary,by="cluster")
  t_tissue_summary$delta <- t_tissue_summary$old_methylation - t_tissue_summary$young_methylation
  t_tissue_summary <- t_tissue_summary[,c("cluster","delta")]
  colnames(t_tissue_summary)[2] <- tissue_label_change(tissue)
  if(nrow(tissue_summary)==0){
    tissue_summary <- t_tissue_summary
  }else{
    tissue_summary <- merge(tissue_summary,t_tissue_summary,by="cluster",all=T)
  }
}
# write.csv(tissue_summary,"data/samples/WGBS/all_tissues_delta_in_H3K9me3_recursion_peaks_mCG_allCG_ratio.csv",row.names = 1)
tissue_summary <- read.csv("data/samples/WGBS/all_tissues_delta_in_H3K9me3_recursion_peaks_mCG_allCG_ratio.csv")
to_plot <- tissue_summary
to_plot$cluster <- factor(to_plot$cluster,levels=c("kmeans1-chrY","kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"))
to_plot <- to_plot[order(to_plot$cluster),]
to_plot <- to_plot[,-1]
rownames(to_plot) <- to_plot$cluster
to_plot<- to_plot[,-1]
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- c(seq(-10, -2, length.out = 40), seq(-1.9, 1.9, length.out = 20), seq(2, 10, length.out = 40))
# tissue_order <-c("Mammary Gland","Cecum","Thymus","Uterus","iWAT","Stomach","Skin","Spleen","Muscle","Bone Marrow","Liver","Ileum","Testis","Cortex","Jejunum","Tongue","Hippocampus","Colon","Bladder",
#                  "Aorta","Cerebellum","Lung","Heart","Kidney","BAT","Ovary","Pancreas")
tissue_order <-c("Lung","Cerebellum","BAT","Muscle","Heart","Aorta","Skin","Kidney","Hippocampus","Cortex","Liver","Tongue","Uterus","Testis","Bladder","Ovary",
                 "Colon","Stomach","Thymus","Cecum","Jejunum","Pancreas","Bone.Marrow","Ileum","Spleen","iWAT","Mammary.Gland")
to_plot <- to_plot[,tissue_order]
pheatmap::pheatmap(to_plot,cluster_rows = F,cluster_cols = F,color = color_palette,breaks = breaks,border_color = "black",filename = "result/Sup_figures/WGBS_delta_in_H3K9me3_peaks.pdf",width = 8,height = 6)

to_plot_long <-to_plot
to_plot_long$cluster <- rownames(to_plot_long)
to_plot_long <- reshape2::melt(to_plot_long)
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(as.character(unique(to_plot_long$variable))))
to_plot_long$cluster <- factor(to_plot_long$cluster,levels=c("kmeans1-chrY","kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"))
ggplot(to_plot_long[which(to_plot_long$variable %in% tissue_order[c(1:10)]),], aes(x = cluster, y =value, fill = cluster)) +  
  geom_boxplot(alpha = 0.7,outliers = F) +  
  geom_point(aes(color = variable), size = 2) + 
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
  )+
  geom_signif(
    comparisons = list(c("Stable", "random out of peak"),
                       c("kmeans2", "random out of peak")),
    textsize = 4,test = "t.test",
    map_signif_level = TRUE,
    y_position = c(0.01, 0.02),
    tip_length = c(1/100,1/100)
  )
to_plot_test <- reshape2::dcast(to_plot_long,variable~cluster)
t.test(to_plot_test$Stable,to_plot_test$`random out of peak`,paired = T,alternative = "less")
to_plot_long <- to_plot_long[!is.na(to_plot_long$value), ]
p_value_data.frame <- data.frame(`random out of peak` = numeric(6), `random whole genome` = numeric(6), stringsAsFactors = FALSE)
rownames(p_value_data.frame) <- c("kmeans1-chrY","kmeans1","kmeans2","kmeans3","kmeans4","Stable")
colnames(p_value_data.frame) <- c("random out of peak","random whole genome")
for(i in c(1:6)){
  for(j in c(1:2)){
    test <- wilcox.test(to_plot_long$value[which(to_plot_long$cluster==rownames(p_value_data.frame)[i])],to_plot_long$value[which(to_plot_long$cluster==colnames(p_value_data.frame)[j])])    
    p_value_data.frame[i,j] <- test$p.value
  }
}


p_value_data_frame <- data.frame()
for(tissue in tissue_order){
  print(tissue)
  t_p_value_data_frame <- data.frame()
  for(condition in c("kmeans1-chrY","kmeans1","kmeans2","kmeans3","kmeans4","Stable")){
    if(condition == "kmeans1-chrY" & tissue %in% c("Ovary","Mammary Gland","Uterus")){
      t_p_value_data_frame_condition <- data.frame(condition=condition,p_value=NA)
    }else{
      test <- wilcox.test(to_plot$value[which(to_plot$cluster==condition & to_plot$variable==tissue)],to_plot$value[which(to_plot$cluster=="random out of peak" & to_plot$variable==tissue)])
      t_p_value_data_frame_condition <- data.frame(condition=condition,p_value=test$p.value)
    }
    t_p_value_data_frame <- rbind(t_p_value_data_frame,t_p_value_data_frame_condition)
  }
  colnames(t_p_value_data_frame)[2] <- tissue
  if(nrow(p_value_data_frame)==0){
    p_value_data_frame <- t_p_value_data_frame
  }else{
    p_value_data_frame <- merge(p_value_data_frame,t_p_value_data_frame,by="condition")
  }
}
mark_significance <- function(p_value) {
  if (is.na(p_value)) {
    return(NA)
  } else if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return(NA)
  }
}
p_value_data_frame_significant <- p_value_data_frame %>%
  mutate(across(-condition, ~ sapply(., mark_significance)))
