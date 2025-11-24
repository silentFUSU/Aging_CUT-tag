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
tissue_summary <- read.csv("data/samples/WGBS/all_tissues_delta_in_H3K9me3_recursion_peaks.csv",row.names = 1)
colnames(tissue_summary)[which(colnames(tissue_summary)=="Mammary.Gland")] <- "Mammary Gland"
colnames(tissue_summary)[which(colnames(tissue_summary)=="Bone.Marrow")] <- "Bone Marrow"
colnames(tissue_summary)[which(colnames(tissue_summary)=="IWAT")] <- "iWAT"

annotation <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation_add_stable.csv",row.names = 1)
colnames(annotation)[1]<-"label"
random_regions1 <- read.table("data/samples/all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY.bed")
random_regions1$cluster <- "random whole genome"
random_regions1$label <- paste0(random_regions1$V1,":",random_regions1$V2,"-",random_regions1$V3)
annotation <- rbind(annotation,random_regions1[,c("label","cluster")])

random_regions2 <- read.table("data/samples/all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY_out_recursion_peaks.bed")
random_regions2$cluster <- "random out of peak"
random_regions2$label <- paste0(random_regions2$V1,":",random_regions2$V2,"-",random_regions2$V3)
annotation <- rbind(annotation,random_regions2[,c("label","cluster")])

to_plot <- merge(tissue_summary,annotation,by="label")
to_plot <- reshape2::melt(to_plot)
# to_plot$value <- abs(to_plot$value)
# to_plot <- to_plot %>%
#   mutate(cluster = ifelse(grepl("chrY", label) & cluster == "kmeans1", "kmeans1-chrY", cluster))

to_plot$cluster <- factor(to_plot$cluster,levels=c("kmeans1-chrY","kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"))
to_plot <- to_plot[!is.na(to_plot$value), ]
ggplot(to_plot, aes(x = cluster, y =value, fill = cluster)) +  
  geom_boxplot(alpha = 0.7,outliers = F) +  
  # scale_fill_manual(values = color) +
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
  ) +
  geom_signif(
    comparisons = list(c("Stable", "random out of peak"),c("Stable", "random whole genome"),
                       c("kmeans2", "random out of peak"),c("kmeans2", "random whole genome"),
                       c("kmeans1", "random out of peak"),c("kmeans1", "random whole genome"),
                       c("kmeans1-chrY", "random out of peak"),c("kmeans1-chrY", "random whole genome"),
                       c("kmeans3", "random out of peak"),c("kmeans3", "random whole genome"),
                       c("kmeans4", "random out of peak"),c("kmeans4", "random whole genome")
                       ),
    textsize = 4,test = "t.test",
    map_signif_level = TRUE,
    y_position = c(0.05, 0.06, 0.07, 0.08, 0.09, 0.10,0.11,0.12,0.13,0.14,0.15,0.16),
    tip_length = c(1/100,1/100,1/100,1/100,1/100,1/100,1/100,1/100,1/100,1/100,1/100,1/100)
  )

median_values <- to_plot %>%
  group_by(cluster, variable) %>%
  summarize(median_value = median(value, na.rm = TRUE), .groups = "drop")

to_plot_median <- reshape2::dcast(median_values,cluster~variable)
rownames(to_plot_median) <- to_plot_median$cluster
to_plot_median <- to_plot_median[,-1]
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- c(seq(-0.1, -0.02, length.out = 40), seq(-0.019, 0.019, length.out = 20), seq(0.02, 0.1, length.out = 40))
# tissue_order <-c("Mammary Gland","Cecum","Thymus","Uterus","iWAT","Stomach","Skin","Spleen","Muscle","Bone Marrow","Liver","Ileum","Testis","Cortex","Jejunum","Tongue","Hippocampus","Colon","Bladder",
#                  "Aorta","Cerebellum","Lung","Heart","Kidney","BAT","Ovary","Pancreas")
tissue_order <-c("Lung","Cerebellum","BAT","Muscle","Heart","Aorta","Skin","Kidney","Hippocampus","Cortex","Liver","Tongue","Uterus","Testis","Bladder","Ovary",
                 "Colon","Stomach","Thymus","Cecum","Jejunum","Pancreas","Bone Marrow","Ileum","Spleen","iWAT","Mammary Gland")
to_plot_median <- to_plot_median[,tissue_order]
pheatmap::pheatmap(to_plot_median,cluster_rows = F,cluster_cols = F,color = color_palette,breaks = breaks,border_color = "black",filename = "result/Sup_figures/WGBS_delta_in_H3K9me3_peaks.pdf",width = 8,height = 6)

to_plot_median_long <-to_plot_median
to_plot_median_long$cluster <- rownames(to_plot_median_long)
to_plot_median_long <- reshape2::melt(to_plot_median_long)
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(as.character(unique(to_plot_median_long$variable))))
to_plot_median_long$cluster <- factor(to_plot_median_long$cluster,levels=c("kmeans1-chrY","kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"))
ggplot(to_plot_median_long[which(to_plot_median_long$variable %in% tissue_order[c(1:10)]),], aes(x = cluster, y =value, fill = cluster)) +  
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
to_plot_median_test <- reshape2::dcast(to_plot_median_long,variable~cluster)
t.test(to_plot_median_test$Stable,to_plot_median_test$`random out of peak`,paired = T,alternative = "less")
to_plot_median_long <- to_plot_median_long[!is.na(to_plot_median_long$value), ]
p_value_data.frame <- data.frame(`random out of peak` = numeric(6), `random whole genome` = numeric(6), stringsAsFactors = FALSE)
rownames(p_value_data.frame) <- c("kmeans1-chrY","kmeans1","kmeans2","kmeans3","kmeans4","Stable")
colnames(p_value_data.frame) <- c("random out of peak","random whole genome")
for(i in c(1:6)){
  for(j in c(1:2)){
    test <- wilcox.test(to_plot_median_long$value[which(to_plot_median_long$cluster==rownames(p_value_data.frame)[i])],to_plot_median_long$value[which(to_plot_median_long$cluster==colnames(p_value_data.frame)[j])])    
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
