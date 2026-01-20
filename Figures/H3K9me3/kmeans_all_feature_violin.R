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
library(DSS)
library(patchwork)

tissues <- c("liver","lung","kidney","ileum","Hip","mammarygland","skin","bonemarrow",
             "jejunum","colon","ovary","CB","BAT","thymus","testis","stomach","heart",
             "muscle","bladder","aorta","tongue","spleen","pancreas","brain",
             "cecum","uterus","iWAT")
color <- setNames(c("#e64b35","#4dbbd5","#00a087","#3c5488","#f39b7f","grey","#636363"),c("kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"))

#H3K9me3 signal
to_plot <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/feature_table_combined_chrY/kmeans_H3K9me3_signal_detail.csv")
to_plot <- to_plot %>%
  group_by(label, cluster) %>%
  summarize(RPKM = mean(mean_RPKM, na.rm = TRUE))
to_plot <- to_plot[which(to_plot$cluster %in% c("kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak")),]
to_plot$cluster <- factor(to_plot$cluster,levels=c("kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"))
to_plot <- as.data.frame(to_plot)
# to_plot <- to_plot %>%
#   group_by(cluster) %>%
#   mutate(
#     Q1 = quantile(RPKM, 0.25),
#     Q3 = quantile(RPKM, 0.75),
#     IQR = Q3 - Q1,
#     lower_bound = Q1 - 1.5 * IQR,
#     upper_bound = Q3 + 1.5 * IQR
#   ) %>%
#   filter(RPKM >= lower_bound & RPKM <= upper_bound) %>%
#   select(-Q1, -Q3, -IQR, -lower_bound, -upper_bound)  


p1 <-ggplot(to_plot, aes(x = cluster, y = RPKM, fill=cluster)) +
  geom_violin(color = "black")+
  # geom_boxplot(color = "black",outliers = F)+
  theme_bw() +
  scale_fill_manual(values =color) +
  theme(
    axis.text.x = element_blank(),   
    axis.title.x = element_blank(),  
    legend.position = "none",        
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    panel.background = element_blank(),
    panel.grid.major =  element_blank(),
    panel.grid.minor =  element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +ylim(0,1)

ggsave("result/figures/H3K9me3_kmeans_H3K9me3_signal.pdf",p1,height = 3,width = 6)

#H3K27me3 signal
to_plot <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/feature_table_combined_chrY/kmeans_H3K27me3_signal_detail.csv")
to_plot <- to_plot %>%
  group_by(label, cluster) %>%
  summarize(RPKM = mean(mean_RPKM, na.rm = TRUE))
to_plot <- to_plot[which(to_plot$cluster %in% c("kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak")),]
to_plot$cluster <- factor(to_plot$cluster,levels=c("kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"))
to_plot <- as.data.frame(to_plot)
# to_plot <- to_plot %>%
#   group_by(cluster) %>%
#   mutate(
#     Q1 = quantile(RPKM, 0.25),
#     Q3 = quantile(RPKM, 0.75),
#     IQR = Q3 - Q1,
#     lower_bound = Q1 - 1.5 * IQR,
#     upper_bound = Q3 + 1.5 * IQR
#   ) %>%
#   filter(RPKM > lower_bound & RPKM < upper_bound) %>%
#   select(-Q1, -Q3, -IQR, -lower_bound, -upper_bound) 
p2 <-ggplot(to_plot, aes(x = cluster, y = RPKM, fill=cluster)) +
  geom_violin(color = "black")+
  # geom_boxplot(color = "black",outliers = F)+
  theme_bw() +
  scale_fill_manual(values =color) +
  theme(
    axis.text.x = element_blank(),   
    axis.title.x = element_blank(),  
    legend.position = "none",        
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +ylim(0,2)
ggsave("result/figures/H3K9me3_kmeans_H3K27me3_signal.pdf",p,height = 3,width = 6)

#DNA methylaion
to_plot <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/feature_table_combined_chrY/kmeans_median_DNA_methylation_detail.csv")
to_plot <- to_plot[,-1]
to_plot_mean <- to_plot[,c(1:28)]
to_plot_mean$mean <- rowMeans(to_plot_mean[,-1],na.rm = T)
to_plot_mean <- to_plot_mean[,c(1,29)]
to_plot <- merge(to_plot_mean,to_plot[,c(1,29)],by="label")
to_plot <- to_plot[which(to_plot$cluster %in% c("kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak")),]

to_plot$cluster <- factor(to_plot$cluster,levels=c("kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"))
p3 <-ggplot(to_plot, aes(x = cluster, y = mean, fill=cluster)) +
  geom_violin(color = "black")+
  # geom_boxplot(color = "black",outliers = F)+
  theme_bw() +
  scale_fill_manual(values =color) +
  theme(
    axis.text.x = element_blank(),   
    axis.title.x = element_blank(),  
    legend.position = "none",        
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +ylim(50,100)
ggsave("result/figures/H3K9me3_kmeans_DNA_methylation.pdf",p,height = 3,width = 6)

#Compartment B
to_plot <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/feature_table_combined_chrY/kmeans_mean_compartmentB_percentage_detail.csv")
to_plot <- to_plot[,-1]
to_plot_mean <- to_plot[,c(1:(ncol(to_plot)-1))]
to_plot_mean$mean <- rowMeans(to_plot_mean[,-1],na.rm = T)
to_plot_mean <- to_plot_mean[,c(1,ncol(to_plot_mean))]
to_plot <- merge(to_plot_mean,to_plot[,c(1,ncol(to_plot))],by="label")
to_plot <- to_plot[which(to_plot$cluster %in% c("kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak")),]
to_plot$cluster <- factor(to_plot$cluster,levels=c("kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"))
p4 <-ggplot(to_plot, aes(x = cluster, y = mean, fill=cluster)) +
  geom_violin(color = "black")+
  # geom_boxplot(color = "black",outliers = F)+
  theme_bw() +
  scale_fill_manual(values =color) +
  theme(
    axis.text.x = element_blank(),   
    axis.title.x = element_blank(),  
    legend.position = "none",        
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +ylim(0,100)
ggsave("result/figures/H3K9me3_kmeans_compartmentB.pdf",p,height = 3,width = 6)

#TE
#LTR
to_plot <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/feature_table_combined_chrY/kmeans_median_LTR_median_TE_length_percentage_detail.csv")
to_plot <- to_plot[which(to_plot$condition %in% c("kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak")),]
to_plot$condition <- factor(to_plot$condition,levels=c("kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"))
p5 <-ggplot(to_plot, aes(x = condition, y = TE_length_percentage, fill=condition)) +
  geom_violin(color = "black")+
  # geom_boxplot(color = "black",outliers = F)+
  theme_bw() +
  scale_fill_manual(values =color) +
  theme(
    axis.text.x = element_blank(),   
    axis.title.x = element_blank(),  
    legend.position = "none",        
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +ylim(0,40)
ggsave("result/figures/H3K9me3_kmeans_LTR.pdf",p,height = 3,width = 6)

to_plot <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/feature_table_combined_chrY/kmeans_median_ERV1_median_TE_length_percentage_detail.csv")
to_plot <- to_plot[which(to_plot$condition %in% c("kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak")),]
to_plot$condition <- factor(to_plot$condition,levels=c("kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"))
p6 <-ggplot(to_plot, aes(x = condition, y = TE_length_percentage, fill=condition)) +
  geom_violin(color = "black")+
  # geom_boxplot(color = "black",outliers = F)+
  theme_bw() +
  scale_fill_manual(values =color) +
  theme(
    axis.text.x = element_blank(),   
    axis.title.x = element_blank(),  
    legend.position = "none",        
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +ylim(0,20)
ggsave("result/figures/H3K9me3_kmeans_ERV1.pdf",p,height = 3,width = 6)

to_plot <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/feature_table_combined_chrY/kmeans_median_ERVK_median_TE_length_percentage_detail.csv")
to_plot <- to_plot[which(to_plot$condition %in% c("kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak")),]
to_plot$condition <- factor(to_plot$condition,levels=c("kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"))
p7 <-ggplot(to_plot, aes(x = condition, y = TE_length_percentage, fill=condition)) +
  geom_violin(color = "black")+
  # geom_boxplot(color = "black",outliers = F)+
  theme_bw() +
  scale_fill_manual(values =color) +
  theme(
    axis.text.x = element_blank(),   
    axis.title.x = element_blank(),  
    legend.position = "none",        
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +ylim(0,30)
ggsave("result/figures/H3K9me3_kmeans_ERVK.pdf",p,height = 3,width = 6)

#phast cons
to_plot <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/feature_table_combined_chrY/kmeans_median_phastCons_score_detail.csv")
to_plot <- to_plot[which(to_plot$cluster %in% c("kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak")),]
to_plot$cluster <- factor(to_plot$cluster,levels=c("kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"))
p8 <-ggplot(to_plot, aes(x = cluster, y = avg_V4, fill=cluster)) +
  geom_violin(color = "black")+
  # geom_boxplot(color = "black",outliers = F)+
  theme_bw() +
  scale_fill_manual(values =color) +
  theme(
    axis.text.x = element_blank(),   
    axis.title.x = element_blank(),  
    legend.position = "none",        
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +ylim(0,0.6)
ggsave("result/figures/H3K9me3_kmeans_phastcons.pdf",p,height = 3,width = 6)

combined_plot <- p1 / p2 / p3 / p4/p5/p6/p7/p8
ggsave("result/figures/H3K9me3_all_features_boxplot.pdf",combined_plot,width = 6,height = 8)
