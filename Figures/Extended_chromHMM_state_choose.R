rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(corrplot)
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver","ileum",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
state <- "15"

chromHMM_compare <- read.table("result/all/ChromHMM/all_tissues_normal_chr/comparedir/compare_to_25_state_models.txt",header = T,row.names = 1)
rownames(chromHMM_compare) <- paste0("State",rownames(chromHMM_compare))
colnames(chromHMM_compare) <- 2:25

color_palette <- colorRampPalette(c("white", "#ffffd2","#e64b35"))(100) 
breaks <- c(seq(0, 0.29, length.out = 40), seq(0.3, 0.6, length.out = 20), seq(0.61, 1, length.out = 40))

pheatmap::pheatmap(chromHMM_compare,cluster_rows = F,cluster_cols = F,color = color_palette,breaks = breaks,filename = "result/Sup_figures/chromHMM_state_choose_heatmap.pdf",width = 8,height = 6)

medians <- apply(chromHMM_compare, 2, median)
to_plot <- data.frame(state=2:25,medians=medians)

p <- ggplot(to_plot, aes(x = state, y = medians)) +
  geom_line() +
  geom_point() + # 如果想在每个数据点上加个点
  labs(x = "Model", y = "Median correlation") +
  theme_bw()+
  ylim(0,1)+
  scale_x_continuous(breaks = 2:25) 
ggsave("result/Sup_figures/chromHMM_state_choose_line_plot.pdf",p,width = 8,height = 6)



