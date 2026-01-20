rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
library(ggrepel)
library(tidyr)
library(dplyr)
library(tidyverse)
library(plotly)
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
load("data/samples/GRN/grn_union_tissue.rdata")
load("data/samples/GRN/grn_union_skin.rdata")
grn_tissue[["skin"]] <- grn_union
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")

### limma
summary <- data.frame()
for(tissue in tissues){
  df <- grn_tissue[[tissue]]
  pagerank <- read.table(paste0("data/samples/GRN/TF_pagerank_limma_remove_zero_row_log/",tissue,"_TF_diff.txt"))
  df <- merge(df,pagerank[,c("logFC","P.Value","adj.P.Val")],by.x="TF",by.y="row.names")
  
  gene <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
  gene <- gene[which(gene$Significant !="Stable"),]
  df <- df[which(df$TF %in% gene$X),]

  colnames(df)[which(colnames(df)%in% c("logFC","P.Value","adj.P.Val"))] <- c("pagerank_logFC","P.Value","adj.P.Val")
  
  if(nrow(df) > 0){
    df$tf_logFC <- log2(df$tf_old/df$tf_young)
    df$gene_logFC <- log2(df$gene_old/df$gene_young)
    df$peak_logFC <- log2(df$peak_old/df$peak_young)
    df$condition <- "other"
    df$condition[which(df$gene_logFC >0 & df$tf_logFC >0 & df$peak_logFC >0)] <- "Up"
    df$condition[which(df$gene_logFC <0 & df$tf_logFC <0 & df$peak_logFC <0)] <- "Down"
    
    TF_counts <- as.data.frame(table(df$TF))
    df <- merge(df,TF_counts,by.x="TF",by.y="Var1")
    df$TF <- paste0(tissue_label_change(tissue),"-",df$TF)

    
    t_summary <- df[,c("TF","gene_logFC","Freq","pagerank_logFC","P.Value","adj.P.Val")]
    summary <- rbind(summary,t_summary)
  }
}

mean_logFC_summary <- summary %>%
  filter(!is.infinite(gene_logFC)) %>%  # 过滤掉无穷大值（如 Inf 和 -Inf）
  group_by(TF,pagerank_logFC,P.Value,adj.P.Val,Freq) %>%
  filter(n() >= 50) %>%  # 仅保留至少有10行数据的TF
  summarise(mean_logFC = mean(gene_logFC, na.rm = TRUE))  # 计算中位数
mean_logFC_summary <- mean_logFC_summary[order(mean_logFC_summary$mean_logFC,decreasing = T),]


to_plot <- mean_logFC_summary
to_plot$TF <- factor(to_plot$TF,levels=mean_logFC_summary$TF)
highest_points <- to_plot[1:10,]
lowest_points <- to_plot[(nrow(to_plot)-9):nrow(to_plot),]
annotate_points <- bind_rows(highest_points, lowest_points)

ggplot(to_plot, aes(x = TF, y = mean_logFC)) +
  geom_point(color = "black") +
  labs(title = "Boxplot of logFC by TF", x = "TF", y = "logFC") +
  theme_bw()+
  theme(
    axis.line.x = element_blank(),  
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank()  # 去掉次要网格线
  )+
  geom_text_repel(data = annotate_points, aes(label = TF), color = "red", size = 4)

to_plot_sig <- to_plot[which(to_plot$P.Value < 0.05),]
to_plot_sig <- to_plot_sig[order(to_plot_sig$pagerank_logFC,decreasing = T),]
highest_points <- to_plot_sig[1:10,]
lowest_points <- to_plot_sig[(nrow(to_plot_sig)-9):nrow(to_plot_sig),]
annotate_points <- bind_rows(highest_points, lowest_points)

ggplot(to_plot[which(to_plot$P.Value <0.05),], aes(x = pagerank_logFC, y = mean_logFC)) +
  geom_point(aes(size = Freq),color = "black",alpha=0.5) +
  labs(x = "pagerank logFC", y = "TF related gene logFC") +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank()  # 去掉次要网格线
  )+
  geom_text_repel(data = annotate_points, aes(label = TF), color = "red", size = 4)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "blue")


to_plot_sig <- to_plot[which(to_plot$adj.P.Val < 0.05),]
to_plot_sig <- to_plot_sig[order(to_plot_sig$pagerank_logFC,decreasing = T),]
highest_points <- to_plot_sig[1:10,]
lowest_points <- to_plot_sig[(nrow(to_plot_sig)-9):nrow(to_plot_sig),]
annotate_points <- bind_rows(highest_points, lowest_points)

ggplot(to_plot[which(to_plot$adj.P.Val <0.05),], aes(x = pagerank_logFC, y = mean_logFC)) +
  geom_point(aes(size = Freq),color = "black",alpha=0.5) +
  labs(x = "pagerank logFC", y = "TF related gene logFC") +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank()  # 去掉次要网格线
  )+
  geom_text_repel(data = annotate_points, aes(label = TF), color = "red", size = 4)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "blue")



