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
summary <- data.frame()
for(tissue in tissues){
  df <- grn_tissue[[tissue]]
  pagerank <- read.table(paste0("data/samples/GRN/TF_pagerank_limma_remove_zero_row_log//",tissue,"_TF_diff.txt"))
  df <- merge(df,pagerank[,c("logFC"),drop=F],by.x="TF",by.y="row.names")
  colnames(df)[which(colnames(df)=="logFC")] <- "tf_logFC"
  df$gene_logFC <- log2(df$gene_old/df$gene_young)
  df$peak_logFC <- log2(df$peak_old/df$peak_young)
  df$condition <- "other"
  df$condition[which(df$gene_logFC >0 & df$tf_logFC >0 & df$peak_logFC >0)] <- "Up"
  df$condition[which(df$gene_logFC <0 & df$tf_logFC <0 & df$peak_logFC <0)] <- "Down"
  
  t_summary <- as.data.frame(table(df$condition[which(df$condition !="other")]))
  t_summary$tissue <- tissue_label_change(tissue)
  summary <- rbind(summary,t_summary)
}
summary$Var1 <- factor(summary$Var1,levels = c("Up","Down"))
p <- ggplot(summary, aes(x = tissue, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  labs(x = "Tissue", y = "Frequency", fill = "Var1") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)  # 如果有太多 tissue 列，可以旋转标签
  )
p
ggsave("result/Sup_figures/TF_cre_gene_same_trend.pdf",p,width = 8,height = 6)




