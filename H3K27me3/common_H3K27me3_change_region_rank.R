rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)

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
common_increase <- read.csv("data/samples/all/H3K27me3/common_increase_10kb_bins.csv")
common_decrease <- read.csv("data/samples/all/H3K27me3/common_decrease_10kb_bins.csv")
common_increase_region <- common_increase$Geneid[which(common_increase$n > 1)]
common_decrease_region <- common_decrease$Geneid[which(common_decrease$n > 1)]

tissues <- c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
             "thymus","skin","bladder","bonemarrow","Hip","heart",
             "muscle","jejunum","uterus","ovary","liver","tongue",
             "cecum","colon","testis","stomach","pancreas","iWAT","ileum")
increase_logFC_rank <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_10kb_bins_diff.csv"))
  df <- df[which(df$Geneid %in% common_increase_region),]
  mean_logFC <- mean(df$LogFC.old.young)
  t_increase_logFC_rank <- data.frame(tissue=tissue_label_change(tissue),mean_logFC=mean_logFC)  
  increase_logFC_rank <- rbind(increase_logFC_rank,t_increase_logFC_rank)
}

color <- read.table("data/samples/30_distinct_color.txt")
color <- color$V1
color <- setNames(color,sort(unique(increase_logFC_rank$tissue)))
increase_logFC_rank <- increase_logFC_rank[order(increase_logFC_rank$mean_logFC),]
increase_logFC_rank$tissue <- factor(increase_logFC_rank$tissue,levels=increase_logFC_rank$tissue)
ggplot(increase_logFC_rank,mapping = aes(x=mean_logFC,y=tissue,fill = tissue))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+ylab("")+xlab("Mean log2(Fold Change)")+ggtitle(paste0("Increase"))+
  theme(text = element_text(size = 18))+ scale_fill_manual(values = color) +guides(fill= guide_legend(title = ""))

increase_var <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_10kb_bins_diff.csv"))
  df <- df[which(df$Geneid %in% common_increase_region),c("Geneid","LogFC.old.young")]
  colnames(df)[2] <- tissue_label_change(tissue)
  if(ncol(increase_var)==0){
    increase_var <- df    
  }else{
    increase_var <- merge(increase_var,df,by="Geneid")
  }
}
rownames(increase_var) <- increase_var$Geneid
increase_var <- increase_var[,-1]
increase_sd <- as.data.frame(apply(increase_var, 1, sd))
colnames(increase_sd) <- "sd"
increase_mean <- apply(increase_var, 1, mean)  
increase_cv <- as.data.frame(increase_sd$`apply(increase_var, 1, sd)` / increase_mean)
colnames(increase_cv) <- "cv"
hist(increase_sd$sd, breaks = 30, main = "Distribution of Row Standard Deviations", xlab = "Standard Deviation")  
hist(increase_cv$cv, breaks = 30, main = "Distribution of log2(Fold Change) Coefficient of Variation", xlab = "Coefficient of Variation")  

decrease_logFC_rank <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_10kb_bins_diff.csv"))
  df <- df[which(df$Geneid %in% common_decrease_region),]
  mean_logFC <- mean(df$LogFC.old.young)
  t_decrease_logFC_rank <- data.frame(tissue=tissue_label_change(tissue),mean_logFC=mean_logFC)  
  decrease_logFC_rank <- rbind(decrease_logFC_rank,t_decrease_logFC_rank)
}

color <- read.table("data/samples/30_distinct_color.txt")
color <- color$V1
color <- setNames(color,sort(unique(decrease_logFC_rank$tissue)))
decrease_logFC_rank <- decrease_logFC_rank[order(decrease_logFC_rank$mean_logFC,decreasing = T),]
decrease_logFC_rank$tissue <- factor(decrease_logFC_rank$tissue,levels=decrease_logFC_rank$tissue)
ggplot(decrease_logFC_rank,mapping = aes(x=mean_logFC,y=tissue,fill = tissue))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+ylab("")+xlab("Mean log2(Fold Change)")+ggtitle(paste0("Decrease"))+
  theme(text = element_text(size = 18))+ scale_fill_manual(values = color) +guides(fill= guide_legend(title = ""))

decrease_var <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_10kb_bins_diff.csv"))
  df <- df[which(df$Geneid %in% common_decrease_region),c("Geneid","LogFC.old.young")]
  colnames(df)[2] <- tissue_label_change(tissue)
  if(ncol(decrease_var)==0){
    decrease_var <- df    
  }else{
    decrease_var <- merge(decrease_var,df,by="Geneid")
  }
}
rownames(decrease_var) <- decrease_var$Geneid
decrease_var <- decrease_var[,-1]
decrease_sd <- as.data.frame(apply(decrease_var, 1, sd))
colnames(decrease_sd) <- "sd"
decrease_mean <- apply(decrease_var, 1, mean)  
decrease_mean <- -decrease_mean
decrease_cv <- as.data.frame(decrease_sd$sd/ decrease_mean)
colnames(decrease_cv) <- "cv"
hist(decrease_sd$sd, breaks = 30, main = "Distribution of Row Standard Deviations", xlab = "Standard Deviation")  
hist(decrease_cv$cv, breaks = 30, main = "Distribution of log2(Fold Change) Coefficient of Variation", xlab = "Coefficient of Variation")  
