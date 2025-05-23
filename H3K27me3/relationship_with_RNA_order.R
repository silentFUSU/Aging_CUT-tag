rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(stringr)
library(dplyr)
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

tissues <- sort(c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
                  "thymus","skin","bladder","bonemarrow","Hip","heart",
                  "muscle","jejunum","uterus","ovary","liver","tongue",
                  "cecum","colon","testis","stomach","pancreas","iWAT","ileum"))
summary <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_gene_TSS_diff_after_remove_batch_effect.csv"))
  df <- df[which(df$Significant != "Stable"),]
  RNA <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
  RNA <- RNA[which(RNA$Significant != "Stable"),]
  colnames(RNA)[1] <- "Geneid"
  df <- merge(df[,c("Geneid","LogFC.old.young")],RNA[,c("Geneid","logFC")],by="Geneid")
  second_quadrant <- df[which(df$LogFC.old.young > 0 & df$logFC < 0),]
  fourth_quadrant <- df[which(df$LogFC.old.young < 0 & df$logFC > 0),]
  t_summary <- data.frame(condition=c("second","fourth"),count=c(nrow(second_quadrant),nrow(fourth_quadrant)))
  colnames(t_summary)[2] <- tissue_label_change(tissue)
  if(nrow(summary)==0){
    summary <-t_summary
  }else{
    summary <- merge(summary,t_summary,by="condition")
  }
}
summary <- reshape2::melt(summary)
summary$value[which(summary$condition=="second")] <- -summary$value[which(summary$condition=="second")]
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(unique(summary$variable)))
summary_aggregate <- summary %>%  
  group_by(variable) %>%  
  summarise(total_abs_value = sum(abs(value)))  
summary_aggregate <- summary_aggregate[order(summary_aggregate$total_abs_value),]
summary$variable <- factor(summary$variable,summary_aggregate$variable)
ggplot(summary, aes(x = variable, y = value, fill = variable,color = condition)) +  
  geom_bar(stat = "identity",aes(alpha = ifelse(condition == "second", 0.8, 1))) + 
  theme_minimal() +  
  xlab(NULL)+
  ylab("Count")+
  scale_y_continuous(labels = abs) +  
  scale_fill_manual(values=color) +
  theme_bw()+
  theme(  
    axis.title.x = element_text(size = 14),     
    axis.title.y = element_text(size = 14),    
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),   
    axis.text.y = element_text(size = 12),    
    plot.title = element_text(size = 16, face = "bold"),
  ) +  
  scale_color_manual(values = c("second" = "black", "fourth" = "white")) 

gene_expression_summary <- data.frame()
for(tissue in tissues){
  RNA <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
  colnames(RNA)[1] <- "Geneid"
  RNA <- RNA[which(RNA$Significant != "Stable"),]
  increase <- RNA[which(RNA$logFC > 0),]
  decrease <- RNA[which(RNA$logFC < 0),]
  t_summary <- data.frame(condition=c("Increase","Decrease"),count=c(nrow(increase),nrow(decrease)))
  colnames(t_summary)[2] <- tissue_label_change(tissue)
  if(nrow(gene_expression_summary)==0){
    gene_expression_summary <-t_summary
  }else{
    gene_expression_summary <- merge(gene_expression_summary,t_summary,by="condition")
  }
}
gene_expression_summary <- reshape2::melt(gene_expression_summary)
gene_expression_summary$value[which(gene_expression_summary$condition=="Decrease")] <- -gene_expression_summary$value[which(gene_expression_summary$condition=="Decrease")]
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(unique(gene_expression_summary$variable)))
gene_expression_summary_aggregate <- gene_expression_summary %>%  
  group_by(variable) %>%  
  summarise(total_abs_value = sum(abs(value)))  
gene_expression_summary_aggregate <- gene_expression_summary_aggregate[order(gene_expression_summary_aggregate$total_abs_value),]
gene_expression_summary$variable <- factor(gene_expression_summary$variable,summary_aggregate$variable)

summary$variable <- factor(summary$variable,summary_aggregate$variable)
ggplot(gene_expression_summary, aes(x = variable, y = value, fill = variable,color = condition)) +  
  geom_bar(stat = "identity",aes(alpha = ifelse(condition == "Decrease", 0.8, 1))) + 
  geom_bar(data = summary,aes(x = variable, y = value), stat = "identity", fill = "black", alpha = 0.5) + 
  theme_minimal() +  
  xlab(NULL)+
  ylab("Count")+
  scale_y_continuous(labels = abs) +  
  scale_fill_manual(values=color) +
  theme_bw()+
  theme(  
    axis.title.x = element_text(size = 14),     
    axis.title.y = element_text(size = 14),    
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),   
    axis.text.y = element_text(size = 12),    
    plot.title = element_text(size = 16, face = "bold"),
  ) +  
  scale_color_manual(values = c("Increase" = "white","Decrease" = "black")) 


second_summary <- summary[which(summary$condition=="second"),]
decrease_gene_expression_summary <- gene_expression_summary[which(gene_expression_summary$condition=="Decrease"),]
second_summary <- merge(second_summary[,c(1:3)],decrease_gene_expression_summary[,c(2,3)], by="variable")
second_summary$percentage <- abs(second_summary$value.x)/abs(second_summary$value.y) *100
second_summary <- second_summary[,c("variable","condition","percentage")]

fourth_summary <- summary[which(summary$condition=="fourth"),]
increase_gene_expression_summary <- gene_expression_summary[which(gene_expression_summary$condition=="Increase"),]
fourth_summary <- merge(fourth_summary[,c(1:3)],increase_gene_expression_summary[,c(2,3)], by="variable")
fourth_summary$percentage <- abs(fourth_summary$value.x)/abs(fourth_summary$value.y) *100
fourth_summary <- fourth_summary[,c("variable","condition","percentage")]

percentage_summary <- rbind(second_summary,fourth_summary)
percentage_summary$variable <- factor(percentage_summary$variable,summary_aggregate$variable)
percentage_summary$percentage[which(percentage_summary$condition=="second")] <- -percentage_summary$percentage[which(percentage_summary$condition=="second")]
percentage_summary <- merge(percentage_summary,chi_label_summary,by.x="variable",by.y="tissue")

ggplot(percentage_summary, aes(x = variable, y = percentage, fill = variable,color=condition)) +  
  geom_bar(stat = "identity",aes(alpha = ifelse(condition == "second", 0.8, 1))) + 
  theme_minimal() +  
  xlab(NULL)+
  ylab("Percentgae")+
  scale_y_continuous(labels = abs, limits = c(-50, 50)) + 
  scale_fill_manual(values=color) +
  theme_bw()+
  geom_text(  
    data = percentage_summary[which(percentage_summary$condition=="fourth"),],
    aes(x = variable, y = percentage,label = label),   
    vjust = 0,
    color = "black" 
  ) +
  theme(  
    axis.title.x = element_text(size = 14),     
    axis.title.y = element_text(size = 14),    
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),   
    axis.text.y = element_text(size = 12),    
    plot.title = element_text(size = 16, face = "bold"),
  ) +  
  scale_color_manual(values = c("second" = "black", "fourth" = "white")) 
