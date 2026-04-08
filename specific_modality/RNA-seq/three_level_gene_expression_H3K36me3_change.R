rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
library(dplyr)
tissue <- "brain"
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
ref <- read.table("/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/TSS/refBed/mm10_refGene.bed")
ref$length <- ref$V3-ref$V2+1
ref <- ref %>%  
  group_by(V5) %>%  
  filter(length == max(length)) %>%  
  ungroup()  
three_levels_gene_expression <- function(tissue){
  dir.create(paste0("data/samples/RNA/",tissue,"/bed/"))
  RNA <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene.csv"))
  RNA <- RNA %>%  
    mutate(expression_level = cut(logCPM,   
                                  breaks = quantile(logCPM, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),   
                                  labels = c("low", "medium", "high"),  
                                  include.lowest = TRUE))  
  conditions <- c("low","medium","high")
  to_plot <- data.frame()
  for(condition in conditions){
    t_to_plot <- RNA[which(RNA$expression_level==condition),]  
    t_to_plot <- as.data.frame(table(t_to_plot$Significant))
    t_to_plot$percent <- t_to_plot$Freq/sum(t_to_plot$Freq)*100
    t_to_plot$condition <- condition
    to_plot <- rbind(to_plot,t_to_plot)
  }
  to_plot$condition <- factor(to_plot$condition,levels=c("low","medium","high"))
  to_plot$Var1 <- factor(to_plot$Var1,levels=c("Stable","Up","Down"))
  ggplot(to_plot, aes(x = condition, y = percent, fill = Var1)) +  
    geom_bar(stat = "identity") +  
    labs(x = "Expression Level", y = "Percentage", fill = "condition") +  
    theme_minimal() +  
    theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    ggtitle(paste0(tissue_label_change(tissue))) +
    labs(x = "", y = "Percentage", fill = "Category")
  
  low <- ref[which(ref$V5 %in% RNA$X[which(RNA$expression_level=="low")]),]
  medium <- ref[which(ref$V5 %in% RNA$X[which(RNA$expression_level=="medium")]),]
  high <- ref[which(ref$V5 %in% RNA$X[which(RNA$expression_level=="high")]),]
  write.table(low,paste0("data/samples/RNA/",tissue,"/bed/low_expression_gene.bed"),append = F,quote = F,row.names = F,col.names = F,sep = "\t")
  write.table(medium,paste0("data/samples/RNA/",tissue,"/bed/medium_expression_gene.bed"),append = F,quote = F,row.names = F,col.names = F,sep = "\t")
  write.table(high,paste0("data/samples/RNA/",tissue,"/bed/high_expression_gene.bed"),append = F,quote = F,row.names = F,col.names = F,sep = "\t")
  
  stable <- ref[which(ref$V5 %in% RNA$X[which(RNA$Significant=="Stable")]),]
  write.table(stable,paste0("data/samples/RNA/",tissue,"/bed/Stable_expression_gene.bed"),append = F,quote = F,row.names = F,col.names = F,sep = "\t")
  Up <- ref[which(ref$V5 %in% RNA$X[which(RNA$Significant=="Up")]),]
  write.table(Up,paste0("data/samples/RNA/",tissue,"/bed/Up_expression_gene.bed"),append = F,quote = F,row.names = F,col.names = F,sep = "\t")
  Down <- ref[which(ref$V5 %in% RNA$X[which(RNA$Significant=="Down")]),]
  write.table(Down,paste0("data/samples/RNA/",tissue,"/bed/Down_expression_gene.bed"),append = F,quote = F,row.names = F,col.names = F,sep = "\t")
  
}

tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","ileum","mammarygland","iWAT")
for(tissue in tissues){
  three_levels_gene_expression(tissue)
}
