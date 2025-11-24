rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(ggrepel)
library(gridExtra)
library(grid) 
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
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
diff_number <-data.frame()
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  diff <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))  
  sig <- data.frame(Var1=c("Up","Stable","Down"),Freq=c(0,0,0))
  t_sig<-as.data.frame(table(diff$Significant))
  sig <- merge(sig, t_sig, by="Var1", all.x=TRUE) 
  sig$Freq.x <- ifelse(is.na(sig$Freq.y), sig$Freq.x, sig$Freq.y)  
  colnames(sig)[2] <- "Freq"
  sig <- sig[, -3]   
  sig$tissue <- tissue
  sig$antibody <- "RNA"
  diff_number<-rbind(diff_number,sig)
}
color <- read.table("data/samples/30_distinct_color.txt")
color <- color$V1
diff_number$tissue_label <- sapply(diff_number$tissue,tissue_label_change)
color <- setNames(color,sort(unique(diff_number$tissue_label)))
conditions <- c("Up","Down")
color <- setNames(c("#4589C8FF","#EE7C7AFF"),c("Down","Up"))


df <- diff_number[which(diff_number$Var1 %in% conditions),]
df$tissue_label <- factor(df$tissue_label,levels=rev(sort(unique(df$tissue_label))))
p <- ggplot(df,mapping = aes(x=Freq,y=tissue_label,fill = Var1))+
  # geom_bar(stat = "identity",position = position_dodge2(),aes(alpha = ifelse(Var1 == "Down", 0.9, 1)),color = "black")+
  geom_bar(stat = "identity",position = position_dodge2())+
  theme_bw()+ylab("")+
  xlab("RNA")+
  theme(  
    text = element_text(size = 10),  
    axis.text.y = element_text(size = 11),  # 改变 y 轴刻度标签的字体大小  
    axis.title.y = element_blank(),  
    axis.ticks.y = element_blank(),  
    legend.position = "none"  
  ) +
  scale_fill_manual(values = color) +
  theme(legend.position = "none") +
  scale_x_continuous(
    limits = c(0,7000), 
    breaks = c(0,7000)
  )
ggsave(paste0("result/figures/all_tissues_DEG_after_remove_batch_effect.pdf"), plot = p, width = 4, height = 6)


p_list <- list()
for(i in c(1:length(conditions))){
  condition <- conditions[i]
  df <- diff_number[which(diff_number$Var1==condition),]
  # df$tissue_label <- factor(df$tissue_label, levels=result$tissue_label)
  if(i == 1){
    p_list[[i]] <-  ggplot(df,mapping = aes(x=Freq,y=tissue_label,fill = tissue_label))+
      geom_bar(stat = "identity", position = position_dodge2())+
      theme_bw()+ylab("")+
      xlab(condition)+
      theme(  
        text = element_text(size = 10),  
        axis.text.y = element_text(size = 11),  # 改变 y 轴刻度标签的字体大小  
        axis.title.y = element_blank(),  
        axis.ticks.y = element_blank(),  
        legend.position = "none"  
      ) +
      scale_fill_manual(values = color) +
      theme(legend.position = "none") + 
      geom_text(aes(label = Freq), position = position_dodge2(width = 0.9), hjust = 0.1, size = 3) +
      xlim(0,7000)
  }else{
    p_list[[i]] <-  ggplot(df,mapping = aes(x=Freq,y=tissue_label,fill = tissue_label))+
      geom_bar(stat = "identity", position = position_dodge2())+
      theme_bw()+ylab("")+
      xlab(condition)+
      theme(  
        text = element_text(size = 10),  
        axis.text.y = element_text(size = 11),  # 改变 y 轴刻度标签的字体大小  
        axis.title.y = element_blank(),  
        axis.ticks.y = element_blank(),  
        legend.position = "none"  
      ) + 
      scale_fill_manual(values = color) +    
      theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),legend.position = "none") +   
      geom_text(aes(label = Freq), position = position_dodge2(width = 0.9), hjust = 0.1, size = 3) +
      xlim(0,7000)
  }
  combined_plot <- arrangeGrob(  
    grobs = p_list,  
    ncol = length(p_list),  
    widths = c(1.4,rep(1,(length(p_list)-1))),
    top = textGrob(paste0("Differential expression genes number"), gp = gpar(fontsize = 15, fontface = "bold"))  
  )  
  ggsave(paste0("result/Sup_figures/gene_expression_count_after_remove_batch_effect.pdf"), plot = combined_plot, width = 6, height = 6)
  grid.draw(combined_plot) 
}
