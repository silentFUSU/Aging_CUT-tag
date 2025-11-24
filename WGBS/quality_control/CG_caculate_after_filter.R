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
search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
# summary <- data.frame(sample = as.character(),
#                       CG = as.numeric(),
#                       tissue = as.character())
summary <- data.frame()
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  t_search_table <- search_table[which(search_table$tissue == tissue),]
  samples <- t_search_table$sample_name
  for(sample in samples){
    df <- fread(paste0("data/samples/WGBS/",tissue,"/bdg/",sample,"_CpG.bdg"),sep = "\t")
    CG <- sum(df$V4)/sum(df$V5)
    t_summary <- data.frame(sample = sample, CG = CG,tissue = tissue)
    summary <- rbind(summary,t_summary)
  }
}
# write.csv(summary,"data/samples/WGBS/CG_manual.csv",row.names = F)

summary <- read.csv("data/samples/WGBS/CG_manual.csv")
summary$tissue <- sapply(summary$tissue,tissue_label_change)
summary <- merge(summary,search_table[,c("sample_name","age")],by.x="sample",by.y="sample_name")
result <- summary %>%
  group_by(tissue, age) %>%
  summarise(mean_CG = mean(CG, na.rm = TRUE), .groups = 'drop')
result$mean_CG <- result$mean_CG *100

result_young <- result[which(result$age=="3M"),]
result_old <- result[which(result$age=="24M"),]
to_plot <- merge(result_young,result_old,by="tissue")
to_plot$delta <- to_plot$mean_CG.y - to_plot$mean_CG.x
to_plot <- to_plot[order(to_plot$delta),]

to_plot_CG <- to_plot[,c(1,3,5)]
colnames(to_plot_CG)[c(2,3)] <- c("young","old")
to_plot_CG <- reshape2::melt(to_plot_CG)
color <- setNames(c("black","gray"),c("young","old"))
p <- ggplot(to_plot_CG, aes(x = tissue, y = value, fill = variable)) +  
  # geom_bar(stat = 'identity', aes(alpha = variable), position = position_dodge(width = 0.9),color="black") +
  geom_bar(stat = 'identity', position = position_dodge(width = 0.9),color="black") +
  theme_bw() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90,hjust = 1),
        text = element_text(size = 14),legend.title = element_blank()) +
  ylab("CG %")+
  ylim(0,100)+
  ggtitle("DNA methylation")
ggsave("result/Sup_figures/WGBS_CG_barplot.pdf",p,height = 6,width = 8)

color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(to_plot$tissue))
to_plot$tissue <- factor(to_plot$tissue,levels=to_plot$tissue)
tissues_order <- to_plot$tissue
to_plot$condition <- "Up"
to_plot$condition[which(to_plot$delta < 0 )] <- "Down"
color <- setNames(c("#f39b7f","#4dbbd5"),c("Up","Down"))
p <- ggplot(to_plot, aes(x = tissue, y = delta, fill = condition)) +  
  geom_bar(stat = 'identity') +   
  theme_bw() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Delta")+
  ylim(-6,6) +
  ggtitle(NULL) +
  guides(fill = FALSE) +
  geom_hline(yintercept = c(-1, 1), color = "black", linetype = "dashed")
ggsave("result/figures/WGBS_delta_barplot.pdf",p,width = 6,height = 4)
to_plot <- to_plot[,c("tissue","delta")]
rownames(to_plot) <- to_plot$tissue
to_plot <- to_plot[,-1,drop=F]
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- c(seq(-5, -0.61, length.out = 40), seq(-0.6, 0.6, length.out = 20), seq(0.61, 5, length.out = 40))
pheatmap::pheatmap(to_plot,cluster_rows = F,cluster_cols = F,breaks = breaks,color = color_palette,filename = "result/figures/WGBS_delta_heatmap.pdf",width = 2,height = 5)

