rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(stringr)
library(dplyr)
library(dbplyr)
library(clusterProfiler)
library(GSVA)
library(enrichplot)
library(MASS)  
library(RANSAC)
options(scipen = 0) 
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


tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
                  "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT","ileum"))
summary_antibody <- data.frame()
region <- "HMD"
for(tissue in tissues){
  tissue_overall_summary <- read.csv("data/samples/WGBS/all/DNA_methylation_change_in_PMD_HMD_overall_summary.csv",row.names = 1)
  tissue_overall_summary <- tissue_overall_summary[which(tissue_overall_summary$V5==region),]
  df <- tissue_overall_summary[which(tissue_overall_summary$tissue==tissue_label_change(tissue)),]
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  young <- mean(df$methylation[which(df$sample %in% c(search_table$sample_name[which(search_table$age=="3M")]))])
  old <- mean(df$methylation[which(df$sample %in% c(search_table$sample_name[which(search_table$age=="24M")]))])
  t_summary_antibody <- data.frame(tissue=tissue_label_change(tissue),delta=old-young)  
  summary_antibody <- rbind(summary_antibody,t_summary_antibody)
}
to_plot <- summary_antibody 
to_plot <- to_plot[order(to_plot$delta),]
if(region=="PMD"){
  # tissue_order <- to_plot$tissue
  tissue_order <- c("Mammary Gland","Cecum","Thymus","Uterus","iWAT","Stomach","Skin","Spleen",
                    "Muscle","Bone Marrow","Liver","Ileum","Testis","Cortex","Jejunum","Tongue",
                    "Hippocampus","Colon","Bladder","Aorta","Cerebellum","Lung","Heart","Kidney","BAT","Ovary","Pancreas")
}
to_plot$condition <- "Up"
to_plot$condition[which(to_plot$delta < 0 )] <- "Down"
to_plot$tissue <- factor(to_plot$tissue,levels=tissue_order)
color <- setNames(c("#f39b7f","#4dbbd5"),c("Up","Down"))

p <- ggplot(to_plot, aes(x = tissue, y = delta, fill = condition)) +  
  geom_bar(stat = 'identity') +   
  theme_bw() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Delta")+
  ylim(-8,8) +
  ggtitle(region) +
  guides(fill = FALSE) +
  geom_hline(yintercept = c(-1, 1), color = "black", linetype = "dashed")
p
