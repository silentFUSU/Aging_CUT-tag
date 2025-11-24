rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(ggsignif)
library(stringr)
options(scipen = 999)  

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
    }
  }
  return(tissue_label)
}
resolution <- "20000"

tissues <- sort(c("brain","CB","kidney","liver","lung","bonemarrow","colon","heart","Hip","mammarygland","stomach","thymus","skin","muscle","cecum","ileum","pancreas","spleen"))
tissue_summary <- data.frame()
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  young_samples <- search_table$sample_name[which(search_table$age=="3M")]
  old_samples <- search_table$sample_name[which(search_table$age=="24M")]
  samples_list <- list(young=young_samples,old=old_samples)
  ages <- c("young","old")
  tad_summary <- data.frame()
  for(age in ages){
    for(i in c(1:length(samples_list[[age]]))){
      sample <- samples_list[[age]][i]
      df <- read.csv(paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",sample,"_",resolution,"_tads.csv"))
      if(nrow(tad_summary)==0){
        tad_summary <- data.frame(sample=sample,tissue=tissue,age=age,counts=nrow(df))
      }else{
        tad_summary <- rbind(tad_summary,data.frame(sample=sample,tissue=tissue,age=age,counts=nrow(df)))
      }
    }
  }
  rounded_average_counts <- tad_summary %>%
    group_by(age) %>%
    summarise(counts = round(mean(counts)))
  rounded_average_counts$tissue <- tissue_label_change(tissue)
  tissue_summary <- rbind(tissue_summary,rounded_average_counts)
}
tissue_summary$age <- factor(tissue_summary$age,levels = c("young","old"))

color <- read.table("data/samples/30_distinct_color.txt")
tissues_label <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
                   "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
tissues_label <- sapply(tissues_label, tissue_label_change)
color <- setNames(color$V1,sort(tissues_label))
color <- color[which(names(color) %in% sapply(tissues, tissue_label_change))]

p <- ggplot(tissue_summary, aes(x = age, y = counts)) +
  geom_boxplot(aes(fill = age), outlier.shape = NA, alpha = 0.5) +  # 绘制填充颜色基于age的箱线图
  geom_point(aes(color = tissue), size = 2) +  # 根据tissue颜色绘制点
  geom_line(aes(group = tissue, color = tissue), size = 0.8, alpha = 0.7) +  # 根据tissue连接点
  scale_color_manual(values = color)+
  geom_signif(comparisons = list(c("young", "old")), map_signif_level = TRUE,test = "t.test") +
  labs(x = "Age",
       y = "Counts",
       fill = "Age",
       color = "Tissue") +
  theme_bw()+
  ggtitle("TAD")+
  ylim(0,5000)
ggsave("result/figures/TAD_number.pdf",p,width = 5,height = 6)  

