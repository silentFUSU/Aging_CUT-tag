rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)
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


tissues <- c("brain","CB","kidney","liver","lung","bonemarrow","colon","heart","Hip","mammarygland","stomach","thymus","muscle","skin","cecum","ileum","pancreas","spleen")

for(tissue in tissues){
  print(tissue)
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  samples <- search_table$sample_name[which(search_table$tissue==tissue)]
  chromosomes <- paste0("chr",c(1:19,"X"))
  resolution <- "50000"
  tissue_summary <- data.frame()
  for(sample in samples){
    print(sample)
    t_interaction_summary <- data.frame()
    compartment <- read.table(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_50000.PC1.txt"))
    compartment <- compartment[,c(1,6)]
    compartment$compartmeent <- ifelse(compartment[, 2] > 0, "A", "B")
    compartment <- compartment[,c(1,3)]
    colnames(compartment) <- c("label","compartment")
    for(chr in chromosomes){
      df <- read.table(paste0("data/samples/HiC/",tissue,"/ob_ex_matrix/",sample,"/",sample,"_",resolution,"_",chr,"_ob_ex_Matrix.txt"),skip = 1)
      rownames(df) <- df$V1
      df <- df[,-c(1)]
      colnames(df) <- c("bin",rownames(df))
      df <- reshape2::melt(df)
      colnames(df) <- c("bin1","bin2","value")
      df <- merge(df,compartment,by.x = "bin1",by.y = "label")
      colnames(df)[which(colnames(df)=="compartment")] <- "compartment_bin1"
      df <- merge(df,compartment,by.x = "bin2",by.y = "label")
      colnames(df)[which(colnames(df)=="compartment")] <- "compartment_bin2"
      t_interaction_summary <- rbind(t_interaction_summary,df)
    }
    t_interaction_summary$condition <- paste0(t_interaction_summary$compartment_bin1,"-",t_interaction_summary$compartment_bin2)
    t_interaction_summary$condition[which(t_interaction_summary$condition=="B-A")] <- "A-B"
    t_interaction_summary$log2_value <- log2(t_interaction_summary$value)
    saveRDS(t_interaction_summary,paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/",sample,"_compartment_50000_interaction_bin_level.rds"))
    interaction_summary <- t_interaction_summary %>%
      group_by(condition) %>%
      summarize(median_log2_value = median(log2_value, na.rm = TRUE))
    interaction_summary$sample <- sample
    tissue_summary <- rbind(tissue_summary,interaction_summary)
  }
  write.csv(tissue_summary,paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/compartment_50000_interaction_median.csv"))
}


summary <- data.frame()

for(tissue in tissues){
  df <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/compartment_50000_interaction_median.csv"),row.names = 1)
  summary <- rbind(summary,df)
  }

search_table <- read.csv("data/samples/all/HiC_search_table.csv")
summary <- merge(summary,search_table,by.x="sample",by.y="sample_name")

summary$condition <- factor(summary$condition,levels = c("A-A","B-B","A-B"))
summary$age <- factor(summary$age, levels=c("3M","24M"))
summary$tissue_label <- sapply(summary$tissue, tissue_label_change)
ggplot(summary,aes(x=condition,y=median_log2_value,color = age,shape=age))+
  geom_boxplot(aes(fill = age), outlier.shape = NA, alpha = 0.5) +
  ggtitle(paste0("Compartment interaction"))+ylim(-1.2,0.7)+
  theme_bw()+xlab("")+ylab("log2(ob/ex)") +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
t.test(summary$median_log2_value[which(summary$age=="24M" & summary$condition=="A-A")],summary$median_log2_value[which(summary$age=="3M" &  summary$condition=="A-A")])
t.test(summary$median_log2_value[which(summary$age=="24M" & summary$condition=="B-B")],summary$median_log2_value[which(summary$age=="3M" &  summary$condition=="B-B")])
t.test(summary$median_log2_value[which(summary$age=="24M" & summary$condition=="A-B")],summary$median_log2_value[which(summary$age=="3M" &  summary$condition=="A-B")])

summary_avg <- summary %>%
  group_by(condition, tissue_label, age) %>%
  summarise(mean_median_log2_value = mean(median_log2_value, na.rm = TRUE))
summary_avg$cond_age <- interaction(summary_avg$condition, summary_avg$age)
summary_avg$cond_age <- factor(summary_avg$cond_age,levels=c("A-A.3M","A-A.24M","B-B.3M","B-B.24M","A-B.3M","A-B.24M"))

color <- read.table("data/samples/30_distinct_color.txt")
tissues_label <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
                   "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
tissues_label <- sapply(tissues_label, tissue_label_change)
color <- setNames(color$V1,sort(tissues_label))
color <- color[which(names(color) %in% sapply(tissues, tissue_label_change))]


p <- ggplot(summary_avg,aes(x=cond_age,y=mean_median_log2_value))+
  geom_boxplot(aes(fill = age), outlier.shape = NA, alpha = 0.5)+
  geom_point(aes(color = tissue_label), size = 2) +
  geom_line(aes(group = interaction(tissue_label, condition), color = tissue_label), alpha = 0.5) +
  scale_color_manual(values = color)+
  ggtitle(paste0("Compartment interaction"))+ylim(-1.2,0.7)+
  theme_bw()+xlab("")+ylab("log2(ob/ex)") +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
ggsave("result/figures/Compartment_interaction.pdf",p,width = 5,height = 6)  
t.test(summary_avg$mean_median_log2_value[which(summary_avg$age=="24M" & summary_avg$condition=="A-A")],summary_avg$mean_median_log2_value[which(summary_avg$age=="3M" &  summary_avg$condition=="A-A")])
t.test(summary_avg$mean_median_log2_value[which(summary_avg$age=="24M" & summary_avg$condition=="B-B")],summary_avg$mean_median_log2_value[which(summary_avg$age=="3M" &  summary_avg$condition=="B-B")])
t.test(summary_avg$mean_median_log2_value[which(summary_avg$age=="24M" & summary_avg$condition=="A-B")],summary_avg$mean_median_log2_value[which(summary_avg$age=="3M" &  summary_avg$condition=="A-B")])


# ggplot(summary,aes(x=condition,y=median_log2_value,color = age,shape=age))+
#   geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.3), size = 2, alpha = 0.7)+
#   ggtitle(paste0("Compartment interaction"))+ylim(-1.2,0.7)+
#   theme_bw()+xlab("")+ylab("log2(ob/ex)") + facet_wrap(~ tissue_label, ncol = 4)+
#   theme(
#     plot.title = element_text(size = 16, hjust = 0.5),
#     axis.title.y = element_text(size = 14),
#     axis.text = element_text(size = 12),
#     strip.text = element_text(size = 14),
#     legend.title = element_text(size = 12),
#     legend.text = element_text(size = 10)
#   )
