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
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
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
tissue <- "lung"

tissues <- c("brain","CB","kidney","liver","lung","bonemarrow","colon","heart","Hip","mammarygland","stomach","thymus","muscle","skin","cecum","ileum","pancreas","spleen")
resolution=50000
min=1000000
max=60000000
for(tissue in tissues){
  print(tissue)
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  samples <- search_table$sample_name[which(search_table$tissue==tissue)]
  chromosomes <- paste0("chr",c(1:19,"X"))
  resolution <- "50000"
  tissue_summary <- data.frame()
  for(sample in samples){
    t_interaction_summary <- readRDS(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/",sample,"_compartment_",resolution,"_interaction_bin_level.rds"))
    t_interaction_summary$chr <- sapply(strsplit(as.character(t_interaction_summary$bin1), "-"), `[`, 1)
    t_interaction_summary$bin1_start <- sapply(strsplit(as.character(t_interaction_summary$bin1), "-"), `[`, 2)
    t_interaction_summary$bin2_start <- sapply(strsplit(as.character(t_interaction_summary$bin2), "-"), `[`, 2)
    t_interaction_summary$bin1_start <- as.numeric(t_interaction_summary$bin1_start)
    t_interaction_summary$bin2_start <- as.numeric(t_interaction_summary$bin2_start)
    # t_interaction_summary <- t_interaction_summary[which(abs(t_interaction_summary$bin1_start - t_interaction_summary$bin2_start) > 1000000 & abs(t_interaction_summary$bin1_start - t_interaction_summary$bin2_start) < 40000000),]
    t_interaction_summary <- t_interaction_summary[which(abs(t_interaction_summary$bin1_start - t_interaction_summary$bin2_start) >= min & abs(t_interaction_summary$bin1_start - t_interaction_summary$bin2_start) <= max),]
    interaction_summary <- t_interaction_summary %>%
      group_by(condition,chr) %>%
      summarize(median_log2_value = median(log2_value, na.rm = TRUE))
    interaction_summary$chr <- factor(interaction_summary$chr, levels = paste0("chr",c(1:19,"X","Y")))
    interaction_summary$sample <- sample
    tissue_summary <- rbind(tissue_summary,interaction_summary)
  }
  write.csv(tissue_summary,paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/compartment_",resolution,"_interaction_median_per_chromosome_filter_distance_",(min/1000000),"mb_",(max/1000000),"mb.csv"))
  tissue_summary <- merge(tissue_summary,search_table,by.x="sample",by.y="sample_name")
  conditions <- c("A-A","A-B","B-B")
  p_list <- list()
  for(condition in conditions){
    to_plot <- tissue_summary[which(tissue_summary$condition==condition),]
    to_plot <- to_plot[order(to_plot$chr),]
    t_search_table <- search_table[which(search_table$tissue==tissue),]
    t_search_table$age <- factor(t_search_table$age,levels=c("3M","24M"))
    t_search_table <- t_search_table[order(t_search_table$age),]
    t_search_table$color <- c("#f38181", "#ff2e63", "#112d4e", "#3f72af")
    color <- setNames(t_search_table$color,t_search_table$sample_name)
    p_list[[condition]] <- ggplot(to_plot, aes(x = chr, y = median_log2_value,  color = sample, group = sample)) +
      geom_line() +
      scale_color_manual(values = color)+
      geom_point() +
      xlab(NULL)+
      ylab("log2(ob/ex)")+
      ggtitle(tissue_label_change(tissue),paste0(condition," interaction"))+
      theme_minimal()+
      theme(
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  }
  combined_plot <- plot_a_list(p_list,no_of_rows = 3,no_of_cols = 1)
  ggsave(paste0("result/HiC/",tissue,"/compartment/",tissue,"_homer_compartment_",resolution,"_interaction_ob_ex_change_per_chromosome_filter_distance_",(min/1000000),"mb_",(max/1000000),"mb.png"),combined_plot,width = 8,height = 10,type="cairo")
}

p_list <- list()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/compartment_",resolution,"_interaction_median_per_chromosome_filter_distance_",(min/1000000),"mb_",(max/1000000),"mb.csv"),row.names = 1)
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  to_plot <- merge(df,search_table,by.x="sample",by.y="sample_name")
  search_table$age <- factor(search_table$age,c("3M","24M"))
  search_table <- search_table[order(search_table$age),]
  search_table$color <- c("#f38181","#ff2e63","#112d4e","#3f72af")
  color <- setNames(search_table$color,search_table$sample_name)
  to_plot$sample <- factor(to_plot$sample,levels = search_table$sample_name)
  to_plot$condition <- factor(to_plot$condition, levels = c("A-A","B-B","A-B"))
  p_list[[tissue]] <- ggplot(to_plot,aes(x=condition,y=median_log2_value,color = sample))+
    geom_boxplot() +
    scale_color_manual(values = color) +
    ggtitle(paste0(tissue_label_change(tissue)," compartment interaction"))+
    theme_bw()+xlab("")+ylab("log2(ob/ex)") +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5),
      axis.title.y = element_text(size = 14),
      axis.text = element_text(size = 12),
      strip.text = element_text(size = 14),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  }
combined_plot <- plot_a_list(p_list,no_of_rows = 4,no_of_cols = 4)
ggsave(paste0("result/HiC/all_tissues_50000_homer_compartment_interaction_median_per_chromosome_filter_distance_",(min/1000000),"mb_",(max/1000000),"mb.png"),combined_plot,width = 25,height = 25,type="cairo")


summary <- data.frame()
resolution <- "50000"
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/compartment_",resolution,"_interaction_median_per_chromosome_filter_distance_",(min/1000000),"mb_",(max/1000000),"mb.csv"),row.names = 1)
  summary <- rbind(summary,df)
}
search_table <- read.csv("data/samples/all/HiC_search_table.csv")
summary <- merge(summary,search_table,by.x="sample",by.y="sample_name")
summary$condition <- factor(summary$condition,levels = c("A-A","B-B","A-B"))
summary$age <- factor(summary$age, levels=c("3M","24M"))
summary$tissue_label <- sapply(summary$tissue, tissue_label_change)
summary <- summary[order(summary$age),]
ggplot(summary,aes(x=condition,y=median_log2_value,color = age,shape=age))+
  geom_boxplot() +
  ggtitle(paste0("Compartment interaction"))+
  theme_bw()+xlab("")+ylab("log2(ob/ex)") + facet_wrap(~ tissue_label, ncol = 4)+
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

conditions <- c("A-A","A-B","B-B")
p_value_summary <- data.frame(
  `A-A` = rep(NA, length(tissues)),
  `A-B` = rep(NA, length(tissues)),
  `B-B` = rep(NA, length(tissues))
)
colnames(p_value_summary) <- c("A-A","A-B","B-B")
rownames(p_value_summary) <- sapply(tissues, tissue_label_change)

for(condition in conditions){
  for(tissue in tissues){
    df <- summary[which(summary$condition==condition & summary$tissue==tissue),]
    test <- wilcox.test(df$median_log2_value[which(df$age=="3M")], df$median_log2_value[which(df$age=="24M")])
    p_value_summary[tissue_label_change(tissue),condition] <- test$p.value
  }
}

mark_significance <- function(p_value) {
  if (is.na(p_value)) {
    return(NA)
  } else if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return(NA)
  }
}
p_value_summary <- p_value_summary %>%
  mutate(
    `A-A` = sapply(`A-A`, mark_significance),
    `A-B` = sapply(`A-B`, mark_significance),
    `B-B` = sapply(`B-B`, mark_significance)
  )

# ggplot(summary,aes(x=condition,y=median_log2_value,color = age,shape=age))+
#   geom_boxplot() +
#   ggtitle(paste0("Compartment interaction"))+
#   theme_bw()+xlab("")+ylab("log2(ob/ex)") + 
#   theme(
#     plot.title = element_text(size = 16, hjust = 0.5),
#     axis.title.y = element_text(size = 14),
#     axis.text = element_text(size = 12),
#     strip.text = element_text(size = 14),
#     legend.title = element_text(size = 12),
#     legend.text = element_text(size = 10)
#   )
# 
# p_value_summary <- data.frame(
#   `A-A` = rep(NA, 1), 
#   `A-B` = rep(NA, 1), 
#   `B-B` = rep(NA, 1)  
# )
# for(condition in conditions){
#   df <- summary[which(summary$condition==condition),]
#   test <- wilcox.test(df$median_log2_value[which(df$age=="3M")], df$median_log2_value[which(df$age=="24M")])
#   p_value_summary[1,condition] <- test$p.value
# }
