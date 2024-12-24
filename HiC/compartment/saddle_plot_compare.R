rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
tissue <- "brain"
resolution <- "50000"
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
tissue <- "brain"
sample <- "WJH-Cerebellum-96"
resolution <- "50000"
saddle_plot <- function(tissue,sample,resolution){
  df <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/cooler/",sample,"interaction_sum_count.csv"))
  left_average <- reshape2::melt(df[1:3,1:3])
  left_average <- mean(left_average$value)
  
  right_average <- reshape2::melt(df[38:40,38:40])
  right_average <- mean(right_average$value)
  
  data_long <- reshape2::melt(df, variable.name = "Column", value.name = "Value")  
  data_long$Row <- rep(1:nrow(df), times = ncol(df))  
  data_long$Row <- factor(data_long$Row,levels = nrow(df):1)
  data_long$Value[which(data_long$Value > 2)] <- 2
  data_long$Value[which(data_long$Value < 0.5)] <- 0.5
  p <- ggplot(data_long, aes(x = Column, y = as.factor(Row), fill = Value)) +  
    geom_tile() +  
    scale_fill_gradient2(low = "#3f72af", mid = "white", high = "#e84545", midpoint = 1, limits = c(0.5, 2)) + # 设置颜色渐变  
    theme_minimal() + # 使用简约主题  
    labs(x = NULL, y = NULL, fill = "Average obs/exp") + # 设置标签  
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    ggtitle(sample) +
    theme(  
      axis.title.x = element_blank(),  # 去除 X 轴标题  
      axis.title.y = element_blank(),  # 去除 Y 轴标题  
      axis.text.x = element_blank(),   # 去除 X 轴文本  
      axis.text.y = element_blank(),   # 去除 Y 轴文本  
      axis.ticks = element_blank()     # 去除刻度线  
    ) +  
    annotate("text", x = 1, y = 40, label = sprintf("%.2f", left_average), hjust = 0, vjust = 1,size=10) +  
    annotate("text", x = 40, y = 1, label = sprintf("%.2f", right_average), hjust = 1, vjust = 0,size=10)  +
    theme(text = element_text(size = 20))
  print(p)
}
search_table <- read.csv("data/samples/all/HiC_search_table.csv")
search_table <- search_table[which(search_table$tissue==tissue),]
samples <- search_table$sample_name
for(sample in samples){
  saddle_plot(tissue,sample,resolution)
}

saddle_compare <- function(tissue,resolution,compare){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  young <- search_table$sample_name[which(search_table$age=="3M")]
  old <- search_table$sample_name[which(search_table$age=="24M")]
  compare_table <- as.data.frame(expand.grid(young = young, old = old))
  for(i in c(1:nrow(compare_table))){
    young_saddle <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/cooler/",compare_table[i,"young"],"interaction_sum_count.csv"))
    old_saddle <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/cooler/",compare_table[i,"old"],"interaction_sum_count.csv"))
    difference_saddle <- old_saddle - young_saddle  
    
    data_long <- reshape2::melt(difference_saddle, variable.name = "Column", value.name = "Value")  
    data_long$Row <- rep(1:nrow(difference_saddle), times = ncol(difference_saddle))  
    data_long$Row <- factor(data_long$Row,levels = nrow(difference_saddle):1)
    data_long$Value[which(data_long$Value > 0.25)] <- 0.25
    data_long$Value[which(data_long$Value < -0.25)] <- -0.25
    p <- ggplot(data_long, aes(x = Column, y = as.factor(Row), fill = Value)) +  
      geom_tile() +  
      scale_fill_gradient2(low = "#3f72af", mid = "white", high = "#e84545", midpoint = 0, limits = c(-0.25, 0.25)) + # 设置颜色渐变  
      theme_minimal() + # 使用简约主题  
      labs(x = NULL, y = NULL, fill = "Difference") + # 设置标签  
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      ggtitle(paste0(compare_table[i,"old"],"-",compare_table[i,"young"])) +
      theme(  
        axis.title.x = element_blank(),  # 去除 X 轴标题  
        axis.title.y = element_blank(),  # 去除 Y 轴标题  
        axis.text.x = element_blank(),   # 去除 X 轴文本  
        axis.text.y = element_blank(),   # 去除 Y 轴文本  
        axis.ticks = element_blank()     # 去除刻度线  
      ) +
      theme(text = element_text(size = 20))
    print(p)
    } 
}
