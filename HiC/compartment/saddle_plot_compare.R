rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(dplyr)
library(stringr)
tissue <- "Hip"
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
# tissue <- "brain"
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}
saddle_plot <- function(tissue,sample,resolution){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$sample_name==sample),]
  if(search_table$age=="3M"){
    age_label <- "Young"
  }else{
    age_label <- "Old"
  }
  df <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/cooler/",sample,"interaction_sum_count_homer.csv"))
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
    ggtitle(paste(sample,age_label)) +
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
  return(p)
}
tissues <- c("kidney","colon","liver","lung","brain","CB")
tissues <- c("bonemarrow","stomach","heart")
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  samples <- search_table$sample_name
  p_list <- list()
  for(sample in samples){
    p_list[[sample]] <- saddle_plot(tissue,sample,resolution)
  }
  combined_plot <- plot_a_list(p_list,2,2)
  ggsave(paste0("result/HiC/",tissue,"/compartment/",tissue,"_compartment_saddle_plot.png"),combined_plot,width = 12,height = 10,type="cairo")
  
}

saddle_compare <- function(tissue,resolution){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  young <- search_table$sample_name[which(search_table$age=="3M")]
  old <- search_table$sample_name[which(search_table$age=="24M")]
  compare_table <- as.data.frame(expand.grid(young = young, old = old))
  p_list <- list()
  for(i in c(1:nrow(compare_table))){
    young_saddle <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/cooler/",compare_table[i,"young"],"interaction_sum_count_homer.csv"))
    old_saddle <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/cooler/",compare_table[i,"old"],"interaction_sum_count_homer.csv"))
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
      ggtitle(paste0(compare_table[i,"old"]," Old -",compare_table[i,"young"]," Young")) +
      theme(  
        axis.title.x = element_blank(),  # 去除 X 轴标题  
        axis.title.y = element_blank(),  # 去除 Y 轴标题  
        axis.text.x = element_blank(),   # 去除 X 轴文本  
        axis.text.y = element_blank(),   # 去除 Y 轴文本  
        axis.ticks = element_blank()     # 去除刻度线  
      ) +
      theme(text = element_text(size = 12))
    print
    p_list[[i]] <- p
  } 
  combined_plot <- plot_a_list(p_list,2,2)
  ggsave(paste0("result/HiC/",tissue,"/compartment/",tissue,"_cooler_compartment_saddle_plot_compare_homer.png"),combined_plot,width = 11,height = 10,type="cairo")
  
  }

for(tissue in tissues){
  saddle_compare(tissue,resolution)
}
tissues <- c("brain","CB", "kidney", "liver", "lung", "bonemarrow", "colon", "heart", "Hip", "mammarygland", "stomach", "thymus")
smooth_3x3 <- function(mat) {  
  nrow_mat <- nrow(mat)  
  ncol_mat <- ncol(mat)  
  padded_mat <- matrix(NA, nrow = nrow_mat + 2, ncol = ncol_mat + 2)  
  padded_mat[2:(nrow_mat + 1), 2:(ncol_mat + 1)] <- as.matrix(mat)
  smoothed_mat <- matrix(0, nrow = nrow_mat, ncol = ncol_mat)  
  for (i in 1:nrow_mat) {  
    for (j in 1:ncol_mat) {  
      sub_matrix <- padded_mat[i:(i + 2), j:(j + 2)]  
      smoothed_mat[i, j] <- mean(c(as.matrix(sub_matrix)),na.rm = T)
    }  
  }  
  
  return(smoothed_mat)  
}  
saddle_compare_average <- function(tissue,resolution){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  young <- search_table$sample_name[which(search_table$age=="3M")]
  old <- search_table$sample_name[which(search_table$age=="24M")]
  for(i in c(1:length(young))){
    if(i == 1){
      young_data <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/cooler/",young[i],"interaction_sum_count_homer.csv"))
    }else{
      df <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/cooler/",young[i],"interaction_sum_count_homer.csv"))
      young_data <- young_data + df
    }
  }
  young_data <- young_data/length(young)
  
  for(i in c(1:length(old))){
    if(i == 1){
      old_data <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/cooler/",old[i],"interaction_sum_count_homer.csv"))
    }else{
      df <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/cooler/",old[i],"interaction_sum_count_homer.csv"))
      old_data <- old_data + df
    }
  }
  old_data <- old_data/length(old)

  difference_saddle <- old_data - young_data
  smoothed_difference_saddle <- smooth_3x3(difference_saddle)  
  difference_saddle <- as.data.frame(smoothed_difference_saddle)
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
    ggtitle(paste0(tissue_label_change(tissue)," Old vs Young")) +
    theme(  
      axis.title.x = element_blank(),  # 去除 X 轴标题  
      axis.title.y = element_blank(),  # 去除 Y 轴标题  
      axis.text.x = element_blank(),   # 去除 X 轴文本  
      axis.text.y = element_blank(),   # 去除 Y 轴文本  
      axis.ticks = element_blank()     # 去除刻度线  
    ) +
    theme(text = element_text(size = 12))
  return(p)
}
tissues <- c("lung","thymus","mammarygland","heart","stomach","liver","kidney","bonemarrow","CB","brain","Hip","colon")
p_list <- list()
for(tissue in tissues){
  p_list[[tissue]] <- saddle_compare_average(tissue,resolution)
}
combined_plot <- plot_a_list(p_list,4,3)
ggsave(paste0("result/HiC/all_tissues_homer_compartment_saddle_plot_smooth.png"),combined_plot,width = 10,height = 13,type="cairo")


saddle_summary <- data.frame()
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  tissue_saddle_summary <- data.frame()
  for(sample in search_table$sample_name){
    df <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/cooler/",sample,"interaction_sum_count_homer.csv"))
    BB_average <- reshape2::melt(df[1:3,1:3])
    BB_average <- mean(BB_average$value)
    
    AA_average <- reshape2::melt(df[38:40,38:40])
    AA_average <- mean(AA_average$value)
    
    AB_average <- reshape2::melt(df[1:3,38:40])  
    AB_average <- mean(AB_average$value)
    t_saddle_summary <- data.frame(sample=sample,
                                   tissue=tissue_label_change(tissue),
                                   age=search_table$age[which(search_table$sample_name==sample)],
                                   BB_average=BB_average,
                                   AA_average=AA_average,
                                   AB_average=AB_average)    

    tissue_saddle_summary <- rbind(tissue_saddle_summary,
                            t_saddle_summary)
  }
  columns_to_center <- c("BB_average", "AA_average", "AB_average")  
  tissue_saddle_summary[,columns_to_center] <- lapply(tissue_saddle_summary[,columns_to_center], function(x) x - mean(x))  
  saddle_summary <- rbind(saddle_summary,tissue_saddle_summary)
}
saddle_summary$age[which(saddle_summary$age=="3M")] <- "Young"
saddle_summary$age[which(saddle_summary$age=="24M")] <- "Old"
saddle_summary$age <- factor(saddle_summary$age,levels=c("Young","Old")) 
rownames(saddle_summary) <- saddle_summary$sample
annotation <- saddle_summary[,c("tissue","age")]

conditions <- c("BB_average","AA_average","AB_average")
for(condition in conditions){
  to_plot <- saddle_summary[,c("sample","tissue","age",condition)]
  to_plot$age <- factor(to_plot$age,levels=c("Young","Old")) 
  to_plot <- to_plot %>%
    arrange(tissue,age)
  to_plot$age <- paste0(to_plot$age,c(1,2))
  to_plot <- to_plot[,-1]
  to_plot <-to_plot %>%  
    pivot_wider(names_from = age, values_from = colnames(to_plot)[3])
  to_plot <- as.data.frame(to_plot)
  rownames(to_plot) <- to_plot$tissue
  to_plot <- to_plot[,-1]
  pheatmap::pheatmap(to_plot,cluster_cols = F,main = condition, breaks = seq(-0.1, 0.1, length.out = 101))
}


