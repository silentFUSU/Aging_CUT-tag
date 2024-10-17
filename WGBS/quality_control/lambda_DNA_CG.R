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
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
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
    }
  }
  return(tissue_label)
}

tissue <- "BAT"
CpG_overview <- function(tissue){
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  search_table <- search_table[which(search_table$tissue == tissue),]
  df_list <- list()
  cpg_sum <- data_frame(CG = as.numeric(),
                        depth = as.numeric(),
                        sample = as.character())
  for(i in c(1:nrow(search_table))){
    df_list[[i]] <- fread(paste0("data/samples/WGBS/",tissue,"/bdg/",search_table$sample_name[i],"_CpG.bdg"),sep = "\t")
    df_list[[i]] <- df_list[[i]][which(df_list[[i]]$V1 == "chrL"),]
    colnames(df_list[[i]])[5] <- "depth"
    df_list[[i]]$percent <- df_list[[i]]$V4/df_list[[i]]$depth*100
    df_list[[i]][, label := paste0(V1, "-", V2, "-", V3)]
    names(df_list)[i] <- search_table$sample_name[i]
    t_cpg_sum <- data_frame(CG=sum(df_list[[i]]$V4), depth=sum(df_list[[i]]$depth),sample=search_table$sample_name[i])
    cpg_sum <- rbind(cpg_sum,t_cpg_sum)
  }
  merge_list <- lapply(names(df_list), function(name) {  
    df_list[[name]] %>%  
      filter(depth > 0) %>%  
      select(label, !!name := percent)
  }) 
  merge_list <- setNames(merge_list, names(df_list))
  merge_df <- Reduce(function(x, y) merge(x, y, by = "label",all=TRUE), merge_list)  
  merge_df_to_plot <- reshape2::melt(merge_df)
  merge_df_to_plot <- merge_df_to_plot %>%  
    filter(!is.na(value))  
  search_table$age[which(search_table$age == "3M")] <- "young"
  search_table$age[which(search_table$age == "24M")] <- "old"
  search_table$age <- factor(search_table$age, levels= c("young","old"))
  search_table <- search_table[order(search_table$age),]
  colnames(search_table)[3] <- "variable"
  variable_order <- paste0(search_table$variable,"-",search_table$mouse_ID,"-",search_table$age)
  merge_df_to_plot <- merge(merge_df_to_plot, search_table, by = "variable" )
  merge_df_to_plot$variable_label <- paste0(merge_df_to_plot$variable, "-", merge_df_to_plot$mouse_ID, "-", merge_df_to_plot$age)
  merge_df_to_plot$variable_label <- factor(merge_df_to_plot$variable_label, levels = variable_order)
  t <- t.test(merge_df_to_plot$value[which(merge_df_to_plot$age=="old")],merge_df_to_plot$value[which(merge_df_to_plot$age=="young")])
  young_cg <- sum(cpg_sum$CG[which(cpg_sum$sample)])
  p <- ggplot(merge_df_to_plot, aes(x = variable_label, y = value, fill= age)) +  
    geom_violin(adjust = 2.5) +          
    scale_fill_brewer(palette = "Pastel1") +
    geom_boxplot(width = 0.1, color = "black", fill = "white", outlier.shape = NA) +  
    theme_minimal()+
    ggtitle(tissue_label_change(tissue))+
    theme(text = element_text(size = 20),legend.position = "none")+
    labs(x = NULL,y = "CpG%") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    annotate("text", x = Inf, y = -Inf, label = paste("p-value =",  format(t$p.value, scientific = TRUE, digits = 3)  ),   
             hjust = 1.1, vjust = -1.1, size = 5, colour = "red") +
    annotate("text", x = Inf, y = Inf, label = paste("old mean =",  round(t$estimate[[1]],2)  ),   
             hjust = 1.1, vjust = 1.1, size = 5, colour = "red") +
    annotate("text", x = -Inf, y = Inf, label = paste("young mean =",  round(t$estimate[[2]],2)  ),   
             hjust = 0, vjust = 1.1, size = 5, colour = "red")
  return(p)
}

lambda_conversion <- function(tissue){
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  search_table <- search_table[which(search_table$tissue == tissue),]
  df_list <- list()
  for(i in c(1:nrow(search_table))){
    df_list[[i]] <- fread(paste0("data/samples/WGBS/",tissue,"/bdg/",search_table$sample_name[i],"_CpG.bdg"),sep = "\t")
    df_list[[i]] <- df_list[[i]][which(df_list[[i]]$V1 == "chrL"),]
    colnames(df_list[[i]])[5] <- "depth"
    df_list[[i]]$percent <- df_list[[i]]$V4/df_list[[i]]$depth*100
    df_list[[i]][, label := paste0(V1, "-", V2, "-", V3)]
    names(df_list)[i] <- search_table$sample_name[i]
    t_cpg_sum <- data_frame(CG=sum(df_list[[i]]$V4), depth=sum(df_list[[i]]$depth),sample=search_table$sample_name[i])
    cpg_sum <<- rbind(cpg_sum,t_cpg_sum)
  }
}
cpg_sum <- data_frame(CG = as.numeric(),
                      depth = as.numeric(),
                      sample = as.character())
tissues <- c("liver","lung","mammarygland","kidney","ileum","Hip","skin","bonemarrow","jejunum","colon","ovary","CB","BAT","thymus","testis")
for(tissue in tissues){
  lambda_conversion(tissue)
}
