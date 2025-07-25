rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(tidyverse)  
library(dplyr) 
library(limma)
library(data.table)
library(edgeR)
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect",axis_titles = "collect",axes = "collect")
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
    }else if(tissue_label == "Mammarygland"){
      tissue_label <- "Mammary Gland"
    }else if(tissue_label == "Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
} 
resolution <- "20000"
tissue <- "lung"
TAD_diff_in_compartment <- function(tissue,resolution){
  TAD <- read.csv(paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",tissue,"_redundant_",resolution,"_TAD_diff.csv"))
  colnames(TAD)[1] <- "Geneid"
  increase <- TAD[which(TAD$Significant=="Up"),"Geneid",drop=F]
  decrease <- TAD[which(TAD$Significant=="Down"),"Geneid",drop=F]
  increase <- increase %>%
    separate(Geneid, into = c("Chr", "Start", "End"), sep = "-")
  increase$Start <- as.numeric(increase$Start)
  increase$End <- as.numeric(increase$End)
  increase <- as.data.table(increase)
  setDT(increase)
  setkey(increase,Chr,Start,End)
  
  decrease <- decrease %>%
    separate(Geneid, into = c("Chr", "Start", "End"), sep = "-")
  decrease$Start <- as.numeric(decrease$Start)
  decrease$End <- as.numeric(decrease$End)
  decrease <- as.data.table(decrease)
  setDT(decrease)
  setkey(decrease,Chr,Start,End)
  
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  samples <- search_table$sample_name
  increase_summary <- data.frame()
  decrease_summary <- data.frame()
  for(sample in samples){
    df <-  read.table(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_50000.PC1.txt"))
    df <- df[which(df$V2 %in% c(paste0("chr",c(1:19,"X")))),]
    df <- df[,c(2:4,6)]
    df$compartmeent <- ifelse(df[, 4] > 0, "A", "B")
    df$V3 <- df$V3 + 1
    df <- as.data.table(df)
    setDT(df)
    setkey(df,V2,V3,V4)
    overlaps <- foverlaps(increase, df, type = "any", nomatch = 0L)  
    t_increase_summary <- as.data.frame(table(overlaps$compartmeent))
    t_increase_summary$sample <- sample
    t_increase_summary$percent <- t_increase_summary$Freq/sum(t_increase_summary$Freq)*100
    overlaps <- foverlaps(decrease, df, type = "any", nomatch = 0L)
    t_decrease_summary <- as.data.frame(table(overlaps$compartmeent))
    t_decrease_summary$sample <- sample
    t_decrease_summary$percent <- t_decrease_summary$Freq/sum(t_decrease_summary$Freq)*100
    
    increase_summary <- rbind(increase_summary,t_increase_summary)  
    decrease_summary <- rbind(decrease_summary,t_decrease_summary)
  }
  
  p1 <- ggplot(increase_summary, aes(x = sample, y = percent, fill = Var1)) +  
    geom_bar(stat = 'identity',colour = "white") +   
    theme_minimal() +   
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 20),legend.title = element_blank()) +
    ylab("Proportion")+
    ggtitle(paste0(tissue_label_change(tissue)),"Increased TAD")
  
  p2 <- ggplot(decrease_summary, aes(x = sample, y = percent, fill = Var1)) +  
    geom_bar(stat = 'identity',colour = "white") +   
    theme_minimal() +   
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 20),legend.title = element_blank()) +
    ylab("Proportion")+
    ggtitle(paste0(tissue_label_change(tissue)),"Decreased TAD")
  combined_plot <- plot_a_list(list(p1,p2),no_of_rows = 2,no_of_cols = 1)
  ggsave(paste0("result/HiC/",tissue,"/differential_analysis/changed_TAD_compartment_proportion.png"),combined_plot,height = 8,width = 5,type="cairo")
}
tissues <- c("brain","CB","kidney","lung","heart","Hip","mammarygland","stomach","thymus")
for(tissue in tissues){
  TAD_diff_in_compartment(tissue,resolution)
}
