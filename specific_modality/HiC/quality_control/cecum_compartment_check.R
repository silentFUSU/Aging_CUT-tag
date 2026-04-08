rm(list=ls())
.libPaths(c("/storrep/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storrep/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyverse)  
check_row <- function(row) {  
  all(row[-1] == row[-1][1])  
}  
options(scipen = 999)  
tissue <- "cecum"
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
    }else if(tissue_label == "Mammarygland"){
      tissue_label <- "Mammary Gland"
    }else if(tissue_label == "Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
} 


df_list <- list(rep1=list(),rep2=list())
search_table <- read.csv("data/samples/all/HiC_search_table.csv")
search_table <- search_table[which(search_table$tissue==tissue),]
search_table$rep <- c("rep1","rep1","rep2","rep2")
rep1_samples <- search_table$sample_name[which(search_table$rep=="rep1")]
rep2_samples <- search_table$sample_name[which(search_table$rep=="rep2")]
samples_list <- list(rep1=rep1_samples,rep2=rep2_samples)
for(rep in c("rep1","rep2")){
  for(i in c(1:length(samples_list[[rep]]))){
    sample <- samples_list[[rep]][i]
    if(sample=="WJH-Liver-103"){
      df <- read.table(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_",resolution,".PC1_recorrect.txt"))
    }else{
      df <- read.table(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_",resolution,".PC1.txt"))
    }
    df <- df[which(df$V2 %in% c(paste0("chr",c(1:19,"X")))),]
    df <- df[,c(1,6)]
    df$compartmeent <- ifelse(df[, 2] > 0, "A", "B")
    colnames(df)[3] <- sample
    df <- df[,c(1,3)]
    df_list[[rep]][[i]] <- df
    names(df_list[[rep]])[i]<-sample
  }
}


rep1_data <- Reduce(function(x, y) merge(x, y, by = "V1"), df_list[["rep1"]])  
rep2_data <- Reduce(function(x, y) merge(x, y, by = "V1"), df_list[["rep2"]])  

rep1_data <- rep1_data[apply(rep1_data, 1, check_row), ] 
rep1_data <- rep1_data[,c(1,2)]
colnames(rep1_data)[2] <- "rep1"

rep2_data <- rep2_data[apply(rep2_data, 1, check_row), ] 
rep2_data <- rep2_data[,c(1,2)]
colnames(rep2_data)[2] <- "rep2"

merged_data <- merge(rep1_data,rep2_data,by="V1")

merged_data$condition <- paste0(merged_data$rep1,"-",merged_data$rep2)
to_plot <- as.data.frame(table(merged_data$condition))
to_plot$Percentrep <- to_plot$Freq/sum(to_plot$Freq)*100
to_plot$Label <- paste0(to_plot$Var1, " (", round(to_plot$Percentrep, 1), "%)")
colors <- read.table("data/samples/7_distinct_color.txt")
colors <-setNames(colors$V1,to_plot$Label)
p <- ggplot(to_plot, aes(x = "", y = Freq, fill = Label)) +  
  geom_bar(width = 1, stat = "identity", color = "white") +  
  scale_fill_manual(values = colors)+
  coord_polar("y", start = 0) +  
  theme_void() + 
  labs(fill = NULL) +  
  ggtitle(paste0(tissue_label_change(tissue))) + 
  theme(legend.position = "right",plot.title = element_text(hjust = 0.5),text = element_text(size = 16))

