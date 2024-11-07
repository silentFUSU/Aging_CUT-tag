rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(stringr)
tissue <- "lung"
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
pairsqc_plot <- function(tissue){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  to_plot <- data.frame()
  for(i in c(1:length(search_table$sample_name))){
    sample <- search_table$sample_name[i]    
    x=read.table(paste0("data/samples/HiC/",tissue,"/4DN_pairs/",sample,"_cis_report/sample.plot_table.out"),sep="\t",stringsAsFactors=F,header=T)
    x <- x[,c("distance","log10prob")]
    colnames(x)[2] <- sample
    if(nrow(to_plot)==0){
      to_plot <- x
    }else{
      to_plot <- merge(to_plot,x,by="distance")
    }
  }
  to_plot <- reshape2::melt(to_plot,id.vars = "distance")
  colnames(to_plot)[2] <- "sample_name"
  to_plot <- merge(to_plot,search_table,by="sample_name")
  to_plot$age <- factor(to_plot$age, levels=c("3M","24M"))
  ylim=range(to_plot$value)
  ggplot(to_plot, aes(x = distance, y = value, color = sample_name)) +  
    geom_line() +  
    labs(x = "distance (10^x)", y = "Contact probability (10^y)") +  
    theme_minimal() + 
    theme(text = element_text(size = 20))+ xlim(3,max(to_plot$distance))+
    ggtitle(paste0(tissue_label_change(tissue)," Contact probability vs Distance"))
}

