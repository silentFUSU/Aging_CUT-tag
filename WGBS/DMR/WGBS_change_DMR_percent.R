rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
library(ggplot2)
set.seed(1)
tissues <- c("liver","lung","mammarygland","kidney","ileum","Hip")
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

DMR_percent <- function(tissues){
  DMR_summary <- data.frame(condition = as.character(),
                            count = as.numeric(),
                            percent = as.numeric(),
                            tissue = as.character())
  for(i in c(1:length(tissues))){
    tissue <- tissues[i]
    DMR <- read.table(paste0("data/samples/WGBS/",tissue,"/DSS_table/",tissue,"_DMR_delta0.txt"),header = T)
    increase <- nrow(DMR[which(DMR$areaStat > 0),])
    decrease <- nrow(DMR[which(DMR$areaStat < 0),])
    tissue_label <- tissue_label_change(tissue)
    t_DMR_percent <- data.frame(condition = c("Increase","Decrease"),
                                count = c(increase,decrease),
                                percent = c(increase/nrow(DMR)*100, decrease/nrow(DMR)*100),
                                tissue = c(tissue_label,tissue_label)) 
    DMR_summary <- rbind(DMR_summary,t_DMR_percent)
  }
  p1 <- ggplot(DMR_summary, aes(x = tissue, y = count, fill = condition)) +  
    geom_bar(stat = 'identity') +   
    theme_minimal() +   
    scale_fill_brewer(palette = "Pastel1") +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 20),legend.title = element_blank()) +
    ylab("Counts")+
    ggtitle("DMR")
  print(p1)
  
  DMR_summary$position <- 100
  DMR_summary$position[which(DMR_summary$condition=="Increase")] <- DMR_summary$percent[which(DMR_summary$condition == "Increase")]
  p2 <- ggplot(DMR_summary, aes(x = tissue, y = percent, fill = condition)) +  
    geom_bar(stat = 'identity') +   
    theme_minimal() +   
    scale_fill_brewer(palette = "Pastel1") +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 20),legend.title = element_blank()) +
    geom_text(data = subset(DMR_summary),   
              aes(label = paste0(round(percent, 2),"%"), y = position),   
              color = "black", size = 5, vjust = 0.5) + 
    ylab("Percent (%)")+
    ggtitle("DMR")
  print(p2)
}
DMR_percent(tissues)

DMR_overlap_heterochromatin_switch_percent <- function(tissues){
  DMR_summary <- data.frame(condition = as.character(),
                            count = as.numeric(),
                            percent = as.numeric(),
                            tissue = as.character())
  for(i in c(1:length(tissues))){
    tissue <- tissues[i]
    increase <- read.table(paste0("data/samples/WGBS/",tissue,"/DSS_table/bed/",tissue,"_DMR_increase_overlap_heterochromatin_switch.bed"))
    decrease <- read.table(paste0("data/samples/WGBS/",tissue,"/DSS_table/bed/",tissue,"_DMR_decrease_overlap_heterochromatin_switch.bed"))
    increase <- nrow(increase)
    decrease <- nrow(decrease)
    tissue_label <- tissue_label_change(tissue)
    t_DMR_percent <- data.frame(condition = c("Increase","Decrease"),
                                count = c(increase,decrease),
                                percent = c(increase/(increase+decrease)*100, decrease/(increase+decrease)*100),
                                tissue = c(tissue_label,tissue_label)) 
    DMR_summary <- rbind(DMR_summary,t_DMR_percent)
  }
  ggplot(DMR_summary, aes(x = tissue, y = count, fill = condition)) +  
    geom_bar(stat = 'identity') +   
    theme_minimal() +   
    scale_fill_brewer(palette = "Pastel1") +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 20),legend.title = element_blank()) +
    ylab("Counts")+
    ggtitle("DMR overlap heterochromtin switch")
  
  DMR_summary$position <- 100
  DMR_summary$position[which(DMR_summary$condition=="Increase")] <- DMR_summary$percent[which(DMR_summary$condition == "Increase")]
  ggplot(DMR_summary, aes(x = tissue, y = percent, fill = condition)) +  
    geom_bar(stat = 'identity') +   
    theme_minimal() +   
    scale_fill_brewer(palette = "Pastel1") +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 20),legend.title = element_blank()) +
    geom_text(data = subset(DMR_summary),   
              aes(label = paste0(round(percent, 2),"%"), y = position),   
              color = "black", size = 5, vjust = 0.5) + 
    ylab("Percent (%)")+
    ggtitle("DMR overlap heterochromtin switch")
  }
