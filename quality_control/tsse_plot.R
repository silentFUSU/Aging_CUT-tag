rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(ggrepel)
# antibody <- "ATAC"
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
      tissue_label <- "IWAT"
    }
  }
  return(tissue_label)
} 
antibodys <- c("ATAC","H3K27ac","H3K4me3")
search_table_cut <- read.csv("data/samples/all/CUTTag_search_table.csv")
search_table_atac <- read.csv("data/samples/all/ATAC_search_table.csv")
search_table <- rbind(search_table_cut,search_table_atac)
for(i in c(1:length(antibodys))){
  antibody <- antibodys[i]
  tsse <- read.table(paste0("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/all/tsse/",antibody,"/tsse.txt"), header = FALSE)
  pattern <- ".*(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  tsse$V2 <- gsub(pattern, "\\1", tsse$V2)
  colnames(tsse)[2] <- "sample_name"
  tsse <- merge(search_table,tsse[,c(2:3)],by="sample_name")
 
  tsse$age <- factor(tsse$age, levels = c("3m","24m"))
  ggplot(tsse,aes(x=tissue,y=V3,color = age))+    
    geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.7)+
    geom_text(aes(label = mouse_ID), position = position_jitter(width = 0.2), vjust = -1, size = 3) +
    scale_fill_brewer(palette="Set3")+
    geom_hline(yintercept = 4, color = "red", linetype = "dashed") + 
    geom_hline(yintercept = 5, color = "blue", linetype = "dashed") + 
    annotate("text", x = Inf, y = 3.5, label = "tsse=4", color = "red", vjust = -0.5, hjust = 1.1, size = 5) +   
    annotate("text", x = Inf, y = 4.5, label = "tsse=5", color = "blue", vjust = -0.5, hjust = 1.1, size = 5) + 
    ggtitle(paste0(antibody," TSSe"))+ylim(0,20)+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")+labs(fill = "", color = "") +ylab("TSSe")
  
  ggsave(paste0("result/all/QC/tsse/",antibody,"_tsse.png"),width=15,height = 10,type="cairo")
  
}

tissues <- c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
                        "thymus","skin","bladder","bonemarrow","Hip","heart",
                        "muscle","jejunum","uterus","ovary","liver","tongue",
                        "cecum","colon","testis","stomach","ileum","pancreas")
search_table_cut <- read.csv("data/samples/all/CUTTag_search_table.csv")
search_table_atac <- read.csv("data/samples/all/ATAC_search_table.csv")
search_table <- rbind(search_table_cut,search_table_atac)
file_paths <- list(  
  ATAC = "/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/all/tsse/ATAC/tsse.txt",  
  H3K27ac = "/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/all/tsse/H3K27ac/tsse.txt",  
  H3K4me3 = "/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/all/tsse/H3K4me3/tsse.txt"  
)  
for(tissue in tissues){
  data_list <- list()  
  for (mark in names(file_paths)) {  
    data_list[[mark]] <- read.table(file_paths[[mark]], header = FALSE)[which(read.table(file_paths[[mark]], header = FALSE)$V1 == tissue),]  
  }  
  tissue_tsse <- do.call(rbind, data_list)  
  pattern <- ".*(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  tissue_tsse$V2 <- gsub(pattern, "\\1", tissue_tsse$V2)
  colnames(tissue_tsse)[2] <- "sample_name"
  tissue_tsse <- merge(tissue_tsse,search_table,by="sample_name")
  tissue_tsse$age <- factor(tissue_tsse$age,levels=c("3m","24m"))
  ggplot(tissue_tsse,aes(x=antibody,y=V3,color = age))+    
    geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.7)+
    geom_text(aes(label = mouse_ID), position = position_jitter(width = 0.2), vjust = -1, size = 3) +
    scale_fill_brewer(palette="Set3")+
    geom_hline(yintercept = 4, color = "red", linetype = "dashed") + 
    geom_hline(yintercept = 5, color = "blue", linetype = "dashed") + 
    annotate("text", x = Inf, y = 3.5, label = "tsse=4", color = "red", vjust = -0.5, hjust = 1.1, size = 5) +   
    annotate("text", x = Inf, y = 4.5, label = "tsse=5", color = "blue", vjust = -0.5, hjust = 1.1, size = 5) + 
    ggtitle(paste0(tissue_label_change(tissue)," TSSe"))+ylim(0,20)+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")+labs(fill = "", color = "") +ylab("TSSe")
  
  }

