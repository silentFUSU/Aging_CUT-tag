rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/","/usr/local/lib64/R/library"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(patchwork)
library(stringr)
library(dplyr)
library(tidyr)
library(UpSetR)

tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas")

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}
bin_size <- function(antibody){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    return ("10kb")
  }else{
    return("1kb")
  }  
  
}
antibody <- "H3K27ac"
jaccard_index<-data.frame(antibody = character(),  
                  jaccard_index = numeric(),
                  tissue1 = character(),
                  tissue2 = character(),
                  condition = character(),
                  stringsAsFactors = FALSE)
antibodys <- c("H3K27ac","H3K4me1","H3K4me3","H3K36me3","H3K27me3","H3K9me3")
for (k in c(1:length(antibodys))){
  antibody <- antibodys[k]
  for (i in c(1:(length(tissues)-1))){
    tissue1 <- tissues[i]
    df <- read.csv(paste0("data/samples/",tissue1,"/",antibody,"/",antibody,"_",bin_size(antibody),"_bins_diff_after_remove_batch_effect.csv")) 
    tissue1_peak_increase <- df$Geneid[df$Significant_bar=="Up"]
    tissue1_peak_decrease <- df$Geneid[df$Significant_bar=="Down"]
    for (j in c((i+1):length(tissues))){
      tissue2 <- tissues[j]
      df <- read.csv(paste0("data/samples/",tissue2,"/",antibody,"/",antibody,"_",bin_size(antibody),"_bins_diff_after_remove_batch_effect.csv")) 
      tissue2_peak_increase <- df$Geneid[df$Significant_bar=="Up"]
      tissue2_peak_decrease <- df$Geneid[df$Significant_bar=="Down"]
      t_jaccard_index<-data.frame(antibody = antibody,  
                                  jaccard_index = jaccard(tissue1_peak_increase,tissue2_peak_increase),
                                  tissue1 = tissue1,
                                  tissue2 = tissue2,
                                  condition = "increase",
                                  stringsAsFactors = FALSE)  
      jaccard_index <- rbind(jaccard_index,t_jaccard_index)
      t_jaccard_index<-data.frame(antibody = antibody,  
                                  jaccard_index = jaccard(tissue1_peak_decrease,tissue2_peak_decrease),
                                  tissue1 = tissue1,
                                  tissue2 = tissue2,
                                  condition = "decrease",
                                  stringsAsFactors = FALSE)  
      jaccard_index <- rbind(jaccard_index,t_jaccard_index)
    }
  }
}

write.csv(jaccard_index,"result/all/jaccard_index.csv")
jaccard_index$condition[which(jaccard_index$condition=="increase")] <- "Increase"
jaccard_index$condition[which(jaccard_index$condition=="decrease")] <- "Decrease"


ggplot(jaccard_index,aes(x=antibody,y=jaccard_index,fill =condition))+
  geom_boxplot()+
  # geom_jitter(shape=16,size=1,position = position_jitter(0.2))+
  scale_fill_brewer(palette="Set3")+
  theme_bw()+theme(text = element_text(size = 18))+ylab("Jaccard Index")+
  xlab("")+labs(fill = "", color = "")

for (i in c(1:length(tissues))){
  tissue <- tissues[i]
  df <- jaccard_index[which(jaccard_index$tissue1==tissue | jaccard_index$tissue2==tissue),]
  df$tissue2[which(df$tissue2==tissue)] <- df$tissue1[which(df$tissue2==tissue)] 
  df$tissue1 <- tissue
  
  ggplot(df[which(df$condition=="Decrease"),],mapping = aes(x=antibody,y=jaccard_index,fill = tissue2))+
    geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+xlab("")+ylab("jaccard index")+ggtitle(paste0(tissue," Decrease"))+
    ggsci::scale_fill_d3()+theme(text = element_text(size = 18))
  ggsave(paste0("result/all/jaccard_index/",tissue,"_Decrease_jaccard_index.png"),width=15,height = 10)
  
  ggplot(df[which(df$condition=="Increase"),],mapping = aes(x=antibody,y=jaccard_index,fill = tissue2))+
    geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+xlab("")+ylab("jaccard index")+ggtitle(paste0(tissue," Increase"))+
    ggsci::scale_fill_d3()+theme(text = element_text(size = 18))
  ggsave(paste0("result/all/jaccard_index/",tissue,"_Increase_jaccard_index.png"),width=15,height = 10)
}
