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
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene

tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas")
bin_size <- "1kb"
antibody <- "H3K27ac"
peaks<-data.frame(Geneid = character(),  
                  Chr = character(),
                  Start = numeric(),
                  End = numeric(),
                  Significant = character(),  
                  tissue = character(),
                  stringsAsFactors = FALSE)  
for (i in c(1:length(tissues))){
  tissue <- tissues[i]
  df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect.csv")) 
  df <- df[which(df$Significant_bar!="Stable"),c("Geneid","Chr","Start","End","Significant_bar")]
  if(nrow(df) >0){
    df$tissue <- tissue
    df <- unique(df)
    peaks <- rbind(peaks,df)
  }
}
increase <- as.data.frame(table(peaks[which(peaks$Significant_bar=="Up"),c("Geneid")]))
decrease <- as.data.frame(table(peaks[which(peaks$Significant_bar=="Down"),c("Geneid")]))
bin_file <- read.table(paste0("/storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_",bin_size,"_bins.bed"))
colnames(bin_file)[4]<-"Geneid"
colnames(increase)[1] <- "Geneid"
increase <- merge(increase,bin_file,by="Geneid")
colnames(decrease)[1] <- "Geneid"
decrease <- merge(decrease,bin_file,by="Geneid")
write.csv(increase,paste0("data/samples/all/",antibody,"/common_increase_",bin_size,"_bins.csv"))
write.csv(decrease,paste0("data/samples/all/",antibody,"/common_decrease_",bin_size,"_bins.csv"))
increase_upset <- list()
decrease_upset <- list()

for (i in c(1:length(tissues))){
  tissue <- tissues[i]
  df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect.csv")) 
  df <- df[which(df$Significant_bar!="Stable"),c("Geneid","Chr","Start","End","Significant_bar")]
  if(nrow(df) >0){
    df$tissue <- tissue
    df <- unique(df)
  }
  increase_upset[[i]] <- df$Geneid[which(df$Significant_bar=="Up")]
  names(increase_upset)[i] <- tissue
  decrease_upset[[i]] <- df$Geneid[which(df$Significant_bar=="Down")]
  names(decrease_upset)[i] <- tissue
}
upset(fromList(increase_upset),  # fromList一个函数，用于将列表转换为与UpSetR兼容的数据形式。
      sets=tissues,nsets = 100,     # 绘制的最大集合个数
      nintersects = 20, #绘制的最大交集个数，NA则全部绘制
      order.by = "degree", # 矩阵中的交点是如何排列的。 "freq"根据交集个数排序，"degree"根据
      keep.order = T, # 保持设置与使用sets参数输入的顺序一致。默认值是FALSE，它根据集合的大小排序。
      mb.ratio = c(0.6,0.4),   # 左侧和上方条形图的比例关系
      text.scale = 2 # 文字标签的大小
)
upset(fromList(decrease_upset),  # fromList一个函数，用于将列表转换为与UpSetR兼容的数据形式。
      sets=tissues,nsets = 100,     # 绘制的最大集合个数
      nintersects = 20, #绘制的最大交集个数，NA则全部绘制
      order.by = "degree", # 矩阵中的交点是如何排列的。 "freq"根据交集个数排序，"degree"根据
      keep.order = T, # 保持设置与使用sets参数输入的顺序一致。默认值是FALSE，它根据集合的大小排序。
      mb.ratio = c(0.6,0.4),   # 左侧和上方条形图的比例关系
      text.scale = 2 # 文字标签的大小
)
