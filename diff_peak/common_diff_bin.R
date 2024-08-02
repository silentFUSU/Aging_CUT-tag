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
library(ChIPseeker)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene

tissues = c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas","Hip","cecum","bonemarrow","ileum")
bin_size <- function(antibody){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    return ("10kb")
  }else{
    return("1kb")
  }  
}
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
  df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size(antibody),"_bins_diff_after_remove_batch_effect.csv")) 
  df <- df[which(df$Significant_bar!="Stable"),c("Geneid","Chr","Start","End","Significant_bar")]
  if(nrow(df) >0){
    df$tissue <- tissue
    df <- unique(df)
    peaks <- rbind(peaks,df)
  }
}
bin_file <- read.table(paste0("/storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_",bin_size(antibody),"_bins.bed"))
colnames(bin_file)[4]<-"Geneid"

increase <- peaks[which(peaks$Significant_bar=="Up"),]
increase_count <- increase %>%   
  count(Geneid)
increase_tissue <- increase %>%   
  group_by(Geneid) %>%   
  summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
increase_count <- merge(increase_count,increase_tissue,by="Geneid")
increase_count <- merge(increase_count,bin_file,by="Geneid")
peak_obj <- GRanges(seqnames = increase_count$V1,   
                    ranges = IRanges(start = increase_count$V2, end = increase_count$V3))
peak_anno <- annotatePeak(peak_obj, tssRegion=c(-3000, 3000),
                          TxDb=txdb, annoDb="org.Mm.eg.db")
peak_anno <- unique(as.data.frame(peak_anno))
increase_count$label <- paste0(increase_count$V1,"-",increase_count$V2,"-",increase_count$V3)
peak_anno$label <- paste0(peak_anno$seqnames,"-",peak_anno$start,"-",peak_anno$end)
increase_count <- merge(increase_count,peak_anno[,6:ncol(peak_anno)],by="label")

decrease <- peaks[which(peaks$Significant_bar=="Up"),]
decrease_count <- decrease %>%   
  count(Geneid)
decrease_tissue <- decrease %>%   
  group_by(Geneid) %>%   
  summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
decrease_count <- merge(decrease_count,decrease_tissue,by="Geneid")
decrease_count <- merge(decrease_count,bin_file,by="Geneid")
peak_obj <- GRanges(seqnames = decrease_count$V1,   
                    ranges = IRanges(start = decrease_count$V2, end = decrease_count$V3))
peak_anno <- annotatePeak(peak_obj, tssRegion=c(-3000, 3000),
                          TxDb=txdb, annoDb="org.Mm.eg.db")
peak_anno <- unique(as.data.frame(peak_anno))
decrease_count$label <- paste0(decrease_count$V1,"-",decrease_count$V2,"-",decrease_count$V3)
peak_anno$label <- paste0(peak_anno$seqnames,"-",peak_anno$start,"-",peak_anno$end)
decrease_count <- merge(decrease_count,peak_anno[,6:ncol(peak_anno)],by="label")



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
