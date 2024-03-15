rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(edgeR)
library(ggplot2)
library(ChIPseeker)
library(EnsDb.Mmusculus.v79)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
library(limma)
library(UpSetR)
antibody = "H3K27ac"
tissue = "colon"
window_size = "1000"
gap_size = "3000"
e_value = "100"
GO_database <- 'org.Mm.eg.db'
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
increase_upset <- list()
decrease_upset <- list()

tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle")
tissues <- c("liver","kidney","spleen","muscle")

genes<-data.frame(symbol = character(),  
                              Significant = character(),  
                              tissue = character(),
                              stringsAsFactors = FALSE)  

for (i in c(1:length(tissues))){
  tissue <- tissues[i]
  df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_merge-W",window_size,"-G",gap_size,"-E",e_value,"_diff_after_remove_batch_effect.csv")) 
  df <- df[which(df$Significant!="Stable"),c("symbol","Significant")]
  df$tissue <- tissue
  df <- unique(df)
  genes <- rbind(genes,df)
}
increase <- as.data.frame(table(genes[which(genes$Significant=="Up"),c("symbol")]))
decrease <- as.data.frame(table(genes[which(genes$Significant=="Down"),c("symbol")]))


for (i in c(1:length(tissues))){
  tissue <- tissues[i]
  df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_merge-W",window_size,"-G",gap_size,"-E",e_value,"_diff_after_remove_batch_effect.csv")) 
  df <- df[which(df$Significant!="Stable"),c("symbol","Significant")]
  df <- unique(df)
  increase_upset[[i]] <- df$symbol[which(df$Significant=="Up")]
  names(increase_upset)[i] <- tissue
  decrease_upset[[i]] <- df$symbol[which(df$Significant=="Down")]
  names(decrease_upset)[i] <- tissue
}

##################################
genes<-data.frame(symbol = character(),  
                  Significant = character(),  
                  tissue = character(),
                  stringsAsFactors = FALSE)  
for (i in c(1:length(tissues))){
  tissue <- tissues[i]
  df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_macs_narrowpeak_diff_after_remove_batch_effect.csv")) 
  df <- df[which(df$Significant!="Stable"),c("symbol","Significant")]
  df$tissue <- tissue
  df <- unique(df)
  genes <- rbind(genes,df)
}
increase <- as.data.frame(table(genes[which(genes$Significant=="Up"),c("symbol")]))
decrease <- as.data.frame(table(genes[which(genes$Significant=="Down"),c("symbol")]))

increase_upset <- list()
decrease_upset <- list()
for (i in c(1:length(tissues))){
  tissue <- tissues[i]
  df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_macs_narrowpeak_diff_after_remove_batch_effect.csv")) 
  df <- df[which(df$Significant!="Stable"),c("symbol","Significant")]
  df <- unique(df)
  increase_upset[[i]] <- df$symbol[which(df$Significant=="Up")]
  names(increase_upset)[i] <- tissue
  decrease_upset[[i]] <- df$symbol[which(df$Significant=="Down")]
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

genelist_up <-bitr(increase$Var1[which(increase$Freq==length(tissues))],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_up_GO <-enrichGO( genelist_up$ENTREZID,#GO富集分析
                           OrgDb = GO_database,
                           keyType = "ENTREZID",#设定读取的gene ID类型
                           ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                           pvalueCutoff = 0.05,#设定p值阈值
                           qvalueCutoff = 0.05,#设定q值阈值
                           readable = T)
barplot(genelist_up_GO,font.size=15,label_format = 100)

genelist_down <-bitr(decrease$Var1[which(decrease$Freq==length(tissues))],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
# genelist_down <-bitr(diff$symbol[which(diff$FDR.old.young<0.05 & diff$LogFC.old.young < -0.5)],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_down_GO <-enrichGO( genelist_down$ENTREZID,#GO富集分析
                             OrgDb = GO_database,
                             keyType = "ENTREZID",#设定读取的gene ID类型
                             ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                             pvalueCutoff = 0.05,#设定p值阈值
                             qvalueCutoff = 0.05,#设定q值阈值
                             readable = T)
barplot(genelist_down_GO,font.size=15,label_format = 100)
