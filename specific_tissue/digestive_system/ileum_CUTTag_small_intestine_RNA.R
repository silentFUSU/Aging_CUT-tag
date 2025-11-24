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
library(gg.gap)
library(scales)
library(colorspace)  
library(UpSetR)
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K4me1","H3K4me3","H3K27ac")

bin_size <- function(antibody){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    return ("10kb")
  }else{
    return("1kb")
  }  
}

GO_database <- 'org.Mm.eg.db'
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene

increase <- list()
decrease <- list()
for (i in c(1:length(antibodys))){
  antibody <- antibodys[i]
  diff<-read.csv(paste0("data/samples/ileum/",antibody,"/",antibody,"_",bin_size(antibody),"_bins_diff_after_remove_batch_effect.csv"))
  diff_up <- diff[which(diff$Significant_bar=="Up"),]
  if(nrow(diff_up) > 0 ){
    peak_obj <- GRanges(seqnames = diff_up$Chr,   
                        ranges = IRanges(start = diff_up$Start, end = diff_up$End))
    peak_anno <- annotatePeak(peak_obj, tssRegion=c(-3000, 3000),
                              TxDb=txdb, annoDb="org.Mm.eg.db")
    peak_anno <- unique(as.data.frame(peak_anno))
    gene_list_up <- unique(peak_anno$SYMBOL)
    increase[[i]] <- gene_list_up
    names(increase)[i] <- antibody
  }
  diff_down <- diff[which(diff$Significant_bar=="Down"),]
  if(nrow(diff_down) > 0){
    peak_obj <- GRanges(seqnames = diff_down$Chr,   
                        ranges = IRanges(start = diff_down$Start, end = diff_down$End))
    peak_anno <- annotatePeak(peak_obj, tssRegion=c(-3000, 3000),
                              TxDb=txdb, annoDb="org.Mm.eg.db")
    peak_anno <- unique(as.data.frame(peak_anno))
    gene_list_down <- unique(peak_anno$SYMBOL)
    decrease[[i]] <- gene_list_down
    names(decrease)[i] <- antibody
  }
}
RNA <- read.csv("data/public_data/GSE132040/Small_Intestine.csv")
increase[[i+1]] <- RNA$gene[which(RNA$FDR.m24.m3<0.05 & RNA$logFC.m24.m3>0)]
names(increase)[[i+1]] <- "RNA"
decrease[[i+1]] <- RNA$gene[which(RNA$FDR.m24.m3<0.05 & RNA$logFC.m24.m3<0)]
names(decrease)[[i+1]] <- "RNA"

increase <- Filter(Negate(is.null), increase)
decrease <- Filter(Negate(is.null), decrease)
dir.create(paste0("result/ileum/relationship_RNA/"))
for(i in c(1:(length(increase)-1))){
  pdf(file=paste0("result/ileum/relationship_RNA/",names(increase)[i],"_increase.pdf"),onefile = F,height = 7,width = 10)
  print(upset(fromList(setNames(list(increase[[i]], increase[[6]]), c(names(increase)[i], "RNA"))),  # fromList一个函数，用于将列表转换为与UpSetR兼容的数据形式。
        nsets = 100,     # 绘制的最大集合个数
        nintersects = NA, #绘制的最大交集个数，NA则全部绘制
        order.by = "degree", # 矩阵中的交点是如何排列的。 "freq"根据交集个数排序，"degree"根据
        keep.order = T, # 保持设置与使用sets参数输入的顺序一致。默认值是FALSE，它根据集合的大小排序。
        mb.ratio = c(0.6,0.4),   # 左侧和上方条形图的比例关系
        text.scale = 2 # 文字标签的大小
  ))
  dev.off()
}
for(i in c(1:(length(decrease)-1))){
  pdf(file=paste0("result/ileum/relationship_RNA/",names(increase)[i],"_decrease.pdf"),onefile = F,height = 7,width = 10)
  print(upset(fromList(setNames(list(decrease[[i]], decrease[[6]]), c(names(decrease)[i], "RNA"))),  # fromList一个函数，用于将列表转换为与UpSetR兼容的数据形式。
              nsets = 100,     # 绘制的最大集合个数
              nintersects = NA, #绘制的最大交集个数，NA则全部绘制
              order.by = "degree", # 矩阵中的交点是如何排列的。 "freq"根据交集个数排序，"degree"根据
              keep.order = T, # 保持设置与使用sets参数输入的顺序一致。默认值是FALSE，它根据集合的大小排序。
              mb.ratio = c(0.6,0.4),   # 左侧和上方条形图的比例关系
              text.scale = 2 # 文字标签的大小
  ))
  dev.off()
}


