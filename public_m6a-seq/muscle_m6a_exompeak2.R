rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/","/usr/local/lib64/R/library"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(stringr)
library(dplyr)
library(ChIPseeker)
library(ggplot2)
library(org.Hs.eg.db)
library(AnnotationHub)
library(TxDb.Mfascicularis.UCSC.macFas5.ncbiRefSeq)
diff <- readRDS("data/public_data/m6a_monkey_CRA005942/muscle_exomePeak2_output/diffPeaks.rds")
diff <- as.data.frame(diff)
out<-read.csv("data/public_data/m6a_monkey_CRA005942/muscle_exomePeak2_output/diffPeaks.csv",row.names = 1)
out$Significant <- ifelse(out$fdr < 0.05 & abs(out$diff.log2FC) >= 0, 
                          ifelse(out$diff.log2FC > 0, "Up", "Down"), "Stable")
colour<- c("blue","red")
ggplot(
  # 数据、映射、颜色
  out, aes(x = diff.log2FC, y = -log10(fdr))) +
  geom_point(aes(color = Significant), size=2) +
  scale_color_manual(values = colour) +
  # scale_color_manual(values = c("blue","grey")) +
  # 注释
  # geom_text_repel(
  #   data = subset(out,`FDR.FA-noFA` < 0.05 & abs(out$logFC) >= 1),
  #   aes(label = Geneid),
  #   size = 5,max.overlaps = 100,
  #   box.padding = unit(0.35, "lines"),
  #   point.padding = unit(0.3, "lines")) +
  # 辅助线
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (p-value)") +
  # xlim(-3,3)+
  # 图例
  theme_bw()+
  theme(text = element_text(size = 20))
hub <- AnnotationHub()
peak_obj <- GRanges(seqnames = out$chr,   
                    ranges = IRanges(start = out$chromStart, end = out$chromEnd))
txdb <- TxDb.Mfascicularis.UCSC.macFas5.ncbiRefSeq
peak_anno <- annotatePeak(peak_obj, tssRegion=c(-3000, 3000),
                          TxDb=txdb,annoDb = 'org.Mmu.eg.db')
plotAnnoBar(peak_anno)
peak_anno<-as.data.frame(peak_anno)
