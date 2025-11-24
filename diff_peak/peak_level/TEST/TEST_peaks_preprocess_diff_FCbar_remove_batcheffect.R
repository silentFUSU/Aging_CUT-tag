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
antibody = "H3K27ac"
tissue = "brain"
window_size = "1000"
gap_size = "3000"
e_value = "100"
peak_macs_FCbar <- function(tissue,antibody){
  df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_macs_narrowpeak_diff_after_remove_batch_effect.csv")) 
  df$Significant_bar <- "Stable"
  df$Significant_bar[which(df$FDR.old.young < 0.05 & (df$old_1/df$young_1 > 1.2) & (df$old_2/df$young_2 > 1.2))] <- "Up"
  df$Significant_bar[which(df$FDR.old.young < 0.05 & (df$old_1/df$young_1 < 0.8) & (df$old_2/df$young_2 < 0.8))] <- "Down"
  colour<- list(c("grey"),c("grey","red"),c("blue","grey","red"))
  ggplot(
    # 数据、映射、颜色
    df, aes(x = LogFC.old.young, y = -log10(FDR.old.young))) +
    geom_point(aes(color = Significant_bar), size=2) +
    scale_color_manual(values = colour[[nrow(as.data.frame(table(df$Significant)))]]) +
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
  ggsave(paste0("result/",tissue,"/diffpeaks/",antibody,"_volcano_plot_FCbar_after_remove_batch_effect.png"),width = 10,height = 10)
  write.csv(df,paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_macs_narrowpeak_diff_FCbar_after_remove_batch_effect.csv"),row.names = F)
}

peak_sicer_FCbar <- function(tissue,antibody,window_size,gap_size,e_value){
  df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_merge-W",window_size,"-G",gap_size,"-E",e_value,"_diff_after_remove_batch_effect.csv")) 
  df$Significant_bar <- "Stable"
  df$Significant_bar[which(df$FDR.old.young < 0.05 & (df$old_1/df$young_1 > 1.2) & (df$old_2/df$young_2 > 1.2))] <- "Up"
  df$Significant_bar[which(df$FDR.old.young < 0.05 & (df$old_1/df$young_1 < 0.8) & (df$old_2/df$young_2 < 0.8))] <- "Down"
  colour<- list(c("grey"),c("grey","red"),c("blue","grey","red"))
  ggplot(
    # 数据、映射、颜色
    df, aes(x = LogFC.old.young, y = -log10(FDR.old.young))) +
    geom_point(aes(color = Significant_bar), size=2) +
    scale_color_manual(values = colour[[nrow(as.data.frame(table(df$Significant)))]]) +
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
  ggsave(paste0("result/",tissue,"/diffpeaks/",antibody,"_volcano_plot_FCbar_after_remove_batch_effect.png"),width = 10,height = 10)
  write.csv(df,paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_merge-W",window_size,"-G",gap_size,"-E",e_value,"_diff_FCbar_after_remove_batch_effect.csv"),row.names = F)
}

tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas")
antibodys <- c("H3K27me3","H3K9me3","H3K36me3")
for (j in c(1:length(tissues))){
  tissue <- tissues[j]
  for (i in c(1:length(antibodys))){
    antibody <- antibodys[i]
    peak_sicer_FCbar(tissue,antibody,window_size,gap_size,e_value)
  }
}
antibodys <- c("H3K27ac","H3K4me3","H3K4me1")
for (j in c(1:length(tissues))){
  tissue <- tissues[j]
  for (i in c(1:length(antibodys))){
    antibody <- antibodys[i]
    peak_macs_FCbar(tissue,antibody)
  }
}
