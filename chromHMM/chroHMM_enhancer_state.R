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
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas","Hip","cecum","bonemarrow","ileum","heart","thymus","stomach","skin") 

state  <- "E15"
target_state <- "E11"
plist<-list()
antibody <- "H3K4me1"
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  young1 <- read.delim(paste0("result/all/ChromHMM/15_until_skin/split_1k/",tissue,"_young1_15_segments_1k.bed"),header = F)
  young2 <- read.delim(paste0("result/all/ChromHMM/15_until_skin/split_1k/",tissue,"_young2_15_segments_1k.bed"),header = F)
  young1_state <- young1[which(young1$V4==state),]
  young2_state <- young2[which(young2$V4==state),]
  old1 <- read.delim(paste0("result/all/ChromHMM/15_until_skin/split_1k/",tissue,"_old1_15_segments_1k.bed"),header = F)
  old2 <- read.delim(paste0("result/all/ChromHMM/15_until_skin/split_1k/",tissue,"_old2_15_segments_1k.bed"),header = F)
  old1_state <- old1[which(old1$V4==target_state),]
  old2_state <- old2[which(old2$V4==target_state),]
  young_state <- intersect(young1_state,young2_state)
  old_state <- intersect(old1_state,old2_state)
  young_state$label <- paste0(young_state$V1,"-",young_state$V2,"-",young_state$V3)
  old_state$label <- paste0(old_state$V1,"-",old_state$V2,"-",old_state$V3)
  
  interest_state <- young_state[which(young_state$label %in% old_state$label),c(1:3,5)]
  df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_1kb_bins_diff_after_remove_batch_effect.csv"))
  df_state <- df[,c(1:4)]
  df_state$label <- paste0(df_state$Chr,"-",df_state$Start,"-",df_state$End)
  df_state_young_enhancer <- df_state[which(df_state$label %in% interest_state$label),]
  
  plist[[i]] <- ggplot()+
    # 数据、映射、颜色
    geom_point(data=df, mapping=aes(`LogFC.old.young`, -log10(`FDR.old.young`)),color = "grey",alpha=0.5) +
    geom_point(data=df[which(df$Geneid %in% df_state_young_enhancer$Geneid),], mapping=aes(`LogFC.old.young`, -log10(`FDR.old.young`)),color = "red",alpha=0.5) +
    # geom_point(data=H3K27me3[which(H3K27me3$H3K9me3_sig=="H3K9me3_decrease" & H3K27me3$FDR.old.young <0.05 & H3K27me3$LogFC.old.young<0),], mapping=aes(`LogFC.old.young`, -log10(`FDR.old.young`)),color = "blue",alpha=0.7)+
    # geom_point(data=H3K27me3[which(H3K27me3$H3K9me3_sig=="H3K9me3_increase" & H3K27me3$FDR.old.young <0.05 & H3K27me3$LogFC.old.young<0),], mapping=aes(`LogFC.old.young`, -log10(`FDR.old.young`)),color = "red",alpha=0.7)+
    # geom_point(data=H3K27me3[which(H3K27me3$H3K9me3_sig=="H3K9me3_increase" & H3K27me3$FDR.old.young <0.05 & H3K27me3$LogFC.old.young>0),], mapping=aes(`LogFC.old.young`, -log10(`FDR.old.young`)),color = "red",alpha=0.7)+
    # geom_point(data=H3K27me3[which(H3K27me3$H3K9me3_sig=="H3K9me3_decrease" & H3K27me3$FDR.old.young <0.05 & H3K27me3$LogFC.old.young>0),], mapping=aes(`LogFC.old.young`, -log10(`FDR.old.young`)),color = "blue",alpha=0.7)+
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
    ggtitle(tissue)+
    # 坐标轴
    labs(x="log2(fold change)",
         y="-log10 (p-value)") +
    # xlim(-3,3)+
    # 图例
    theme_bw()+
    theme(text = element_text(size = 30))+
    annotate("text", x = min(df$LogFC.old.young), y = max(-log10(df$FDR.old.young)), label = nrow(df[which(df$Geneid %in% df_state_young_enhancer$Geneid & df$LogFC.old.young<0),]), vjust = 5, hjust = 0,colour="blue",size=8)+
    annotate("text", x = max(df$LogFC.old.young), y = max(-log10(df$FDR.old.young)), label = nrow(df[which(df$Geneid %in% df_state_young_enhancer$Geneid & df$LogFC.old.young>0),]), vjust = 5, hjust = 1.5,colour="red",size=8)
  
}

combined_plot <- plot_a_list(plist, 4, 5)
ggsave(paste0("result/all/ChromHMM/15_until_skin/state_transfer/",antibody,"_",state,"_to_",target_state,".png"),combined_plot,width = 32.5,height = 26,type="cairo")

for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  young1 <- read.delim(paste0("result/all/ChromHMM/15_until_skin/split_1k/",tissue,"_young1_15_segments_1k.bed"),header = F)
  young2 <- read.delim(paste0("result/all/ChromHMM/15_until_skin/split_1k/",tissue,"_young2_15_segments_1k.bed"),header = F)
  young1_state <- young1[which(young1$V4==state),]
  young2_state <- young2[which(young2$V4==state),]
  old1 <- read.delim(paste0("result/all/ChromHMM/15_until_skin/split_1k/",tissue,"_old1_15_segments_1k.bed"),header = F)
  old2 <- read.delim(paste0("result/all/ChromHMM/15_until_skin/split_1k/",tissue,"_old2_15_segments_1k.bed"),header = F)
  # old1_state <- old1[which(old1$V4==state),]
  # old2_state <- old2[which(old2$V4==state),]
  young_state <- intersect(young1_state,young2_state)
  old1_state <- old1[which(paste0(old1$V1,"-",old1$V2,"-",old1$V3) %in% paste0(young_state$V1,"-",young_state$V2,"-",young_state$V3)),]
  old2_state <- old2[which(paste0(old2$V1,"-",old2$V2,"-",old2$V3) %in% paste0(young_state$V1,"-",young_state$V2,"-",young_state$V3)),]
  old_state <- intersect(old1_state,old2_state)
  keep_state <- old_state[which(old_state$V4 != "E11"),]
  df <- read.csv(paste0("data/samples/",tissue,"/H3K27ac/H3K27ac_1kb_bins_diff_after_remove_batch_effect.csv"))
  df_state <- df[,c(1:4)]
  df_state$label <- paste0(df_state$Chr,"-",df_state$Start,"-",df_state$End)
  keep_state$label <- paste0(keep_state$V1,"-",keep_state$V2,"-",keep_state$V3)
  df_state_young_enhancer <- df_state[which(df_state$label %in% keep_state$label),]
  
  plist[[i]] <- ggplot()+
    # 数据、映射、颜色
    geom_point(data=df, mapping=aes(`LogFC.old.young`, -log10(`FDR.old.young`)),color = "grey",alpha=0.5) +
    geom_point(data=df[which(df$Geneid %in% df_state_young_enhancer$Geneid),], mapping=aes(`LogFC.old.young`, -log10(`FDR.old.young`)),color = "red",alpha=0.5) +
    # geom_point(data=H3K27me3[which(H3K27me3$H3K9me3_sig=="H3K9me3_decrease" & H3K27me3$FDR.old.young <0.05 & H3K27me3$LogFC.old.young<0),], mapping=aes(`LogFC.old.young`, -log10(`FDR.old.young`)),color = "blue",alpha=0.7)+
    # geom_point(data=H3K27me3[which(H3K27me3$H3K9me3_sig=="H3K9me3_increase" & H3K27me3$FDR.old.young <0.05 & H3K27me3$LogFC.old.young<0),], mapping=aes(`LogFC.old.young`, -log10(`FDR.old.young`)),color = "red",alpha=0.7)+
    # geom_point(data=H3K27me3[which(H3K27me3$H3K9me3_sig=="H3K9me3_increase" & H3K27me3$FDR.old.young <0.05 & H3K27me3$LogFC.old.young>0),], mapping=aes(`LogFC.old.young`, -log10(`FDR.old.young`)),color = "red",alpha=0.7)+
    # geom_point(data=H3K27me3[which(H3K27me3$H3K9me3_sig=="H3K9me3_decrease" & H3K27me3$FDR.old.young <0.05 & H3K27me3$LogFC.old.young>0),], mapping=aes(`LogFC.old.young`, -log10(`FDR.old.young`)),color = "blue",alpha=0.7)+
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
    ggtitle(tissue)+
    # 坐标轴
    labs(x="log2(fold change)",
         y="-log10 (p-value)") +
    # xlim(-3,3)+
    # 图例
    theme_bw()+
    theme(text = element_text(size = 30))
}
combined_plot <- plot_a_list(plist, 4, 5)
ggsave(paste0("result/all/ChromHMM/15_until_skin/state_transfer/H3K27ac_transfer_",state,".png"),combined_plot,width = 22.5,height = 18,type="cairo")
