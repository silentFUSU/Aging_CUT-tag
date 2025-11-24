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
library(tidyr)
library(patchwork)
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
tissues = c("brain","Hip","liver","testis","kidney",
            "lung","bonemarrow","muscle","thymus","heart",
            "spleen","colon","cecum","ileum","pancreas")
antibody <- "H3K9me3"
p_list <- list()
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  out <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_10kb_bins_diff_after_remove_batch_effect.csv"))
  colour<- list(c("grey"),c("grey","red"),c("blue","grey","red"))
  if(tissue=="brain") tissue<-"FC"
  p_list[[i]]<-ggplot(
    # 数据、映射、颜色
    out, aes(x = `LogFC.old.young`, y = -log10(`FDR.old.young`))) +
    geom_point(aes(color = Significant), size=2,show.legend = FALSE) +
    scale_color_manual(values = colour[[nrow(as.data.frame(table(out$Significant)))]]) +
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
    annotate("text", x = min(out$LogFC.old.young), y = max(-log10(out$FDR.old.young)), label = nrow(out[which(out$Significant_bar=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$LogFC.old.young), y = max(-log10(out$FDR.old.young)), label = nrow(out[which(out$Significant_bar=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
}

combined_plot <- plot_a_list(p_list, 3, 5)
ggsave(paste0("result/all/diff/",antibody,"/all_tissue_bins_valcono_plot.png"),combined_plot,height = 30,width = 30,type="cairo")
