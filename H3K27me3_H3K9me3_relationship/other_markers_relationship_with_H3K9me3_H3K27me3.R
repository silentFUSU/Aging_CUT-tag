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
library(patchwork)
bin_size <- "1kb"
tissues <- c("brain","Hip","lung","colon","spleen",
             "thymus","muscle","bonemarrow","heart","kidney",
             "skin","liver","stomach","testis","cecum")
antibodys <- c("H3K27ac","H3K4me3","H3K4me1")
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,axis_titles = "collect")
}

for(antibody in antibodys){
  p_list <- list()
  j=1
  for(tissue in tissues){
    out <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_1kb_bins_diff_after_remove_batch_effect.csv"))
    colour<- list(c("grey"),c("grey","red"),c("blue","grey","red"))
    quadrants<-c("first","second","third","fourth")
    quadrants_peaks <-list()
    for (i in c(1:length(quadrants))){
      quadrant <- quadrants[i]
      if(file.info(paste0("data/samples/",tissue,"/",antibody,"/bed/H3K27me3_H3K9me3_all_significant_",bin_size,"_",quadrant,"_quadrant_unique_sort.bed"))$size >0){
        peaks<- read.delim(paste0("data/samples/",tissue,"/",antibody,"/bed/H3K27me3_H3K9me3_all_significant_",bin_size,"_",quadrant,"_quadrant_unique_sort.bed"),header = F)
        quadrants_peaks[[i]] <- peaks$V4
        names(quadrants_peaks)[i] <- quadrant
      }
    }
    if(!is.null(quadrants_peaks[[1]])){
      p1<-ggplot() +
        geom_point(out, mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)), color = "grey",alpha=0.2) +
        geom_point(data=out[which(out$Geneid %in%quadrants_peaks[["first"]]),], mapping=aes(x = LogFC.old.young,y = -log10(FDR.old.young)),color = "#00b8a9",alpha=0.7)+
        geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
        # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
        geom_hline(yintercept = -log10(0.05),lty=4,col="red",lwd=0.8) +
        # 坐标轴
        labs(x="log2(fold change)",
             y="-log10 (FDR)") +
        # xlim(-3,3)+
        # 图例
        theme_bw()+
        theme(text = element_text(size = 20))
    }else{
      p1<-ggplot() +
        geom_point(out, mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)), color = "grey",alpha=0.2) +
        geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
        # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
        geom_hline(yintercept = -log10(0.05),lty=4,col="red",lwd=0.8) +
        # 坐标轴
        labs(x="log2(fold change)",
             y="-log10 (FDR)") +
        # xlim(-3,3)+
        # 图例
        theme_bw()+
        theme(text = element_text(size = 20))
    }
    if(!is.null(quadrants_peaks[[2]])){
      p2<-ggplot() +
        geom_point(out, mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)), color = "grey",alpha=0.2) +
        geom_point(data=out[which(out$Geneid %in%quadrants_peaks[["second"]]),], mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)),color ="#f6416c",alpha=0.7)+
        geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
        # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
        geom_hline(yintercept = -log10(0.05),lty=4,col="red",lwd=0.8) +
        # 坐标轴
        labs(x="log2(fold change)",
             y="-log10 (FDR)") +
        # xlim(-3,3)+
        # 图例
        theme_bw()+
        theme(text = element_text(size = 20))
    }else{
      p2<-ggplot() +
        geom_point(out, mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)), color = "grey",alpha=0.2) +
        geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
        # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
        geom_hline(yintercept = -log10(0.05),lty=4,col="red",lwd=0.8) +
        # 坐标轴
        labs(x="log2(fold change)",
             y="-log10 (FDR)") +
        # xlim(-3,3)+
        # 图例
        theme_bw()+
        theme(text = element_text(size = 20))
    }
    if(!is.null(quadrants_peaks[[3]])){
      p3<-ggplot() +
        geom_point(out, mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)), color = "grey",alpha=0.2) +
        geom_point(data=out[which(out$Geneid %in%quadrants_peaks[["third"]]),], mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)),color ="#ffde7d",alpha=0.7)+
        geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
        # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
        geom_hline(yintercept = -log10(0.05),lty=4,col="red",lwd=0.8) +
        # 坐标轴
        labs(x="log2(fold change)",
             y="-log10 (FDR)") +
        # xlim(-3,3)+
        # 图例
        theme_bw()+
        theme(text = element_text(size = 20))
    }else{
      p3<-ggplot() +
        geom_point(out, mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)), color = "grey",alpha=0.2) +
        geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
        # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
        geom_hline(yintercept = -log10(0.05),lty=4,col="red",lwd=0.8) +
        # 坐标轴
        labs(x="log2(fold change)",
             y="-log10 (FDR)") +
        # xlim(-3,3)+
        # 图例
        theme_bw()+
        theme(text = element_text(size = 20))
    }
    if(!is.null(quadrants_peaks[[3]])){
      p4<-ggplot() +
        geom_point(out, mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)), color = "grey",alpha=0.2) +
        geom_point(data=out[which(out$Geneid %in%quadrants_peaks[["fourth"]]),], mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)),color ="#48466d",alpha=0.7)+
        geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
        # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
        geom_hline(yintercept = -log10(0.05),lty=4,col="red",lwd=0.8) +
        # 坐标轴
        labs(x="log2(fold change)",
             y="-log10 (FDR)") +
        # xlim(-3,3)+
        # 图例
        theme_bw()+
        theme(text = element_text(size = 20))
    }else{
      p4<-ggplot() +
        geom_point(out, mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)), color = "grey",alpha=0.2) +
        geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
        # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
        geom_hline(yintercept = -log10(0.05),lty=4,col="red",lwd=0.8) +
        # 坐标轴
        labs(x="log2(fold change)",
             y="-log10 (FDR)") +
        # xlim(-3,3)+
        # 图例s
        theme_bw()+
        theme(text = element_text(size = 20))
    }
    if(tissue == "brain") tissue <- "FC"
    p2 <- p2+ggtitle(tissue)
    p_list[[j]] <- p2
    names(p_list)[j] <- tissue
    j=j+1
  #   p <- (p2+p1)/(p3+p4)
  #   dir.create(paste0("result/",tissue,"/",antibody))
  #   ggsave(paste0("result/",tissue,"/",antibody,"/diff_",antibody,"_with_H3K27me3_and_H3K9me3_relationship.png"),p,width = 15,height=15,type="cairo")
  }
   combined_plot <- plot_a_list(p_list,3,5)
   ggsave(paste0("result/all/diff/",antibody,"/",antibody,"_relationship_skin_stomach_with_H3K27me3_H3K9me3.png"),combined_plot,width = 15,height = 10,type="cairo")
}

