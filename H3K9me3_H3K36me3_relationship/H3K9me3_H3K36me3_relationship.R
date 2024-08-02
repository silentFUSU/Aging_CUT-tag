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
tissues <- c("bladder","thymus","muscle","liver")
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
bin_size="10kb"
plist<-list()
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  # peak_preprocess_bin_level(tissue,"H3K36me3","10kb")
  # peak_preprocess_bin_level(tissue,"H3K9me3","10kb")
  H3K36me3<-read.csv(paste0("data/samples/",tissue,"/H3K36me3/H3K36me3_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
  H3K9me3<-read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
  logFC<-merge(H3K36me3[,c("Geneid","LogFC.old.young","FDR.old.young")],H3K9me3[,c("Geneid","LogFC.old.young","FDR.old.young")],by="Geneid")
  colnames(logFC)[2:5]<-c("logFC_H3K36me3","FDR_H3K36me3","logFC_H3K9me3","FDR_H3K9me3")
  logFC<-merge(logFC,H3K36me3[,c("Geneid","Chr","Start","End")],by="Geneid")
  # logFC <- logFC[which(logFC$FDR_H3K36me3<0.05 & logFC$FDR_H3K9me3<0.05),]
  # cor(logFC$logFC_H3K36me3,logFC$logFC_H3K9me3,method = c("pearson"),use = "complete.obs")
  logFC_sig <- logFC[which(logFC$FDR_H3K36me3<0.05 & logFC$FDR_H3K9me3<0.05),]
  if(tissue=="brain") tissue<-"FC"
  plist[[i]] <- ggplot() +  
    geom_point(data=logFC, mapping=aes(logFC_H3K9me3, logFC_H3K36me3),color = "grey",alpha=0.5) +  
    geom_point(data=logFC_sig[which(logFC_sig$logFC_H3K36me3>0 & logFC_sig$logFC_H3K9me3<0),], mapping=aes(logFC_H3K9me3, logFC_H3K36me3),color = "#f6416c")+
    geom_point(data=logFC_sig[which(logFC_sig$logFC_H3K36me3>0 & logFC_sig$logFC_H3K9me3>0),], mapping=aes(logFC_H3K9me3, logFC_H3K36me3),color = "#00b8a9")+
    geom_point(data=logFC_sig[which(logFC_sig$logFC_H3K36me3<0 & logFC_sig$logFC_H3K9me3<0),], mapping=aes(logFC_H3K9me3, logFC_H3K36me3),color = "#ffde7d")+
    geom_point(data=logFC_sig[which(logFC_sig$logFC_H3K36me3<0 & logFC_sig$logFC_H3K9me3>0),], mapping=aes(logFC_H3K9me3, logFC_H3K36me3),color = "#48466d")+
    geom_hline(yintercept = 0, color = "red") +  
    geom_vline(xintercept = 0, color = "red") +
    ggtitle(tissue)+
    coord_cartesian(xlim = c(-2, 2), ylim = c(-10, 10))+
    labs(x="log2(old/young) H3K9me3",
         y="log2(old/young) H3K36me3") +
    theme_minimal() + theme(text = element_text(size = 20)) +
    annotate("text",label = paste0(nrow(logFC_sig[which(logFC_sig$logFC_H3K36me3>0 & logFC_sig$logFC_H3K9me3>0),])),x=1, y=10,colour="#00b8a9",size=5)+
    annotate("text",label = paste0(nrow(logFC_sig[which(logFC_sig$logFC_H3K36me3<0 & logFC_sig$logFC_H3K9me3<0),])),x=-1, y=-10,colour="#ff9a00",size=5)+
    annotate("text",label = paste0(nrow(logFC_sig[which(logFC_sig$logFC_H3K36me3>0 & logFC_sig$logFC_H3K9me3<0),])),x=-1, y=10,colour="#f6416c",size=5)+
    annotate("text",label = paste0(nrow(logFC_sig[which(logFC_sig$logFC_H3K36me3<0 & logFC_sig$logFC_H3K9me3>0),])),x=1, y=-10,colour="#48466d",size=5)
    dir.create(paste0("result/",tissue,"/H3K36me3_H3K9me3_relationship/"))
    ggsave(paste0("result/",tissue,"/H3K36me3_H3K9me3_relationship/all_intersect_",bin_size,"bins.png"),width = 8,height = 8,type="cairo")
    dir.create(paste0("data/samples/",tissue,"/H3K36me3_H3K9me3_intersect/"))
    write.table(logFC_sig[which(logFC_sig$logFC_H3K36me3>0 & logFC_sig$logFC_H3K9me3>0),c("Chr","Start","End","Geneid")],
                file=paste0("data/samples/",tissue,"/H3K36me3_H3K9me3_intersect/",bin_size,"_all_significant_first_quadrant.bed"),
                sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(logFC_sig[which(logFC_sig$logFC_H3K36me3>0 & logFC_sig$logFC_H3K9me3<0),c("Chr","Start","End","Geneid")],
                file=paste0("data/samples/",tissue,"/H3K36me3_H3K9me3_intersect/",bin_size,"_all_significant_second_quadrant.bed"),
                sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(logFC_sig[which(logFC_sig$logFC_H3K36me3<0 & logFC_sig$logFC_H3K9me3<0),c("Chr","Start","End","Geneid")],
                file=paste0("data/samples/",tissue,"/H3K36me3_H3K9me3_intersect/",bin_size,"_all_significant_third_quadrant.bed"),
                sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(logFC_sig[which(logFC_sig$logFC_H3K36me3<0 & logFC_sig$logFC_H3K9me3>0),c("Chr","Start","End","Geneid")],
                file=paste0("data/samples/",tissue,"/H3K36me3_H3K9me3_intersect/",bin_size,"_all_significant_fourth_quadrant.bed"),
                sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}
combined_plot <- plot_a_list(plist, 2, 2)
# dir.create("result/all/H3K36me3_H3K9me3")
ggsave("result/all/H3K36me3_H3K9me3/all_tissue_H3K36me3_H3K9me3.png",combined_plot,width = 10,height = 10,type="cairo")

bin_size <- "1kb"
antibodys <- c("H3K27ac","H3K4me3","H3K4me1")
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
      if(file.info(paste0("data/samples/",tissue,"/",antibody,"/bed/H3K36me3_H3K9me3_all_significant_",bin_size,"_",quadrant,"_quadrant_unique_sort.bed"))$size >0){
        peaks<- read.delim(paste0("data/samples/",tissue,"/",antibody,"/bed/H3K36me3_H3K9me3_all_significant_",bin_size,"_",quadrant,"_quadrant_unique_sort.bed"),header = F)
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
        theme(text = element_text(size = 20))+
        annotate("text", x = min(out$LogFC.old.young), y = max(-log10(out$FDR.old.young)), label = nrow(out[which(out$Geneid %in%quadrants_peaks[["third"]] & out$LogFC.old.young<0),]), vjust = 5, hjust = 0,colour="blue",size=8)+
        annotate("text", x = max(out$LogFC.old.young), y = max(-log10(out$FDR.old.young)), label = nrow(out[which(out$Geneid %in%quadrants_peaks[["third"]] & out$LogFC.old.young>0),]), vjust = 5, hjust = 1.5,colour="red",size=8)
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
    p3 <- p3+ggtitle(tissue)
    p_list[[j]] <- p3
    names(p_list)[j] <- tissue
    j=j+1
    #   p <- (p2+p1)/(p3+p4)
    #   dir.create(paste0("result/",tissue,"/",antibody))
    #   ggsave(paste0("result/",tissue,"/",antibody,"/diff_",antibody,"_with_H3K36me3_and_H3K9me3_relationship.png"),p,width = 15,height=15,type="cairo")
  }
  combined_plot <- plot_a_list(p_list,2,2)
  ggsave(paste0("result/all/diff/",antibody,"/",antibody,"_relationship_with_H3K36me3_H3K9me3.png"),combined_plot,width = 10,height = 10,type="cairo")
}

tissues <- c("bladder","thymus","muscle","liver")
peaks_list <-list()
GO_database <- 'org.Mm.eg.db'
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
quadrant <- c("third")
bin_size <- "10kb"
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  peak <- read.table(paste0("data/samples/",tissue,"/H3K36me3_H3K9me3_intersect/",bin_size,"_all_significant_",quadrant,"_quadrant.bed"))
  peak_obj <- GRanges(seqnames = peak$V1,   
                      ranges = IRanges(start = peak$V2, end = peak$V3))
  peak_anno <- annotatePeak(peak_obj, tssRegion=c(-3000, 3000),
                            TxDb=txdb, annoDb="org.Mm.eg.db")
  peaks_list[[i]]<-peak_anno
  names(peaks_list)[i] <- tissue
  peak_anno <- as.data.frame(peak_anno)
  promoter <- peak_anno[which(str_detect(peak_anno$annotation,"Promoter")),]
  genebody <- peak_anno[which(str_detect(peak_anno$annotation,"Exon") | str_detect(peak_anno$annotation,"Intron")),]
  write.table(promoter[,1:3],
              file=paste0("data/samples/",tissue,"/H3K36me3_H3K9me3_intersect/",bin_size,"_all_significant_third_quadrant_promoter.bed"),
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(genebody[,1:3],
              file=paste0("data/samples/",tissue,"/H3K36me3_H3K9me3_intersect/",bin_size,"_all_significant_third_quadrant_genebody.bed"),
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}
peaks_list[[1]]@annoStat$Feature <- as.vector(peaks_list[[1]]@annoStat$Feature)
peaks_list[[1]]@annoStat$Feature <- factor(peaks_list[[1]]@annoStat$Feature,levels = c("Promoter (<=1kb)","Promoter (1-2kb)","Promoter (2-3kb)","5' UTR","3' UTR"  ,"1st Exon" ,"Other Exon","1st Intron"  , "Other Intron" ,"Distal Intergenic","Downstream (<=300)" ))
plotAnnoBar(peaks_list)

bin_size <- "1kb"
antibodys <- c("H3K27ac","H3K4me3","H3K4me1")
position <- "genebody"
for(antibody in antibodys){
  p_list <- list()
  j=1
  for(tissue in tissues){
    out <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_1kb_bins_diff_after_remove_batch_effect.csv"))
    colour<- list(c("grey"),c("grey","red"),c("blue","grey","red"))
    quadrant <- "third"
    if(file.info(paste0("data/samples/",tissue,"/",antibody,"/bed/H3K36me3_H3K9me3_all_significant_",bin_size,"_",quadrant,"_quadrant_",position,"_unique_sort.bed"))$size >0){
        peaks<- read.delim(paste0("data/samples/",tissue,"/",antibody,"/bed/H3K36me3_H3K9me3_all_significant_",bin_size,"_",quadrant,"_quadrant_",position,"_unique_sort.bed"),header = F)
        quadrants_peaks <- peaks$V4
      }
    if(!is.null(quadrants_peaks[[3]])){
      p3<-ggplot() +
        geom_point(out, mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)), color = "grey",alpha=0.2) +
        geom_point(data=out[which(out$Geneid %in%quadrants_peaks),], mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)),color ="#ffde7d",alpha=0.7)+
        geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
        # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
        geom_hline(yintercept = -log10(0.05),lty=4,col="red",lwd=0.8) +
        # 坐标轴
        labs(x="log2(fold change)",
             y="-log10 (FDR)") +
        # xlim(-3,3)+
        # 图例
        theme_bw()+
        theme(text = element_text(size = 20))+
        annotate("text", x = min(out$LogFC.old.young), y = max(-log10(out$FDR.old.young)), label = nrow(out[which(out$Geneid %in%quadrants_peaks & out$LogFC.old.young<0),]), vjust = 5, hjust = 0,colour="blue",size=8)+
        annotate("text", x = max(out$LogFC.old.young), y = max(-log10(out$FDR.old.young)), label = nrow(out[which(out$Geneid %in%quadrants_peaks & out$LogFC.old.young>0),]), vjust = 5, hjust = 1.5,colour="red",size=8)
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

    if(tissue == "brain") tissue <- "FC"
    p3 <- p3+ggtitle(tissue)
    p_list[[j]] <- p3
    names(p_list)[j] <- tissue
    j=j+1
    #   p <- (p2+p1)/(p3+p4)
    #   dir.create(paste0("result/",tissue,"/",antibody))
    #   ggsave(paste0("result/",tissue,"/",antibody,"/diff_",antibody,"_with_H3K36me3_and_H3K9me3_relationship.png"),p,width = 15,height=15,type="cairo")
  }
  combined_plot <- plot_a_list(p_list,2,2)
  ggsave(paste0("result/all/diff/",antibody,"/",antibody,"_relationship_with_H3K36me3_H3K9me3_in_",position,".png"),combined_plot,width = 10,height = 10,type="cairo")
}

antibody1="H3K4me1"
antibody2="H3K27ac"
position <- "genebody"
p_list <- list()
j=1
for(tissue in tissues){
  out <- read.csv(paste0("data/samples/",tissue,"/",antibody1,"/",antibody1,"_1kb_bins_diff_after_remove_batch_effect.csv"))
  out2 <- read.csv(paste0("data/samples/",tissue,"/",antibody2,"/",antibody2,"_1kb_bins_diff_after_remove_batch_effect.csv"))
  colour<- list(c("grey"),c("grey","red"),c("blue","grey","red"))
  quadrant <- "third"
  if(file.info(paste0("data/samples/",tissue,"/",antibody,"/bed/H3K36me3_H3K9me3_all_significant_",bin_size,"_",quadrant,"_quadrant_",position,"_unique_sort.bed"))$size >0){
    peaks<- read.delim(paste0("data/samples/",tissue,"/",antibody,"/bed/H3K36me3_H3K9me3_all_significant_",bin_size,"_",quadrant,"_quadrant_",position,"_unique_sort.bed"),header = F)
    quadrants_peaks <- peaks$V4
  }
  quadrants_peaks <- intersect(quadrants_peaks,out2$Geneid[which(out2$Significant_bar=="Up")])
  if(!is.null(quadrants_peaks)){
    p3<-ggplot() +
      geom_point(out, mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)), color = "grey",alpha=0.2) +
      geom_point(data=out[which(out$Geneid %in%quadrants_peaks),], mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)),color ="#ffde7d",alpha=0.7)+
      geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
      # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
      geom_hline(yintercept = -log10(0.05),lty=4,col="red",lwd=0.8) +
      # 坐标轴
      labs(x="log2(fold change)",
           y="-log10 (FDR)") +
      # xlim(-3,3)+
      # 图例
      theme_bw()+
      theme(text = element_text(size = 20))+
      annotate("text", x = min(out$LogFC.old.young), y = max(-log10(out$FDR.old.young)), label = nrow(out[which(out$Geneid %in%quadrants_peaks & out$LogFC.old.young<0),]), vjust = 5, hjust = 0,colour="blue",size=8)+
      annotate("text", x = max(out$LogFC.old.young), y = max(-log10(out$FDR.old.young)), label = nrow(out[which(out$Geneid %in%quadrants_peaks & out$LogFC.old.young>0),]), vjust = 5, hjust = 1.5,colour="red",size=8)
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
  
  if(tissue == "brain") tissue <- "FC"
  p3 <- p3+ggtitle(tissue)
  p_list[[j]] <- p3
  names(p_list)[j] <- tissue
  j=j+1
  #   p <- (p2+p1)/(p3+p4)
  #   dir.create(paste0("result/",tissue,"/",antibody))
  #   ggsave(paste0("result/",tissue,"/",antibody,"/diff_",antibody,"_with_H3K36me3_and_H3K9me3_relationship.png"),p,width = 15,height=15,type="cairo")
}
combined_plot <- plot_a_list(p_list,2,2)



tissues <- c("muscle","liver","spleen","colon","cecum")
mm10_10kb <- read.table("~/ref_data/mm10_10kb_bins.bed")
GO_database <- 'org.Mm.eg.db'
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
annotate_list <- list()
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  H3K9me3 <-  read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  H3K36me3 <- read.csv(paste0("data/samples/",tissue,"/H3K36me3/H3K36me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  H3K9me3_increase_bar <- unique(H3K9me3$Geneid[which(H3K9me3$Significant_bar=="Up")])
  H3K9me3_decrease_bar <- unique(H3K9me3$Geneid[which(H3K9me3$Significant_bar=="Down")])
  H3K36me3_increase_bar <- unique(H3K36me3$Geneid[which(H3K36me3$Significant_bar=="Up")])
  H3K36me3_decrease_bar <- unique(H3K36me3$Geneid[which(H3K36me3$Significant_bar=="Down")])
  
  H3K9me3_H3K36me3_down <- intersect(H3K9me3_decrease_bar,H3K36me3_decrease_bar)
  H3K9me3_H3K36me3_up <- intersect(H3K9me3_increase_bar,H3K36me3_increase_bar)
  
  H3K9me3_H3K36me3_down <- mm10_10kb[which(mm10_10kb$V4 %in% H3K9me3_H3K36me3_down),]
  H3K9me3_H3K36me3_up <- mm10_10kb[which(mm10_10kb$V4 %in% H3K9me3_H3K36me3_up),]
  
  if(nrow(H3K9me3_H3K36me3_down)>0){
    H3K9me3_H3K36me3_down_grange <- GRanges(seqnames = H3K9me3_H3K36me3_down$V1,   
                                            ranges = IRanges(start = H3K9me3_H3K36me3_down$V2, end = H3K9me3_H3K36me3_down$V3))
    H3K9me3_H3K36me3_down_peak_anno <- annotatePeak(H3K9me3_H3K36me3_down_grange, tssRegion=c(-3000, 3000),
                                                    TxDb=txdb, annoDb="org.Mm.eg.db")
  }else{
    H3K9me3_H3K36me3_down_peak_anno <- NULL
  }

  if(nrow(H3K9me3_H3K36me3_up)>0){
  H3K9me3_H3K36me3_up_grange <- GRanges(seqnames = H3K9me3_H3K36me3_up$V1,   
                                         ranges = IRanges(start = H3K9me3_H3K36me3_up$V2, end = H3K9me3_H3K36me3_up$V3))
  H3K9me3_H3K36me3_up_peak_anno <- annotatePeak(H3K9me3_H3K36me3_up_grange, tssRegion=c(-3000, 3000),
                                                 TxDb=txdb, annoDb="org.Mm.eg.db")
  }else{
    H3K9me3_H3K36me3_up_peak_anno <- NULL
  }
  annotate_list[[i*2-1]] <- H3K9me3_H3K36me3_up_peak_anno
  annotate_list[[i*2]] <- H3K9me3_H3K36me3_down_peak_anno 
  names(annotate_list)[[i*2-1]] <- paste0(tissue,"_up")
  names(annotate_list)[[i*2]] <- paste0(tissue,"_down")

}
annotate_list <- Filter(Negate(is.null), annotate_list)
plotAnnoBar(annotate_list,title = "H3K9me3 and H3K36me3")
number_annotate_list <- vector()
for(i in c(1:length(annotate_list))){
  number_annotate_list[i] <- annotate_list[[i]]@peakNum
}




for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  H3K9me3 <-  read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  H3K36me3 <- read.csv(paste0("data/samples/",tissue,"/H3K36me3/H3K36me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  H3K9me3_increase_bar <- unique(H3K9me3$Geneid[which(H3K9me3$Significant_bar=="Up")])
  H3K9me3_decrease_bar <- unique(H3K9me3$Geneid[which(H3K9me3$Significant_bar=="Down")])
  H3K36me3_increase_bar <- unique(H3K36me3$Geneid[which(H3K36me3$Significant_bar=="Up")])
  H3K36me3_decrease_bar <- unique(H3K36me3$Geneid[which(H3K36me3$Significant_bar=="Down")])
  
  H3K9me3_H3K36me3_down <- intersect(H3K9me3_decrease_bar,H3K36me3_decrease_bar)
  H3K9me3_H3K36me3_up <- intersect(H3K9me3_increase_bar,H3K36me3_increase_bar)
  
  H3K9me3_H3K36me3_down <- mm10_10kb[which(mm10_10kb$V4 %in% H3K9me3_H3K36me3_down),]
  H3K9me3_H3K36me3_up <- mm10_10kb[which(mm10_10kb$V4 %in% H3K9me3_H3K36me3_up),]
  H3K9me3_H3K36me3_down_grange <- GRanges(seqnames = H3K9me3_H3K36me3_down$V1,   
                                          ranges = IRanges(start = H3K9me3_H3K36me3_down$V2, end = H3K9me3_H3K36me3_down$V3))
  H3K9me3_H3K36me3_down_peak_anno <- annotatePeak(H3K9me3_H3K36me3_down_grange, tssRegion=c(-3000, 3000),
                                                  TxDb=txdb, annoDb="org.Mm.eg.db")
  H3K9me3_H3K36me3_up_grange <- GRanges(seqnames = H3K9me3_H3K36me3_up$V1,   
                                        ranges = IRanges(start = H3K9me3_H3K36me3_up$V2, end = H3K9me3_H3K36me3_up$V3))
  H3K9me3_H3K36me3_up_peak_anno <- annotatePeak(H3K9me3_H3K36me3_up_grange, tssRegion=c(-3000, 3000),
                                                TxDb=txdb, annoDb="org.Mm.eg.db")
  
  H3K9me3_H3K36me3_up_peak_anno <- unique(as.data.frame(H3K9me3_H3K36me3_up_peak_anno))

  genelist_up <- bitr(H3K9me3_H3K36me3_up_peak_anno$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  genelist_up_GO <-enrichGO( genelist_up$ENTREZID,#GO富集分析
                               OrgDb = GO_database,
                               keyType = "ENTREZID",#设定读取的gene ID类型
                               ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                               pvalueCutoff = 0.05,#设定p值阈值
                               qvalueCutoff = 0.05,#设定q值阈值
                               readable = T)
  barplot(genelist_up_GO,font.size=15,label_format = 100)
  H3K9me3_H3K36me3_down_peak_anno <- unique(as.data.frame(H3K9me3_H3K36me3_down_peak_anno))
  
  genelist_down <- bitr(H3K9me3_H3K36me3_down_peak_anno$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  genelist_down_GO <-enrichGO( genelist_down$ENTREZID,#GO富集分析
                               OrgDb = GO_database,
                               keyType = "ENTREZID",#设定读取的gene ID类型
                               ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                               pvalueCutoff = 0.05,#设定p值阈值
                               qvalueCutoff = 0.05,#设定q值阈值
                               readable = T)
  barplot(genelist_down_GO,font.size=15,label_format = 100)
}