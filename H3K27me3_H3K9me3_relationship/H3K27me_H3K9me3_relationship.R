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
library(dplyr)
library(reshape2)

tissues <- c("spleen","testis","colon","kidney","lung","liver","muscle","Hip","cecum","bonemarrow","brain")
tissues <- c("bonemarrow")
bin_size="10kb"
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  # peak_preprocess_bin_level(tissue,"H3K27me3","10kb")
  # peak_preprocess_bin_level(tissue,"H3K9me3","10kb")
  H3K27me3<-read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
  H3K9me3<-read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
  logFC<-merge(H3K27me3[,c("Geneid","LogFC.old.young","FDR.old.young")],H3K9me3[,c("Geneid","LogFC.old.young","FDR.old.young")],by="Geneid")
  colnames(logFC)[2:5]<-c("logFC_H3K27me3","FDR_H3K27me3","logFC_H3K9me3","FDR_H3K9me3")
  logFC<-merge(logFC,H3K27me3[,c("Geneid","Chr","Start","End")],by="Geneid")
  # logFC <- logFC[which(logFC$FDR_H3K27me3<0.05 & logFC$FDR_H3K9me3<0.05),]
  # cor(logFC$logFC_H3K27me3,logFC$logFC_H3K9me3,method = c("pearson"),use = "complete.obs")
  logFC_sig <- logFC[which(logFC$FDR_H3K27me3<0.05 & logFC$FDR_H3K9me3<0.05),]
  ggplot() +  
    geom_point(data=logFC, mapping=aes(logFC_H3K9me3, logFC_H3K27me3),color = "grey",alpha=0.5) +  
    geom_point(data=logFC_sig[which(logFC_sig$logFC_H3K27me3>0 & logFC_sig$logFC_H3K9me3<0),], mapping=aes(logFC_H3K9me3, logFC_H3K27me3),color = "#f6416c")+
    geom_point(data=logFC_sig[which(logFC_sig$logFC_H3K27me3>0 & logFC_sig$logFC_H3K9me3>0),], mapping=aes(logFC_H3K9me3, logFC_H3K27me3),color = "#00b8a9")+
    geom_point(data=logFC_sig[which(logFC_sig$logFC_H3K27me3<0 & logFC_sig$logFC_H3K9me3<0),], mapping=aes(logFC_H3K9me3, logFC_H3K27me3),color = "#ffde7d")+
    geom_point(data=logFC_sig[which(logFC_sig$logFC_H3K27me3<0 & logFC_sig$logFC_H3K9me3>0),], mapping=aes(logFC_H3K9me3, logFC_H3K27me3),color = "#48466d")+
    geom_hline(yintercept = 0, color = "red") +  
    geom_vline(xintercept = 0, color = "red") +
    coord_cartesian(xlim = c(-2, 2), ylim = c(-10, 10))+
    labs(x="log2(old/young) H3K9me3",
         y="log2(old/young) H3K27me3") +
    theme_minimal() + theme(text = element_text(size = 20)) +
    annotate("text",label = paste0(nrow(logFC_sig[which(logFC_sig$logFC_H3K27me3>0 & logFC_sig$logFC_H3K9me3>0),])),x=1, y=10,colour="#00b8a9",size=5)+
    annotate("text",label = paste0(nrow(logFC_sig[which(logFC_sig$logFC_H3K27me3<0 & logFC_sig$logFC_H3K9me3<0),])),x=-1, y=-10,colour="#ff9a00",size=5)+
    annotate("text",label = paste0(nrow(logFC_sig[which(logFC_sig$logFC_H3K27me3>0 & logFC_sig$logFC_H3K9me3<0),])),x=-1, y=10,colour="#f6416c",size=5)+
    annotate("text",label = paste0(nrow(logFC_sig[which(logFC_sig$logFC_H3K27me3<0 & logFC_sig$logFC_H3K9me3>0),])),x=1, y=-10,colour="#48466d",size=5)
  dir.create(paste0("result/",tissue,"/H3K27me3_H3K9me3_relationship/"))
  ggsave(paste0("result/",tissue,"/H3K27me3_H3K9me3_relationship/all_intersect_",bin_size,"bins.png"),width = 8,height = 8)
  dir.create(paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/")) 
  write.table(logFC_sig[which(logFC_sig$logFC_H3K27me3>0 & logFC_sig$logFC_H3K9me3>0),c("Chr","Start","End","Geneid")], 
              file=paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/",bin_size,"_all_significant_first_quadrant.bed"), 
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(logFC_sig[which(logFC_sig$logFC_H3K27me3>0 & logFC_sig$logFC_H3K9me3<0),c("Chr","Start","End","Geneid")], 
              file=paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/",bin_size,"_all_significant_second_quadrant.bed"), 
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(logFC_sig[which(logFC_sig$logFC_H3K27me3<0 & logFC_sig$logFC_H3K9me3<0),c("Chr","Start","End","Geneid")], 
              file=paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/",bin_size,"_all_significant_third_quadrant.bed"), 
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(logFC_sig[which(logFC_sig$logFC_H3K27me3<0 & logFC_sig$logFC_H3K9me3>0),c("Chr","Start","End","Geneid")], 
              file=paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/",bin_size,"_all_significant_fourth_quadrant.bed"), 
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}
peaks_list <-list(second =list(),fourth=list())
GO_database <- 'org.Mm.eg.db'
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
quadrants <- c("second","fourth")
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(quadrants))){
    quadrant <- quadrants[j]
    
    peak <- read.table(paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/",bin_size,"_all_significant_",quadrant,"_quadrant.bed"))
    peak_obj <- GRanges(seqnames = peak$V1,   
                        ranges = IRanges(start = peak$V2, end = peak$V3))
    peak_anno <- annotatePeak(peak_obj, tssRegion=c(-3000, 3000),
                              TxDb=txdb, annoDb="org.Mm.eg.db")
    peaks_list[[j]][[i]]<-peak_anno
    names(peaks_list[[j]])[i] <- tissue
  }  
}

plotAnnoBar(peaks_list[[1]])
plotAnnoBar(peaks_list[[2]])
# args <- commandArgs(trailingOnly=TRUE)
# arg_list <- list()
# for (arg in args) {
#   split_arg <- strsplit(arg, "=")[[1]]
#   arg_name <- split_arg[1]
#   arg_value <- split_arg[2]
#   arg_list[[arg_name]] <- arg_value
# }
# tissue <- arg_list[["tissue"]]
# window_size = "1000"
# gap_size = "3000"
# e_value = "100"
# txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
# common_peak_diff_in_another_modification<- function(tissue,antibody,window_size,gap_size,e_value){
#   if(antibody=="H3K27me3"){
#     another_antibody <- "H3K9me3"
#   }else{
#     another_antibody <- "H3K27me3"
#   }
#   
#   increase_peak_of_another <- readPeakFile(paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/",another_antibody,"-W1000-G3000-E100_up_in_",antibody,"_unique_sort.bed"))
#   decrease_peak_of_another <- readPeakFile(paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/",another_antibody,"-W1000-G3000-E100_down_in_",antibody,"_unique_sort.bed"))
#   increase_peakAnno_of_another <- annotatePeak(increase_peak_of_another, tssRegion=c(-3000, 3000),
#                                                TxDb=txdb, annoDb="org.Mm.eg.db")
#   decrease_peakAnno_of_another <- annotatePeak(decrease_peak_of_another, tssRegion=c(-3000, 3000),
#                                                TxDb=txdb, annoDb="org.Mm.eg.db")
#   
#   increase_peakAnno_of_another <- unique(as.data.frame(increase_peakAnno_of_another))
#   colnames(increase_peakAnno_of_another)[6]<-"Geneid"
#   decrease_peakAnno_of_another <- unique(as.data.frame(decrease_peakAnno_of_another))
#   colnames(decrease_peakAnno_of_another)[6]<-"Geneid" 
#   
#   df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_merge-W",window_size,"-G",gap_size,"-E",e_value,"_diff_after_remove_batch_effect.csv")) 
#   df$relationship_between_H3K27me3_H3K9me3 <-"Stable"
#   df$relationship_between_H3K27me3_H3K9me3[which(df$Geneid %in% increase_peakAnno_of_another$Geneid)] <-paste0("increase in ",another_antibody)
#   df$relationship_between_H3K27me3_H3K9me3[which(df$Geneid %in% decrease_peakAnno_of_another$Geneid)] <-paste0("decrease in ",another_antibody)
#   df$relationship_between_H3K27me3_H3K9me3[which(df$Significant=="Stable")]<-"Stable"
#   df$relationship_between_H3K27me3_H3K9me3 <- factor(df$relationship_between_H3K27me3_H3K9me3,levels = c(paste0("decrease in ",another_antibody),"Stable",paste0("increase in ",another_antibody)))
#   df_increase_peakAnno_of_another <- df[which(df$relationship_between_H3K27me3_H3K9me3==paste0("increase in ",another_antibody)),]
#   df_decrease_peakAnno_of_another <- df[which(df$relationship_between_H3K27me3_H3K9me3==paste0("decrease in ",another_antibody)),]
#   colour<- list(c("grey"),c("grey","red"),c("blue","grey","red"))
#   ggplot() +
#     geom_point(data=df, mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)), size=2, color="grey") +
#     geom_point(data=df_increase_peakAnno_of_another, mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)), size=2, color="red") +
#     geom_point(data=df_decrease_peakAnno_of_another, mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)), size=2, color="blue") +
#     # scale_color_manual(values = c("blue","grey")) +
#     # 注释
#     # geom_text_repel(
#     #   data = subset(out,`FDR.FA-noFA` < 0.05 & abs(out$logFC) >= 1),
#     #   aes(label = Geneid),
#     #   size = 5,max.overlaps = 100,
#     #   box.padding = unit(0.35, "lines"),
#     #   point.padding = unit(0.3, "lines")) +
#     # 辅助线
#     geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
#     # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
#     geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
#     # 坐标轴
#     labs(x="log2(fold change)",
#          y="-log10 (p-value)") +
#     # xlim(-3,3)+
#     # 图例
#     theme_bw()+labs(colour = "")+
#     theme(text = element_text(size = 20))
#   ggsave(paste0("result/",tissue,"/H3K27me3_H3K9me3_relationship/",another_antibody,"_diff_peak_in_",antibody,"_volcano.png"),width=10,height = 10)
#   write.csv(df,paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_merge-W",window_size,"-G",gap_size,"-E",e_value,"_diff_after_remove_batch_effect_relationship_with_",another_antibody,".csv"))
# }
# 
# common_peak_diff_in_another_modification(tissue,"H3K9me3",window_size,gap_size,e_value)
# common_peak_diff_in_another_modification(tissue,"H3K27me3",window_size,gap_size,e_value)

