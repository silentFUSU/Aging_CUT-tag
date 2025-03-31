rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
library(limma)
library(reshape2)

tissues <- c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
             "thymus","skin","bladder","bonemarrow","Hip","heart",
             "muscle","jejunum","uterus","ovary","liver","tongue",
             "cecum","colon","testis","stomach","pancreas","iWAT","ileum")

plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
tissue_label_change <- function(tissue){
  if(tissue=="brain"){
    tissue_label <- "Cortex"
  }else if(tissue == "Hip"){
    tissue_label <- "Hippocampus"
  }else if(tissue == "CB"){
    tissue_label <- "Cerebellum"
  }else{
    tissue_label <- str_to_title(tissue)
    if(tissue_label == "Bonemarrow"){
      tissue_label <- "Bone Marrow"
    }else if(tissue_label == "Bat"){
      tissue_label <- "BAT"
    }else if(tissue_label == "Mammarygland"){
      tissue_label <- "Mammary Gland"
    }else if(tissue_label == "Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
} 

bin_size="10kb"
plist<-list()
sort_table <- data.frame(tissue=as.character(),count=as.character())
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  H3K27me3<-read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
  H3K9me3<-read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
  logFC<-merge(H3K27me3[,c("Geneid","LogFC.old.young","FDR.old.young")],H3K9me3[,c("Geneid","LogFC.old.young","FDR.old.young")],by="Geneid")
  colnames(logFC)[2:5]<-c("logFC_H3K27me3","FDR_H3K27me3","logFC_H3K9me3","FDR_H3K9me3")
  logFC<-merge(logFC,H3K27me3[,c("Geneid","Chr","Start","End")],by="Geneid")
  logFC$both_sig <- "Stable"
  logFC$both_sig[which(logFC$FDR_H3K27me3<0.05 & logFC$FDR_H3K9me3<0.05)] <- "Both_Significant"
  logFC$quadrant <- "First"
  logFC$quadrant[which(logFC$logFC_H3K27me3>0 & logFC$logFC_H3K9me3<0)] <-"Second"
  logFC$quadrant[which(logFC$logFC_H3K27me3<0 & logFC$logFC_H3K9me3<0)] <-"Third"
  logFC$quadrant[which(logFC$logFC_H3K27me3<0 & logFC$logFC_H3K9me3>0)] <-"Fourth"
  # logFC <- logFC[which(logFC$Chr %in% paste0("chr",c(1:19,"X"))),]
  # H3K9me3_peaks <- read.table(paste0("data/samples/",tissue,"/H3K9me3/bed/H3K9me3_10kb_in_young_merge-W1000-G3000-E100.bed"))
  # logFC <- logFC[-which(logFC$Geneid %in% H3K9me3_peaks$V4),]
  # H3K27me3_peaks <- read.table(paste0("data/samples/",tissue,"/H3K27me3/bed/H3K27me3_10kb_in_old_merge-W1000-G3000-E100.bed"))
  # logFC <- logFC[which(logFC$Geneid %in% H3K27me3_peaks$V4),]
  # logFC <- logFC[which(logFC$Geneid %in% c(H3K27me3_peaks$V4,H3K9me3_peaks$V4)),]
  # dir.create(paste0("data/samples/all/diff_table/H3K27me3_H3K9me3/"))
  # write.csv(logFC,paste0("data/samples/all/diff_table/H3K27me3_H3K9me3/",tissue,"_H3K27me3_H3K9me3_10kb_bins_diff_after_remove_batch_effect.csv"),row.names = F)
  logFC_sig <- logFC[which(logFC$FDR_H3K27me3<0.05 & logFC$FDR_H3K9me3<0.05),]
  plist[[i]] <- ggplot() +  
    geom_point(data=logFC, mapping=aes(logFC_H3K9me3, logFC_H3K27me3),color = "grey",alpha=0.5) +  
    geom_point(data=logFC_sig[which(logFC_sig$logFC_H3K27me3>0 & logFC_sig$logFC_H3K9me3<0),], mapping=aes(logFC_H3K9me3, logFC_H3K27me3),color = "#f6416c")+
    geom_point(data=logFC_sig[which(logFC_sig$logFC_H3K27me3>0 & logFC_sig$logFC_H3K9me3>0),], mapping=aes(logFC_H3K9me3, logFC_H3K27me3),color = "#00b8a9")+
    geom_point(data=logFC_sig[which(logFC_sig$logFC_H3K27me3<0 & logFC_sig$logFC_H3K9me3<0),], mapping=aes(logFC_H3K9me3, logFC_H3K27me3),color = "#ffde7d")+
    geom_point(data=logFC_sig[which(logFC_sig$logFC_H3K27me3<0 & logFC_sig$logFC_H3K9me3>0),], mapping=aes(logFC_H3K9me3, logFC_H3K27me3),color = "#48466d")+
    geom_hline(yintercept = 0, color = "red") +  
    geom_vline(xintercept = 0, color = "red") +
    ggtitle(tissue_label_change(tissue))+
    coord_cartesian(xlim = c(-2, 2), ylim = c(-10, 10))+
    labs(x="log2(old/young) H3K9me3",
         y="log2(old/young) H3K27me3") +
    theme_minimal() + theme(text = element_text(size = 20)) +
    annotate("text",label = paste0(nrow(logFC_sig[which(logFC_sig$logFC_H3K27me3>0 & logFC_sig$logFC_H3K9me3>0),])),x=1, y=10,colour="#00b8a9",size=5)+
    annotate("text",label = paste0(nrow(logFC_sig[which(logFC_sig$logFC_H3K27me3<0 & logFC_sig$logFC_H3K9me3<0),])),x=-1, y=-10,colour="#ff9a00",size=5)+
    annotate("text",label = paste0(nrow(logFC_sig[which(logFC_sig$logFC_H3K27me3>0 & logFC_sig$logFC_H3K9me3<0),])),x=-1, y=10,colour="#f6416c",size=5)+
    annotate("text",label = paste0(nrow(logFC_sig[which(logFC_sig$logFC_H3K27me3<0 & logFC_sig$logFC_H3K9me3>0),])),x=1, y=-10,colour="#48466d",size=5)
    t_sort <- data.frame(tissue=tissue,count=nrow(logFC_sig[which(logFC_sig$logFC_H3K27me3>0 & logFC_sig$logFC_H3K9me3<0),]))
    sort_table <- rbind(sort_table,t_sort)
    names(plist)[i] <- tissue
    dir.create(paste0("result/",tissue,"/H3K27me3_H3K9me3_relationship/"))
    ggsave(paste0("result/",tissue,"/H3K27me3_H3K9me3_relationship/all_intersect_",bin_size,"bins_after_remove_batch_effect.png"),width = 8,height = 8,type="cairo")
    dir.create(paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/"))
    write.table(logFC_sig[which(logFC_sig$logFC_H3K27me3>0 & logFC_sig$logFC_H3K9me3>0),c("Chr","Start","End","Geneid")],
                file=paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/",bin_size,"_all_significant_first_quadrant_after_remove_batch_effect.bed"),
                sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(logFC_sig[which(logFC_sig$logFC_H3K27me3>0 & logFC_sig$logFC_H3K9me3<0),c("Chr","Start","End","Geneid")],
                file=paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/",bin_size,"_all_significant_second_quadrant_after_remove_batch_effect.bed"),
                sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(logFC_sig[which(logFC_sig$logFC_H3K27me3<0 & logFC_sig$logFC_H3K9me3<0),c("Chr","Start","End","Geneid")],
                file=paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/",bin_size,"_all_significant_third_quadrant_after_remove_batch_effect.bed"),
                sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(logFC_sig[which(logFC_sig$logFC_H3K27me3<0 & logFC_sig$logFC_H3K9me3>0),c("Chr","Start","End","Geneid")],
                file=paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/",bin_size,"_all_significant_fourth_quadrant_after_remove_batch_effect.bed"),
                sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}
sort_table <- sort_table[order(sort_table$count,decreasing = TRUE),]
plot_list <- plist[sort_table$tissue]  
combined_plot <- plot_a_list(plot_list, 4, 7)
# dir.create("result/all/H3K27me3_H3K9me3")
ggsave("result/all/H3K27me3_H3K9me3/all_tissue_H3K27me3_H3K9me3_after_remove_batch_effect.png",combined_plot,width = 35,height = 20,type="cairo")


