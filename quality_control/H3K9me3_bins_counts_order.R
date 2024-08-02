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
library(ggsci)
library(RColorBrewer) 
library(patchwork)
tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas","Hip","cecum","bonemarrow","ileum","heart","thymus","stomach","skin","aorta","tongue","bladder")
antibody <- "H3K9me3"
diff_peak_number <- read.csv("data/samples/all/diff_bins_number.csv")
diff_peak_number <- diff_peak_number[which(diff_peak_number$antibody == antibody),]
diff_peak_number$tissue[which(diff_peak_number$tissue=="brain")] <- "brain_FC"
diff_peak_number$tissue[which(diff_peak_number$tissue=="Hip")] <- "brain_Hip"
tissues[which(tissues=="brain")] <- "brain_FC"
tissues[which(tissues=="Hip")] <- "brain_Hip"

color <- pal_d3("category20") 
color(length(tissues))
colors <- setNames(color(length(tissues)),tissues) 


decrease <- diff_peak_number[which(diff_peak_number$Var1=="Down"),]
decrease <-  arrange(decrease, Freq)  
decrease$tissue <- factor(decrease$tissue,levels=decrease$tissue)
decrease$tissue_legend <- factor(decrease$tissue,levels=rev(decrease$tissue))
decrease_p <- ggplot(decrease,mapping = aes(x=Freq,y=tissue,fill = tissue_legend))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+ylab("")+xlab("Bins")+ggtitle(paste0("Decrease"))+
  theme(text = element_text(size = 18))+ scale_fill_manual(values = setNames(color(length(tissues)),tissues))+ xlim(0,(max(decrease$Freq+1000))) +guides(fill= guide_legend(title = ""))

increase <- diff_peak_number[which(diff_peak_number$Var1=="Up"),]
increase <-  arrange(increase, Freq)  
increase$tissue <- factor(increase$tissue,levels=increase$tissue)
increase$tissue_legend <- factor(increase$tissue,levels=rev(increase$tissue))
increase_p<-ggplot(increase,mapping = aes(x=Freq,y=tissue,fill = tissue_legend))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+ylab("")+xlab("Bins")+ggtitle(paste0("Increase"))+
  theme(text = element_text(size = 18))+ scale_fill_manual(values = setNames(color(length(tissues)),tissues)) +xlim(0,(max(decrease$Freq+1000))) +guides(fill= guide_legend(title = ""))

combined_plot <- increase_p + decrease_p 


tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas","Hip","cecum","bonemarrow","ileum","heart","thymus","stomach","skin","aorta","tongue","bladder")
antibody <- "H3K9me3"
diff_peak_number <- read.csv("data/samples/all/diff_bins_overlap_peaks_number.csv")
diff_peak_number <- diff_peak_number[which(diff_peak_number$antibody == antibody),]
diff_peak_number$tissue[which(diff_peak_number$tissue=="brain")] <- "brain_FC"
diff_peak_number$tissue[which(diff_peak_number$tissue=="Hip")] <- "brain_Hip"
tissues[which(tissues=="brain")] <- "brain_FC"
tissues[which(tissues=="Hip")] <- "brain_Hip"

color <- pal_d3("category20") 
color(length(tissues))
colors <- setNames(color(length(tissues)),tissues) 


decrease <- diff_peak_number[which(diff_peak_number$Var1=="Down"),]
decrease <-  arrange(decrease, Freq)  
decrease$tissue <- factor(decrease$tissue,levels=decrease$tissue)
decrease$tissue_legend <- factor(decrease$tissue,levels=rev(decrease$tissue))
decrease_p <- ggplot(decrease,mapping = aes(x=Freq,y=tissue,fill = tissue_legend))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+ylab("")+xlab("Bins overlap with peaks")+ggtitle(paste0("Decrease"))+
  theme(text = element_text(size = 18))+ scale_fill_manual(values = setNames(color(length(tissues)),tissues))  + xlim(0,(max(decrease$Freq+1000)))+guides(fill= guide_legend(title = ""))

increase <- diff_peak_number[which(diff_peak_number$Var1=="Up"),]
increase <-  arrange(increase, Freq)  
increase$tissue <- factor(increase$tissue,levels=increase$tissue)
increase$tissue_legend <- factor(increase$tissue,levels=rev(increase$tissue))
increase_p<-ggplot(increase,mapping = aes(x=Freq,y=tissue,fill = tissue_legend))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+ylab("")+xlab("Bins overlap with peaks")+ggtitle(paste0("Increase"))+
  theme(text = element_text(size = 18))+ scale_fill_manual(values = setNames(color(length(tissues)),tissues))+ xlim(0,(max(decrease$Freq+1000)))+guides(fill= guide_legend(title = ""))
combined_plot <- increase_p + decrease_p 

