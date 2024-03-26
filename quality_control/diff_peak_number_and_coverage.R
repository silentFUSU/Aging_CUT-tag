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
tissues <-  c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas","Hip","cecum","bonemarrow")
antibodys <- c("H3K27me3","H3K9me3","H3K36me3")
diff_peak_number <-data.frame(Var1 = character(),  
                         Freq = numeric(),  
                         tissue = character(),  
                         antibody = character(),  
                         stringsAsFactors = FALSE)  
diff_peak_percent <-data.frame(Var1 = character(),  
                              Freq = numeric(),  
                              tissue = character(),  
                              antibody = character(),  
                              stringsAsFactors = FALSE)  
diff_peak_coverage <-data.frame(Var1 = character(),  
                               Freq = numeric(),  
                               tissue = character(),  
                               antibody = character(),  
                               stringsAsFactors = FALSE)  
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(antibodys))){
    antibody<-antibodys[j]
    diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_10kb_bins_diff_after_remove_batch_effect.csv"))
    sig<-as.data.frame(table(diff$Significant_bar))
    sig$tissue <- tissue
    sig$antibody <- antibody
    diff_peak_number<-rbind(diff_peak_number,sig)
    t_sum=sum(sig$Freq)
    sig$Freq=sig$Freq/t_sum
    diff_peak_percent <- rbind(diff_peak_percent,sig)
    sig$coverage<-0
    for (k in c(1:nrow(sig))){
      sig$coverage[k]<-sum(diff$Length[which(diff$Significant_bar == sig$Var1[k])])
    }
    diff_peak_coverage <- rbind(diff_peak_coverage,sig)
  }
}
antibodys <- c("H3K27ac","H3K4me3","H3K4me1")
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(antibodys))){
    antibody<-antibodys[j]
    diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_1kb_bins_diff_after_remove_batch_effect.csv"))
    sig<-as.data.frame(table(diff$Significant_bar))
    sig$tissue <- tissue
    sig$antibody <- antibody
    diff_peak_number<-rbind(diff_peak_number,sig)
    t_sum=sum(sig$Freq)
    sig$Freq=sig$Freq/t_sum
    diff_peak_percent <- rbind(diff_peak_percent,sig)
    sig$coverage<-0
    for (k in c(1:nrow(sig))){
      sig$coverage[k]<-sum(diff$Length[which(diff$Significant_bar == sig$Var1[k])])
    }
    diff_peak_coverage <- rbind(diff_peak_coverage,sig)
  }
}
write.csv(diff_peak_number,"data/samples/all/diff_peak_number.csv",row.names = F)
write.csv(diff_peak_percent,"data/samples/all/diff_peak_percent.csv",row.names = F)
write.csv(diff_peak_coverage,"data/samples/all/diff_peak_coverage.csv",row.names = F)
ggplot(diff_peak_number[which(diff_peak_number$Var1 !="Stable"),],aes(x=antibody,y=Freq,fill =Var1))+
  geom_boxplot()+
  # geom_jitter(shape=16,size=1,position = position_jitter(0.2))+
  scale_fill_brewer(palette="Set3")+
  theme_bw()+theme(text = element_text(size = 18))+ylab("Peaks")+
  xlab("")+labs(fill = "", color = "")



diff_peak_percent$Freq <- diff_peak_percent$Freq *100

ggplot(diff_peak_percent[which(diff_peak_percent$Var1 !="Stable"),],aes(x=antibody,y=Freq,fill =Var1))+
  geom_boxplot()+
  # geom_jitter(shape=16,size=1,position = position_jitter(0.2))+
  scale_fill_brewer(palette="Set3")+
  theme_bw()+theme(text = element_text(size = 18))+ylab("Percent (%)")+
  xlab("")+labs(fill = "", color = "")

ggplot(diff_peak_coverage[which(diff_peak_coverage$Var1 !="Stable" & diff_peak_coverage$antibody %in%c("H3K27me3","H3K9me3","H3K36me3")),],aes(x=antibody,y=coverage,fill =Var1))+
  geom_boxplot()+
  # geom_jitter(shape=16,size=1,position = position_jitter(0.2))+
  scale_fill_brewer(palette="Set3")+
  theme_bw()+theme(text = element_text(size = 18))+ylab("peak length (bp)")+
  xlab("")+labs(fill = "", color = "")

ggplot(diff_peak_coverage[which(diff_peak_coverage$Var1 !="Stable" & diff_peak_coverage$antibody %in%c("H3K27ac","H3K4me3","H3K4me1")),],aes(x=antibody,y=coverage,fill =Var1))+
  geom_boxplot()+
  # geom_jitter(shape=16,size=1,position = position_jitter(0.2))+
  scale_fill_brewer(palette="Set3")+
  theme_bw()+theme(text = element_text(size = 18))+ylab("peak length (bp)")+
  xlab("")+labs(fill = "", color = "")

custom_colors <- c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#800000", "#008000", "#000080", "#808000", "#800080", "#008080") 
diff_peak_coverage$tissue[which(diff_peak_coverage$tissue=="brain")] <- "brain_FC"
diff_peak_coverage$tissue[which(diff_peak_coverage$tissue=="Hip")] <- "brain_Hip"

diff_peak_number$tissue[which(diff_peak_number$tissue=="brain")] <- "brain_FC"
diff_peak_number$tissue[which(diff_peak_number$tissue=="Hip")] <- "brain_Hip"

diff_peak_percent$tissue[which(diff_peak_percent$tissue=="brain")] <- "brain_FC"
diff_peak_percent$tissue[which(diff_peak_percent$tissue=="Hip")] <- "brain_Hip"

ggplot(diff_peak_number[which(diff_peak_number$Var1 =="Down"),],mapping = aes(x=antibody,y=Freq,fill = tissue))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+xlab("")+ylab("Peaks")+ggtitle(paste0("Decrease"))+
  theme(text = element_text(size = 18))+scale_fill_manual(values =custom_colors)

ggplot(diff_peak_number[which(diff_peak_number$Var1 =="Up"),],mapping = aes(x=antibody,y=Freq,fill = tissue))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+xlab("")+ylab("Peaks")+ggtitle(paste0("Increase"))+
  theme(text = element_text(size = 18))+scale_fill_manual(values =custom_colors)



ggplot(diff_peak_percent[which(diff_peak_percent$Var1 =="Down"),],mapping = aes(x=antibody,y=Freq,fill = tissue))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+xlab("")+ylab("Percent")+ggtitle(paste0("Decrease"))+
  theme(text = element_text(size = 18))+scale_fill_manual(values =custom_colors)

ggplot(diff_peak_percent[which(diff_peak_percent$Var1 =="Up"),],mapping = aes(x=antibody,y=Freq,fill = tissue))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+xlab("")+ylab("Percent")+ggtitle(paste0("Increase"))+
  theme(text = element_text(size = 18))+scale_fill_manual(values =custom_colors)

ggplot(diff_peak_coverage[which(diff_peak_coverage$Var1 =="Down" & diff_peak_coverage$antibody %in% c("H3K27ac","H3K4me1","H3K4me3")),],mapping = aes(x=antibody,y=coverage,fill = tissue))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+xlab("")+ylab("Peak length")+ggtitle(paste0("Decrease"))+
  ggsci::scale_fill_d3()+theme(text = element_text(size = 18))
ggplot(diff_peak_coverage[which(diff_peak_coverage$Var1 =="Down" & diff_peak_coverage$antibody %in% c("H3K27me3","H3K9me3","H3K36me3")),],mapping = aes(x=antibody,y=coverage,fill = tissue))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+xlab("")+ylab("Peak length")+ggtitle(paste0("Decrease"))+
  ggsci::scale_fill_d3()+theme(text = element_text(size = 18))
# gg.gap(plot=p,
#        segments=c(7000000,10000000),
#        ylim=c(0,400000000),
#        tick_width = c(2000000,100000000),)

ggplot(diff_peak_coverage[which(diff_peak_coverage$Var1 =="Up" & diff_peak_coverage$antibody %in% c("H3K27ac","H3K4me1","H3K4me3")),],mapping = aes(x=antibody,y=coverage,fill = tissue))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+xlab("")+ylab("Percent of peaks(%)")+ggtitle(paste0("Increase"))+
  ggsci::scale_fill_d3()+theme(text = element_text(size = 18)) +scale_y_continuous(labels = scientific_format())

ggplot(diff_peak_coverage[which(diff_peak_coverage$Var1 =="Up" & diff_peak_coverage$antibody %in% c("H3K27me3","H3K9me3","H3K36me3")),],mapping = aes(x=antibody,y=coverage,fill = tissue))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+xlab("")+ylab("Peak length")+ggtitle(paste0("Increase"))+
  ggsci::scale_fill_d3()+theme(text = element_text(size = 18))

# 
# gg.gap(plot=p,
#        segments=c(10000000,50000000),
#        ylim=c(0,160000000),
#        tick_width = c(2000000,20000000))

peak_coverage <- data.frame(tissue = character(),  
                            antibody = character(),  
                            stringsAsFactors = FALSE,
                            coverage = numeric()) 
antibodys <- c("H3K27me3","H3K9me3","H3K36me3")
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(antibodys))){
    antibody<-antibodys[j]
    diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_merge-W1000-G3000-E100_diff_after_remove_batch_effect.csv"))
    df <- data.frame(tissue = tissue,  
                     antibody = antibody,  
                     stringsAsFactors = FALSE,
                     coverage = sum(diff$Length)) 
    peak_coverage <- rbind(peak_coverage,df)
  }
}
antibodys <- c("H3K27ac","H3K4me3","H3K4me1")
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(antibodys))){
    antibody<-antibodys[j]
    diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_macs_narrowpeak_diff_after_remove_batch_effect.csv"))
    df <- data.frame(tissue = tissue,  
                     antibody = antibody,  
                     stringsAsFactors = FALSE,
                     coverage = sum(diff$Length)) 
    peak_coverage <- rbind(peak_coverage,df)
  }
}
ggplot(peak_coverage,mapping = aes(x=antibody,y=coverage,fill = tissue))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+xlab("")+ylab("peak length")+ggtitle(paste0("all peaks"))+
  ggsci::scale_fill_d3()+theme(text = element_text(size = 18))
