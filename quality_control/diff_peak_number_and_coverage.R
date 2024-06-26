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
tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas","Hip","cecum","bonemarrow","ileum","heart","thymus","stomach","skin","aorta","tongue","bladder")
# tissues <- c("Hip","testis", "colon", "kidney", "lung", "spleen", "muscle", "pancreas","cecum","bonemarrow","ileum","heart","thymus")
# antibodys <- c("ATAC")
############## all bins ########################
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

mm10_10k <- read.delim("~/ref_data/mm10_10kb_bins.bed")
mm10_1k <- read.delim("~/ref_data/mm10_1kb_bins.bed")
antibodys <- c("H3K27me3","H3K9me3","H3K36me3")
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(antibodys))){
    antibody<-antibodys[j]
    diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_10kb_bins_diff_after_remove_batch_effect.csv"))
    sig <- data.frame(Var1=c("Up","Stable","Down"),Freq=c(0,0,0))
    t_sig<-as.data.frame(table(diff$Significant_bar))
    sig <- merge(sig, t_sig, by="Var1", all.x=TRUE) 
    sig$Freq.x <- ifelse(is.na(sig$Freq.y), sig$Freq.x, sig$Freq.y)  
    colnames(sig)[2] <- "Freq"
    sig <- sig[, -3] 
    sig$tissue <- tissue
    sig$antibody <- antibody
    diff_peak_number<-rbind(diff_peak_number,sig)
    t_sum=sum(sig$Freq)
    sig$Freq=sig$Freq/nrow(mm10_10k)
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
    sig <- data.frame(Var1=c("Up","Stable","Down"),Freq=c(0,0,0))
    t_sig<-as.data.frame(table(diff$Significant_bar))
    sig <- merge(sig, t_sig, by="Var1", all.x=TRUE) 
    sig$Freq.x <- ifelse(is.na(sig$Freq.y), sig$Freq.x, sig$Freq.y)  
    colnames(sig)[2] <- "Freq"
    sig <- sig[, -3] 
    sig$tissue <- tissue
    sig$antibody <- antibody
    diff_peak_number<-rbind(diff_peak_number,sig)
    t_sum=sum(sig$Freq)
    sig$Freq=sig$Freq/nrow(mm10_1k)
    diff_peak_percent <- rbind(diff_peak_percent,sig)
    sig$coverage<-0
    for (k in c(1:nrow(sig))){
      sig$coverage[k]<-sum(diff$Length[which(diff$Significant_bar == sig$Var1[k])])
    }
    diff_peak_coverage <- rbind(diff_peak_coverage,sig)
  }
}
write.csv(diff_peak_number,"data/samples/all/diff_bins_number.csv",row.names = F)
write.csv(diff_peak_percent,"data/samples/all/diff_bins_percent.csv",row.names = F)
write.csv(diff_peak_coverage,"data/samples/all/diff_bins_coverage.csv",row.names = F)

############### bins overlap with peaks #######################
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
antibodys <- c("H3K27me3","H3K9me3","H3K36me3")
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(antibodys))){
    antibody<-antibodys[j]
    diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_10kb_bins_overlap_old_only_peaks_diff_after_remove_batch_effect.csv"))
    sig <- data.frame(Var1=c("Up","Stable","Down"),Freq=c(0,0,0))
    t_sig<-as.data.frame(table(diff$Significant_bar))
    sig <- merge(sig, t_sig, by="Var1", all.x=TRUE) 
    sig$Freq.x <- ifelse(is.na(sig$Freq.y), sig$Freq.x, sig$Freq.y)  
    colnames(sig)[2] <- "Freq"
    sig <- sig[, -3] 
    sig$tissue <- tissue
    sig$antibody <- antibody
    diff_peak_number<-rbind(diff_peak_number,sig)
    t_sum=sum(sig$Freq)
    sig$Freq=sig$Freq/nrow(mm10_10k)
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
    diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_1kb_bins_overlap_peaks_diff_after_remove_batch_effect.csv"))
    sig <- data.frame(Var1=c("Up","Stable","Down"),Freq=c(0,0,0))
    t_sig<-as.data.frame(table(diff$Significant_bar))
    sig <- merge(sig, t_sig, by="Var1", all.x=TRUE) 
    sig$Freq.x <- ifelse(is.na(sig$Freq.y), sig$Freq.x, sig$Freq.y)  
    colnames(sig)[2] <- "Freq"
    sig <- sig[, -3] 
    sig$tissue <- tissue
    sig$antibody <- antibody
    diff_peak_number<-rbind(diff_peak_number,sig)
    t_sum=sum(sig$Freq)
    sig$Freq=sig$Freq/nrow(mm10_1k)
    diff_peak_percent <- rbind(diff_peak_percent,sig)
    sig$coverage<-0
    for (k in c(1:nrow(sig))){
      sig$coverage[k]<-sum(diff$Length[which(diff$Significant_bar == sig$Var1[k])])
    }
    diff_peak_coverage <- rbind(diff_peak_coverage,sig)
  }
}
write.csv(diff_peak_number,"data/samples/all/diff_bins_overlap_peaks_number.csv",row.names = F)
write.csv(diff_peak_percent,"data/samples/all/diff_bins_overlap_peaks_percent.csv",row.names = F)
write.csv(diff_peak_coverage,"data/samples/all/diff_bins_overlap_peaks_coverage.csv",row.names = F)


##################### PLOT ######################
custom_colors <- c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#800000", "#008000", "#000080", "#808000", "#800080", "#008080", "#C0C0C0", "#808080", "#FFA500", "#FF1493", "#00FF00")  
diff_peak_coverage$tissue[which(diff_peak_coverage$tissue=="brain")] <- "brain_FC"
diff_peak_coverage$tissue[which(diff_peak_coverage$tissue=="Hip")] <- "brain_Hip"

diff_peak_number$tissue[which(diff_peak_number$tissue=="brain")] <- "brain_FC"
diff_peak_number$tissue[which(diff_peak_number$tissue=="Hip")] <- "brain_Hip"

diff_peak_percent$tissue[which(diff_peak_percent$tissue=="brain")] <- "brain_FC"
diff_peak_percent$tissue[which(diff_peak_percent$tissue=="Hip")] <- "brain_Hip"

ggplot(diff_peak_number[which(diff_peak_number$Var1 =="Down"),],mapping = aes(x=antibody,y=Freq,fill = tissue))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+xlab("")+ylab("Bins")+ggtitle(paste0("Decrease"))+
  theme(text = element_text(size = 18))+ scale_fill_d3("category20")

ggplot(diff_peak_number[which(diff_peak_number$Var1 =="Up"),],mapping = aes(x=antibody,y=Freq,fill = tissue))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+xlab("")+ylab("Bins")+ggtitle(paste0("Increase"))+
  theme(text = element_text(size = 18))+scale_fill_d3("category20")



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


diff_peak_percent <- read.csv("data/samples/all/diff_bins_overlap_peaks_percent.csv")
diff_peak_percent$tissue[which(diff_peak_percent$tissue=="brain")] <- "brain_FC"
diff_peak_percent$tissue[which(diff_peak_percent$tissue=="Hip")] <- "brain_Hip"
H3K27me3 <- diff_peak_percent[which(diff_peak_percent$antibody=="H3K27me3"),]
H3K27me3 <- H3K27me3[-which(H3K27me3$Var1=="Stable"),]
H3K27me3$Var1 <- factor(H3K27me3$Var1,levels=c("Stable","Up","Down"))
ggplot(H3K27me3, aes( x = tissue, weight = Freq, fill = Var1))+
  geom_bar( position = "stack")+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(text = element_text(size = 15))+labs(x = NULL, y =  "Percent", fill = NULL) 
H3K9me3 <- diff_peak_percent[which(diff_peak_percent$antibody=="H3K9me3"),]
H3K9me3 <- H3K9me3[-which(H3K9me3$Var1=="Stable"),]
H3K9me3$Var1 <- factor(H3K9me3$Var1,levels=c("Stable","Up","Down"))

ggplot(H3K9me3, aes( x = tissue, weight = Freq, fill = Var1))+
  geom_bar( position = "stack")+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(text = element_text(size = 15))+labs(x = NULL, y = "Percent", fill = NULL) 









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


diff_peak_percent$group <- "other tissues"
diff_peak_percent$group[which(diff_peak_percent$tissue %in% c("cecum","colon","ileum"))] <- "intestine"

ggplot(diff_peak_percent[which(diff_peak_percent$Var1 =="Up" & diff_peak_percent$antibody %in% c("H3K27me3","H3K9me3")),],aes(x=antibody,y=Freq,fill =group))+
  geom_boxplot()+
  # geom_jitter(shape=16,size=1,position = position_jitter(0.2))+
  scale_fill_brewer(palette="Set3")+
  theme_bw()+theme(text = element_text(size = 18))+ylab("Percent")+
  xlab("")+labs(fill = "", color = "",title = "heterochromatin marker increase")

ggplot(diff_peak_percent[which(diff_peak_percent$Var1 =="Down" & diff_peak_percent$antibody %in% c("H3K27me3","H3K9me3")),],aes(x=antibody,y=Freq,fill =group))+
  geom_boxplot()+
  # geom_jitter(shape=16,size=1,position = position_jitter(0.2))+
  scale_fill_brewer(palette="Set3")+
  theme_bw()+theme(text = element_text(size = 18))+ylab("Percent")+
  xlab("")+labs(fill = "", color = "",title = "heterochromatin marker decrease")
wilcox.test(diff_peak_percent$Freq[which(diff_peak_percent$Var1 =="Down" & diff_peak_percent$antibody %in% c("H3K27me3") & diff_peak_percent$group == "intestine")],
            diff_peak_percent$Freq[which(diff_peak_percent$Var1 =="Down" & diff_peak_percent$antibody %in% c("H3K27me3") & diff_peak_percent$group == "other tissues")])
wilcox.test(diff_peak_percent$Freq[which(diff_peak_percent$Var1 =="Down" & diff_peak_percent$antibody %in% c("H3K9me3") & diff_peak_percent$group == "intestine")],
            diff_peak_percent$Freq[which(diff_peak_percent$Var1 =="Down" & diff_peak_percent$antibody %in% c("H3K9me3") & diff_peak_percent$group == "other tissues")])

wilcox.test(diff_peak_percent$Freq[which(diff_peak_percent$Var1 =="Up" & diff_peak_percent$antibody %in% c("H3K27me3") & diff_peak_percent$group == "intestine")],
            diff_peak_percent$Freq[which(diff_peak_percent$Var1 =="Up" & diff_peak_percent$antibody %in% c("H3K27me3") & diff_peak_percent$group == "other tissues")])
wilcox.test(diff_peak_percent$Freq[which(diff_peak_percent$Var1 =="Up" & diff_peak_percent$antibody %in% c("H3K9me3") & diff_peak_percent$group == "intestine")],
            diff_peak_percent$Freq[which(diff_peak_percent$Var1 =="Up" & diff_peak_percent$antibody %in% c("H3K9me3") & diff_peak_percent$group == "other tissues")])


diff_peak_number<- read.csv("data/samples/all/diff_bins_number.csv")
diff_peak_number$group <- "other tissues"
diff_peak_number$group[which(diff_peak_number$tissue %in% c("cecum","colon","ileum"))] <- "intestine"

ggplot(diff_peak_number[which(diff_peak_number$Var1 =="Up"),],aes(x=antibody,y=Freq,fill =group))+
  geom_boxplot()+
  # geom_jitter(shape=16,size=1,position = position_jitter(0.2))+
  scale_fill_brewer(palette="Set3")+
  theme_bw()+theme(text = element_text(size = 18))+ylab("Bins")+
  xlab("")+labs(fill = "", color = "",title = "Increase")

ggplot(diff_peak_number[which(diff_peak_number$Var1 =="Down"),],aes(x=antibody,y=Freq,fill =group))+
  geom_boxplot()+
  # geom_jitter(shape=16,size=1,position = position_jitter(0.2))+
  scale_fill_brewer(palette="Set3")+
  theme_bw()+theme(text = element_text(size = 18))+ylab("Bins")+
  xlab("")+labs(fill = "", color = "",title = "Decrease")


wilcox.test(diff_peak_number$Freq[which(diff_peak_number$Var1 =="Up" & diff_peak_number$antibody %in% c("H3K27me3") & diff_peak_number$group == "intestine")],
            diff_peak_number$Freq[which(diff_peak_number$Var1 =="Up" & diff_peak_number$antibody %in% c("H3K27me3") & diff_peak_number$group == "other tissues")])
wilcox.test(diff_peak_number$Freq[which(diff_peak_number$Var1 =="Up" & diff_peak_number$antibody %in% c("H3K9me3") & diff_peak_number$group == "intestine")],
            diff_peak_number$Freq[which(diff_peak_number$Var1 =="Up" & diff_peak_number$antibody %in% c("H3K9me3") & diff_peak_number$group == "other tissues")])
wilcox.test(diff_peak_number$Freq[which(diff_peak_number$Var1 =="Up" & diff_peak_number$antibody %in% c("H3K36me3") & diff_peak_number$group == "intestine")],
            diff_peak_number$Freq[which(diff_peak_number$Var1 =="Up" & diff_peak_number$antibody %in% c("H3K36me3") & diff_peak_number$group == "other tissues")])
wilcox.test(diff_peak_number$Freq[which(diff_peak_number$Var1 =="Up" & diff_peak_number$antibody %in% c("H3K4me3") & diff_peak_number$group == "intestine")],
            diff_peak_number$Freq[which(diff_peak_number$Var1 =="Up" & diff_peak_number$antibody %in% c("H3K4me3") & diff_peak_number$group == "other tissues")])
wilcox.test(diff_peak_number$Freq[which(diff_peak_number$Var1 =="Up" & diff_peak_number$antibody %in% c("H3K4me1") & diff_peak_number$group == "intestine")],
            diff_peak_number$Freq[which(diff_peak_number$Var1 =="Up" & diff_peak_number$antibody %in% c("H3K4me1") & diff_peak_number$group == "other tissues")])
wilcox.test(diff_peak_number$Freq[which(diff_peak_number$Var1 =="Up" & diff_peak_number$antibody %in% c("H3K27ac") & diff_peak_number$group == "intestine")],
            diff_peak_number$Freq[which(diff_peak_number$Var1 =="Up" & diff_peak_number$antibody %in% c("H3K27ac") & diff_peak_number$group == "other tissues")])

wilcox.test(diff_peak_number$Freq[which(diff_peak_number$Var1 =="Down" & diff_peak_number$antibody %in% c("H3K27me3") & diff_peak_number$group == "intestine")],
            diff_peak_number$Freq[which(diff_peak_number$Var1 =="Down" & diff_peak_number$antibody %in% c("H3K27me3") & diff_peak_number$group == "other tissues")])
wilcox.test(diff_peak_number$Freq[which(diff_peak_number$Var1 =="Down" & diff_peak_number$antibody %in% c("H3K9me3") & diff_peak_number$group == "intestine")],
            diff_peak_number$Freq[which(diff_peak_number$Var1 =="Down" & diff_peak_number$antibody %in% c("H3K9me3") & diff_peak_number$group == "other tissues")])
wilcox.test(diff_peak_number$Freq[which(diff_peak_number$Var1 =="Down" & diff_peak_number$antibody %in% c("H3K36me3") & diff_peak_number$group == "intestine")],
            diff_peak_number$Freq[which(diff_peak_number$Var1 =="Down" & diff_peak_number$antibody %in% c("H3K36me3") & diff_peak_number$group == "other tissues")])
wilcox.test(diff_peak_number$Freq[which(diff_peak_number$Var1 =="Down" & diff_peak_number$antibody %in% c("H3K4me3") & diff_peak_number$group == "intestine")],
            diff_peak_number$Freq[which(diff_peak_number$Var1 =="Down" & diff_peak_number$antibody %in% c("H3K4me3") & diff_peak_number$group == "other tissues")])
wilcox.test(diff_peak_number$Freq[which(diff_peak_number$Var1 =="Down" & diff_peak_number$antibody %in% c("H3K4me1") & diff_peak_number$group == "intestine")],
            diff_peak_number$Freq[which(diff_peak_number$Var1 =="Down" & diff_peak_number$antibody %in% c("H3K4me1") & diff_peak_number$group == "other tissues")])
wilcox.test(diff_peak_number$Freq[which(diff_peak_number$Var1 =="Down" & diff_peak_number$antibody %in% c("H3K27ac") & diff_peak_number$group == "intestine")],
            diff_peak_number$Freq[which(diff_peak_number$Var1 =="Down" & diff_peak_number$antibody %in% c("H3K27ac") & diff_peak_number$group == "other tissues")])

