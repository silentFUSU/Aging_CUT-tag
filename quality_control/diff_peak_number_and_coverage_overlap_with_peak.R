rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(ggrepel)
library(gridExtra)
library(ggpubr)
library(grid)
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
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
# diff_peak_number <- read.csv("data/samples/all/diff_bins_number.csv")
# diff_peak_coverage <- read.csv("data/samples/all/diff_bins_coverage.csv")
# diff_peak_percent <- read.csv("data/samples/all/diff_bins_percent.csv")
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(antibodys))){
    antibody<-antibodys[j]
    diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_10kb_bins_diff_overlap_young_old_peaks.csv"))
    sig <- data.frame(Var1=c("Up","Stable","Down"),Freq=c(0,0,0))
    t_sig<-as.data.frame(table(diff$Significant))
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
      sig$coverage[k]<-sum(diff$Length[which(diff$Significant == sig$Var1[k])])
    }
    diff_peak_coverage <- rbind(diff_peak_coverage,sig)
  }
}
antibodys <- c("H3K27ac","H3K4me3","H3K4me1")
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(antibodys))){
    antibody<-antibodys[j]
    diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_1kb_bins_diff_overlap_young_old_peaks.csv"))
    sig <- data.frame(Var1=c("Up","Stable","Down"),Freq=c(0,0,0))
    t_sig<-as.data.frame(table(diff$Significant))
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
      sig$coverage[k]<-sum(diff$Length[which(diff$Significant == sig$Var1[k])])
    }
    diff_peak_coverage <- rbind(diff_peak_coverage,sig)
  }
}
write.csv(diff_peak_number,"data/samples/all/diff_bins_overlap_peaks_number.csv",row.names = F)
write.csv(diff_peak_percent,"data/samples/all/diff_bins_overlap_peaks_percent.csv",row.names = F)
write.csv(diff_peak_coverage,"data/samples/all/diff_bins_overlap_peaks_coverage.csv",row.names = F)

##################### PLOT ######################
color <- read.table("data/samples/30_distinct_color.txt")
color <- color$V1
tissues[which(tissues=="brain")] <- "Cortex"
tissues[which(tissues=="Hip")] <- "Hippocampus"
tissues[which(tissues=="CB")] <- "Cerebellum"
tissues <- str_to_title(tissues)
tissues[which(tissues=="Bonemarrow")] <- "Bone Marrow"
tissues[which(tissues=="Bat")] <- "BAT"
tissues[which(tissues=="Mammarygland")] <- "Mammary Gland"
tissues[which(tissues=="Iwat")] <- "iWAT"
color <- setNames(color,tissues)
diff_peak_number$tissue[which(diff_peak_number$tissue=="brain")] <- "Cortex"
diff_peak_number$tissue[which(diff_peak_number$tissue=="Hip")] <- "Hippocampus"
diff_peak_number$tissue[which(diff_peak_number$tissue=="CB")] <- "Cerebellum"
diff_peak_number$tissue <- str_to_title(diff_peak_number$tissue)
diff_peak_number$tissue[which(diff_peak_number$tissue=="Bonemarrow")] <- "Bone Marrow"
diff_peak_number$tissue[which(diff_peak_number$tissue=="Bat")] <- "BAT"
diff_peak_number$tissue[which(diff_peak_number$tissue=="Mammarygland")] <- "Mammary Gland"
diff_peak_number$tissue[which(diff_peak_number$tissue=="Iwat")] <- "iWAT"

p_list <- list()
conditions <- c("Down","Up")
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me3","H3K4me1")
for(condition in conditions){
  for(i in c(1:length(antibodys))){
    df <- diff_peak_number[which(diff_peak_number$antibody==antibodys[i] & diff_peak_number$Var1==condition),]
    # df <- arrange(df, Freq)  
    df$tissue <- factor(df$tissue,levels=sort(tissues))
    if(i == 1){
      p_list[[i]] <-  ggplot(df,mapping = aes(x=Freq,y=tissue,fill = tissue))+
        geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+ylab("")+xlab(antibodys[[i]])+
        theme(text = element_text(size = 13))+ scale_fill_manual(values = color) +theme(legend.position = "none") + 
        geom_text(aes(label = Freq), position = position_dodge2(width = 0.9), hjust = 0.4, size = 5) + xlim(0,80000)
    }else{
      p_list[[i]] <-  ggplot(df,mapping = aes(x=Freq,y=tissue,fill = tissue))+
        geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+ylab("")+xlab(antibodys[[i]])+
        theme(text = element_text(size = 13))+ scale_fill_manual(values = color) +    
        theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),legend.position = "none") +   
        geom_text(aes(label = Freq), position = position_dodge2(width = 0.9), hjust = 0.4, size = 5) + xlim(0,80000)
    }
  }
  combined_plot <- arrangeGrob(  
    grobs = p_list,  
    ncol = length(p_list),  
    widths = c(1.7,rep(1,(length(p_list)-1))),
    top = textGrob(paste0(condition," bin number in peak"), gp = gpar(fontsize = 15, fontface = "bold"))  
  )  
  grid.draw(combined_plot) 
  ggsave(paste0("result/all/diff/all_tissues_",condition,"_bin_overlap_peaks_number.png"), plot = combined_plot, width = 20, height = 6,type="cairo")  
}

