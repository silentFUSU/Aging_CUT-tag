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
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")

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
mm10_10k <- read.delim("~/ref_data/mm10_10kb_bins.bed")
mm10_1k <- read.delim("~/ref_data/mm10_1kb_bins.bed")
antibodys <- c("ATAC")
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(antibodys))){
    antibody<-antibodys[j]
    diff <- read.csv(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_1kb_bins_diff.csv"))
    sig <- data.frame(Var1=c("Up","Stable","Down"),Freq=c(0,0,0))
    t_sig<-as.data.frame(table(diff$Significant))
    sig <- merge(sig, t_sig, by="Var1", all.x=TRUE) 
    sig$Freq.x <- ifelse(is.na(sig$Freq.y), sig$Freq.x, sig$Freq.y)  
    colnames(sig)[2] <- "Freq"
    sig <- sig[, -3] 
    sig$tissue <- tissue_label_change(tissue)
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

conditions <- c("Down","Up")
for(i in c(1:length(conditions))){
  condition <- conditions[i]
  df <- diff_peak_number[which(diff_peak_number$Var1==condition),]  
  color <- read.table("data/samples/30_distinct_color.txt")
  color <- setNames(color$V1,sort(unique(df$tissue)))
  p <- ggplot(df,mapping = aes(x=Freq,y=tissue,fill = tissue))+
    geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+ylab("")+xlab("ATAC Bin")+
    theme(text = element_text(size = 13))+ scale_fill_manual(values = color) +theme(legend.position = "none") + 
    geom_text(aes(label = Freq), position = position_dodge2(width = 0.9), hjust = 0.4, size = 5) + xlim(0,27000)+ggtitle(paste0(condition," ATAC Bin number"))
  print(p)
}
