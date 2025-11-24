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
library(grid)  
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
    }else if(tissue_label=="Mammarygland"){
      tissue_label <- "Mammary Gland"
    }else if(tissue_label=="Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
}
antibody <- "H3K9me3"
tissue <- "Hip"
diff_peak_number <-data.frame()
diff <- read.csv(paste0("data/public_data/Hippocampus_aging/Hip_young_old_merge-W1000-G3000-E100.csv"))
sig <- data.frame(Var1=c("Up","Stable","Down"),Freq=c(0,0,0))
t_sig<-as.data.frame(table(diff$Significant))
sig <- merge(sig, t_sig, by="Var1", all.x=TRUE) 
sig$Freq.x <- ifelse(is.na(sig$Freq.y), sig$Freq.x, sig$Freq.y)  
colnames(sig)[2] <- "Freq"
sig <- sig[, -3] 

sig$tissue <- tissue
sig$antibody <- antibody
diff_peak_number<-rbind(diff_peak_number,sig)

color <- setNames(c("#fc5185","#00adb5"),c("Up","Down"))

t_diff_peak_number <- diff_peak_number[which(diff_peak_number$tissue==tissue & diff_peak_number$antibody==antibody),]
t_diff_peak_number$Var1 <- factor(t_diff_peak_number$Var1,levels=c("Up","Down","Stable"))
ggplot(t_diff_peak_number[which(t_diff_peak_number$Var1!="Stable"),],mapping = aes(x=Var1,y=Freq,fill =Var1))+
  geom_bar(stat = "identity", position = position_dodge2())+
  theme_bw()+ylab("")+
  xlab(tissue_label_change(tissue))+
  theme(  
    text = element_text(size = 14),  
    legend.position = "none"  
  ) + scale_fill_manual(values = color) +    
  geom_text(aes(label = Freq), position = position_dodge2(width = 0.9), size = 5)

