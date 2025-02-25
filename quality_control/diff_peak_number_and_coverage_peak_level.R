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
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
# antibodys <- c("ATAC")
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
diff_peak_number <-data.frame()
mm10_10k <- read.delim("~/ref_data/mm10_10kb_bins.bed")
mm10_1k <- read.delim("~/ref_data/mm10_1kb_bins.bed")
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me3","H3K4me1")
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(antibodys))){
    antibody<-antibodys[j]
    if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
      diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_young_old_merge-W1000-G3000-E100_diff_after_remove_batch_effect.csv"))
    }else{
      diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_macs_young_old_narrowpeak_diff_after_remove_batch_effect.csv"))
    }
    sig <- data.frame(Var1=c("Up","Stable","Down"),Freq=c(0,0,0))
    t_sig<-as.data.frame(table(diff$Significant))
    sig <- merge(sig, t_sig, by="Var1", all.x=TRUE) 
    sig$Freq.x <- ifelse(is.na(sig$Freq.y), sig$Freq.x, sig$Freq.y)  
    colnames(sig)[2] <- "Freq"
    sig <- sig[, -3] 
    
    sig$tissue <- tissue
    sig$antibody <- antibody
    diff_peak_number<-rbind(diff_peak_number,sig)
  }
}

color <- read.table("data/samples/30_distinct_color.txt")
color <- color$V1
diff_peak_number <- diff_peak_number %>%  
  mutate(tissue_label = sapply(tissue, tissue_label_change))  
color <- setNames(color,sort(unique(diff_peak_number$tissue_label)))
p_list <- list()
conditions <- c("Down","Up")
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me3","H3K4me1")
for(condition in conditions){
  for(i in c(1:length(antibodys))){
    df <- diff_peak_number[which(diff_peak_number$antibody==antibodys[i] & diff_peak_number$Var1==condition),]
    # df <- arrange(df, Freq)  
    df$tissue <- factor(df$tissue,levels=sort(tissues))
    if(i == 1){
      p_list[[i]] <-  ggplot(df,mapping = aes(x=Freq,y=tissue_label,fill = tissue_label))+
        geom_bar(stat = "identity", position = position_dodge2())+
        theme_bw()+ylab("")+
        xlab(antibodys[[i]])+
        theme(  
          text = element_text(size = 10),  
          axis.text.y = element_text(size = 11),  # 改变 y 轴刻度标签的字体大小  
          axis.title.y = element_blank(),  
          axis.ticks.y = element_blank(),  
          legend.position = "none"  
        ) +
        scale_fill_manual(values = color) +
        theme(legend.position = "none") + 
        geom_text(aes(label = Freq), position = position_dodge2(width = 0.9), hjust = 0.1, size = 3) + xlim(0,40000)
    }else{
      p_list[[i]] <-  ggplot(df,mapping = aes(x=Freq,y=tissue_label,fill = tissue_label))+
        geom_bar(stat = "identity", position = position_dodge2())+
        theme_bw()+ylab("")+
        xlab(antibodys[[i]])+
        theme(  
          text = element_text(size = 10),  
          axis.text.y = element_text(size = 11),  # 改变 y 轴刻度标签的字体大小  
          axis.title.y = element_blank(),  
          axis.ticks.y = element_blank(),  
          legend.position = "none"  
        ) + scale_fill_manual(values = color) +    
        theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),legend.position = "none") +   
        geom_text(aes(label = Freq), position = position_dodge2(width = 0.9), hjust = 0.1, size = 3) + xlim(0,40000)
    }
  }
  combined_plot <- arrangeGrob(  
    grobs = p_list,  
    ncol = length(p_list),  
    widths = c(1.4,rep(1,(length(p_list)-1))),
    top = textGrob(paste0(condition," peak number"), gp = gpar(fontsize = 15, fontface = "bold"))  
  )  
  grid.draw(combined_plot) 
  ggsave(paste0("result/all/diff/all_tissues_",condition,"_peak_number_remove_batch_effect.png"), plot = combined_plot, width = 18, height = 6,type="cairo")  
}
