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
options(scipen = 999)
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
tissue_label_change <- function(tissue){
  if(tissue=="brain"){
    tissue_label <- "Frontal Cortex"
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

diff_bin_number <- data.frame()
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me3","H3K4me1")
for(tissue in tissues){
  for(antibody in antibodys){
    if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
      bin_size <- "10kb"
    }else{
      bin_size <- "1kb"
    }
    diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
    sig <- data.frame(Var1=c("Up","Stable","Down"),Freq=c(0,0,0))
    t_sig<-as.data.frame(table(diff$Significant))
    sig <- merge(sig, t_sig, by="Var1", all.x=TRUE) 
    sig$Freq.x <- ifelse(is.na(sig$Freq.y), sig$Freq.x, sig$Freq.y)  
    colnames(sig)[2] <- "Freq"
    sig <- sig[, -3]
    sig$tissue <- tissue_label_change(tissue)
    sig$antibody <- antibody
    diff_bin_number<-rbind(diff_bin_number,sig)
    }
}
### ATAC
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  df <- read.csv(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_macs_young_old_narrowpeak_summits_spm3_diff_after_remove_batch_effect.csv"))
  up <- df[which(df$LogFC.old.young > 0 & df$FDR.old.young < 0.05),]
  down <- df[which(df$LogFC.old.young < 0 & df$FDR.old.young < 0.05),]
  sig <- data.frame(Var1=c("Up","Down"),Freq=c(nrow(up),nrow(down)))
  colnames(sig)[2] <- "Freq"
  sig$tissue <- tissue_label_change(tissue)
  sig$antibody <- "ATAC"
  diff_bin_number<-rbind(diff_bin_number,sig)
}
### WGBS
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  df <- read.table(paste0("data/samples/WGBS/",tissue,"/DSS_table/",tissue,"_DMR_delta01.txt"),header = T)
  sig <- data.frame(Var1=c("Up","Down"),Freq=c(nrow(df[which(df$areaStat>0),]),nrow(df[which(df$areaStat<0),])))
  colnames(sig)[2] <- "Freq"
  sig$tissue <- tissue_label_change(tissue)
  sig$antibody <- "WGBS"
  diff_bin_number<-rbind(diff_bin_number,sig)
}
### RNA
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  diff <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))  
  sig <- data.frame(Var1=c("Up","Stable","Down"),Freq=c(0,0,0))
  t_sig<-as.data.frame(table(diff$Significant))
  sig <- merge(sig, t_sig, by="Var1", all.x=TRUE) 
  sig$Freq.x <- ifelse(is.na(sig$Freq.y), sig$Freq.x, sig$Freq.y)  
  colnames(sig)[2] <- "Freq"
  sig <- sig[, -3]   
  sig$tissue <- tissue_label_change(tissue)
  sig$antibody <- "RNA"
  diff_bin_number<-rbind(diff_bin_number,sig)
}

conditions <- c("Up","Down")
diff_bin_number_rank <- diff_bin_number[which(diff_bin_number$Var1 %in% conditions),]
diff_bin_number_rank  <- diff_bin_number_rank  %>%
  group_by(tissue) %>%
  summarise(total_Freq = sum(Freq, na.rm = TRUE))
diff_bin_number_rank <- diff_bin_number_rank[order(diff_bin_number_rank$total_Freq, decreasing = T),]
tissues_order <- c("Ovary","Mammary Gland","Thymus","Cerebellum","Muscle","BAT","Uterus",  
                  "Skin","Frontal Cortex","Hippocampus","Liver","Lung","Bladder","Aorta",       
                  "iWAT","Heart","Bone Marrow","Spleen","Kidney",     
                  "Tongue","Testis","Stomach","Ileum","Jejunum","Colon","Cecum","Pancreas")
diff_bin_number$tissue <- factor(diff_bin_number$tissue,levels=rev(tissues_order))
color <- setNames(c("#3c5488","#e64b35"),c("Down","Up"))

antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me3","H3K4me1","ATAC","WGBS","RNA")
p_list <- list()
for(i in c(1:length(antibodys))){
  if(antibodys[i]=="ATAC"){
    xmax=12000
  }else if(antibodys[i]=="WGBS"){
    xmax=120000
  }else if(antibodys[i]=="RNA"){
    xmax=7000
  }else{
    xmax=100000
  }
  df <- diff_bin_number[which(diff_bin_number$antibody==antibodys[i] & diff_bin_number$Var1 %in% conditions),]
  # df$tissue <- factor(df$tissue,levels = rev(sort(unique(df$tissue))))
  if(i == 1){
    p_list[[i]] <-  ggplot(df,mapping = aes(x=Freq,y=tissue,fill = Var1))+
      # geom_bar(stat = "identity",position = position_dodge2(),aes(alpha = ifelse(Var1 == "Down", 0.9, 1)),color = "black")+
      geom_bar(stat = "identity",position = position_dodge2())+
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
      scale_x_continuous(
        limits = c(0, xmax), 
        breaks = c(0, xmax)
      )
  }else{
    p_list[[i]] <-  ggplot(df,mapping = aes(x=Freq,y=tissue,fill = Var1))+
      # geom_bar(stat = "identity",position = position_dodge2(),aes(alpha = ifelse(Var1 == "Down", 0.9, 1)),color = "black")+
      geom_bar(stat = "identity",position = position_dodge2())+
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
      theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),legend.position = "none") +   
      scale_x_continuous(
        limits = c(0, xmax), 
        breaks = c(0, xmax)
      )
  }
}
combined_plot <- arrangeGrob(  
  grobs = p_list,  
  ncol = length(p_list),  
  widths = c(1.4,rep(1,(length(p_list)-1)))
)  
grid.draw(combined_plot) 
ggsave(paste0("result/figures/all_tissues_after_remove_batch_effect.pdf"), plot = combined_plot, width = 24, height = 12)


