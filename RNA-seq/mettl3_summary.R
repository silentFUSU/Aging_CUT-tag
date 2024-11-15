rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)

tissues <- c("aorta","BAT","bladder","bonemarrow","FC","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
tissue_label_change <- function(tissue){
  if(tissue=="FC"){
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
mettl14_summary <- data.frame()
mettl3_summary <- data.frame()
Bap1_summary <- data.frame()
setdb1_summary <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene.csv"))
  mettl14 <- df[which(df$X=="Mettl14"),]
  mettl14 <- mettl14[,c("logFC","fdr","Significant")]
  mettl14$tissue <- tissue_label_change(tissue)
  mettl3 <- df[which(df$X=="Mettl3"),]
  mettl3 <- mettl3[,c("logFC","fdr","Significant")]
  mettl3$tissue <- tissue_label_change(tissue)
  
  Bap1 <- df[which(df$X=="Bap1"),]
  Bap1 <- Bap1[,c("logFC","fdr","Significant")]
  Bap1$tissue <- tissue_label_change(tissue)
  setdb1 <- df[which(df$X=="Setdb1"),]
  setdb1 <- setdb1[,c("logFC","fdr","Significant")]
  setdb1$tissue <- tissue_label_change(tissue)
  if(nrow(mettl14_summary)==0){
    mettl14_summary <- mettl14
  }else{
    mettl14_summary <- rbind(mettl14_summary,mettl14)
  }
  
  if(nrow(mettl3_summary)==0){
    mettl3_summary <- mettl3
  }else{
    mettl3_summary <- rbind(mettl3_summary,mettl3)
  }
  
  if(nrow(Bap1_summary)==0){
    Bap1_summary <- Bap1
  }else{
    Bap1_summary <- rbind(Bap1_summary,Bap1)
  }
  if(nrow(setdb1_summary)==0){
    setdb1_summary <- setdb1
  }else{
    setdb1_summary <- rbind(setdb1_summary,setdb1)
  }
  
}

ggplot(mettl3_summary,mapping = aes(x=logFC,y=tissue,fill = Significant))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+ylab("")+
  theme(text = element_text(size = 13))+ 
  ggtitle("Mettl3 RNA") +
  geom_text(data = mettl3_summary[which(mettl3_summary$logFC<0),],aes(label = round(logFC,3)), position = position_dodge2(width = 0.9), hjust = 0.5, size = 5) +
  geom_text(data = mettl3_summary[which(mettl3_summary$logFC>0),],aes(label = round(logFC,3)), position = position_dodge2(width = 0.9), hjust = 0.5, size = 5)

ggplot(mettl14_summary,mapping = aes(x=logFC,y=tissue,fill = Significant))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+ylab("")+
  theme(text = element_text(size = 13))+ 
  ggtitle("Mettl14 RNA") +
  geom_text(data = mettl14_summary[which(mettl14_summary$logFC<0),],aes(label = round(logFC,3)), position = position_dodge2(width = 0.9), hjust = 0.5, size = 5) +
  geom_text(data = mettl14_summary[which(mettl14_summary$logFC>0),],aes(label = round(logFC,3)), position = position_dodge2(width = 0.9), hjust = 0.5, size = 5)


ggplot(setdb1_summary,mapping = aes(x=logFC,y=tissue,fill = Significant))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+ylab("")+
  theme(text = element_text(size = 13))+ 
  ggtitle("setdb1 RNA") +
  geom_text(data = setdb1_summary[which(setdb1_summary$logFC<0),],aes(label = round(logFC,3)), position = position_dodge2(width = 0.9), hjust = 0.5, size = 5) +
  geom_text(data = setdb1_summary[which(setdb1_summary$logFC>0),],aes(label = round(logFC,3)), position = position_dodge2(width = 0.9), hjust = 0.5, size = 5)



ggplot(Bap1_summary,mapping = aes(x=logFC,y=tissue,fill = Significant))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+ylab("")+
  theme(text = element_text(size = 13))+ 
  ggtitle("BAP1 RNA") +
  geom_text(data = Bap1_summary[which(Bap1_summary$logFC<0),],aes(label = round(logFC,3)), position = position_dodge2(width = 0.9), hjust = 0.5, size = 5) +
  geom_text(data = Bap1_summary[which(Bap1_summary$logFC>0),],aes(label = round(logFC,3)), position = position_dodge2(width = 0.9), hjust = 0.5, size = 5)
