rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(corrplot)
library(data.table)
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
tissue <- "lung"
state_num <- "15"

tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver","ileum",
                  "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT"))

summary_df <- data.frame()
for(tissue in tissues){
  file_dir <- paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/")  
  files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_young[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
  file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
  chromHMM_young <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
  chromHMM_young <- chromHMM_young[which(chromHMM_young$V1 == "chrY"),]
 
  t_summary_df <- as.data.frame(table(chromHMM_young$V4))
  t_summary_df$percent <- t_summary_df$Freq/sum(t_summary_df$Freq)*100
  colnames(t_summary_df)[3] <- tissue_label_change(tissue) 
  if(nrow(summary_df)==0){
    summary_df <- t_summary_df[,c(1,3)]
  }else{
    summary_df <- merge(summary_df,t_summary_df[,c(1,3)],by="Var1",all=T)
  }
}
to_plot<- reshape2::melt(summary_df)
to_plot$Var1 <- factor(to_plot$Var1,levels=paste0("E",1:15))
to_plot$variable <- factor(to_plot$variable,levels=c("Kidney","Muscle","Skin","Liver","Aorta","Testis","Cortex","Tongue","Uterus","Bladder","Stomach","Heart","Hippocampus","Cerebellum","BAT","Lung","Mammary Gland","Pancreas","Bone Marrow","iWAT","Cecum","Colon","Jejunum","Spleen","Thymus","Ileum","Ovary"))
color <- read.table("data/samples/20_distinct_color.txt")
color <- setNames(color$V1,paste0("E",1:15))
p<-ggplot(to_plot, aes(x = variable, y = value, fill = Var1)) +  
  geom_bar(stat = 'identity',color="white") +   
  theme_minimal() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Proportion")







