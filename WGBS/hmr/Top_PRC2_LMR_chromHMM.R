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
tissue <- "colon"
state_num <- "15"
summary_df <- data.frame()
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
                  "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT","ileum")) 
for(tissue in tissues){
  regions <- read.table(paste0("data/samples/WGBS/",tissue,"/hmr/PRC2_LMR_top1000.bed"))
  regions <- as.data.table(regions)
  setDT(regions)
  setkey(regions,V1,V2,V3)
  file_dir <- paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/split_1k/")  
  files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_young[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
  file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
  chromHMM_young <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
  chromHMM_young <- chromHMM_young[which(chromHMM_young$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  chromHMM_young$V2 <- chromHMM_young$V2+1
  chromHMM_young <- data.table(chromHMM_young)
  setDT(chromHMM_young) 
  setkey(chromHMM_young, V1, V2, V3) 
  
  overlaps_young <- foverlaps(regions , chromHMM_young, type = "any", nomatch = 0L)  
  t_summary_df <- as.data.frame(table(overlaps_young$V4))
  t_summary_df$percent <- t_summary_df$Freq/sum(t_summary_df$Freq)*100
  colnames(t_summary_df)[3] <- tissue_label_change(tissue) 
  if(nrow(summary_df)==0){
    summary_df <- t_summary_df[,c(1,3)]
  }else{
    summary_df <- merge(summary_df,t_summary_df[,c(1,3)],by="Var1",all=T)
  }
}
dictionary <- list("E1"=1, "E2"=2, "E3"=3,
                   "E4"=4, "E5"=5, "E6"=7,
                   "E7"=8, "E8"=6, "E9"=9,
                   "E10"=10,"E11"=11,"E12"=15,
                   "E13"=12,"E14"=13,"E15"=14)
keys <- names(dictionary)
values <- unlist(dictionary)
to_plot <- summary_df

to_plot$Var1 <- values[match(to_plot$Var1, keys)]
to_plot <- to_plot[order(to_plot$Var1),]
to_plot[is.na(to_plot)] <- 0

to_plot$Var1 <- paste0("state",to_plot$Var1)
to_plot <- reshape2::melt(to_plot)

to_plot$Var1 <- factor(to_plot$Var1,levels=paste0("state",1:15))
color <- read.table("data/samples/20_distinct_color.txt")
color <- setNames(color$V1,paste0("state",1:15))
ggplot(to_plot, aes(x = variable, y = value, fill = Var1)) +  
  geom_bar(stat = "identity") +  
  xlab(NULL)+
  ylab(NULL)+
  scale_fill_manual(values = color) +  
  theme(text = element_text(size = 20),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(fill = guide_legend(title = "State"))  
