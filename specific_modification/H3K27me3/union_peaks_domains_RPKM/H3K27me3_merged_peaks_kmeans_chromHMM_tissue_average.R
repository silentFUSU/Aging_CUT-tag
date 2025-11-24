rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(corrplot)
library(data.table)
tissue <- "lung"
state_num <- 15

kmeans_df <- read.csv("data/samples/all/H3K27me3/peaks_merged/kmeans_annotation_larger_27_tissues_RPKM.csv")
colnames(kmeans_df)[1] <- "label"
split_chr <- strsplit(as.character(kmeans_df$label), ":")  
chr_column <- sapply(split_chr, `[[`, 1)  
split_start_end <- strsplit(sapply(split_chr, `[[`, 2), "-")  
start_column <- sapply(split_start_end, `[[`, 1)  
end_column <- sapply(split_start_end, `[[`, 2)  
kmeans_df <- data.frame(chr = chr_column,start = start_column, end = end_column, cluster=kmeans_df$cluster)
kmeans_df$start <- as.numeric(kmeans_df$start)
kmeans_df$end <- as.numeric(kmeans_df$end)
kmeans_df <- as.data.table(kmeans_df)
setDT(kmeans_df)
setkey(kmeans_df,chr,start,end)
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver","ileum",
                  "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT"))
kmeans_summary <- data.frame()
for(kmean in c(1:3)){
  summary_df <- data.frame()
  for(tissue in tissues){
    file_dir <- paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/")  
    files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_old[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
    file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
    chromHMM_old <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
    if(tissue %in% c("mammarygland","uterus","ovary")){
      chromHMM_old <- chromHMM_old[which(chromHMM_old$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
    }else{
      chromHMM_old <- chromHMM_old[which(chromHMM_old$V1 %in% paste0("chr",c(1:19,"X"))),]
    }
    
    chromHMM_old$V2 <- chromHMM_old$V2+1
    chromHMM_old <- data.table(chromHMM_old)
    setDT(chromHMM_old) 
    setkey(chromHMM_old, V1, V2, V3) 
    if(tissue %in% c("mammarygland","uterus","ovary")){
      overlaps_old <- foverlaps(kmeans_df[cluster == kmean & chr %in% paste0("chr", c(1:19, "X"))] , chromHMM_old, type = "any", nomatch = 0L)  
    }else{
      overlaps_old <- foverlaps(kmeans_df[cluster == kmean] , chromHMM_old, type = "any", nomatch = 0L)  
    }
    t_summary_df <- as.data.frame(table(overlaps_old$V4))
    t_summary_df$percent <- t_summary_df$Freq/sum(t_summary_df$Freq)*100
    colnames(t_summary_df)[3] <- tissue
    if(nrow(summary_df)==0){
      summary_df <- t_summary_df[,c(1,3)]
    }else{
      summary_df <- merge(summary_df,t_summary_df[,c(1,3)],by="Var1")
    }
  }
  summary_df$mean <- rowMeans(summary_df[,-1])
  summary_df <- summary_df[,c("Var1","mean")]
  colnames(summary_df) <- c("State",paste0(kmean))
  if(nrow(kmeans_summary)==0){
    kmeans_summary <- summary_df
  }else{
    kmeans_summary <- merge(kmeans_summary,summary_df,by="State")
  }
}
to_plot<- reshape2::melt(kmeans_summary)
to_plot$State <- sub("^E", "state", to_plot$State)
to_plot$State <- factor(to_plot$State,levels=paste0("state",1:state_num))
to_plot$variable <- paste0("kmeans",to_plot$variable)
color <- read.table("data/samples/20_distinct_color.txt")
color <- setNames(color$V1,paste0("state",1:state_num))
p<-ggplot(to_plot, aes(x = variable, y = value, fill = State)) +  
  geom_bar(stat = 'identity',color="black") +   
  theme_minimal() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Proportion")
