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
state_num <- "11"
condition <- "peak"
if(condition == "domain"){
  kmeans_df <- read.csv("data/samples/all/H3K27me3/edd_domain_merged/kmeans_annotation_RPKM.csv")
}else{
  kmeans_df <- read.csv("data/samples/all/H3K27me3/peaks_merged/kmeans_annotation_larger_27_tissues_RPKM.csv")
}
split_chr <- strsplit(as.character(kmeans_df$X), ":")  
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
for(kmean in c(1:3)){
  summary_df <- data.frame()
  for(tissue in tissues){
    file_dir <- paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/")  
    files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_young[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
    file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
    chromHMM_young <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
    chromHMM_young <- chromHMM_young[which(chromHMM_young$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
    chromHMM_young$V2 <- chromHMM_young$V2+1
    chromHMM_young <- data.table(chromHMM_young)
    setDT(chromHMM_young) 
    setkey(chromHMM_young, V1, V2, V3) 
    
    overlaps_young <- foverlaps(kmeans_df[cluster == kmean] , chromHMM_young, type = "any", nomatch = 0L)  
    t_summary_df <- as.data.frame(table(overlaps_young$V4))
    t_summary_df$percent <- t_summary_df$Freq/sum(t_summary_df$Freq)*100
    colnames(t_summary_df)[3] <- tissue_label_change(tissue) 
    if(nrow(summary_df)==0){
      summary_df <- t_summary_df[,c(1,3)]
    }else{
      summary_df <- merge(summary_df,t_summary_df[,c(1,3)],by="Var1")
    }
  }
  to_plot<- reshape2::melt(summary_df)
  to_plot$Var1 <- factor(to_plot$Var1,levels=paste0("E",1:11))
  to_plot$variable <- factor(to_plot$variable,levels=c("Uterus","Ovary","Thymus","Liver","Bladder","Kidney","Pancreas","Stomach","Mammary Gland","Skin","Cecum","Ileum","Tongue","iWAT","Testis",
                                                       "Jejunum","Spleen","Cerebellum","Lung","Bone Marrow","Colon","Heart","Muscle","Aorta","BAT","Cortex","Hippocampus"))
  color <- read.table("data/samples/20_distinct_color.txt")
  color <- setNames(color$V1,paste0("E",1:11))
  p<-ggplot(to_plot, aes(x = variable, y = value, fill = Var1)) +  
    geom_bar(stat = 'identity',color="white") +   
    theme_minimal() +   
    scale_fill_manual(values = color) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 20),legend.title = element_blank()) +
    ylab("Proportion")
  if(condition=="domain"){
    ggsave(paste0("result/all/diff/H3K27me3/recursion_domain_kmeans/chromHMM/RPKM_logFC_kmeans",kmean,"_chromHMM_annotation_young.png"),p,width = 20,height = 10,type="cairo",create.dir = T)
  }else{
    ggsave(paste0("result/all/diff/H3K27me3/peaks_merged/chromHMM/RPKM_logFC_kmeans",kmean,"_chromHMM_annotation_young.png"),p,width = 20,height = 10,type="cairo",create.dir = T)
  }
  
}

for(kmean in c(1:3)){
  summary_df <- data.frame()
  for(tissue in tissues){
    file_dir <- paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/")  
    files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_old[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
    file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
    chromHMM_old <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
    chromHMM_old <- chromHMM_old[which(chromHMM_old$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
    chromHMM_old$V2 <- chromHMM_old$V2+1
    chromHMM_old <- data.table(chromHMM_old)
    setDT(chromHMM_old) 
    setkey(chromHMM_old, V1, V2, V3) 
    
    overlaps_old <- foverlaps(kmeans_df[cluster == kmean] , chromHMM_old, type = "any", nomatch = 0L)  
    t_summary_df <- as.data.frame(table(overlaps_old$V4))
    t_summary_df$percent <- t_summary_df$Freq/sum(t_summary_df$Freq)*100
    colnames(t_summary_df)[3] <- tissue_label_change(tissue) 
    if(nrow(summary_df)==0){
      summary_df <- t_summary_df[,c(1,3)]
    }else{
      summary_df <- merge(summary_df,t_summary_df[,c(1,3)],by="Var1")
    }
  }
  to_plot<- reshape2::melt(summary_df)
  to_plot$Var1 <- factor(to_plot$Var1,levels=paste0("E",1:11))
  to_plot$variable <- factor(to_plot$variable,levels=c("Uterus","Ovary","Thymus","Liver","Bladder","Kidney","Pancreas","Stomach","Mammary Gland","Skin","Cecum","Ileum","Tongue","iWAT","Testis",
                                                       "Jejunum","Spleen","Cerebellum","Lung","Bone Marrow","Colon","Heart","Muscle","Aorta","BAT","Cortex","Hippocampus"))
  color <- read.table("data/samples/20_distinct_color.txt")
  color <- setNames(color$V1,paste0("E",1:11))
  p<-ggplot(to_plot, aes(x = variable, y = value, fill = Var1)) +  
    geom_bar(stat = 'identity',color="white") +   
    theme_minimal() +   
    scale_fill_manual(values = color) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 20),legend.title = element_blank()) +
    ylab("Proportion")
  if(condition=="domain"){
    ggsave(paste0("result/all/diff/H3K27me3/recursion_domain_kmeans/chromHMM/RPKM_logFC_kmeans",kmean,"_chromHMM_annotation_old.png"),p,width = 20,height = 10,type="cairo",create.dir = T)
  }else{
    ggsave(paste0("result/all/diff/H3K27me3/peaks_merged/chromHMM/RPKM_logFC_kmeans",kmean,"_chromHMM_annotation_old.png"),p,width = 20,height = 10,type="cairo",create.dir = T)
  }
  
}
