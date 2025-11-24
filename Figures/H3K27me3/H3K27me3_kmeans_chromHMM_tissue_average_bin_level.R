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
conditions <- c("Up","Down")
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver","ileum",
                  "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT"))
for(condition in conditions){
  summary_df_young <- data.frame()
  summary_df_old <- data.frame()
  for(tissue in tissues){
    df <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_10kb_bins_diff_after_remove_batch_effect.csv"))
    df <- df[which(df$Significant==condition),c("Chr","Start","End","Significant")]
    df$Start <- df$Start +1 
    df <- as.data.table(df)
    setDT(df)
    setkey(df,Chr,Start,End)
    
    file_dir <- paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/split_1k/")  
    files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_young[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
    file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
    chromHMM_young <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
    chromHMM_young <- chromHMM_young[which(chromHMM_young$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
    chromHMM_young$V2 <- chromHMM_young$V2+1
    chromHMM_young <- data.table(chromHMM_young)
    setDT(chromHMM_young) 
    setkey(chromHMM_young, V1, V2, V3) 
    
    overlaps_young <- foverlaps(df, chromHMM_young, type = "any", nomatch = 0L)  
    t_summary_df_young <- as.data.frame(table(overlaps_young$V4))
    t_summary_df_young$percent <- t_summary_df_young$Freq/sum(t_summary_df_young$Freq)*100
    colnames(t_summary_df_young)[3] <- tissue
    
    
    file_dir <- paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/split_1k/")  
    files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_old[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
    file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
    chromHMM_old <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
    chromHMM_old <- chromHMM_old[which(chromHMM_old$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
    chromHMM_old$V2 <- chromHMM_old$V2+1
    chromHMM_old <- data.table(chromHMM_old)
    setDT(chromHMM_old) 
    setkey(chromHMM_old, V1, V2, V3) 
    
    overlaps_old <- foverlaps(df, chromHMM_old, type = "any", nomatch = 0L)  
    t_summary_df_old <- as.data.frame(table(overlaps_old$V4))
    t_summary_df_old$percent <- t_summary_df_old$Freq/sum(t_summary_df_old$Freq)*100
    colnames(t_summary_df_old)[3] <- tissue
    
    if(nrow(summary_df_young)==0){
      summary_df_young <- t_summary_df_young[,c(1,3)]
    }else{
      summary_df_young <- merge(summary_df_young,t_summary_df_young[,c(1,3)],by="Var1")
    }
    
    if(nrow(summary_df_old)==0){
      summary_df_old <- t_summary_df_old[,c(1,3)]
    }else{
      summary_df_old <- merge(summary_df_old,t_summary_df_old[,c(1,3)],by="Var1")
    }
  }
  
  summary_df_young$mean <- rowMeans(summary_df_young[,-1])
  summary_df_young$age <- "young"
  
  summary_df_old$mean <- rowMeans(summary_df_old[,-1])
  summary_df_old$age <- "old"
  
  summary_df_young <- summary_df_young[,c("Var1","mean","age")]
  summary_df_old <- summary_df_old[,c("Var1","mean","age")]
  to_plot <- rbind(summary_df_young,summary_df_old)
  write.csv(to_plot,paste0("tmp_H3K27me3_",condition,"_chromHMM.csv"))
}
condition <- "Up"
to_plot <- read.csv(paste0("tmp_H3K27me3_",condition,"_chromHMM.csv"),row.names = 1)
dictionary <- list("E1"=1, "E2"=2, "E3"=3,
                   "E4"=4, "E5"=5, "E6"=7,
                   "E7"=8, "E8"=6, "E9"=9,
                   "E10"=10,"E11"=11,"E12"=15,
                   "E13"=12,"E14"=13,"E15"=14)
keys <- names(dictionary)
values <- unlist(dictionary)
to_plot$Var1 <- values[match(to_plot$Var1, keys)]
to_plot$Var1 <- paste0("E",to_plot$Var1)
to_plot$Var1 <- sub("^E", "state", to_plot$Var1)
to_plot$Var1 <- factor(to_plot$Var1,levels=paste0("state",1:state_num))
color <- read.table("data/samples/20_distinct_color.txt")
color <- setNames(color$V1,paste0("state",1:state_num))
to_plot$age <- factor(to_plot$age,levels=c("young","old"))
p <- ggplot(to_plot, aes(x = age, y = mean, fill = Var1)) +  
  geom_bar(stat = 'identity',color="black") +   
  theme_minimal() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Proportion")
ggsave("result/Sup_figures/H3K27me3_up_bin_chromHMM_state.pdf",p,width = 6,height = 8)

