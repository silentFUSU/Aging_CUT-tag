rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(edgeR)
library(ggalluvial)  
library(tools)
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}

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
state_num <- 15
age="young"
tissues <-  sort(c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
                   "thymus","skin","bladder","bonemarrow","Hip","heart",
                   "muscle","jejunum","uterus","ovary","liver","tongue",
                   "cecum","colon","testis","stomach","pancreas","iWAT","ileum"))
age_percentage_summary <- data.frame()
tissue_summary_list <- list()
for(age in c("young","old")){
  tissue_summary <- data.frame()
  for(tissue in tissues){
    peaks <- read.table(paste0("data/samples/ATAC/",tissue,"/ATAC/bed/ATAC_macs_young_old_narrowpeak_summits_spm3.bed"))
    peaks <- peaks[,c(1:3)]
    # peaks$V2 <- peaks$V2 + 1
    peaks <- as.data.table(peaks)
    setDT(peaks)
    setkey(peaks,V1,V2,V3)
    
    file_dir <- paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/split_1k/")  
    files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_",age,"[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
    file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
    chromHMM <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
    chromHMM <- chromHMM[which(chromHMM$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
    if(tissue %in% c("mammarygland","ovary","uterus")){
      chromHMM <- chromHMM[which(chromHMM$V1 %in% paste0("chr",c(1:19,"X"))),]
    }
    chromHMM$V2 <- chromHMM$V2+1
    chromHMM <- data.table(chromHMM)
    
    setDT(chromHMM) 
    setkey(chromHMM, V1, V2, V3) 
    overlaps<- foverlaps(peaks, chromHMM, type = "any", nomatch = 0L)  
    
    overlaps$label <- paste(overlaps$V1,overlaps$V2,overlaps$V3,sep = "-")
    overlaps <- as.data.frame(overlaps)
    result <- overlaps %>%  
      group_by(V4) %>%  
      summarize(freq = n())  
    
    result$percent <- result$freq / sum(result$freq) *100
    result <- result[,c(1,3)]
    colnames(result)[2] <- tissue_label_change(tissue)
    
    if(nrow(tissue_summary)==0){
      tissue_summary <- result
    }else{
      tissue_summary <- merge(tissue_summary,result,by="V4",all=T)
    }
  }
  tissue_summary_list[[age]] <- tissue_summary
  t_age_summary <- reshape2::melt(tissue_summary)
  t_age_summary <- t_age_summary %>%
    group_by(V4) %>%
    summarise(mean_percent = mean(value, na.rm = TRUE))
  t_age_summary$age <- age
  age_percentage_summary <- rbind(age_percentage_summary,t_age_summary)
  }

to_plot <- reshape2::melt(tissue_summary_list[["old"]])
to_plot$V4 <- factor(to_plot$V4,levels = paste0("E",1:15))
to_plot$variable <- factor(to_plot$variable,levels=sort(unique(to_plot$variable)))
color <- read.table("data/samples/20_distinct_color.txt")
color <- setNames(color$V1,paste0("E",1:15))
ggplot(to_plot, aes(x = variable, y = value, fill = V4)) +  
  geom_bar(stat = 'identity',color="white") +   
  theme_minimal() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Proportion")

to_plot <- as.data.frame(age_percentage_summary)
dictionary <- list("E1"=1, "E2"=2, "E3"=3,
                   "E4"=4, "E5"=5, "E6"=7,
                   "E7"=8, "E8"=6, "E9"=9,
                   "E10"=10,"E11"=11,"E12"=15,
                   "E13"=12,"E14"=13,"E15"=14)
keys <- names(dictionary)
values <- unlist(dictionary)
to_plot$V4 <- values[match(to_plot$V4, keys)]
to_plot$V4 <- paste0("E",to_plot$V4)
to_plot$V4 <- factor(to_plot$V4, levels = paste0("E",1:15))
to_plot$age <- factor(to_plot$age,levels=c("young","old"))

p <- ggplot(to_plot, aes(x = age, y = mean_percent, fill = V4)) +
  geom_bar(stat = 'identity',color="white") +
  theme_bw() +
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Proportion")
ggsave("result/Sup_figures/ATAC_peak_chromHMM_annotation.pdf",p,height = 6,width = 4)


conditions <- c("up","down")
tissue_condition_summary <- list(up=list(young=data.frame(),old=data.frame()),down=list(young=data.frame(),old=data.frame()))
for(condition in conditions){
  for(age in c("young","old")){
    for(tissue in tissues){
      df <- read.csv(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_macs_young_old_narrowpeak_summits_spm3_diff_after_remove_batch_effect.csv"))
      if(condition == "up"){
        peaks <- df[which(df$Significant=="Up"),c("Chr","Start","End")]
      }else{
        peaks <- df[which(df$Significant=="Down"),c("Chr","Start","End")]
      }
      if(nrow(peaks) >= 20){
        colnames(peaks) <- c("V1","V2","V3")
        peaks <- as.data.table(peaks)
        setDT(peaks)
        setkey(peaks,V1,V2,V3)
        file_dir <- paste0("result/all/ChromHMM/all_tissues_previous/",state_num,"_all_tissues/split_1k/")  
        files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_",age,"[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
        file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
        chromHMM <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
        chromHMM <- chromHMM[which(chromHMM$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
        if(tissue %in% c("mammarygland","ovary","uterus")){
          chromHMM <- chromHMM[which(chromHMM$V1 %in% paste0("chr",c(1:19,"X"))),]
        }
        chromHMM$V2 <- chromHMM$V2+1
        chromHMM <- data.table(chromHMM)
        
        setDT(chromHMM) 
        setkey(chromHMM, V1, V2, V3) 
        overlaps<- foverlaps(peaks, chromHMM, type = "any", nomatch = 0L)  
        
        overlaps$label <- paste(overlaps$V1,overlaps$V2,overlaps$V3,sep = "-")
        overlaps <- as.data.frame(overlaps)
        result <- overlaps %>%  
          group_by(V4) %>%  
          summarize(freq = n())  
        
        result$percent <- result$freq / sum(result$freq) *100
        result <- result[,c(1,3)]
        colnames(result)[2] <- tissue_label_change(tissue)
        if(nrow(tissue_condition_summary[[condition]][[age]]) > 0){
          tissue_condition_summary[[condition]][[age]] <- merge(tissue_condition_summary[[condition]][[age]],result,by="V4",all=T)
        }else{
          tissue_condition_summary[[condition]][[age]] <- result
        }
      }
    }
  }
}

for(condition in conditions){
  for(age in c("young","old")){
    to_plot <- tissue_condition_summary[[condition]][[age]]
    to_plot[is.na(to_plot)] <- 0
    to_plot <- reshape2::melt(to_plot)
    to_plot$V4 <- factor(to_plot$V4,levels = paste0("E",1:15))
    color <- read.table("data/samples/20_distinct_color.txt")
    color <- setNames(color$V1,paste0("E",1:15))
    ggplot(to_plot, aes(x = variable, y = value, fill = V4)) +  
      geom_bar(stat = 'identity',color="white") +   
      theme_minimal() +   
      scale_fill_manual(values = color) +
      theme(axis.title.x = element_blank(), 
            axis.text.x = element_text(angle = 45, hjust = 1),
            text = element_text(size = 20),legend.title = element_blank()) +
      ggtitle(paste0(toTitleCase(condition)," ",age))+
      ylab("Proportion")
  }
}

for(condition in conditions){
  to_plot <- data.frame()
  for(age in c("young","old")){
    t_to_plot <- tissue_condition_summary[[condition]][[age]]
    t_to_plot[is.na(t_to_plot)] <- 0
    t_to_plot <- reshape2::melt(t_to_plot)
    t_to_plot <- t_to_plot %>%
      group_by(V4) %>%
      summarise(mean_percent = mean(value, na.rm = TRUE))
    t_to_plot$age <- age
    to_plot <- rbind(to_plot,t_to_plot)
  }
}
to_plot$V4 <- factor(to_plot$V4, levels = paste0("E",1:15))
to_plot$age <- factor(to_plot$age,levels = c("young","old"))
ggplot(to_plot, aes(x = age, y = mean_percent, fill = V4)) +
  geom_bar(stat = 'identity',color="white") +
  theme_minimal() +
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ggtitle(toTitleCase(condition))+
  ylab("Proportion")
