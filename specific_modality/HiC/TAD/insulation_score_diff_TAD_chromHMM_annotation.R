rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggalluvial)  
library(data.table)
library(tidyverse)
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
    }else if(tissue_label == "Mammarygland"){
      tissue_label <- "Mammary Gland"
    }else if(tissue_label == "Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
} 
tissues <- c("brain","CB","kidney","liver","lung","bonemarrow","colon","heart","Hip","mammarygland","stomach","thymus")
state_num <- 11
for(tissue in tissues){
  out <- read.csv(paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",tissue,"_redundant_20000_TAD_diff.csv"))
  out <- out %>%  
    separate(X, into = c("chr", "start", "end"), sep = "-", convert = TRUE)  
  up <- out[which(out$Significant=="Up"),c("chr","start","end")]
  down <- out[which(out$Significant=="Down"),c("chr","start","end")]
  
  file_dir <- paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/")  
  files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_young[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
  file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
  chromHMM_young <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
  chromHMM_young <- chromHMM_young[which(chromHMM_young$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  chromHMM_young$V2 <- chromHMM_young$V2+1
  chromHMM_young <- data.table(chromHMM_young)
  setDT(chromHMM_young) 
  setkey(chromHMM_young, V1, V2, V3) 
  
  file_dir <- paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/")  
  files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_old[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
  file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
  chromHMM_old <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
  chromHMM_old <- chromHMM_old[which(chromHMM_old$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  chromHMM_old$V2 <- chromHMM_old$V2+1
  chromHMM_old <- data.table(chromHMM_old)
  setDT(chromHMM_old) 
  setkey(chromHMM_old, V1, V2, V3) 
  
  if(nrow(up) > 20){
    up <- as.data.table(up)
    setDT(up)
    setkey(up,chr,start,end)
    overlaps_young <- foverlaps(up, chromHMM_young, type = "any", nomatch = 0L)  
    overlaps_old <- foverlaps(up, chromHMM_old, type = "any", nomatch = 0L)  
    
    overlaps_young$label <- paste(overlaps_young$chr,overlaps_young$V2,overlaps_young$V3,sep = "-")
    overlaps_old$label <- paste(overlaps_old$chr,overlaps_old$V2,overlaps_old$V3,sep = "-")
    overlaps <- merge(overlaps_young[,c("V4","label")],overlaps_old[,c("V4","label")],by="label")
    overlaps <- as.data.frame(overlaps)
    colnames(overlaps)[2:3] <- c("Young_state","Old_state") 
    to_plot <- overlaps %>%  
      group_by(Young_state, Old_state) %>%  
      summarize(freq = n())  
    to_plot$Young_state <- factor(to_plot$Young_state,levels=paste0("E",1:state_num))
    to_plot$Old_state <- factor(to_plot$Old_state,levels=paste0("E",1:state_num))
    color <- read.table("data/samples/20_distinct_color.txt")
    color <- setNames(color$V1,paste0("E",1:11))
    ggplot(to_plot, aes(axis1 = Young_state, axis2 = Old_state, y = freq)) +  
      geom_alluvium(aes(fill = Young_state)) +  
      geom_stratum() +  
      geom_text(stat = "stratum", aes(label = after_stat(stratum))) +  
      theme_minimal() +  
      theme(legend.position = "none")+
      labs(y = "Count", x = "State Transition", 
           fill = "Young State")+
      ggtitle(paste0(tissue_label_change(tissue)," increased TAD regions"),"Sankey Plot of State Transitions")
  }
  if(nrow(down) > 20){
    down <- as.data.table(down)
    setDT(down)
    setkey(down,chr,start,end)
    overlaps_young <- foverlaps(down, chromHMM_young, type = "any", nomatch = 0L)  
    overlaps_old <- foverlaps(down, chromHMM_old, type = "any", nomatch = 0L)  
    
    overlaps_young$label <- paste(overlaps_young$chr,overlaps_young$V2,overlaps_young$V3,sep = "-")
    overlaps_old$label <- paste(overlaps_old$chr,overlaps_old$V2,overlaps_old$V3,sep = "-")
    overlaps <- merge(overlaps_young[,c("V4","label")],overlaps_old[,c("V4","label")],by="label")
    overlaps <- as.data.frame(overlaps)
    colnames(overlaps)[2:3] <- c("Young_state","Old_state") 
    to_plot <- overlaps %>%  
      group_by(Young_state, Old_state) %>%  
      summarize(freq = n())  
    to_plot$Young_state <- factor(to_plot$Young_state,levels=paste0("E",1:state_num))
    to_plot$Old_state <- factor(to_plot$Old_state,levels=paste0("E",1:state_num))
    color <- read.table("data/samples/20_distinct_color.txt")
    color <- setNames(color$V1,paste0("E",1:11))
    ggplot(to_plot, aes(axis1 = Young_state, axis2 = Old_state, y = freq)) +  
      geom_alluvium(aes(fill = Young_state)) +  
      geom_stratum() +  
      geom_text(stat = "stratum", aes(label = after_stat(stratum))) +  
      theme_minimal() +  
      theme(legend.position = "none")+
      labs(y = "Count", x = "State Transition", 
           fill = "Young State")+
      
      ggtitle(paste0(tissue_label_change(tissue)," decreased TAD regions"),"Sankey Plot of State Transitions")
  }
}

tissues <- c("brain","CB","kidney","liver","lung","bonemarrow","colon","heart","Hip","mammarygland","stomach","thymus")
state_num <- 11
up_summary <- data.frame()
down_summary <- data.frame()
for(tissue in tissues){
  out <- read.csv(paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",tissue,"_redundant_20000_TAD_diff.csv"))
  out <- out %>%  
    separate(X, into = c("chr", "start", "end"), sep = "-", convert = TRUE)  
  up <- out[which(out$Significant=="Up"),c("chr","start","end")]
  down <- out[which(out$Significant=="Down"),c("chr","start","end")]
  
  file_dir <- paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/")  
  files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_young[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
  file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
  chromHMM_young <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
  chromHMM_young <- chromHMM_young[which(chromHMM_young$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  chromHMM_young$V2 <- chromHMM_young$V2+1
  chromHMM_young <- data.table(chromHMM_young)
  setDT(chromHMM_young) 
  setkey(chromHMM_young, V1, V2, V3) 
  if(nrow(up) > 20){
    up <- as.data.table(up)
    setDT(up)
    setkey(up,chr,start,end)
    overlaps_young <- foverlaps(up, chromHMM_young, type = "any", nomatch = 0L)  
    overlaps_young_summary <- as.data.frame(table(overlaps_young$V4))
    overlaps_young_summary$percentage <- overlaps_young_summary$Freq/sum(overlaps_young_summary$Freq) * 100
    overlaps_young_summary$tissue <- tissue_label_change(tissue)
    overlaps_young_summary <- overlaps_young_summary[,c("tissue","Var1","percentage")]  
    up_summary <- rbind(up_summary,overlaps_young_summary)
    }
  if(nrow(down) > 20){
    down <- as.data.table(down)
    setDT(down)
    setkey(down,chr,start,end)
    overlaps_young <- foverlaps(down, chromHMM_young, type = "any", nomatch = 0L)  
    overlaps_young_summary <- as.data.frame(table(overlaps_young$V4))
    overlaps_young_summary$percentage <- overlaps_young_summary$Freq/sum(overlaps_young_summary$Freq) * 100
    overlaps_young_summary$tissue <- tissue_label_change(tissue)
    overlaps_young_summary <- overlaps_young_summary[,c("tissue","Var1","percentage")]  
    down_summary <- rbind(down_summary,overlaps_young_summary)
  }
}

up_summary$Var1 <- factor(up_summary$Var1,levels=paste0("E",c(1:state_num)))
down_summary$Var1 <- factor(down_summary$Var1,levels=paste0("E",c(1:state_num)))

color <- read.table("data/samples/20_distinct_color.txt")
color <- setNames(color$V1,paste0("E",1:state_num))

ggplot(up_summary, aes(x = tissue, y = percentage, fill = Var1)) +  
  geom_bar(stat = 'identity',color="white") +   
  theme_minimal() +   
  ggtitle("Increase TAD")+
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Proportion")

ggplot(down_summary, aes(x = tissue, y = percentage, fill = Var1)) +  
  geom_bar(stat = 'identity',color="white") +   
  theme_minimal() +   
  ggtitle("Decrease TAD")+
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Proportion")

##### chromHMM in old state
tissues <- c("brain","CB","kidney","liver","lung","bonemarrow","colon","heart","Hip","mammarygland","stomach","thymus")
state_num <- 11
up_summary <- data.frame()
down_summary <- data.frame()
for(tissue in tissues){
  out <- read.csv(paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",tissue,"_redundant_20000_TAD_diff.csv"))
  out <- out %>%  
    separate(X, into = c("chr", "start", "end"), sep = "-", convert = TRUE)  
  up <- out[which(out$Significant=="Up"),c("chr","start","end")]
  down <- out[which(out$Significant=="Down"),c("chr","start","end")]
  
  file_dir <- paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/")  
  files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_old[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
  file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
  chromHMM_old <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
  chromHMM_old <- chromHMM_old[which(chromHMM_old$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  chromHMM_old$V2 <- chromHMM_old$V2+1
  chromHMM_old <- data.table(chromHMM_old)
  setDT(chromHMM_old) 
  setkey(chromHMM_old, V1, V2, V3) 
  if(nrow(up) > 20){
    up <- as.data.table(up)
    setDT(up)
    setkey(up,chr,start,end)
    overlaps_old <- foverlaps(up, chromHMM_old, type = "any", nomatch = 0L)  
    overlaps_old_summary <- as.data.frame(table(overlaps_old$V4))
    overlaps_old_summary$percentage <- overlaps_old_summary$Freq/sum(overlaps_old_summary$Freq) * 100
    overlaps_old_summary$tissue <- tissue_label_change(tissue)
    overlaps_old_summary <- overlaps_old_summary[,c("tissue","Var1","percentage")]  
    up_summary <- rbind(up_summary,overlaps_old_summary)
  }
  if(nrow(down) > 20){
    down <- as.data.table(down)
    setDT(down)
    setkey(down,chr,start,end)
    overlaps_old <- foverlaps(down, chromHMM_old, type = "any", nomatch = 0L)  
    overlaps_old_summary <- as.data.frame(table(overlaps_old$V4))
    overlaps_old_summary$percentage <- overlaps_old_summary$Freq/sum(overlaps_old_summary$Freq) * 100
    overlaps_old_summary$tissue <- tissue_label_change(tissue)
    overlaps_old_summary <- overlaps_old_summary[,c("tissue","Var1","percentage")]  
    down_summary <- rbind(down_summary,overlaps_old_summary)
  }
}

up_summary$Var1 <- factor(up_summary$Var1,levels=paste0("E",c(1:state_num)))
down_summary$Var1 <- factor(down_summary$Var1,levels=paste0("E",c(1:state_num)))

color <- read.table("data/samples/20_distinct_color.txt")
color <- setNames(color$V1,paste0("E",1:state_num))

ggplot(up_summary, aes(x = tissue, y = percentage, fill = Var1)) +  
  geom_bar(stat = 'identity',color="white") +   
  theme_minimal() +   
  ggtitle("Increase TAD")+
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Proportion")

ggplot(down_summary, aes(x = tissue, y = percentage, fill = Var1)) +  
  geom_bar(stat = 'identity',color="white") +   
  theme_minimal() +   
  ggtitle("Decrease TAD")+
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Proportion")

