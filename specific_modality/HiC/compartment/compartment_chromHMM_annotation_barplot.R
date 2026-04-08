rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggalluvial)  
library(data.table)
tissue <- "mammarygland"
state_num <- 11
resolution <- "50000"
check_row <- function(row) {  
  all(row[-1] == row[-1][1])  
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
    }else if(tissue_label == "Mammarygland"){
      tissue_label <- "Mammary Gland"
    }else if(tissue_label == "Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
} 
compartment_chromHMM <- function(tissue,state_num,resolution){
  ages <- c("young","old")
  df_list <- list(young=list(),old=list())
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  young_samples <- search_table$sample_name[which(search_table$age=="3M")]
  old_samples <- search_table$sample_name[which(search_table$age=="24M")]
  samples_list <- list(young=young_samples,old=old_samples)
  
  for(age in ages){
    for(i in c(1:length(samples_list[[age]]))){
      sample <- samples_list[[age]][i]
      if(sample=="WJH-Liver-103"){
        df <- read.table(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_",resolution,".PC1_recorrect.txt"))
      }else{
        df <- read.table(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_",resolution,".PC1.txt"))
      }
      df <- df[which(df$V2 %in% c(paste0("chr",c(1:19,"X","Y")))),]
      df <- df[,c(2:4,6)]
      df$compartmeent <- ifelse(df[, 4] > 0, "A", "B")
      colnames(df)[5] <- sample
      df <- df[,c(1:3,5)]
      df$label <- paste(df$V2,df$V3,df$V4,sep = "-")
      df_list[[age]][[i]] <- df
      names(df_list[[age]])[i]<-sample
    }
  }
  
  young_data <- Reduce(function(x, y) merge(x, y, by = "label"), df_list[["young"]])  
  old_data <- Reduce(function(x, y) merge(x, y, by = "label"), df_list[["old"]])   
  young_data <- young_data[,c(1:5,9)]
  old_data <- old_data[,c(1:5,9)]
  colnames(young_data)[2:4] <- c("Chr","Start","End")
  colnames(old_data)[2:4] <- c("Chr","Start","End")
  young_data$Start <- young_data$Start +1
  old_data$Start <- old_data$Start+1
  
  young_data <- young_data[young_data[,5]==young_data[,6], ] 
  old_data <- old_data[old_data[,5]==old_data[,6], ] 
  
  young_data <- young_data[,c(2:5)]
  old_data <- old_data[,c(2:5)]
  young_data <- as.data.table(young_data)
  old_data <- as.data.table(old_data)
  
  setDT(young_data)
  setDT(old_data)
  setkey(young_data,Chr,Start,End)
  setkey(old_data,Chr,Start,End)
  
  colnames(young_data)[4] <- "young_compartment"
  colnames(old_data)[4] <- "old_compartment"
  
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
  
  young_overlap <- foverlaps(young_data, chromHMM_young, type = "any", nomatch = 0L)  
  old_overlap <- foverlaps(old_data, chromHMM_old, type = "any", nomatch = 0L)  
  overlap_list <- list(young=young_overlap,old=old_overlap)
  p_list <- list()
  for(age in ages){
    df <- overlap_list[[age]]
    A_compartment <- df[which(df[,7] == "A"),]
    B_compartment <- df[which(df[,7] == "B"),]
    A_compartment <- as.data.frame(table(A_compartment$V4))
    A_compartment$percent <- A_compartment$Freq/sum(A_compartment$Freq) * 100
    B_compartment <- as.data.frame(table(B_compartment$V4))
    B_compartment$percent <- B_compartment$Freq/sum(B_compartment$Freq) * 100
    A_compartment$compartment <- "CompartmentA"
    B_compartment$compartment <- "CompartmentB"
    to_plot <- rbind(A_compartment,B_compartment)
    to_plot$Var1 <- factor(to_plot$Var1, levels=paste0("E",1:14))
    color <- read.table("data/samples/20_distinct_color.txt")
    color <- setNames(color$V1,paste0("E",1:11))
    p_list[[age]] <- ggplot(to_plot, aes(x = compartment, y = percent, fill = Var1)) +  
      geom_bar(stat = 'identity',color="white") +   
      theme_minimal() +   
      scale_fill_manual(values = color) +
      theme(axis.title.x = element_blank(), 
            axis.text.x = element_text(angle = 45, hjust = 1),
            text = element_text(size = 20),legend.title = element_blank()) +
      ylab("Proportion")+
      ggtitle(paste0(tissue_label_change(tissue)," ",age))
  
  }
  p <- p_list[["young"]] + p_list[["old"]]
  print(p)
  ggsave(paste0("result/HiC/",tissue,"/compartment/",tissue,"_homer_compartment_chromHMM_state_annotation.png"),p,width = 10,height = 6,type="cairo")
}
tissues <- c("brain","CB", "kidney", "liver", "lung", "bonemarrow", "colon", "heart", "Hip", "mammarygland", "stomach", "thymus")
for(tissue in tissues){
  compartment_chromHMM(tissue,state_num, resolution)
}

compartment_change_chromHMM_annotation <- function(tissues,state_num,resolution){
  summary <- data.frame()
  for(tissue in tissues){
    compartment_change <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/compartment_change_",resolution,".csv"))
    compartment_change_A_B <- compartment_change[which(compartment_change$condition=="A-B"),c(1:3)]
    compartment_change_A_B <- as.data.table(compartment_change_A_B)
    setDT(compartment_change_A_B)
    setkey(compartment_change_A_B,chr,start,end)
    compartment_change_B_A <- compartment_change[which(compartment_change$condition=="B-A"),c(1:3)]
    setDT(compartment_change_B_A)
    setkey(compartment_change_B_A,chr,start,end)
    }
  }
