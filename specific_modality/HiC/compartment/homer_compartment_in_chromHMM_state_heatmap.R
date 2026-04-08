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
options(bitmapType = "cairo")  
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
state_num <- 14
resolution <- "50000"
compartment_in_chromHMM_state <- function(tissue,state_num,resolution){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue == tissue),]
  search_table$age <- factor(search_table$age, c("3M","24M"))
  search_table <- search_table[order(search_table$age),]
  
  file_dir <- paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/")  
  files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_young[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
  file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
  chromHMM_young <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
  chromHMM_young <- chromHMM_young[which(chromHMM_young$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  chromHMM_young$V2 <- chromHMM_young$V2+1
  chromHMM_young <- data.table(chromHMM_young)
  
  file_dir <- paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/")  
  files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_old[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
  file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
  chromHMM_old <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
  chromHMM_old <- chromHMM_old[which(chromHMM_old$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  chromHMM_old$V2 <- chromHMM_old$V2+1
  chromHMM_old <- data.table(chromHMM_old)
  
  for(i in c(1:length(search_table$sample_name))){
    sample <- search_table$sample_name[i]
    df <- read.table(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_",resolution,".PC1.txt"))
    df <- df[which(df$V2 %in% c(paste0("chr",c(1:19,"X","Y")))),]
    df <- df[,c(2:4,6)]
    df$compartment <- ifelse(df[, 4] > 0, "A", "B")
    df <- df[,c(1:3,5)]
    df$V3 <- df$V3+1
    setDT(df)
    setkey(df, V2, V3, V4) 
    if(search_table$age[i] == "3M"){
      chromHMM <- chromHMM_young
    }else{
      chromHMM <- chromHMM_old
    }
    setDT(chromHMM)  
    setkey(chromHMM, V1, V2, V3) 
    overlaps <- foverlaps(chromHMM,df, type = "any", nomatch = 0L)  
    A_compartment <- overlaps[which(overlaps$compartment=="A")]
    A_compartment <- as.data.frame(table(A_compartment$i.V4))
    A_compartment$percent <- A_compartment$Freq/sum(A_compartment$Freq)*100
    colnames(A_compartment)[3] <- "CompartmentA"
    
    B_compartment <- overlaps[which(overlaps$compartment=="B")]
    B_compartment <- as.data.frame(table(B_compartment$i.V4))
    B_compartment$percent <- B_compartment$Freq/sum(B_compartment$Freq)*100
    colnames(B_compartment)[3] <- "CompartmentB"
    
    compartment <- merge(A_compartment[,c(1,3)],B_compartment[,c(1,3)],by="Var1")
    compartment$Var1 <- factor(compartment$Var1,levels=c(paste0("E",1:state_num)))
    compartment <- compartment[order(compartment$Var1),]
    rownames(compartment) <- compartment$Var1
    compartment <- compartment[,-1]
    color_palette <- colorRampPalette(c("white", "blue"))(101) 
    pheatmap::pheatmap(compartment,cluster_rows = F,cluster_cols = F,breaks = seq(0, 100, length.out=101),color = color_palette,main = paste(tissue,sample,search_table$age[i]))
    }
}

compartment_change_in_chromHMM_state <- function(tissue,state_num,resolution){
  file_dir <- paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/")  
  files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_young[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
  file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
  chromHMM_young <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
  chromHMM_young <- chromHMM_young[which(chromHMM_young$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  chromHMM_young$V2 <- chromHMM_young$V2+1
  chromHMM_young <- data.table(chromHMM_young)
  chromHMM <- chromHMM_young
  setDT(chromHMM)  
  setkey(chromHMM, V1, V2, V3) 
  compartment <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/compartment_change_",resolution,".csv"))
  compartment <- compartment[which(compartment$condition %in% c("A-B","B-A")),]
  
  compartment <- as.data.table(compartment)
  setDT(compartment)  
  setkey(compartment, chr, start, end) 
  overlaps <- foverlaps(chromHMM, compartment, type = "any", nomatch = 0L)  
  A_B_compartment <- overlaps[which(overlaps$condition=="A-B")]
  A_B_compartment <- as.data.frame(table(A_B_compartment$V4))
  A_B_compartment$percent <- A_B_compartment$Freq/sum(A_B_compartment$Freq)*100
  colnames(A_B_compartment)[3] <- "A-B"
  
  B_A_compartment <- overlaps[which(overlaps$condition=="B-A")]
  B_A_compartment <- as.data.frame(table(B_A_compartment$V4))
  B_A_compartment$percent <- B_A_compartment$Freq/sum(B_A_compartment$Freq)*100
  colnames(B_A_compartment)[3] <- "B-A"
  
  compartment <- merge(A_B_compartment[,c(1,3)],B_A_compartment[,c(1,3)],by="Var1",all = T)
  compartment$Var1 <- factor(compartment$Var1,levels=c(paste0("E",1:state_num)))
  compartment <- compartment[order(compartment$Var1),]
  rownames(compartment) <- compartment$Var1
  compartment <- compartment[,-1]
  color_palette <- colorRampPalette(c("white", "blue"))(101) 
  pheatmap::pheatmap(compartment,cluster_rows = F,cluster_cols = F,breaks = seq(0, 100, length.out=101),filename = paste0("result/HiC/",tissue,"/compartment/",tissue,"_homer_compartment_change_in_chromHMM_state.png"),width = 5,height = 6,color = color_palette,main = paste(tissue_label_change(tissue),"homer compartment change"),display_numbers = T)
}
tissues <- c("kidney","colon","liver","lung","CB","brain")
for(tissue in tissues){
  compartment_change_in_chromHMM_state(tissue, state_num, resolution)
}
