rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(data.table)
library(tidyr)
tissue <- "lung"
resolution <- "50000"
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

contact_matrix_combine <- function(tissue, resolution){
  chromosomes <- paste0("chr",c(1:19,"X","Y"))
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  for(i in c(2:length(search_table$sample_name))){
    sample <- search_table$sample_name[i]
    compartment <- read.table(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_",resolution,".PC1.txt"))
    compartment <- compartment[which(compartment$V2 %in% c(paste0("chr",c(1:19,"X","Y")))),]
    compartment <- compartment[,c(1,6)]
    compartment$compartmeent <- ifelse(compartment[, 2] > 0, "A", "B")
    colnames(compartment)[3] <- "compartment"
    compartment <- as.data.table(compartment)
    setkey(compartment,V1)
    chr_list <- list()
    for(chr in chromosomes){
      df <- read.table(paste0("data/samples/HiC/",tissue,"/ob_ex_matrix/",sample,"/",sample,"_",resolution,"_",chr,"_ob_ex_Matrix.txt"),sep = "\t")
      colnames(df) <- df[1,]
      df <- df[-1,]
      rownames(df) <- df[,1]
      df <- df[,-c(1:2)]
      upper_tri_indices <- which(row(df) <= col(df), arr.ind=TRUE)  
      chr_df <- data.frame(  
        bin1 = rownames(df)[upper_tri_indices[, 1]],   
        bin2 = colnames(df)[upper_tri_indices[, 2]],  
        Values = df[upper_tri_indices]  
      )  
      chr_df <- chr_df %>%  
        separate(bin1, into = c("bin1.chr", "bin1.start"), sep = "-", remove = FALSE) %>%  
        separate(bin2, into = c("bin2.chr", "bin2.start"), sep = "-", remove = FALSE)  
      chr_df <- chr_df[which( abs(as.numeric(chr_df$bin2.start) - as.numeric(chr_df$bin1.start)) > 1000000 ),]
      chr_df <- chr_df[which(chr_df$bin1 %in% compartment$V1 & chr_df$bin2 %in% compartment$V1),]
      chr_df <- chr_df[,c(1,4,7)]
      chr_df <- as.data.table(chr_df)  
      chr_df <- chr_df[compartment, on = .(bin1 = V1), bin1.compartment := i.compartment]
      chr_df <- chr_df[compartment, on = .(bin2 = V1), bin2.compartment := i.compartment]
      chr_df <- as.data.frame(chr_df)
      chr_list[[chr]] <- chr_df
      }
    ob_ex <- Reduce(function(x, y) rbind(x, y),chr_list)
    write.csv(ob_ex,paste0("data/samples/HiC/",tissue,"/ob_ex_matrix/",sample,"/",sample,"_",resolution,"_ob_ex_Matrix.csv"))
  }
}

compartment_contact <- function(tissue,resolution){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  contact_list <- list()
  for(i in c(1:length(search_table$sample_name))){
    sample <- search_table$sample_name[i]
    df <- fread(paste0("data/samples/HiC/",tissue,"/ob_ex_matrix/",sample,"/",sample,"_",resolution,"_ob_ex_Matrix.csv"))
    df[, condition := paste0(bin1.compartment, "-", bin2.compartment)] 
    sum_by_compartment <- as.data.frame(df[, .(TotalValues = sum(Values, na.rm = TRUE)), by = condition])
    sum_by_compartment <- data.frame(condition=c("A-A","B-B","A-B"),value=c(sum_by_compartment[which(sum_by_compartment$condition=="A-A"),2],
                                                                              sum_by_compartment[which(sum_by_compartment$condition=="B-B"),2],
                                                                              sum_by_compartment[which(sum_by_compartment$condition=="A-B"),2]+sum_by_compartment[which(sum_by_compartment$condition=="B-A"),2]))
    sum_by_compartment$percent <- sum_by_compartment$value/sum(sum_by_compartment$value)*100
  
    colnames(sum_by_compartment)[3] <- sample  
    contact_list[[i]] <- sum_by_compartment[,c(1,3)]
    }
  to_plot <- Reduce(function(x, y) merge(x, y, by = "condition"), contact_list) 
  
  to_plot_melt <- reshape2::melt(to_plot)
  to_plot_melt$condition <- factor(to_plot_melt$condition,c("A-A","B-B","A-B"))
  to_plot_melt$position <- to_plot_melt$value
  to_plot_melt$position[which(to_plot_melt$condition=="A-A")] <- 100
  to_plot_melt$position[which(to_plot_melt$condition=="B-B")] <- 100- to_plot_melt$value[which(to_plot_melt$condition=="A-A")]
  ggplot(to_plot_melt, aes(x = variable, y = value, fill = condition)) +  
    geom_bar(width = 1, stat = "identity", color = "white") +
    theme_minimal()+
    labs(fill = NULL) + 
    xlab(NULL)+
    ylab("Proportion")+
    ggtitle(paste0(tissue_label_change(tissue))) +   
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right",
          plot.title = element_text(hjust = 0.5),
          text = element_text(size = 16))+
    geom_text(data =to_plot_melt,   
              aes(label = paste0(round(to_plot_melt$value,1),"%"), y = position),   
              color = "black", size = 5, vjust = 0.5)
  
  
  to_plot_heatmap <- to_plot
  to_plot_heatmap$condition <- factor(to_plot_heatmap$condition, levels=c("A-A","B-B","A-B"))
  rownames(to_plot_heatmap) <- to_plot$condition  
  to_plot_heatmap <- to_plot_heatmap[order(to_plot_heatmap$condition),]
  
  to_plot_heatmap <- to_plot_heatmap[,-1]
  
  annotation <- search_table[,c("sample_name","age")]
  annotation$age <- factor(annotation$age,c("3M","24M"))
  rownames(annotation) <- annotation$sample_name
  annotation <- annotation[,-1,drop=F]
  pheatmap::pheatmap(to_plot_heatmap,main = tissue_label_change(tissue),cluster_cols = F,cluster_rows=F, annotation_col = annotation, scale = "row")
  
  
  }


