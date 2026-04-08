rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)
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

compare_compartment_contact <- function(tissue,resolution){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  young_samples <- search_table$sample_name[which(search_table$age=="3M")]
  old_samples <- search_table$sample_name[which(search_table$age=="24M")]
  samples_list <- list(young=young_samples,old=old_samples)
  combinations <- as.data.frame(expand.grid(young = samples_list$young, old = samples_list$old)) 
  compartment_list <- list()
  for(i in c(1:length(search_table$sample_name))){
    sample <- search_table$sample_name[i]
    df <- read.table(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_",resolution,".PC1.txt"))
    df <- df[which(df$V2 %in% c(paste0("chr",c(1:19,"X","Y")))),]
    df <- df[,c(2:4,6)]
    df$compartment <- ifelse(df[, 4] > 0, "A", "B")
    compartment_list[[i]] <- df
    names(compartment_list)[i]<-sample
  }
  
  contact_list <- list()
  for(i in c(1:length(search_table$sample_name))){
    sample <- search_table$sample_name[i]
    bed <- read.table(paste0("data/samples/HiC/",tissue,"/raw_matrix/",sample,"_",resolution,"_abs.bed"))
    contact <- fread(paste0("data/samples/HiC/",tissue,"/ice_matrix/",sample,"_",resolution,"_iced.matrix"))
    contact <- contact[abs(V1 - V2) > 20] 
    bed$label <- paste0(bed$V1,"-",bed$V2,"-",bed$V3)
    compartment_list[[sample]]$label <- paste0(compartment_list[[sample]]$V2,"-",compartment_list[[sample]]$V3,"-",compartment_list[[sample]]$V4)
    bed <- merge(bed[,c("V1","label","V4")],compartment_list[[sample]][,c("label","compartment")],by="label")
    bed <- as.data.table(bed)
    setkey(bed,V4)
    contact <- contact[bed, on = .(V1 = V4), V1.compartment := i.compartment]
    contact <- contact[bed, on = .(V2 = V4), V2.compartment := i.compartment]
    contact <- contact[!is.na(V1.compartment) & !is.na(V2.compartment)]  
    contact <- contact[bed, on = .(V1 = V4), V1.chr := i.V1]
    contact <- contact[bed, on = .(V2 = V4), V2.chr := i.V1]
    contact <- contact[V1.chr==V2.chr] 
    contact[, compartment_pair := paste0(V1.compartment, "-", V2.compartment)]  
    contact[, contact_pair := paste0(V1, "-", V2)] 
    contact <- as.data.frame(contact)
    contact <- contact[,c(3,8,9)]
    colnames(contact)[1] <- sample
    contact_list[[i]] <- contact
    names(contact_list)[i] <- sample 
  }
  compare_list <- list()
  for(i in c(1:nrow(combinations))){
    young <- as.character(combinations[i,"young"])
    old <- as.character(combinations[i,"old"])
    df <- merge(contact_list[[young]], contact_list[[old]], by="contact_pair")
    df <- df[which(df$compartment_pair.x == df$compartment_pair.y),]
    df <- df[,c(1,2,4,5)]
    df$log2FC <- log2(df[,3]/df[,2])
    df <- df[,c(1,4,5)]
    colnames(df)[2] <- "compartment_pair"
    compare_list[[paste0(old,"_",young)]] <- df
  }
  # saveRDS(compare_list,paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/compartment_contact_compare_",resolution,".rds"))
  to_plot <- list()
  for(i in c(1:length(compare_list))){
    df <- compare_list[[i]]
    df <- df %>%  
      mutate(compartment_pair = ifelse(compartment_pair == "B-A", "A-B", compartment_pair)) 
    result <- df %>%  
      group_by(compartment_pair) %>%  
      summarise(mean_log2FC = mean(log2FC, na.rm = TRUE))
    colnames(result)[2] <- names(compare_list)[i]
    to_plot[[i]] <- result
  }
  to_plot <- Reduce(function(x, y) merge(x, y, by = "compartment_pair"), to_plot)
  rownames(to_plot) <- to_plot$compartment_pair
  to_plot$compartment_pair <- factor(to_plot$compartment_pair,levels=c("A-A","B-B","A-B"))
  to_plot <- to_plot[sort(to_plot$compartment_pair),]
  to_plot <- to_plot[,-1]
  pheatmap::pheatmap(to_plot,main = tissue_label_change(tissue),cluster_cols = F,cluster_rows=F,breaks = seq(-0.5, 0.5, length.out=101),display_numbers = T,fontsize_number = 15)
}

compare_compartment_contact_ob_ex <- function(tissue,resolution){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  young_samples <- search_table$sample_name[which(search_table$age=="3M")]
  old_samples <- search_table$sample_name[which(search_table$age=="24M")]
  samples_list <- list(young=young_samples,old=old_samples)
  combinations <- as.data.frame(expand.grid(young = samples_list$young, old = samples_list$old)) 
  compartment_list <- list()
  
  contact_list <- list()
  for(i in c(1:length(search_table$sample_name))){
    sample <- search_table$sample_name[i]
    df <- fread(paste0("data/samples/HiC/",tissue,"/ob_ex_matrix/",sample,"/",sample,"_",resolution,"_ob_ex_Matrix.csv"))
    df[, compartment_pair := paste0(bin1.compartment, "-", bin2.compartment)] 
    df[, contact_pair := paste0(bin1,"_",bin2)]
    df <- df[, .(contact_pair, compartment_pair, Values)]
    df <- as.data.frame(df)
    contact_list[[i]] <- df
    names(contact_list)[i] <- sample 
  }
  
  compare_list <- list()
  for(i in c(1:nrow(combinations))){
    young <- as.character(combinations[i,"young"])
    old <- as.character(combinations[i,"old"])
    df <- merge(contact_list[[young]], contact_list[[old]], by="contact_pair")
    df <- df[which(df$compartment_pair.x == df$compartment_pair.y),]
    df <- df[,c(1,3,4,5)]
    df$log2FC <- log2(df[,4]/df[,2])
    df <- df[,c(1,3,5)]
    colnames(df)[2] <- "compartment_pair"
    compare_list[[paste0(old,"_",young)]] <- df
  }
  # saveRDS(compare_list,paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/compartment_contact_compare_ob_ex_",resolution,".rds"))
  # compare_list <-readRDS(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/compartment_contact_compare_ob_ex_",resolution,".rds"))
  to_plot <- list()
  for(i in c(1:length(compare_list))){
    df <- compare_list[[i]]
    df <- df %>%  
      mutate(compartment_pair = ifelse(compartment_pair == "B-A", "A-B", compartment_pair)) 
    df <- as.data.frame(df)
    df <- df %>%  
      filter(!is.na(log2FC) & is.finite(log2FC))
    result <- df %>%  
      group_by(compartment_pair) %>%  
      summarise(mean_log2FC = mean(log2FC, na.rm = TRUE))
    colnames(result)[2] <- names(compare_list)[i]
    to_plot[[i]] <- result
  }
  to_plot <- Reduce(function(x, y) merge(x, y, by = "compartment_pair"), to_plot)
  rownames(to_plot) <- to_plot$compartment_pair
  to_plot$compartment_pair <- factor(to_plot$compartment_pair,levels=c("A-A","B-B","A-B"))
  to_plot <- to_plot[sort(to_plot$compartment_pair),]
  to_plot <- to_plot[,-1]
  pheatmap::pheatmap(to_plot,main = tissue_label_change(tissue),cluster_cols = F,cluster_rows=F,breaks = seq(-0.3, 0.3, length.out=101),display_numbers = T,fontsize_number = 15)
  }

