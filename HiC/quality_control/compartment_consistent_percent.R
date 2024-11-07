rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(ggrepel)
library(ggVennDiagram)  
library(ggplot2)  

# 关闭全局科学计数法显示  
options(scipen = 999)  
tissue <- "lung"


compartment_consistent_percent <- function(tissue){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  samples <- search_table$sample_name
  df_list <- list()
  for(i in c(1:length(samples))){
    df <- read.table(paste0("data/samples/HiC/lung/fanc_compartment/first_eigenvector/fanc_",samples[i],"_1mb.compartment.bed"))
    df <- df[which(df$V1 %in% c(paste0("chr",c(1:19,"X","Y")))),]
    df_new <- data_frame()
    for(j in c(1:nrow(df))){
      if(df[j,3]-df[j,2]+1 > 1000000){
          V1=df[j,1]  
          V4=df[j,4]
          V5=df[j,5]
          V6=df[j,6]
          for(start in seq(df[j,2],df[j,3], by = 1000000)){
            t_df <- data.frame(V1=V1,V2=start,V3=min(start+999999,df[j,3]),V4=V4,V5=V5,V6=V6)
            df_new <- rbind(df_new,t_df)
          }
      }else{
          df_new <- rbind(df_new,df[j,])
        }
    }
    df_list[[i]] <- df_new
    names(df_list)[i] <- samples[i]
    colnames(df_list[[i]])[5] <- samples[i]
    df_list[[i]]$label <- paste0(df_list[[i]]$V1,"-",df_list[[i]]$V2,"-",df_list[[i]]$V3)
    df_list[[i]]$state <- paste0(df_list[[i]]$V1,"-",df_list[[i]]$V2,"-",df_list[[i]]$V3,"-",df_list[[i]]$V4)
  }
  if(tissue == "lung"){
    df_list[["WJH-106-Lung"]][which(df_list[["WJH-106-Lung"]]$V1=="chr19"),5] <- -df_list[["WJH-106-Lung"]][which(df_list[["WJH-106-Lung"]]$V1=="chr19"),5]
    df_list[["WJH-106-Lung"]] <- df_list[["WJH-106-Lung"]] %>%  
      mutate(V4 = ifelse(V1 == "chr19",  
                         case_when(  
                           V4 == "A" ~ "B",  
                           V4 == "B" ~ "A",  
                           TRUE ~ V4  
                         ),  
                         V4)) 
    df_list[["WJH-106-Lung"]]$state <-  paste0(df_list[["WJH-106-Lung"]]$V1,"-",df_list[["WJH-106-Lung"]]$V2,"-",df_list[["WJH-106-Lung"]]$V3,"-",df_list[["WJH-106-Lung"]]$V4)
  }
  extracted_columns <- lapply(df_list, function(df) df[, c(5, 7)])  
  merged_data <- Reduce(function(x, y) merge(x, y, by = "label", all = TRUE), extracted_columns)  
  merged_data <- na.omit(merged_data)
  rownames(search_table) <- search_table$sample_name
  annotation <- search_table[,c("age"),drop=F]
  annotation$age <- factor(annotation$age,levels=c("3M","24M"))
  pheatmap::pheatmap(cor(merged_data[,-1]),annotation_row = annotation,breaks = seq(0, 1, length.out=101))
  young <- search_table$sample_name[which(search_table$age=="3M")]
  old <- search_table$sample_name[which(search_table$age=="24M")]
  
  young_table <- df_list[young]
  old_table <- df_list[old]
  young_venn_data <- lapply(young_table, function(x) x$state)
  
  ggVennDiagram(young_venn_data,   
                label = "count",   
                category.names = young) +  
    scale_fill_gradient(low = "white", high = "#0073C2FF")

  old_venn_data <- lapply(old_table, function(x) x$state)
  
  ggVennDiagram(old_venn_data,   
                label = "count",   
                category.names = old) +  
    scale_fill_gradient(low = "white", high = "#0073C2FF")
  
  
  
  }

