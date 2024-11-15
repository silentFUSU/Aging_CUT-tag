rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(diffHic)
library(edgeR)
library(data.table)
library(GenomicRanges) 
library(rtracklayer)  
library(Matrix)
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

pad_matrix <- function(m, nrow_target, ncol_target) {  
  current_nrow <- nrow(m)  
  current_ncol <- ncol(m)  
  if (current_nrow < nrow_target) {  
    m <- rbind(m, matrix(0, nrow = nrow_target - current_nrow, ncol = current_ncol))  
  }  
  if (current_ncol < ncol_target) {  
    m <- cbind(m, matrix(0, nrow = nrow(m), ncol = ncol_target - current_ncol))  
  }  
  return(m)  
}  

diffHiC_input_convert <- function(tissue,resolution){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  
  matrix_list <- list()
  for(sample in search_table$sample_name){
    matrix <- fread(paste0("data/samples/HiC/",tissue,"/raw_matrix/",sample,"_",resolution,".matrix"))
    sparse_matrix <- sparseMatrix(matrix$V1,matrix$V2,x = matrix$V3)
    matrix_list[[sample]] <- sparse_matrix
   }
  matrix_list <- lapply(matrix_list, pad_matrix, 273121, 273121)
  
  contact_list <- list()
  for(sample in search_table$sample_name){
    bin <- import.bed(paste0("data/samples/HiC/",tissue,"/raw_matrix/",sample,"_",resolution,"_abs.bed"))
    contact <- ContactMatrix(matrix_list[[sample]],bin,bin)
    contact_list[[sample]] <- contact
  }
  
  to.keep <- as.matrix(contact_list[[1]]) != 0  
  for (i in 2:length(contact_list)) {  
    to.keep <- to.keep | (as.matrix(contact_list[[i]]) != 0)  
  }  
  
  contact_keep_list <- list()
  for(sample in search_table$sample_name){
    iset <- deflate(contact_list[[sample]], extract=to.keep)    
    contact_keep_list[[sample]] <- iset
  }
  
  data <-  Reduce(cbind, contact_keep_list)  
  interactions(data) <- as(interactions(data), "ReverseStrictGInteractions")
}

tissue <- "lung"
resolution <- "10000"
data <- diffHiC_input_convert(tissue,resolution)
rm(list=setdiff(ls(), "data"))  


keep <- aveLogCPM(asDGEList(data)) > 0

