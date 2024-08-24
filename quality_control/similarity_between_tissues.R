rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(patchwork)
library(edgeR)
library(corrplot) 
antibody <- "H3K9me3"
method <-"pearson"
options(bitmapType = "cairo")  
search_table_cut <- read.csv("data/samples/all/CUTTag_search_table.csv")
search_table_atac <- read.csv("data/samples/all/ATAC_search_table.csv")
search_table <- rbind(search_table_cut,search_table_atac)
correlation_clustering <- function(antibody,method){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    bin_size <- "10kb"
  }else{
    bin_size <- "1kb"
  }
  df <- read.csv(paste0("data/samples/all/",antibody,"/merge-",bin_size,"_bins.counts"),sep = "\t",skip = 1)
  sample_info <- search_table
  rownames(sample_info) <- sample_info$sample_name
  counts<-df[,c(7:ncol(df))]
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  if(antibody == "ATAC"){
    tissues <- sub(".*ATAC\\.([^.]*)\\.ATAC.*", "\\1", colnames(counts))
  }else{
    tissues <- sub(".*samples\\.([^.]*)\\.H3.*", "\\1", colnames(counts))
  }
  colnames(counts) <- gsub(pattern, "\\1", colnames(counts))
  counts <- counts[,which(colnames(counts) %in% sample_info$sample_name)]
  sample_info <- sample_info[which(sample_info$sample_name %in% colnames(counts)),]
  sample_info$sample_name <- factor(sample_info$sample_name,levels=colnames(counts))
  sample_info <- sample_info[order(sample_info$sample_name),]
  y= DGEList(counts=counts)
  keep = which(rowSums(cpm(y)>1)>=10)
  y = y[keep,]
  logCPMs <- as.data.frame(cpm(y, log = TRUE))
  batch <- c()
  for(i in c(1:length(unique(tissues)))){
    if(tissues[i]=="brain"){
      batch <- c(batch,c("batch1","batch2","batch1","batch2"))
    }
      batch <- c(batch,rep(c(rep("batch1", 2), rep("batch2", 2)), 1))
  }
  
  
  for(i in c(1:length(unique(sample_info$tissue)))){
    if(i == 1){
      t_logCPMs <- logCPMs[,c(((i-1)*4+1):(i*4))]
      t_batch <- batch[c(((i-1)*4+1):(i*4))]
      logCPMs_corrected <- limma::removeBatchEffect(t_logCPMs, batch = t_batch)
    }else{
      t_logCPMs <- logCPMs[,c(((i-1)*4+1):(i*4))]
      t_batch <- batch[c(((i-1)*4+1):(i*4))]
      t_logCPMs_corrected <- limma::removeBatchEffect(t_logCPMs, batch = t_batch)
      logCPMs_corrected <- cbind(logCPMs_corrected,t_logCPMs_corrected)
    }
  }
  logCPMs_corrected  <- as.data.frame(logCPMs_corrected)
  if(antibody == "ATAC"){ 
    logCPMs_corrected[,c(1:4)] <- logCPMs[,c(1:4)]
  }
  logCPMs <- logCPMs_corrected
  sample_correlation <- cor(logCPMs, method = method)  
  
  
  annotation <- data.frame(tissues=sample_info$tissue)
  annotation$tissues <- factor(annotation$tissues)
  rownames(annotation) <- sample_info$sample_name
  annotation$sample_name <- rownames(annotation)
  annotation <- merge(annotation,sample_info[,c("sample_name","age")],by="sample_name")
  rownames(annotation) <- annotation$sample_name
  annotation <- annotation[,-1]
  
  colors <- read.table("data/samples/30_distinct_color.txt")
  colors <- colors$V1[1:24]
  names(colors) <- levels(annotation$tissues)
  annotation_colors <- list(tissues=colors)
  pheatmap::pheatmap(sample_correlation,show_rownames = F,show_colnames = F,
                     annotation_row = annotation,annotation_colors = annotation_colors,main = antibody,
                     filename = paste0("result/all/clustering/",antibody,"_",method,".png"),type="png",width = 8,height = 7,border_color = NA)
  
}
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me3","H3K4me1")
for(antibody in antibodys){
  correlation_clustering(antibody,"pearson")
}
