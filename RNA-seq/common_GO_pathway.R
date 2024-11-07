rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/","/usr/local/lib64/R/library"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(edgeR)
library(ggplot2)
library(clusterProfiler)
library(stringr)
tissues <-  c("BAT","mammarygland","CB","lung","kidney","aorta","FC","spleen",
              "thymus","skin","bladder","bonemarrow","Hip","heart",
              "muscle","jejunum","uterus","ovary","liver","tongue",
              "cecum","colon","testis","stomach","pancreas","iWAT","ileum")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
tissue_label_change <- function(tissue){
  if(tissue=="FC"){
    tissue_label <- "Cortex"
  }else if(tissue == "Hip"){
    tissue_label <- "Hippocampus"
  }else if(tissue == "CB"){
    tissue_label <- "Cerebellum"
  }else{
    tissue_label <- str_to_title(tissue)
    if(tissue_label == "Bonemarrow"){
      tissue_label <- "Bone Marrow"
    }else if(tissue_label == "Mammarygland"){
      tissue_label <- "Mammary gland"
    }else if(tissue_label == "Iwat"){
      tissue_label <- "iWAT"
    }else if(tissue_label == "Bat"){
      tissue_label <- "BAT"
    }
  }
  return(tissue_label)
} 

common_GO_pathway <- function(tissues,condition){
  GO_list <- vector()
  GO_table <- data.frame()
  top <- 5
  for(tissue in tissues){
    if(file.exists(paste0("result/RNA/GO/table/",tissue_label_change(tissue),"_",condition,"_gene_GO.csv"))){
      table <- read.csv(paste0("result/RNA/GO/table/",tissue_label_change(tissue),"_",condition,"_gene_GO.csv"),row.names = 1)
      t_table <- table[,c("ID","Description","p.adjust")]
      colnames(t_table)[3] <- tissue_label_change(tissue)
      t_GO_list <- t_table[which(t_table[,3]<0.05),"Description"]
      t_GO_list <- t_GO_list[1:min(length(t_GO_list),top)]
      GO_list <- union(t_GO_list,GO_list)
      if(nrow(GO_table)==0){
        GO_table <- t_table
      }else{
        GO_table <- merge(GO_table,t_table[,c(2,3)],by="Description",all=T)
      }
    }
  }
  to_plot <- GO_table[,-2]
  to_plot <- to_plot[which(to_plot$Description %in% GO_list),]
  rownames(to_plot) <- to_plot$Description
  to_plot <- to_plot[,-1]
  to_plot <- -log10(to_plot)
  
  pheatmap::pheatmap(to_plot,cluster_rows = F,cluster_cols = F,main = paste0(condition," GO pathway"),breaks = seq(-log10(0.05), 30, length.out=101),na_col = "white")
  }


  