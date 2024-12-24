rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(edgeR)
library(MASS) 
library(gridExtra)
library(dplyr)
library(stringr)
tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas","Hip","cecum","bonemarrow","ileum","heart","thymus","stomach","skin","aorta","tongue","bladder","CB","jejunum","uterus","ovary")
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","ATAC","H3K27ac","H3K4me1","H3K4me3")
bin_size <- function(antibody){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    return("10kb")
  }else{
    return("1kb")
  }
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
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}

per_tissue_PCA <- function(tissue,antibodys){
  p_list <- list()  
  for(i in c(1:length(antibodys))){
    antibody <- antibodys[i]
    if (antibody == "ATAC"){
      search_table <- read.csv("data/samples/all/ATAC_search_table.csv")
      data_path <- "data/samples/ATAC/"
    }else{
      search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff.csv")
      data_path <- "data/samples/"
    }
    tab = read.delim(paste0(data_path,tissue,"/",antibody,"/",antibody,"_",bin_size(antibody),"_bins.counts"),skip=1)  
    pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
    colnames(tab)[7:length(tab)] <- gsub(pattern, "\\1", colnames(tab)[7:length(tab)])
    counts <- tab[7:length(tab)]
    if(tissue == "iWAT"){
      counts <- counts[, !grepl("^CKJ", names(counts))]  
    }
    search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
    # search_table <- search_table[which(search_table$mouse_ID != "110"),]
    counts <- counts[,which(colnames(counts) %in% search_table$sample_name)]
    search_table$sample_name <- factor(search_table$sample_name, levels = colnames(counts))
    search_table <- search_table[order(search_table$sample_name),]
    age <- search_table$age
    y= DGEList(counts=counts,group = age)
    keep = which(rowSums(cpm(y)>1)>=2)
    y = y[keep,]
    logCPMs <- cpm(y, log = TRUE)
    pca <- prcomp(t(logCPMs))
    to_plot <- data.frame(pca$x, age = paste0(y$samples$group))
    to_plot$rownames <- rownames(to_plot)
    percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
    use.pcs <- c(1,2)
    labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
    
    table <- search_table[which(search_table$sample_name  %in% to_plot$rownames),c(3:5)]
    to_plot$rownames <- paste0(table$sample_name,"-",table$mouse_ID,"-",table$age)
    to_plot$age <- factor(to_plot$age, levels = c("3m","24m"))
    p_list[[i]] <-
      ggplot(to_plot, aes(x=PC1, y=PC2, color=age)) + 
      geom_point(size=5) +theme_bw()+
      xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
      geom_text_repel(  
        data = to_plot,  
        aes(x = PC1, y = PC2, label = rownames, color = age),  
        size = 5,  
        box.padding = unit(0.35, "lines"),  
        point.padding = unit(0.3, "lines")  
      ) +
      ggtitle(paste(tissue_label_change(tissue), antibody))
  }
  combined_plot <- plot_a_list(p_list,2,4)
  ggsave(paste0("result/all/pca/per_tissue_plot/",tissue_label_change(tissue),"_all_Histone_modification_PCA.png"),combined_plot,width = 24,height = 10,type="cairo")
}

per_tissue_PCA_remove_batcheffect <- function(tissue,antibodys){
  p_list <- list()  
  for(i in c(1:length(antibodys))){
    antibody <- antibodys[i]
    if (antibody == "ATAC"){
      search_table <- read.csv("data/samples/all/ATAC_search_table.csv")
      data_path <- "data/samples/ATAC/"
    }else{
      search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff.csv")
      data_path <- "data/samples/"
    }
    tab = read.delim(paste0(data_path,tissue,"/",antibody,"/",antibody,"_",bin_size(antibody),"_bins.counts"),skip=1)  
    pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
    colnames(tab)[7:length(tab)] <- gsub(pattern, "\\1", colnames(tab)[7:length(tab)])
    counts <- tab[7:length(tab)]
    if(antibody == "ATAC"){
      counts <- counts[,c(1:4)]
    }
    search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
    # search_table <- search_table[which(search_table$mouse_ID != "110"),]
    counts <- counts[,which(colnames(counts) %in% search_table$sample_name)]
    search_table$sample_name <- factor(search_table$sample_name, levels = colnames(counts))
    search_table <- search_table[order(search_table$sample_name),]
    age <- search_table$age
    y= DGEList(counts=counts,group = age)
    keep = which(rowSums(cpm(y)>1)>=2)
    y = y[keep,]
    logCPMs <- cpm(y, log = TRUE)
    # batch=c("batch1","batch1","batch2","batch2","batch3","batch3","batch4","batch4")
    batch=c("batch1","batch1","batch2","batch2")
    logCPMs_corrected <- limma::removeBatchEffect(logCPMs, batch = batch)
    pca <- prcomp(t(logCPMs_corrected))
    to_plot <- data.frame(pca$x, age = paste0(y$samples$group))
    to_plot$rownames <- rownames(to_plot)
    percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
    use.pcs <- c(1,2)
    labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
    table <- search_table[which(search_table$sample_name  %in% to_plot$rownames),c(3:5)]
    to_plot$rownames <- paste0(table$sample_name,"-",table$mouse_ID,"-",table$age)
    
    p_list[[i]] <-
      ggplot(to_plot, aes(x=PC1, y=PC2, color=age)) + 
      geom_point(size=5) +theme_bw()+
      xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
      geom_text_repel(  
        data = to_plot,  
        aes(x = PC1, y = PC2, label = rownames, color = age),  
        size = 5,  
        box.padding = unit(0.35, "lines"),  
        point.padding = unit(0.3, "lines")  
      ) +
      ggtitle(paste(tissue_label_change(tissue), antibody))
  }
  combined_plot <- plot_a_list(p_list,2,4)
  ggsave(paste0("result/all/pca/per_tissue_plot_remove_batch_effect/",tissue_label_change(tissue),"_all_Histone_modification_PCA.png"),combined_plot,width = 24,height = 10,type="cairo")
}
tissues <- c("iWAT","BAT","mammarygland")
for(tissue in tissues){
  per_tissue_PCA(tissue,antibodys)
}
for(tissue in tissues){
  per_tissue_PCA_remove_batcheffect(tissue,antibodys)
}
