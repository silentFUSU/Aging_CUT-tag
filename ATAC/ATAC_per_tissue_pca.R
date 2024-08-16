rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(patchwork)
library(edgeR)
library(MASS) 
tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas","Hip","cecum","bonemarrow","ileum","heart","thymus","stomach","skin","aorta","tongue","bladder","CB","jejunum","uterus","ovary")
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
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
    }
  }
  return(tissue_label)
} 
search_table <- read.csv("data/samples/all/ATAC_search_table.csv")
per_tissue_PCA <- function(tissue){
  tab = read.delim(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_1kb_bins.counts"),skip=1)  
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames <- colnames(tab)[7:length(tab)]
  new_colnames <- gsub(pattern, "\\1", colnames)
  colnames(tab)[7:length(tab)] <- new_colnames
  sorted_index <- order(new_colnames)
  order_colnames <- new_colnames[sorted_index] 
  counts <- tab[,order_colnames] 
  if(tissue=="ovary"){
    age <- c("old","young","old","young")
  }else if(tissue=="brain"){
    age <- c("young","young","old","old")
  }else{
    age <- c("young","old","young","old")
  }
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
  table <- search_table[which(search_table$sample_name %in% to_plot$rownames),]
  to_plot$rownames <- paste0(table$sample_name,"-",table$mouse_ID,"-",table$age)
  p <-  ggplot(to_plot, aes(x=PC1, y=PC2, color=age)) + 
    geom_point(size=5) +theme_bw()+
    xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
    geom_text_repel(  
      data = to_plot,  
      aes(x = PC1, y = PC2, label = rownames, color = age),  
      size = 5,  
      box.padding = unit(0.35, "lines"),  
      point.padding = unit(0.3, "lines")  
    ) +
    ggtitle(paste(tissue_label_change(tissue)))
  return(p)
}
p_list <- list()  
for( i in c(1:length(tissues))){
  tissue <- tissues[i]
  p <- per_tissue_PCA(tissue)
  p_list[[i]] <- p
  }
combined_plot <- plot_a_list(p_list,4,6)
ggsave(paste0("result/all/pca/per_tissue_plot/All_tissues_ATAC_PCA.png"),combined_plot,width = 35,height = 20,type="cairo")


per_tissue_PCA_remove_batch_effect <- function(tissue){
  tab = read.delim(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_1kb_bins.counts"),skip=1)  
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames <- colnames(tab)[7:length(tab)]
  new_colnames <- gsub(pattern, "\\1", colnames)
  colnames(tab)[7:length(tab)] <- new_colnames
  sorted_index <- order(new_colnames)
  order_colnames <- new_colnames[sorted_index] 
  counts <- tab[,order_colnames] 
  if(tissue=="ovary"){
    age <- c("old","young","old","young")
  }else if(tissue=="brain"){
    age <- c("young","young","old","old")
  }else{
    age <- c("young","old","young","old")
  }
  y= DGEList(counts=counts,group = age)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  logCPMs <- cpm(y, log = TRUE)
  if(tissue != "brain"){
    batch=c("batch1","batch1","batch2","batch2")
    logCPMs_corrected <- limma::removeBatchEffect(logCPMs, batch = batch)
  }else{
    logCPMs_corrected <- logCPMs
  }
  pca <- prcomp(t(logCPMs_corrected))
  to_plot <- data.frame(pca$x, age = paste0(y$samples$group))
  to_plot$rownames <- rownames(to_plot)
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
  use.pcs <- c(1,2)
  labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
  table <- search_table[which(search_table$sample_name %in% to_plot$rownames),]
  to_plot$rownames <- paste0(table$sample_name,"-",table$mouse_ID,"-",table$age)
  p <-  ggplot(to_plot, aes(x=PC1, y=PC2, color=age)) + 
    geom_point(size=5) +theme_bw()+
    xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
    geom_text_repel(  
      data = to_plot,  
      aes(x = PC1, y = PC2, label = rownames, color = age),  
      size = 5,  
      box.padding = unit(0.35, "lines"),  
      point.padding = unit(0.3, "lines")  
    ) +
    ggtitle(paste(tissue_label_change(tissue)))
  return(p)
}
p_list <- list()  
for( i in c(1:length(tissues))){
  tissue <- tissues[i]
  p <- per_tissue_PCA_remove_batch_effect(tissue)
  p_list[[i]] <- p
}
combined_plot <- plot_a_list(p_list,4,6)
ggsave(paste0("result/all/pca/per_tissue_plot_remove_batch_effect/All_tissues_ATAC_PCA.png"),combined_plot,width = 35,height = 20,type="cairo")
