rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
library(limma)
library(ggrepel)
antibody <- "H3K27me3"
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
bin_size <- function(antibody){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    return("10kb")
  }else{
    return("1kb")
  }
}
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}
per_markers_all_tissue_pca <- function(antibody){
  tab = read.delim(paste0("data/samples/all/",antibody,"/merge-",bin_size(antibody),"_bins.counts"),row.names = 1,skip=1)    
  tab <- tab[-which(tab$Chr=="chrY"),]
  counts <- tab[,c(6:ncol(tab))]
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(counts) <- gsub(pattern, "\\1", colnames(counts))
  if(antibody=="ATAC"){
    search_table <- read.csv("data/samples/all/ATAC_search_table.csv")
  }else{
    search_table <- read.csv("data/samples/all/CUTTag_search_table.csv")
  }
  search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
  y= DGEList(counts=counts)
  keep = which(rowSums(cpm(y)>1)>=5)
  y = y[keep,]
  logCPMs <- cpm(y, log = TRUE)
  pca <- prcomp(t(logCPMs))
  to_plot <- data.frame(pca$x)
  to_plot$sample_name <- rownames(to_plot)
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
  use.pcs <- c(1,2)
  labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
  to_plot <- merge(to_plot,search_table,by="sample_name")
  to_plot$sample_name <- paste0(to_plot$sample_name,"-",to_plot$mouse_ID,"-",to_plot$age)
  to_plot$age <- factor(to_plot$age, levels = c("3m","24m"))
  tissue <- sort(unique(to_plot$tissue))
  colours <- read.table("data/samples/30_distinct_color.txt")
  colours <- setNames(colours$V1,tissue)
  p <- ggplot(to_plot, aes(x=PC1, y=PC2, color=tissue,shape=age)) + 
    scale_color_manual(values = colours) +
    geom_point(size=5) +
    theme_bw()+
    xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
    ggtitle(antibody)
  return(p)
}
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","ATAC","H3K27ac","H3K4me1","H3K4me3")
p_list <- list()
for(i in c(1:length(antibodys))){
  p_list[[i]] <- per_markers_all_tissue_pca(antibodys[i])
}
combined_plot <- plot_a_list(p_list,no_of_rows = 2,no_of_cols = 4)
ggsave("result/all/pca/all_tissues_plot/per_markers_all_tissues_PCA.png",width = 25,height = 10,type="cairo")

tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
window_size <- "5000"
gap_size <- "10000"
per_markers_all_tissue_pca_peak_level(tissues,antibody){
  search_table <- read.csv(paste0("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv"))
  counts_summary <- data.frame()
  for(tissue in tissues){
    tab <- read.table(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_young_old_merge-W",window_size,"-G",gap_size,"-E100_bedtools.counts"),header = T)
    tab <- tab[which(tab$Length>100000 & tab$Chr !="chrY"),]
    pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
    counts <- tab[,c(7:ncol(tab))]
    colnames(counts) <- gsub(pattern, "\\1", colnames(counts))
    counts <- counts[,which(colnames(counts) %in% search_table$sample_name)]
    counts$Geneid <- tab$Geneid
    if(ncol(counts_summary)==0){
      counts_summary <- counts
    }else{
      counts_summary <- merge(counts_summary,counts,by="Geneid")
    }
  }
  counts_summary <- counts_summary[,-1]
  y <-  DGEList(counts=counts_summary)
  logCPMs <- cpm(y, log = TRUE)
  pca <- prcomp(t(logCPMs))
  to_plot <- data.frame(pca$x)
  to_plot$sample_name <- rownames(to_plot)
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
  use.pcs <- c(1,2)
  labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
  to_plot$sample_name <- rownames(to_plot)
  to_plot <- merge(to_plot,search_table,by="sample_name")
  to_plot$tissue <- sapply(to_plot$tissue, tissue_label_change)  
  to_plot$age <- factor(to_plot$age, levels = c("3m","24m"))
  tissue <- sort(unique(to_plot$tissue))
  colours <- read.table("data/samples/30_distinct_color.txt")
  colours <- setNames(colours$V1,tissue)
  p <- ggplot(to_plot, aes(x=PC1, y=PC2, color=tissue,shape=age)) + 
    scale_color_manual(values = colours) +
    geom_point(size=5) +
    theme_bw()+
    xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
    ggtitle(antibody)
  }