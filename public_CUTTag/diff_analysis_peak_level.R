rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(ggrepel)
antibody <- "H3K9me3"
tissue <- "cellline"
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
      tissue_label <- "IWAT"
    }
  }
  return(tissue_label)
}
window_size=5000
gap_size=10000
peak_preprocess_peak_level <- function(tissue,antibody){
  search_table <- read.csv("data/public_data/werner_syndrom_aging/CUTTag_search_table_used_in_diff_batch.csv")
  # if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
  #   bin_size <- "10kb"
  # }else{
  #   bin_size <- "1kb"
  # }
  tab = read.delim(paste0("data/public_data/werner_syndrom_aging/H3K9me3_young_old_merge-W",window_size,"-G",gap_size,"-E100.counts"),skip=1)
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  pattern <- ".*bam\\.(SRR[0-9]+).*"
  colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
  # counts <- counts[,-5]
  # colnames(counts) <- c("JC_R1_H3K9me3","JC_R2_H3K9me3","VC_R1_H3K9me3","VC_R2_H3K9me3")
  # counts <- counts[,c(1,2,5,6)]
  search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
  counts <- counts[,which(colnames(counts) %in% search_table$sample_name)]
  search_table$sample_name <- factor(search_table$sample_name, levels = colnames(counts))
  search_table <- search_table[order(search_table$sample_name),]
  age <- search_table$age
  mouse_ID <- search_table$mouse_ID
  # age[which(age=="3m")] <- "young"
  # age[which(age=="24m")] <- "old"
  colnames(counts) <- paste0(colnames(counts),"-",age,"-",mouse_ID)
  y= DGEList(counts=counts,group=age)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  y$samples$year <- age
  y$samples$year <- factor(y$samples$year,c("young","old"))
  y <- calcNormFactors(y)
  design <- model.matrix(~year, y$samples)
  y<-estimateCommonDisp(y)
  y<-estimateGLMTagwiseDisp(y,design)
  fit_tag = glmFit(y,design)
  lrt = glmLRT(fit_tag, coef = which(colnames(design) == "yearold"))
  tab<-tab[keep,]
  
  out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
              "PValue.old-young"=lrt$table$PValue,"FDR.old-young"= p.adjust(lrt$table$PValue,method="BH"),
              "LogFC.old-young"=lrt$table$logFC)
  
  out$Significant <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$`LogFC.old-young`) >= log2(1.2), 
                            ifelse(out$`LogFC.old-young` > log2(1.2), "Up", "Down"), "Stable")
  out_sort <- out[order(out$`FDR.old-young`),]
  # dir.create(paste0("data/samples/all/diff_table/",antibody))
  write.csv(out_sort,paste0("data/public_data/werner_syndrom_aging/H3K9me3-W",window_size,"-G",gap_size,"_peak_diff.csv"),row.names = F)
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  p <- ggplot(
    out, aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`))) +
    geom_point(aes(color = Significant, size = Length),alpha=0.2) +
    scale_color_manual(values = colour) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (p-value)") +
    theme_bw()+
    theme(text = element_text(size = 20))+
    ggtitle(paste0(tissue_label_change(tissue)," ",antibody))+
    annotate("text", x = min(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  out <- out[which(out$Length > 100000),]
  p <- ggplot(
    out, aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`))) +
    geom_point(aes(color = Significant, size = Length),alpha=0.2) +
    scale_color_manual(values = colour) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (p-value)") +
    theme_bw()+
    theme(text = element_text(size = 20))+
    ggtitle(paste0(tissue_label_change(tissue)," ",antibody))+
    annotate("text", x = min(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  return(p)
}


MA_plot <- function(tissue,antibody){
  p_list <- list()
  df <- read.csv(paste0("data/public_data/werner_syndrom_aging/H3K9me3-W",window_size,"-G",gap_size,"_peak_diff.csv"))
  df <- df[,c("Length","logCPM","LogFC.old.young","Significant")]  
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  ggplot(
    df, aes(x = `logCPM`, y = `LogFC.old.young`)) +
    geom_point(aes(color = Significant, size = Length),alpha=0.2) +
    scale_color_manual(values = colour) +
    labs(x="Log2(CPM)",
         y="Log2(Fold Change)") +
    theme_bw()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ggtitle(paste0(tissue_label_change(tissue)," ",antibody))+
    annotate("text", x = max(df$logCPM), y = min(df$LogFC.old.young), label = nrow(df[which(df$Significant=="Down"),]), vjust = 0, hjust = 1,colour="blue",size=5)+
    annotate("text", x = max(df$logCPM), y = max(df$LogFC.old.young), label = nrow(df[which(df$Significant=="Up"),]), vjust = 1, hjust = 1,colour="red",size=5)
  df <- df[which(df$Length > 100000),]
  ggplot(
    df, aes(x = `logCPM`, y = `LogFC.old.young`)) +
    geom_point(aes(color = Significant, size = Length),alpha=0.2) +
    scale_color_manual(values = colour) +
    labs(x="Log2(CPM)",
         y="Log2(Fold Change)") +
    theme_bw()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ggtitle(paste0(tissue_label_change(tissue)," ",antibody))+
    annotate("text", x = max(df$logCPM), y = min(df$LogFC.old.young), label = nrow(df[which(df$Significant=="Down"),]), vjust = 0, hjust = 1,colour="blue",size=5)+
    annotate("text", x = max(df$logCPM), y = max(df$LogFC.old.young), label = nrow(df[which(df$Significant=="Up"),]), vjust = 1, hjust = 1,colour="red",size=5)
  return(p_list)
}