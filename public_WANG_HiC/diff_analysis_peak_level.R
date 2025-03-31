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
window_size=5000
gap_size=10000
peak_preprocess_peak_level <- function(antibody){
  tab = read.delim(paste0("data/public_data/cellular_aging_GSE133292/Chip_seq/H3K9me3_young_old_merge-W",window_size,"-G",gap_size,"-E100.counts"),skip=1)
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  pattern <- ".*bam\\.(SRR[0-9]+).*"
  colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
  
  
  age <- c("G","G","DS","DS")
  colnames(counts) <- paste0(colnames(counts),"-",age)
  y= DGEList(counts=counts,group=age)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  y$samples$year <- age
  y$samples$year <- factor(y$samples$year,c("G","DS"))
  y <- calcNormFactors(y)
  design <- model.matrix(~year, y$samples)
  y<-estimateCommonDisp(y)
  y<-estimateGLMTagwiseDisp(y,design)
  fit_tag = glmFit(y,design)
  lrt = glmLRT(fit_tag, coef = 2)
  tab<-tab[keep,]
  
  out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
              "PValue.old-young"=lrt$table$PValue,"FDR.old-young"= p.adjust(lrt$table$PValue,method="BH"),
              "LogFC.old-young"=lrt$table$logFC)
  
  out$Significant <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$`LogFC.old-young`) >= log2(1.2), 
                            ifelse(out$`LogFC.old-young` > log2(1.2), "Up", "Down"), "Stable")
  out_sort <- out[order(out$`FDR.old-young`),]
  write.csv(out, paste0("data/public_data/cellular_aging_GSE133292/Chip_seq/H3K9me3-W",window_size,"-G",gap_size,"_peak_diff.csv"),row.names = F)
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  # if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
  #   peaks <- read.delim(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_10kb_in_young_old_merge-W1000-G3000-E100.bed"),header = F)
  # }else{
  #   peaks <- read.delim(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_1kb_in_young_old_merge_macs_narrowpeak.bed"),header = F)
  # }
  # out <- out[which(out$Geneid %in% peaks$V4),]
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
    ggtitle("H3K9me3 DS vs G")+
    annotate("text", x = min(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  print(p)
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
    ggtitle("H3K9me3 DS vs G")+
    annotate("text", x = min(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  print(p)
  }

MA_plot <- function(){
  df <- read.csv(paste0("data/public_data/cellular_aging_GSE133292/Chip_seq/H3K9me3-W",window_size,"-G",gap_size,"_peak_diff.csv"))
  df <- df[,c("Geneid","Length","logCPM","LogFC.old.young","Significant")]  
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  ggplot(
    df, aes(x = `logCPM`, y = `LogFC.old.young`)) +
    geom_point(aes(color = Significant, size = Length),alpha=0.2) +
    scale_color_manual(values = colour) +
    labs(x="Log2(CPM)",
         y="Log2(Fold Change)") +
    theme_bw()+
    theme(text = element_text(size = 20))+
    ggtitle("H3K9me3 DS vs G")+
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
    theme(text = element_text(size = 20))+
    ggtitle("H3K9me3 DS vs G")+
    annotate("text", x = max(df$logCPM), y = min(df$LogFC.old.young), label = nrow(df[which(df$Significant=="Down"),]), vjust = 0, hjust = 1,colour="blue",size=5)+
    annotate("text", x = max(df$logCPM), y = max(df$LogFC.old.young), label = nrow(df[which(df$Significant=="Up"),]), vjust = 1, hjust = 1,colour="red",size=5)
  }
