rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyr)
library(edgeR)
library(stringr)
library(dplyr)
library(ggrepel)
library(gridExtra)
options(bitmapType="cairo")  
antibody <- "H3K27me3"
tissue <- "kidney"
batch_or_not_scatter_plot <- function(tissue,antibody){
  if (!file.exists(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_10kb_bins_diff_after_remove_batch_effect.csv"))) {  
    return(0)  
  }  
  batch <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_10kb_bins_diff_after_remove_batch_effect.csv"))
  not_batch <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_10kb_bins_diff.csv"))
  colnames(batch)[which(colnames(batch)=="LogFC.old.young")] <- "remove_Batch" 
  colnames(not_batch)[which(colnames(not_batch)=="LogFC.old.young")] <- "not_remove_Batch"
  df <- merge(batch[,c("Geneid","remove_Batch")],not_batch[,c("Geneid","not_remove_Batch")])
  png(filename = paste0("result/",tissue,"/batch_or_not_logFC_scatterplot.png"), width = 1500, height = 1600, res = 300)  
  par(cex.lab = 1.5, cex.axis = 1.2,cex.main=1.5)  
  smoothScatter(df[,3] ~ df[,2],xlab = colnames(df)[2],ylab = colnames(df)[3],main = paste0(tissue," logFC"))
  abline(a = 0, b = 1, col = "red", lty = 2)  
  dev.off()  
  
  batch <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_10kb_bins_diff_after_remove_batch_effect.csv"))
  not_batch <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_10kb_bins_diff.csv"))
  colnames(batch)[which(colnames(batch)=="FDR.old.young")] <- "remove_Batch" 
  colnames(not_batch)[which(colnames(not_batch)=="FDR.old.young")] <- "not_remove_Batch"
  df <- merge(batch[,c("Geneid","remove_Batch")],not_batch[,c("Geneid","not_remove_Batch")])
  df[,c(2:3)] <- -log10(df[,c(2:3)])
  cols <- df[, 2:3]  
  finite_values <- as.vector(cols[is.finite(as.matrix(cols))])  
  png(filename = paste0("result/",tissue,"/batch_or_not_fdr_scatterplot.png"), width = 1500, height = 1600, res = 300)  
  par(cex.lab = 1.5, cex.axis = 1.2)  
  smoothScatter(df[,3] ~ df[,2],xlab = colnames(df)[2],ylab = colnames(df)[3],main = paste0(tissue,"-log10(fdr)"),xlim = c(0,max(finite_values)),ylim =  c(0,max(finite_values)))
  abline(a = 0, b = 1, col = "red", lty = 2)  
  dev.off()
}
tissues <- c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
             "thymus","skin","bladder","bonemarrow","Hip","heart",
             "muscle","jejunum","uterus","ovary","liver","tongue",
             "cecum","colon","testis","stomach","pancreas","iWAT")
for(tissue in tissues){
  batch_or_not_scatter_plot(tissue,"H3K9me3")
  
}

two_batch_compare <- function(tissue,antibody){
  search_table <- read.csv("data/samples/all/CUTTag_search_table.csv")
  if(antibody %in% c("H3K9me3","H3K27me3","H3K36me3")){
    bin_size <- "10kb"
  }else{
    bin_size <- "1kb"
  }
  tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins.counts"),skip=1)
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
  search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
  search_table$sample_name <- factor(search_table$sample_name, levels = colnames(counts))
  search_table <- search_table[order(search_table$sample_name),]
  age <- search_table$age
  mouse_ID <- search_table$mouse_ID
  age[which(age=="3m")] <- "young"
  age[which(age=="24m")] <- "old"
  colnames(counts) <- paste0(colnames(counts),"-",age,"-",mouse_ID)
  y= DGEList(counts=counts,group=age)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  cpm <- as.data.frame(cpm(y,log = T))
  cpm$geneid <- rownames(cpm)
  to_plot <- reshape2::melt(cpm)
  to_plot$batch <- "batch1"
  search_table$age[which(search_table$age=="3m")] <- "young"
  search_table$age[which(search_table$age=="24m")] <- "old"
  to_plot$batch[which(to_plot$variable %in% paste0(search_table$sample_name[3:4],"-",search_table$age[3:4],"-",search_table$mouse_ID[3:4]))] <- "batch2"
  t <- t.test(to_plot$value[which(to_plot$batch=="batch2")],to_plot$value[which(to_plot$batch=="batch1")])
  ggplot(to_plot, aes(x = variable, y = value, fill= batch)) +  
    scale_fill_brewer(palette = "Pastel1") +
    geom_boxplot() +  
    theme_minimal()+
    xlab(NULL)+
    ylab("log(CPM)")+
    theme(text = element_text(size = 20),axis.text.x = element_text(angle = 75, hjust = 1))+
    annotate("text", x = Inf, y = -Inf, label = paste("p-value =",  format(t$p.value, scientific = TRUE, digits = 3)  ),   
             hjust = 1.1, vjust = -1.1, size = 5, colour = "red") +
    annotate("text", x = Inf, y = Inf, label = paste("batch2 mean =",  round(t$estimate[[1]],2)),   
             hjust = 1.1, vjust = 1.1, size = 5, colour = "red") +
    annotate("text", x = -Inf, y = Inf, label = paste("batch1 mean =",  round(t$estimate[[2]],2)),   
             hjust = 0, vjust = 1.1, size = 5, colour = "red")
  
  peak <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
  t <- t.test(to_plot$value[which(to_plot$batch=="batch2" & to_plot$geneid %in% peak$Geneid[which(peak$Significant=="Up")])],to_plot$value[which(to_plot$batch=="batch1" &  to_plot$geneid %in% peak$Geneid[which(peak$Significant=="Up")])])
  ggplot(to_plot[which(to_plot$geneid %in% peak$Geneid[which(peak$Significant=="Up")]),], aes(x = variable, y = value, fill= batch)) +  
    scale_fill_brewer(palette = "Pastel1") +
    geom_boxplot() +  
    theme_minimal()+
    xlab(NULL)+
    ylab("log(CPM)")+
    theme(text = element_text(size = 20),axis.text.x = element_text(angle = 75, hjust = 1))+
    annotate("text", x = Inf, y = -Inf, label = paste("p-value =",  format(t$p.value, scientific = TRUE, digits = 3)  ),   
             hjust = 1.1, vjust = -1.1, size = 5, colour = "red") +
    annotate("text", x = Inf, y = Inf, label = paste("batch2 mean =",  round(t$estimate[[1]],2)),   
             hjust = 1.1, vjust = 1.1, size = 5, colour = "red") +
    annotate("text", x = -Inf, y = Inf, label = paste("batch1 mean =",  round(t$estimate[[2]],2)),   
             hjust = 0, vjust = 1.1, size = 5, colour = "red")
  
  t <- t.test(to_plot$value[which(to_plot$batch=="batch2" & to_plot$geneid %in% peak$Geneid[which(peak$Significant=="Down")])],to_plot$value[which(to_plot$batch=="batch1" &  to_plot$geneid %in% peak$Geneid[which(peak$Significant=="Down")])])
  ggplot(to_plot[which(to_plot$geneid %in% peak$Geneid[which(peak$Significant=="Down")]),], aes(x = variable, y = value, fill= batch)) +  
    scale_fill_brewer(palette = "Pastel1") +
    geom_boxplot() +  
    theme_minimal()+
    xlab(NULL)+
    ylab("log(CPM)")+
    theme(text = element_text(size = 20),axis.text.x = element_text(angle = 75, hjust = 1))+
    annotate("text", x = Inf, y = -Inf, label = paste("p-value =",  format(t$p.value, scientific = TRUE, digits = 3)  ),   
             hjust = 1.1, vjust = -1.1, size = 5, colour = "red") +
    annotate("text", x = Inf, y = Inf, label = paste("batch2 mean =",  round(t$estimate[[1]],2)),   
             hjust = 1.1, vjust = 1.1, size = 5, colour = "red") +
    annotate("text", x = -Inf, y = Inf, label = paste("batch1 mean =",  round(t$estimate[[2]],2)),   
             hjust = 0, vjust = 1.1, size = 5, colour = "red")
  }

batch <- batch[-which(batch$Significant =="Stable"),]
not_batch <- not_batch[-which(not_batch$Significant == "Stable"),]

df <- batch[-which(batch$Geneid %in% not_batch$Geneid),]
df2 <- not_batch[-which(not_batch$Geneid %in% batch$Geneid),]
df3 <-batch[which(batch$Geneid %in% not_batch$Geneid),]

