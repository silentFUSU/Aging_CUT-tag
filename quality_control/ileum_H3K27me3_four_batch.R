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
tab = read.delim(paste0("data/samples/ileum/H3K27me3/H3K27me3_10kb_bins.counts"),skip=1)  
rownames(tab) <- tab$Geneid
counts <- tab[,c(7:ncol(tab))]
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
colnames(counts) <- gsub(pattern, "\\1", colnames(counts))
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
to_plot$age <- c("young","old","young","old","young","old","young","old")
to_plot$age <- factor(to_plot$age, levels = c("young","old"))
ggplot(to_plot, aes(x=PC1, y=PC2, color=age,shape=age)) + 
  geom_point(size=5) +
  theme_bw()+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
  geom_text_repel(  
    data = to_plot,  
    aes(x = PC1, y = PC2, label = sample_name, color = age),  
    size = 5,  
    box.padding = unit(0.35, "lines"),  
    point.padding = unit(0.3, "lines")  
  ) +
  ggtitle("ileum")

batch=c("batch1","batch1","batch2","batch2","batch3","batch3","batch4","batch4")
logCPMs_corrected <- limma::removeBatchEffect(logCPMs, batch = batch)
pca <- prcomp(t(logCPMs_corrected))
to_plot <- data.frame(pca$x)
to_plot$sample_name <- rownames(to_plot)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
to_plot$age <- c("young","old","young","old","young","old","young","old")
to_plot$age <- factor(to_plot$age, levels = c("young","old"))
ggplot(to_plot, aes(x=PC1, y=PC2, color=age,shape=age)) + 
  geom_point(size=5) +
  theme_bw()+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
  geom_text_repel(  
    data = to_plot,  
    aes(x = PC1, y = PC2, label = sample_name, color = age),  
    size = 5,  
    box.padding = unit(0.35, "lines"),  
    point.padding = unit(0.3, "lines")  
  ) +
  ggtitle("ileum")
