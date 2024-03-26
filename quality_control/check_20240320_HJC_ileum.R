rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(edgeR)
library(ggplot2)
library(ChIPseeker)
library(EnsDb.Mmusculus.v79)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
library(tidyr)
library(MASS) 
tab = read.delim("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240320_HJC_all_ileum/merge-1kb_bins.counts",skip=1)

counts = tab[,c(7:ncol(tab))]
colnames(counts) <-sub('.*\\.HJC_(\\d+)_S.*', 'HJC\\1', colnames(counts))
counts <- counts[,-which(colnames(counts) %in% c("HJC105","HJC106","HJC107","HJC108"))]
logCPMs <- cpm(y, log = T)
y= DGEList(counts=counts)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
y <- calcNormFactors(y)
logCPMs <- as.data.frame(cpm(y, log = T))


pca <- prcomp(t(logCPMs))
to_plot <- data.frame(pca$x, y$samples)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
info <- read.csv("data/raw_data/20240320_HJC_all_ileum/20240230_HJC_ileum.csv",row.names = 1,header = F)
to_plot <- merge(to_plot,info,by="row.names")
rownames(to_plot) <- to_plot$Row.names
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
ggplot(to_plot, aes(x=PC1, y=PC2,color=V4)) + 
  geom_point(size=5) +theme_bw()+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+    
  geom_text_repel(
    aes(label = rownames(to_plot)),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))
