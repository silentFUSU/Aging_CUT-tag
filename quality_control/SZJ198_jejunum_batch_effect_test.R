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
library(limma)
tab = read.delim("data/samples/jejunum/all_test_10kb_bins.counts",skip=1)
colnames <- colnames(tab)[7:length(tab)]
new_colnames <- str_extract(colnames, "SZJ\\d+")
new_colnames[5] <- "SZJ198_NuoHe"
new_colnames[6] <- "SZJ198_JiangBei"
counts <- tab[7:length(tab)]
colnames(counts) <- new_colnames
y= DGEList(counts=counts)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
logCPMs <- cpm(y, log = TRUE)
pca <- prcomp(t(logCPMs))
to_plot <- data.frame(pca$x)
to_plot$rownames <- rownames(to_plot)

percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
ggplot(to_plot, aes(x=PC1, y=PC2)) + 
  geom_point(size=5) +theme_bw()+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+    
  geom_text_repel(
    aes(label = rownames),
    size = 10,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))
