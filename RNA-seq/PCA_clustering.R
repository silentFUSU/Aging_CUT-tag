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
library(biomaRt) 
library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
tab_ourdata <- read.table("data/samples/RNA/combined-chrM.counts",header = T)
# tab_ourdata <- read.table("data/raw_data/20240430_LLX/combined-chrM.counts",header = T)
rownames(tab_ourdata) <- tab_ourdata$Geneid
tab_ourdata <- tab_ourdata[,-1]
colnames <- colnames(tab_ourdata)[6:length(tab_ourdata)]
pattern <- ".*bam\\.(LLX[0-9]+).*"
new_colnames <- gsub(pattern, "\\1", colnames)
colnames(tab_ourdata)[6:length(tab_ourdata)] <- new_colnames
sorted_index <- order(new_colnames)
order_colnames <- new_colnames[sorted_index] 
counts <- tab_ourdata[,order_colnames] 

group <- read.csv("data/samples/RNA/sample_tissue_info.csv",sep = '\t')
group <- group[order(group$SampleID),]
group <- group$TissueName
age <- rep(c("3m","24m"),12)
counts[] <- lapply(counts, as.numeric)  
y= DGEList(counts=counts,group = group)
y$age <- age
# y$samples$group[which(rownames(y$samples)=="LLX501")] <- "colon"
# y$samples$group[which(rownames(y$samples)=="LLX505")] <- "cecum"
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
logCPMs <- cpm(y, log = TRUE)
pca <- prcomp(t(logCPMs))
to_plot <- data.frame(pca$x, group = paste0(y$samples$group,"_",y$age))
to_plot$rownames <- rownames(to_plot)

percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
ggplot(to_plot, aes(x=PC1, y=PC2, color=group)) + 
  geom_point(size=5) +theme_bw()+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+    
  geom_text_repel(
    aes(label = group),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))
