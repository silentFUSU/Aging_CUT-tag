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
library(Polychrome)
tab_ourdata <- read.table("data/samples/RNA/all_tissues_combined-chrM.counts",header = T)
tab_ourdata <- tab_ourdata[!grepl("chrY", tab_ourdata$Chr), ] 
new_tissues <- c("Ovary")
rownames(tab_ourdata) <- tab_ourdata$Geneid
tab_ourdata <- tab_ourdata[,-1]
colnames <- colnames(tab_ourdata)[6:length(tab_ourdata)]
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|HM[0-9]+).*"
new_colnames <- gsub(pattern, "\\1", colnames)
colnames(tab_ourdata)[6:length(tab_ourdata)] <- new_colnames
sorted_index <- order(new_colnames)
order_colnames <- new_colnames[sorted_index] 
counts <- tab_ourdata[,order_colnames] 

search_table <- read.csv("data/samples/all/RNA_search_table.csv")
counts[] <- lapply(counts, as.numeric)  
y= DGEList(counts=counts)
# y$samples$group[which(rownames(y$samples)=="LLX501")] <- "colon"
# y$samples$group[which(rownames(y$samples)=="LLX505")] <- "cecum"
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
logCPMs <- as.data.frame(cpm(y, log = TRUE))
pca <- prcomp(t(logCPMs))
to_plot <- data.frame(pca$x)
to_plot$sample_name <- rownames(to_plot)
to_plot <- merge(to_plot,search_table,by="sample_name")
to_plot$age <- factor(to_plot$age,levels = c("3m","24m"))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))

to_plot$tissue[which(to_plot$tissue=="Bat")] <-"BAT"
to_plot$tissue[which(to_plot$tissue=="Iwat")] <-"iWAT"
to_plot$tissue[which(to_plot$tissue=="Mammary gland")] <-"Mammary Gland"
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(unique(to_plot$tissue)))
ggplot(to_plot, aes(x=PC1, y=PC2, color=tissue, shape=age)) + 
  geom_point(size=5) +theme_bw()+
  scale_color_manual(values = color) +
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+ 
  ggtitle("RNA")+
  geom_text_repel(  
    data = subset(to_plot, to_plot$tissue %in% new_tissues),  
    aes(x = PC1, y = PC2, label = sample_name, color = tissue),  
    size = 5,  
    box.padding = unit(0.35, "lines"),  
    point.padding = unit(0.3, "lines")  
  ) +  
  geom_text_repel(  
    data = subset(to_plot, !to_plot$tissue %in% new_tissues),  
    aes(x = PC1, y = PC2, label = tissue, color = tissue),   
    size = 5,  
    box.padding = unit(0.35, "lines"),  
    point.padding = unit(0.3, "lines")  
  ) 

p <- ggplot(to_plot, aes(x=PC1, y=PC2, color=tissue, shape=age)) + 
  geom_point(size=5) +theme_bw()+
  scale_color_manual(values = color) +
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+ 
  ggtitle("RNA")
ggsave("result/Sup_figures/RNA_all_tissues_PCA.pdf",p,width = 6,height = 3)
