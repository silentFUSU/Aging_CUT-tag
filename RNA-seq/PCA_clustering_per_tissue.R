rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
tissue<-"skin"
PCA_per_tissue <- function(tissue){
  tab = read.delim(paste0("data/samples/RNA/",tissue,"/combined-chrM.nodup.counts"),skip=1)
  rownames(tab) <- tab$Geneid
  tab <- tab[,-1]
  colnames <- colnames(tab)[6:length(tab)]
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+).*"
  new_colnames <- gsub(pattern, "\\1", colnames)
  colnames(tab)[6:length(tab)] <- new_colnames
  sorted_index <- order(new_colnames)
  order_colnames <- new_colnames[sorted_index] 
  counts <- tab[,order_colnames]   
  group <- read.csv("data/samples/RNA/sample_tissue_info.csv",sep = ',')
  age <- group[which(group$SampleID %in% colnames(counts)),"Age"]

  y= DGEList(counts=counts,group=age)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  
  logCPMs <- cpm(y, log = TRUE)
  pca <- prcomp(t(logCPMs))
  to_plot <- data.frame(pca$x, age = paste0(y$samples$group))
  to_plot$rownames <- rownames(to_plot)
  
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
  use.pcs <- c(1,2)
  labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
  
  ggplot(to_plot, aes(x=PC1, y=PC2, color=age)) + 
    geom_point(size=5) +theme_bw()+
    xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
    geom_text_repel(  
      data = to_plot,  
      aes(x = PC1, y = PC2, label = rownames, color = age),  
      size = 5,  
      box.padding = unit(0.35, "lines"),  
      point.padding = unit(0.3, "lines")  
    ) +
    ggtitle(tissue)
  
}