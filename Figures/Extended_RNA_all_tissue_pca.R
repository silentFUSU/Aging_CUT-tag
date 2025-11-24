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
tab_ourdata <- read.table("data/samples/RNA/all_tissues_combined-chrM.counts",header = T)
tab_ourdata <- tab_ourdata[!grepl("chrY", tab_ourdata$Chr), ] 
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

keep = which(rowSums(cpm(y)>1)>=5)
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
to_plot$tissue <- sapply(to_plot$tissue, tissue_label_change)

color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(unique(to_plot$tissue)))
p<- ggplot(to_plot, aes(x=PC1, y=PC2, color=tissue, shape=age)) + 
  geom_point(size=5) +theme_bw()+
  scale_color_manual(values = color) +
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+ 
  ggtitle("RNA")
ggsave("result/Sup_figures/RNA_all_tissues_PCA.pdf",p,width = 10,height = 6)
