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

set.seed(1)
umap_result <- umap::umap(t(logCPMs))
umap_df <- data.frame(UMAP1 = umap_result$layout[, 1], UMAP2 = umap_result$layout[, 2], names = colnames(logCPMs))
umap_df <- merge(umap_df,search_table[,c("sample_name","tissue","age")],by.x="names",by.y="sample_name")
umap_df$tissue <- sapply(umap_df$tissue,tissue_label_change)

to_plot <- umap_df
tissue <- sort(unique(to_plot$tissue))
colours <- read.table("data/samples/30_distinct_color.txt")
colours <- setNames(colours$V1,tissue)
p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2,color=tissue,shape=age)) + 
  scale_color_manual(values = colours) +
  geom_point(size=5) +
  theme_bw()+
  theme(text = element_text(size = 20))+
  ggtitle("RNA")
ggsave("result/Sup_figures/RNA_all_tissues_UMAP.pdf",p,width = 10,height = 6)
