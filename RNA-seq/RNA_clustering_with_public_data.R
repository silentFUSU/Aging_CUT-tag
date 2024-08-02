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
tab = read.csv("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/GSE132040/GSE132040_190214_A00111_0269_AHH3J3DSXX_190214_A00111_0270_BHHMFWDSXX.csv")
tab <- tab[-c(54353:54357),]
tab <- as.data.frame(t(tab))
colnames(tab) <- tab[1,]
tab <- tab[-1,]
new_row_name <- sub("\\..*", "", rownames(tab)) 
tab$Sample.name <- new_row_name
table = read.delim("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/GSE132040/GSE132040_MACA_Bulk_metadata.csv",sep = ',')
tab <- merge(tab,table[,c(1,3,5,7)],by="Sample.name")
tab <- tab[!grepl("^NA", tab$source.name), ] 
rownames(tab) <- tab$Sample.name
tab <- tab[,-1]
tab <- tab[which(tab$characteristics..sex=="m"),]
tab <- tab[which(tab$characteristics..age %in% c("3","24")),]
count <- as.data.frame(t(tab[,-c(54353:54355)])) 

tab_ourdata <- read.table("data/samples/RNA/combined-chrM.counts",header = T)
# tab_ourdata <- read.table("data/raw_data/20240430_LLX/combined-chrM.counts",header = T)
rownames(tab_ourdata) <- tab_ourdata$Geneid
tab_ourdata <- tab_ourdata[,-1]
colnames <- colnames(tab_ourdata)[6:ncol(tab_ourdata)]
pattern <- ".*bam\\.(LLX[0-9]+).*"
new_colnames <- gsub(pattern, "\\1", colnames)
sorted_index <- order(new_colnames)
order_colnames <- new_colnames[sorted_index] 
colnames(tab_ourdata)[6:ncol(tab_ourdata)] <- new_colnames
counts <- tab_ourdata[,order_colnames] 


count_merge <- merge(count, counts,by = 'row.names')
rownames(count_merge) <- count_merge$Row.names
count_merge <- count_merge[,-1]
count_merge[] <- lapply(count_merge, as.numeric) 

pbmc <- CreateSeuratObject(count=count_merge)
group_info <- tab[,c(54353:54355)]
group_info$source.name <- sub("_(\\w+)_\\d*|_(\\w+)$", "\\1",  group_info$source.name)  
group_info$source.name <- paste0("TMR-",group_info$source.name)
group <- read.csv("data/samples/RNA/sample_tissue_info.csv",sep = '\t')
group <- group[order(group$SampleID),]
group <- group$TissueName
age <- rep(c("3m","24m"),ncol(counts)/2)


new_row <- data.frame(source.name=group,characteristics..age=age,characteristics..sex="m")
# rownames(new_row) <- paste0(new_row$source.name,"-",new_row$characteristics..age,"-",rep(c("rep1","rep2"),7))
group_info <- rbind(group_info,new_row)

pbmc$group <- group_info$source.name
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:30)
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.5,verbose = FALSE)
# pbmc <- subset(pbmc, cells = rownames(pbmc@meta.data)[which(pbmc$group %in% c("colon","cecum","jejunum","ileum","SmallIntestine"))])
pbmc$rownames <- rownames(pbmc@meta.data)
DimPlot(pbmc, label = T, pt.size = 1.5, label.size = 3,group.by = "group",repel = T) 
DimPlot(pbmc, label = T, pt.size = 1.5, label.size = 3,group.by = "rownames",repel = T) 

ggsave("result/RNA/test_clustering_public_data_umap.png",width = 30,height = 30,type="cairo")





count_merge[] <- lapply(count_merge, as.numeric)  
y= DGEList(counts=count_merge,group = group_info$source.name)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
logCPMs <- cpm(y, log = TRUE)
pca <- prcomp(t(logCPMs))
to_plot <- data.frame(pca$x, group_info$source.name)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
ggplot(to_plot, aes(x=PC1, y=PC2, color=group_info.source.name)) + 
  geom_point(size=10) +theme_bw()+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+    
  geom_text_repel(
    aes(label = group_info.source.name),
    size = 10,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))
ggsave("result/RNA/test_clustering_public_data.png",width = 30,height = 30,type="cairo")

