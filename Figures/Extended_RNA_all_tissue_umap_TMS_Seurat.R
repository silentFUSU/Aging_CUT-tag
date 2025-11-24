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

tab_ourdata <- read.table("data/samples/RNA/all_tissues_combined-chrM.counts",header = T)
rownames(tab_ourdata) <- tab_ourdata$Geneid
tab_ourdata <- tab_ourdata[,-1]
colnames <- colnames(tab_ourdata)[6:ncol(tab_ourdata)]
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|HM[0-9]+).*"
new_colnames <- gsub(pattern, "\\1", colnames)
sorted_index <- order(new_colnames)
order_colnames <- new_colnames[sorted_index] 
colnames(tab_ourdata)[6:ncol(tab_ourdata)] <- new_colnames
counts <- tab_ourdata[,order_colnames] 


count_merge <- merge(count, counts,by = 'row.names')
rownames(count_merge) <- count_merge$Row.names
count_merge <- count_merge[,-1]
count_merge[] <- lapply(count_merge, as.numeric) 


group_info <- tab[,c(54353:54355)]
group_info$source.name <- sub("_(\\w+)_\\d*|_(\\w+)$", "\\1",  group_info$source.name)  
group_info$source.name <- paste0("TMS-",group_info$source.name)

pbmc <- CreateSeuratObject(count=count_merge)
group <- read.csv("data/samples/all/RNA_search_table.csv",sep = ',')
rownames(group) <- group$sample_name
group <- group[which(group$sample_name %in% rownames(pbmc@meta.data)),]
group <- group[order(group$sample_name),]
tissue <- group$tissue
age <- group$age

new_row <- data.frame(source.name=tissue,characteristics..age=age,characteristics..sex="m")
group_info <- rbind(group_info,new_row)
group_info$characteristics..age[which(group_info$characteristics..age %in% c("3", "3m"))] <- "young"
group_info$characteristics..age[which(group_info$characteristics..age %in% c("24","24m"))] <- "old"

pbmc$group <- group_info$source.name
pbmc$age <- group_info$characteristics..age

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:30)
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.5,verbose = FALSE)
pbmc$rownames <- rownames(pbmc@meta.data)
DimPlot(pbmc,label = T, pt.size = 1.5, label.size = 3,group.by = "group",repel = T) +ggtitle("RNA wtih TMS data")

pbmc$tissue <-gsub("^TMS-", "", pbmc$group)
pbmc$condition <- ifelse(startsWith(pbmc$group, "TMS"), "TMS", "Own")
search_table <- read.csv("data/samples/all/RNA_search_table.csv")
search_table <- search_table[-which(search_table$tissue =="MEF"),]
color <- read.table("data/samples/30_distinct_color.txt")
tissues <- sort(unique(search_table$tissue))
tissues <- c(tissues,sort(unique(pbmc$tissue[which(! pbmc$tissue %in% tissues)])))
color <- setNames(c(color$V1,"#FF5733","#33FF57","#5733FF","#FF33A1","#33A1FF","#A1FF33"),tissues)
condition_shape <- c("TMS" = 0, "Own" = 16)

umap_df <- as.data.frame(pbmc@reductions$umap@cell.embeddings)
feature <- pbmc@meta.data

umap_df <- merge(umap_df,feature,by="row.names")
umap_df$shape_combo <- interaction(umap_df$age, umap_df$condition, sep = "_")
condition_shape <- c("young_Own" = 19, "old_Own" = 17, "young_TMS" = 21, "old_TMS" = 24)
unique_tissue_df <- umap_df[!duplicated(umap_df$group), ]

p <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2,color=tissue,shape=shape_combo, alpha = condition)) + 
  scale_color_manual(values = color) +
  scale_shape_manual(values = condition_shape) +
  scale_alpha_manual(values = c("TMS" = 0.5, "Own" = 1)) + 
  geom_point(size=3) +
  theme_bw()+
  geom_text_repel(data = unique_tissue_df, aes(label = group),
                  size = 5,color="black", hjust = 0.5, vjust = -1) +
  theme(text = element_text(size = 20))+
  ggtitle("RNA with TMS data")
  
ggsave("result/Sup_figures/RNA_with_TMS_all_tissues_UMAP.pdf",p,width = 15,height = 12)
