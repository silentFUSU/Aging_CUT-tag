rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(Seurat)
df <-  read.delim(paste0("data/samples/CB/merge-1kb_bins.counts"),skip=1)
counts = df[,c(7:ncol(df))]
rownames(counts)= df$Geneid
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
search_table <- read.table("data/samples/all/combined_analysis_enhancer/cellmarkfiletable.txt")
search_table <- search_table[which(search_table$V1==paste0("CB","_young1")),]
search_table$V3 <- sub("\\.bed$", "", search_table$V3)  

counts <- counts[,search_table$V3]
colnames(counts) <- search_table$V2
counts <- cbind(df[,1:6],counts)
rownames(counts) <- paste0(counts$Chr,":",counts$Start,"-",counts$End)
chromHMM <- read.delim(paste0("result/all/ChromHMM/all_tissues/14_all_tissues/split_1k/CB_young1_14_segments_1k.bed"),header = F)
rownames(chromHMM) <- paste0(chromHMM$V1,":",chromHMM$V2,"-",chromHMM$V3)
chromHMM <- chromHMM[,4,drop=F]

counts <- counts[,-c(1:6)]
counts <- counts[c(1:10000),]
pbmc <- CreateSeuratObject(count=t(counts))
pbmc$group <- chromHMM
pbmc$group <- factor(pbmc$group, levels=paste0("E",c(1:14)))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:2)
pbmc <- FindNeighbors(pbmc, dims = 1:2)
pbmc <- FindClusters(pbmc, resolution = 0.5,verbose = FALSE)
# pbmc <- subset(pbmc, cells = rownames(pbmc@meta.data)[which(pbmc$group %in% c("colon","cecum","jejunum","ileum","SmallIntestine"))])
# pbmc$rownames <- rownames(pbmc@meta.data)
p<-DimPlot(pbmc,label = T, pt.size = 1.5, label.size = 3,group.by = "group",repel = T) +ggtitle(NULL)
ggsave("result/CB/CB_young1_chromHMM_umap.png",p,width=9,height=8)