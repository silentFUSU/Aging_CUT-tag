rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(Polychrome)
library(Seurat)
tissue<-"pancreas"
tissue_label_change <- function(tissue){
  if(tissue=="FC"){
    tissue_label <- "Cortex"
  }else if(tissue == "Hip"){
    tissue_label <- "Hippocampus"
  }else if(tissue == "CB"){
    tissue_label <- "Cerebellum"
  }else{
    tissue_label <- str_to_title(tissue)
    if(tissue_label == "Bonemarrow"){
      tissue_label <- "Bone Marrow"
    }
  }
  return(tissue_label)
} 
search_table <- read.csv("data/samples/all/RNA_search_table.csv")
tab = read.delim(paste0("data/raw_data/pancreas_RNA_test/combined-chrM.nodup.counts"),skip=1)
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
group <- group[which(group$SampleID %in% colnames(counts)),]
group <- group[order(group$SampleID),]
age <- group[which(group$SampleID %in% colnames(counts)),"Age"]

y= DGEList(counts=counts,group=age)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]

logCPMs <- cpm(y, log = TRUE)
logCPMs_corrected <- limma::removeBatchEffect(logCPMs, batch = c("batch1","batch1","batch1","batch2","batch2","batch2","batch2"))
pca <- prcomp(t(logCPMs_corrected))
to_plot <- data.frame(pca$x, age = paste0(y$samples$group))
to_plot$rownames <- rownames(to_plot)
table <- search_table[which(search_table$sample_name %in% to_plot$rownames),]
to_plot$rownames <- paste0(table$sample_name,"-",table$mouse_ID,"-",table$age)
to_plot$age <- factor(to_plot$age,levels=c("3m","24m"))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))

p <- ggplot(to_plot, aes(x=PC1, y=PC2, color=age)) + 
  geom_point(size=5) +theme_bw()+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
  geom_text_repel(  
    data = to_plot,  
    aes(x = PC1, y = PC2, label = rownames, color = age),  
    size = 5,  
    box.padding = unit(0.35, "lines"),  
    point.padding = unit(0.3, "lines")  
  ) +
  ggtitle(tissue_label_change(tissue))
print(p)

all = read.delim(paste0("data/samples/RNA/combined-chrM.nodup.counts"),skip=1)
tab = read.delim(paste0("data/raw_data/pancreas_RNA_test/combined-chrM.nodup.counts"),skip=1)
tab <- merge(tab[,c(1,7:9)],all[,c(1,7:ncol(all))],by="Geneid")
rownames(tab) <- tab$Geneid
tab <- tab[,-1]

colnames <- colnames(tab)
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+).*"
new_colnames <- gsub(pattern, "\\1", colnames)
colnames(tab)<- new_colnames
sorted_index <- order(new_colnames)
order_colnames <- new_colnames[sorted_index] 
counts <- tab[,order_colnames]   
group <- read.csv("data/samples/RNA/sample_tissue_info.csv",sep = ',')
group <- group[which(group$SampleID %in% colnames(counts)),]
group <- group[order(group$SampleID),]
age <- group[which(group$SampleID %in% colnames(counts)),"Age"]

y= DGEList(counts=counts,group=age)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
logCPMs <- cpm(y, log = TRUE)
pca <- prcomp(t(logCPMs))
to_plot <- data.frame(pca$x)
to_plot$sample_name <- rownames(to_plot)
to_plot <- merge(to_plot,search_table,by="sample_name")
to_plot$age <- factor(to_plot$age,levels = c("3m","24m"))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
color <-alphabet.colors(26)
color <- setNames(color,unique(to_plot$tissue))
new_tissues <- "Pancreas"
ggplot(to_plot, aes(x=PC1, y=PC2, color=tissue, shape=age)) + 
  geom_point(size=5) +theme_bw()+
  scale_color_manual(values = color) +
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+    
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

count_merge <- merge(count, counts,by = 'row.names')
rownames(count_merge) <- count_merge$Row.names
count_merge <- count_merge[,-1]
count_merge[] <- lapply(count_merge, as.numeric) 

pbmc <- CreateSeuratObject(count=count_merge)
group_info <- tab[,c(54353:54355)]
group_info$source.name <- sub("_(\\w+)_\\d*|_(\\w+)$", "\\1",  group_info$source.name)  
group_info$source.name <- paste0("TMS-",group_info$source.name)
group <- read.csv("data/samples/RNA/sample_tissue_info.csv",sep = ',')
group <- group[which(group$SampleID %in% rownames(pbmc@meta.data)),]
group <- group[order(group$SampleID),]
tissue <- group$TissueName

age <- group$Age

new_row <- data.frame(source.name=tissue,characteristics..age=age,characteristics..sex="m")
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
DimPlot(pbmc, label = T,cells.highlight = c("LLX539","LLX540","LLX541"), pt.size = 1.5, label.size = 3,group.by = "group",repel = T) + ggtitle("Tissue")
