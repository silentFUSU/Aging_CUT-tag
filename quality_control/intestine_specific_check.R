rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(Seurat)
library(edgeR)
antibody <- "H3K27ac"
search_table <- read.csv("data/samples/all/CUTTag_search_table.csv")
if(antibody %in% c("H3K9me3","H3K36me3","H3K27me3")){
  bin_size <- "10kb"
}else{
  bin_size <- "1kb"
}
tab <- read.delim(paste0("data/samples/intestine/",antibody,"_",bin_size,"_bins.counts"),skip=1)
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
colnames <- colnames(tab)[7:length(tab)]
new_colnames <- gsub(pattern, "\\1", colnames)
colnames(tab)[7:length(tab)] <- new_colnames
sorted_index <- order(new_colnames)
order_colnames <- new_colnames[sorted_index] 
counts <- tab[,order_colnames] 
y= DGEList(counts=counts)
keep = which(rowSums(cpm(y)>1)>=4)
y = y[keep,]
logCPMs <- cpm(y, log = TRUE)

pca <- prcomp(t(logCPMs))
to_plot <- data.frame(pca$x)
to_plot$sample_name <- rownames(to_plot)
to_plot <- merge(to_plot, search_table, by="sample_name")
to_plot$label <- paste0(to_plot$sample_name,"-",to_plot$mouse_ID,"-",to_plot$age)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
color <- read.table("data/samples/7_distinct_color.txt")
color <- setNames(color$V1,unique(to_plot$tissue))
to_plot$age <- factor(to_plot$age, levels = c("3m","24m"))
ggplot(to_plot, aes(x=PC1, y=PC2, color=tissue,shape=age)) + 
  geom_point(size=5) +theme_bw()+ scale_color_manual(values = color)+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
  geom_text_repel(  
    data = to_plot,  
    aes(x = PC1, y = PC2, label = label, color = tissue),  
    size = 3,  
    box.padding = unit(0.35, "lines"),  
    point.padding = unit(0.3, "lines")  
  ) +
  ggtitle(paste0("Intestine tissues ", antibody))


correlation_matrix <- cor(logCPMs,method = "spearman")
search_table <- read.csv("data/samples/all/CUTTag_search_table.csv")
search_table <- search_table[which(search_table$sample_name %in% colnames(correlation_matrix)),]
annotation <- search_table[which(search_table$sample_name %in% colnames(correlation_matrix)),]
rownames(annotation) <- paste0(search_table$sample_name,"-",search_table$mouse_ID,"-",search_table$age)
annotation <- annotation[,c("tissue","age")]
annotation$age <- factor(annotation$age,c("3m","24m"))
color <- read.table("data/samples/7_distinct_color.txt")
antibody_color <- setNames(color$V1,unique(annotation$tissue))
age_color <- setNames(c("red","blue"),c("3m","24m"))
annotation_colors <- list(age=age_color,antibody=antibody_color)
search_table$sample_name <- factor(search_table$sample_name,levels = c(colnames(correlation_matrix)))
search_table <- search_table[order(search_table$sample_name),]
colnames(correlation_matrix) <- paste0(search_table$sample_name,"-",search_table$mouse_ID,"-",search_table$age)
rownames(correlation_matrix) <- paste0(search_table$sample_name,"-",search_table$mouse_ID,"-",search_table$age)
pheatmap::pheatmap(correlation_matrix,annotation = annotation, 
                   annotation_colors = annotation_colors, main = paste0("Intestine tissues ", antibody))


pbmc <- CreateSeuratObject(count=counts)
annotation <- search_table[which(search_table$sample_name %in% colnames(counts)),]
annotation <- annotation[,c("tissue","age")]
pbmc <- AddMetaData(pbmc, metadata = annotation)  
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 10)
pbmc <- RunUMAP(pbmc, dims = 1:10,n.neighbors = 4)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5,verbose = FALSE)
pbmc$label <- paste0(rownames(pbmc@meta.data),"-", pbmc$mouse_ID,"-",pbmc$age)


antibody <- "H3K27ac"
if(antibody %in% c("H3K9me3","H3K36me3","H3K27me3")){
  bin_size <- "10kb"
}else{
  bin_size <- "1kb"
}
tab <- read.delim(paste0("data/samples/intestine/",antibody,"_cecum_colon_",bin_size,"_bins.counts"),skip=1)
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
colnames <- colnames(tab)[7:length(tab)]
new_colnames <- gsub(pattern, "\\1", colnames)
colnames(tab)[7:length(tab)] <- new_colnames
sorted_index <- order(new_colnames)
order_colnames <- new_colnames[sorted_index] 
rownames(tab) <- tab$Geneid
counts <- tab[,order_colnames] 
y= DGEList(counts=counts)
keep = which(rowSums(cpm(y)>1)>=4)
y = y[keep,]
logCPMs <- cpm(y, log = TRUE)

young <-counts[,c(1,3,5,7)]
y= DGEList(counts=young,group=c("cecum","cecum","colon","colon"))
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
y <- calcNormFactors(y)
design <- model.matrix(~group, y$samples)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)
lrt = glmLRT(fit_tag, coef = 2)
tab<-tab[keep,]

out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
            "PValue.old-young"=lrt$table$PValue,"FDR.old-young"= p.adjust(lrt$table$PValue,method="BH"),
            "LogFC.old-young"=lrt$table$logFC)
out$Significant <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$`LogFC.old-young`) >= 0, 
                          ifelse(out$`LogFC.old-young` > 0, "Up", "Down"), "Stable")
diff_bin <- out$Geneid[which(out$Significant!="Stable")]
logCPMs <- logCPMs[which(rownames(logCPMs) %in% diff_bin),]
search_table <- read.csv("data/samples/all/CUTTag_search_table.csv")
search_table <- search_table[which(search_table$sample_name %in% colnames(logCPMs)),]
search_table$sample_name <- factor(search_table$sample_name,levels = c(colnames(logCPMs)))
search_table <- search_table[order(search_table$sample_name),]
colnames(logCPMs) <- paste0(search_table$sample_name,"-",search_table$mouse_ID,"-",search_table$age)
pheatmap::pheatmap(logCPMs,show_rownames = F,scale = "row",main = paste0(antibody," Areas of difference in young samples"))
