rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(ggrepel)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(clusterProfiler)
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
tab <-  read.delim(paste0("data/samples/RNA/skin/combined-chrM.counts"),skip=1)
tab2 <- read.delim(paste0("data/raw_data_transposon2/20250601_DYQ_RNA/combined-chrM.counts"),skip = 1)
tab3 <- read.delim(paste0("data/raw_data_transposon2/20250622_HM_RNA/combined-chrM.counts"),skip=1)
tab <- merge(tab,tab2[,c(1,7:ncol(tab2))],by="Geneid")
tab <- merge(tab,tab3[,c(1,7:ncol(tab3))],by="Geneid")
rownames(tab) <- tab$Geneid
tab <- tab[,-1]
colnames <- colnames(tab)[6:length(tab)]
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|DYQ[0-9]+|HM[0-9]+).*"
new_colnames <- gsub(pattern, "\\1", colnames)
colnames(tab)[6:length(tab)] <- new_colnames
sorted_index <- order(new_colnames)
order_colnames <- new_colnames[sorted_index] 
counts <- tab[,order_colnames]   
# counts <- counts[,-which(colnames(counts)%in%c("DYQ184","DYQ185"))]
age <- c("young","old","young","young","old","old","young","young","old","old")
y= DGEList(counts=counts,group=age)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]

logCPMs <- cpm(y, log = TRUE)
pca <- prcomp(t(logCPMs))
to_plot <- data.frame(pca$x)
to_plot$sample_name <- rownames(to_plot)
# to_plot <- merge(to_plot,search_table,by="sample_name")
to_plot$age <- age
to_plot$rownames <- paste0(to_plot$sample_name,"-",to_plot$age)
to_plot$age <- factor(to_plot$age,levels=c("young","old"))
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
  ggtitle("Skin")
print(p)

diff <- read.csv("data/samples/RNA/skin/diff_expression_gene_change_filter_bar.csv")
diff <- diff[which(diff$Significant != "Stable"),]
genes <- diff$X
to_plot <- as.data.frame(logCPMs[which(rownames(logCPMs) %in% genes),])
annotation <- data.frame(sample=colnames(counts),age= c("young","old","young","young","old","old","young","young","old","old"))
rownames(annotation) <- annotation$sample
annotation <- annotation[,-1,drop =F]
pheatmap::pheatmap(to_plot,show_rownames = F,scale = "row",annotation_col = annotation)

# increase_GO <- read.csv("result/RNA/GO/table/Skin_increased_gene_GO.csv")
# increase_GO$geneID[1]
# genes <- strsplit(increase_GO$geneID[1], "/")[[1]]
# to_plot <- as.data.frame(logCPMs[which(rownames(logCPMs) %in% genes),])
# pheatmap::pheatmap(to_plot,show_rownames = T,scale = "row")
# 
decrease_GO <- read.csv("result/RNA/GO/table/Skin_decreased_gene_GO.csv")
# decrease_GO$geneID[1]
genes <- strsplit(decrease_GO$geneID[17], "/")[[1]]
to_plot <- as.data.frame(logCPMs[which(rownames(logCPMs) %in% genes),])
pheatmap::pheatmap(to_plot,show_rownames = T,scale = "row",annotation_col = annotation)



logCPMs_corrected <- limma::removeBatchEffect(logCPMs, batch = c("batch1","batch1","batch2","batch2","batch2","batch2","batch3","batch3","batch3","batch3"))
pca <- prcomp(t(logCPMs_corrected))
to_plot <- data.frame(pca$x)
to_plot$sample_name <- rownames(to_plot)
# to_plot <- merge(to_plot,search_table,by="sample_name")
to_plot$age <- age
to_plot$rownames <- paste0(to_plot$sample_name,"-",to_plot$age)
to_plot$age <- factor(to_plot$age,levels=c("young","old"))
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
  ggtitle("Skin")
print(p)
increase_GO <- read.csv("result/RNA/GO/table/Skin_increased_gene_GO.csv")
increase_GO$geneID[1]
genes <- strsplit(increase_GO$geneID[1], "/")[[1]]
to_plot <- as.data.frame(logCPMs_corrected[which(rownames(logCPMs_corrected) %in% genes),])
pheatmap::pheatmap(to_plot,show_rownames = T,scale = "row")

decrease_GO <- read.csv("result/RNA/GO/table/Skin_decreased_gene_GO.csv")
# decrease_GO$geneID[1]
genes <- strsplit(decrease_GO$geneID[17], "/")[[1]]
to_plot <- as.data.frame(logCPMs_corrected[which(rownames(logCPMs_corrected) %in% genes),])
pheatmap::pheatmap(to_plot,show_rownames = T,scale = "row",annotation_col = annotation)



diff <- read.csv("data/samples/RNA/skin/diff_expression_gene_change_filter_bar.csv")
diff <- diff[which(diff$Significant != "Stable"),]
genes <- diff$X
to_plot <- as.data.frame(logCPMs_corrected[which(rownames(logCPMs_corrected) %in% genes),])

pheatmap::pheatmap(to_plot,show_rownames = F,scale = "row")

# old <- counts[,c("DYQ188","LLX636","LLX637")]
# group <- c("WT","Mut","Mut")

old <- counts[,c("DYQ187","LLX634","LLX635")]
group <- c("batch2","batch1","batch1")

y= DGEList(counts=old,group=group)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
# y$samples$group <- factor(y$samples$group, levels=c("WT","Mut"))
y$samples$group <- factor(y$samples$group, levels=c("batch1","batch2"))
design <- model.matrix(~group, y$samples)
y <- calcNormFactors(y)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)
lrt = glmLRT(fit_tag, coef = 2)
tab<-tab[keep,]
out = cbind(cpm(y),lrt$table, "fdr"=p.adjust(lrt$table$PValue,method="BH"))
out$Significant <- ifelse(out$fdr< 0.05 & abs(out$logFC) >= 0, 
                          ifelse(out$logFC > 0, "Up", "Down"), "Stable")
colour=setNames(c("blue","grey","red"),c("Down","Stable","Up"))
ggplot(
  # 数据、映射、颜色
  out, aes(x = logFC, y = -log10(fdr))) +
  geom_point(aes(color = Significant), size=2) +
  scale_color_manual(values = colour) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (fdr)") +
  theme_bw()+
  theme(text = element_text(size = 20))+
  ggtitle("Skin Mut vs WT")+
  annotate("text", x = min(out$logFC), y = max(-log10(out$fdr)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
  annotate("text", x = max(out$logFC), y = max(-log10(out$fdr)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
genelist_up <- bitr(rownames(out)[which(out$Significant=="Up")],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_up_GO <- enrichGO( genelist_up$ENTREZID,#GO富集分析
                            OrgDb = GO_database,
                            keyType = "ENTREZID",#设定读取的gene ID类型
                            ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                            pvalueCutoff = 0.05,#设定p值阈值
                            qvalueCutoff = 0.05,#设定q值阈值
                            readable = T)
# barplot(genelist_up_GO,title = paste0("Mut vs WT Increased gene GO pathway"),label_format = 50)
barplot(genelist_up_GO,title = paste0("batch2 vs batch1 Increased gene GO pathway"),label_format = 50)
genelist_down <- bitr(rownames(out)[which(out$Significant=="Down")],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_down_GO <- enrichGO( genelist_down$ENTREZID,#GO富集分析
                            OrgDb = GO_database,
                            keyType = "ENTREZID",#设定读取的gene ID类型
                            ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                            pvalueCutoff = 0.05,#设定p值阈值
                            qvalueCutoff = 0.05,#设定q值阈值
                            readable = T)
# barplot(genelist_down_GO,title = paste0("Mut vs WT Decreased gene GO pathway"),label_format = 50)
barplot(genelist_down_GO,title = paste0("batch2 vs batch1 Decreased gene GO pathway"),label_format = 50)



young<- counts[,c("DYQ182","LLX634","LLX635")]
group <- c("New","Previous","Previous") 

y= DGEList(counts=young,group=group)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
y$samples$group <- factor(y$samples$group, levels=c("New","Previous"))
design <- model.matrix(~group, y$samples)
y <- calcNormFactors(y)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)
lrt = glmLRT(fit_tag, coef = 2)
tab<-tab[keep,]
out = cbind(cpm(y),lrt$table, "fdr"=p.adjust(lrt$table$PValue,method="BH"))
out$Significant <- ifelse(out$fdr< 0.05 & abs(out$logFC) >= 0, 
                          ifelse(out$logFC > 0, "Up", "Down"), "Stable")
colour=setNames(c("blue","grey","red"),c("Down","Stable","Up"))
ggplot(
  # 数据、映射、颜色
  out, aes(x = logFC, y = -log10(fdr))) +
  geom_point(aes(color = Significant), size=2) +
  scale_color_manual(values = colour) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (fdr)") +
  theme_bw()+
  theme(text = element_text(size = 20))+
  ggtitle("Skin Previous vs new")+
  annotate("text", x = min(out$logFC), y = max(-log10(out$fdr)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
  annotate("text", x = max(out$logFC), y = max(-log10(out$fdr)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)



df <- counts[,c("DYQ187","DYQ188")]
group <- c("young","old")

# old <- counts[,c("DYQ187","LLX634","LLX635")]
# group <- c("batch2","batch1","batch1") 

y= DGEList(counts=df,group=group)
keep = which(rowSums(cpm(y)>0)>=2)
y = y[keep,]
y$samples$group <- factor(y$samples$group, levels=c("young","old"))
# y$samples$group <- factor(y$samples$group, levels=c("batch1","batch2"))
design <- model.matrix(~group, y$samples)
y <- calcNormFactors(y)
bcv <- 0.1
lrt <- exactTest(y,dispersion = bcv^2, pair = c("young","old"))
tab<-tab[keep,]
out = cbind(cpm(y),lrt$table, "fdr"=p.adjust(lrt$table$PValue,method="BH"))
out$Significant <- ifelse(out$fdr< 0.05 & abs(out$logFC) >= 0, 
                          ifelse(out$logFC > 0, "Up", "Down"), "Stable")
colour=setNames(c("blue","grey","red"),c("Down","Stable","Up"))
ggplot(
  # 数据、映射、颜色
  out, aes(x = logFC, y = -log10(fdr))) +
  geom_point(aes(color = Significant), size=2) +
  scale_color_manual(values = colour) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (fdr)") +
  theme_bw()+
  theme(text = element_text(size = 20))+
  ggtitle("Skin old vs young")+
  annotate("text", x = min(out$logFC), y = max(-log10(out$fdr)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
  annotate("text", x = max(out$logFC), y = max(-log10(out$fdr)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
genelist_up <- bitr(rownames(out)[which(out$Significant=="Up")],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_up_GO <- enrichGO( genelist_up$ENTREZID,#GO富集分析
                            OrgDb = GO_database,
                            keyType = "ENTREZID",#设定读取的gene ID类型
                            ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                            pvalueCutoff = 0.05,#设定p值阈值
                            qvalueCutoff = 0.05,#设定q值阈值
                            readable = T)
barplot(genelist_up_GO,title = paste0("old vs young Increased gene GO pathway"),label_format = 50)
genelist_down <- bitr(rownames(out)[which(out$Significant=="Down")],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_down_GO <- enrichGO( genelist_down$ENTREZID,#GO富集分析
                              OrgDb = GO_database,
                              keyType = "ENTREZID",#设定读取的gene ID类型
                              ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                              pvalueCutoff = 0.05,#设定p值阈值
                              qvalueCutoff = 0.05,#设定q值阈值
                              readable = T)
barplot(genelist_down_GO,title = paste0("old vs young Decreased gene GO pathway"),label_format = 50)
