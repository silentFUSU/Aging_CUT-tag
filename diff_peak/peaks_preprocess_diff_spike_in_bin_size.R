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
library(limma)
GO_database <- 'org.Mm.eg.db'
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
tab <- read.delim(paste0("data/raw_data/LLY001-LLY004/1kb_bins.counts"),skip=1)
tab <- read.delim(paste0("data/raw_data/LLY001-LLY004/LLY001-LLY004_macs_narrowpeak.counts"),skip=1)
counts = tab[,c(7:10)]
colnames(counts) = c("test_1","test_2","control_1","control_2")
group =c("test","test","control","control")

design <- model.matrix(~ 0 + group)
lvls = levels(factor(group))
colnames(design) = lvls
len = length(lvls)
myContrasts = c("test-control")
contrast.matrix = eval(as.call(c(as.symbol("makeContrasts"),as.list(myContrasts),levels=list(design))))

### use spike in 
spike_in_factor <- read.table("data/raw_data/LLY001-LLY004/all_sample.qc.txt",header = T)
spike_in_factor <- spike_in_factor$SpikeIn.
y= DGEList(counts=counts,group=group)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]

norm.factors <- spike_in_factor / y$samples$lib.size
norm.factors <- norm.factors / prod(norm.factors)^(1/length(norm.factors))
y$samples$norm.factors <- norm.factors
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)
lrt = glmLRT(fit_tag,contrast = contrast.matrix)
tab<-tab[keep,]
out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
            "PValue.test-control"=lrt$table$PValue,"FDR.test-control"= p.adjust(lrt$table$PValue,method="BH"),
            "LogFC.test-control"=lrt$table$logFC)
out$Significant <- ifelse(out$`FDR.test-control` < 0.05 & abs(out$`LogFC.test-control`) >= 0, 
                          ifelse(out$`LogFC.test-control` > 0, "Up", "Down"), "Stable")
colour<- list(c("grey"),c("grey","red"),c("blue","grey","red"))
ggplot(
  # 数据、映射、颜色
  out, aes(x = `LogFC.test-control`, y = -log10(`FDR.test-control`))) +
  geom_point(aes(color = Significant), size=2) +
  scale_color_manual(values = colour[[nrow(as.data.frame(table(out$Significant)))]]) +
  # scale_color_manual(values = c("blue","grey")) +
  # 注释
  # geom_text_repel(
  #   data = subset(out,`FDR.FA-noFA` < 0.05 & abs(out$logFC) >= 1),
  #   aes(label = Geneid),
  #   size = 5,max.overlaps = 100,
  #   box.padding = unit(0.35, "lines"),
  #   point.padding = unit(0.3, "lines")) +
  # 辅助线
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (p-value)") +
  # xlim(-3,3)+
  # 图例
  theme_bw()+
  theme(text = element_text(size = 20))+
  annotate("text", x = min(out$`LogFC.test-control`), y = max(-log10(out$`FDR.test-control`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
  annotate("text", x = max(out$`LogFC.test-control`), y = max(-log10(out$`FDR.test-control`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)



#use calcnromFactors
y= DGEList(counts=counts,group=group)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
y <- calcNormFactors(y)
calcNormFactor <- y$samples$norm.factors
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)
lrt = glmLRT(fit_tag,contrast = contrast.matrix)
tab<-tab[keep,]
out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
            "PValue.test-control"=lrt$table$PValue,"FDR.test-control"= p.adjust(lrt$table$PValue,method="BH"),
            "LogFC.test-control"=lrt$table$logFC)
out$Significant <- ifelse(out$`FDR.test-control` < 0.05 & abs(out$`LogFC.test-control`) >= 0, 
                          ifelse(out$`LogFC.test-control` > 0, "Up", "Down"), "Stable")
colour<- list(c("grey"),c("grey","red"),c("blue","grey","red"))
ggplot(
  # 数据、映射、颜色
  out, aes(x = `LogFC.test-control`, y = -log10(`FDR.test-control`))) +
  geom_point(aes(color = Significant), size=2) +
  scale_color_manual(values = colour[[nrow(as.data.frame(table(out$Significant)))]]) +
  # scale_color_manual(values = c("blue","grey")) +
  # 注释
  # geom_text_repel(
  #   data = subset(out,`FDR.FA-noFA` < 0.05 & abs(out$logFC) >= 1),
  #   aes(label = Geneid),
  #   size = 5,max.overlaps = 100,
  #   box.padding = unit(0.35, "lines"),
  #   point.padding = unit(0.3, "lines")) +
  # 辅助线
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (p-value)") +
  # xlim(-3,3)+
  # 图例
  theme_bw()+
  theme(text = element_text(size = 20))+
  annotate("text", x = min(out$`LogFC.test-control`), y = max(-log10(out$`FDR.test-control`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
  annotate("text", x = max(out$`LogFC.test-control`), y = max(-log10(out$`FDR.test-control`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  increase <- out[which(out$Significant=="Up"),]
  peak_obj <- GRanges(seqnames = increase$Chr,   
                      ranges = IRanges(start = increase$Start, end = increase$End))
  peak_anno <- annotatePeak(peak_obj, tssRegion=c(-3000, 3000),
                            TxDb=txdb, annoDb="org.Mm.eg.db")
  peak_anno <- unique(as.data.frame(peak_anno))
  genelist_up <- bitr(peak_anno$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  genelist_up_GO <- enrichGO( genelist_up$ENTREZID,#GO富集分析
                              OrgDb = GO_database,
                              keyType = "ENTREZID",#设定读取的gene ID类型
                              ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                              pvalueCutoff = 0.05,#设定p值阈值
                              qvalueCutoff = 0.05,#设定q值阈值
                              readable = T)
  barplot(genelist_up_GO)
  result <- genelist_up_GO@result
  write.csv(result,"data/raw_data/LLY001-LLY004/increase_GO.csv")
  decrease <- out[which(out$Significant=="Down"),]
  peak_obj <- GRanges(seqnames = decrease$Chr,   
                      ranges = IRanges(start = decrease$Start, end = decrease$End))
  peak_anno <- annotatePeak(peak_obj, tssRegion=c(-3000, 3000),
                            TxDb=txdb, annoDb="org.Mm.eg.db")
  peak_anno <- unique(as.data.frame(peak_anno))
  genelist_down <- bitr(peak_anno$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  genelist_down_GO <- enrichGO( genelist_down$ENTREZID,#GO富集分析
                              OrgDb = GO_database,
                              keyType = "ENTREZID",#设定读取的gene ID类型
                              ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                              pvalueCutoff = 0.05,#设定p值阈值
                              qvalueCutoff = 0.05,#设定q值阈值
                              readable = T)
  barplot(genelist_down_GO)
  result <- genelist_down_GO@result
  write.csv(result,"data/raw_data/LLY001-LLY004/decrease_GO.csv")
### use yuxiao 1/scalefactor
spike_in_factor <- read.table("data/raw_data/LLY001-LLY004/all_sample.qc.txt",header = T)
spike_in_factor <- spike_in_factor$ScaleFactor
norm.factors <- 1/spike_in_factor
y= DGEList(counts=counts,group=group)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
y$samples$norm.factors <- norm.factors
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)
lrt = glmLRT(fit_tag,contrast = contrast.matrix)
tab<-tab[keep,]
out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
            "PValue.test-control"=lrt$table$PValue,"FDR.test-control"= p.adjust(lrt$table$PValue,method="BH"),
            "LogFC.test-control"=lrt$table$logFC)
out$Significant <- ifelse(out$`FDR.test-control` < 0.05 & abs(out$`LogFC.test-control`) >= 0, 
                          ifelse(out$`LogFC.test-control` > 0, "Up", "Down"), "Stable")
colour<- list(c("grey"),c("grey","red"),c("blue","grey","red"))
ggplot(
  # 数据、映射、颜色
  out, aes(x = `LogFC.test-control`, y = -log10(`FDR.test-control`))) +
  geom_point(aes(color = Significant), size=2) +
  scale_color_manual(values = colour[[nrow(as.data.frame(table(out$Significant)))]]) +
  # scale_color_manual(values = c("blue","grey")) +
  # 注释
  # geom_text_repel(
  #   data = subset(out,`FDR.FA-noFA` < 0.05 & abs(out$logFC) >= 1),
  #   aes(label = Geneid),
  #   size = 5,max.overlaps = 100,
  #   box.padding = unit(0.35, "lines"),
  #   point.padding = unit(0.3, "lines")) +
  # 辅助线
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (p-value)") +
  # xlim(-3,3)+
  # 图例
  theme_bw()+
  theme(text = element_text(size = 20))


simplify_ratio <- function(numbers) {  
  gcd_value <- function(a, b) {  
    if (b == 0) return(a)  
    return(gcd_value(b, a %% b))  
  }  
  
  gcd <- Reduce(gcd_value, numbers)  
  simplified <- numbers / gcd  
  
  return(simplified)  
} 
simplify_ratio(norm.factors)
simplify_ratio(c(0.691,1.196,1.209,1))
