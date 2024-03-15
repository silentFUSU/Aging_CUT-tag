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
antibody = "H3K27ac"
tissue = "testis"
tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_macs_narrowpeak.counts"),skip=1)
counts = tab[,c(7:10)]
rownames(counts)= tab$Geneid
colnames(counts) = c("young_1","old_1","young_2","old_2")
group =c("batch1","batch1","batch2","batch2")
design = model.matrix(~0+group)
lvls = levels(factor(group))
colnames(design) = lvls
len = length(lvls)
myContrasts = c("batch2-batch1")
contrast.matrix = eval(as.call(c(as.symbol("makeContrasts"),as.list(myContrasts),levels=list(design))))
counts_cpm <- as.data.frame(cpm(counts))
# keep = which(rowSums(counts_cpm>1)>=2)
# counts_cpm_normalize <- counts_cpm
# counts_cpm_normalize[,1:2] <- counts_cpm[, 1:2] / rowSums(counts_cpm[, 1:2])
# counts_cpm_normalize[,3:4] <- counts_cpm[, 3:4] / rowSums(counts_cpm[, 3:4])
# counts_cpm_normalize = counts_cpm_normalize[keep,]
# y<-estimateCommonDisp(y)
# y<-estimateGLMTagwiseDisp(y,design)
# fit_tag = glmFit(y,design)

y= DGEList(counts=counts,group=group)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
y =  calcNormFactors(y)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)
lrt = glmLRT(fit_tag, contrast = contrast.matrix)
qBH = p.adjust(lrt$table$PValue,method="BH")
tab<-tab[keep,]
out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion))
for (i in 1:length(myContrasts)){
  lrt = glmLRT(fit_tag, contrast = contrast.matrix[,i])
  out[,paste0("PValue.",myContrasts[i])] = lrt$table$PValue
  out[,paste0("FDR.",myContrasts[i])] = p.adjust(lrt$table$PValue,method="BH")
  out[,paste0("LogFC.",myContrasts[i])] = lrt$table$logFC
}

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
# peak <- readPeakFile(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_macs2_mergecounts.bed"))
peak <- readPeakFile(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_macs_narrowpeak.bed"))
peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnno<-unique(as.data.frame(peakAnno))

out$annotation <- peakAnno$annotation[which(peakAnno$V4 %in% out$Geneid)]
out$ens_id <- peakAnno$geneId[which(peakAnno$V4 %in% out$Geneid)]
out$symbol <- peakAnno$SYMBOL[which(peakAnno$V4 %in% out$Geneid)]

out$Significant <- ifelse(out$`FDR.batch2-batch1` < 0.05 & abs(out$`LogFC.batch2-batch1`) >= 0, 
                          ifelse(out$`LogFC.batch2-batch1` > 0, "Up", "Down"), "Stable")
colour<- list(c("grey"),c("grey","red"),c("blue","grey","red"))
ggplot(
  # 数据、映射、颜色
  out, aes(x = `LogFC.batch2-batch1`, y = -log10(`FDR.batch2-batch1`))) +
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
