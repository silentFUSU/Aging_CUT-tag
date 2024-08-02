rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(DiffBind)
library(tidyverse)
library(ggplot2)
library(ChIPseeker)
library(EnsDb.Mmusculus.v79)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
tissue = "brain"
############### plan A ###################
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
peak <- readPeakFile(paste0("data/samples/",tissue,"/H3K27ac/bed/H3K27ac_macs2_mergecounts.bed"))
# peak <- readPeakFile(paste0("data/samples/",tissue,"/",antibody,"/bed/W1000-G3000-W5000-G10000-intersect.bed"))
peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnno<-unique(as.data.frame(peakAnno))
peakAnno<-peakAnno[-which(str_detect(peakAnno$annotation,"Promoter")),]
write.table(peakAnno[,c(1:3,6)], file=paste0("data/samples/",tissue,"/H3K27ac/bed/H3K27ac_macs2_mergecounts_without_promoter.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) 

tab = read.delim(paste0("data/samples/",tissue,"/H3K4me1_H3K4me3_merge.counts"),skip=1)

counts = tab[,c(7:14)]
rownames(counts)= tab$Geneid
colnames(counts) = c("me1_1","me1_2","me1_3","me1_4","me3_1","me3_2","me3_3","me3_4")
group =c("me1","me1","me1","me1","me3","me3","me3","me3")
design = model.matrix(~0+group)
lvls = levels(factor(group))
colnames(design) = lvls
len = length(lvls)
myContrasts = c("me1-me3")
contrast.matrix = eval(as.call(c(as.symbol("makeContrasts"),as.list(myContrasts),levels=list(design))))
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
out = cbind(tab[,1:6],cpm(y),lrt$table,qBH)
out[,c(7:14)] = sweep(out[,c(7:14)],1,out$Length,'/')*1000
for (i in 1:length(myContrasts)){
  lrt = glmLRT(fit_tag, contrast = contrast.matrix[,i])
  out[,paste0("PValue.",myContrasts[i])] = lrt$table$PValue
  out[,paste0("FDR.",myContrasts[i])] = p.adjust(lrt$table$PValue,method="BH")
}
out$Significant <- ifelse(out$`FDR.me1-me3` < 0.05 & abs(out$logFC) >= 0, 
                          ifelse(out$logFC > 0, "Up", "Down"), "Stable")
colour<- list(c("grey"),c("grey","red"),c("blue","grey","red"))
ggplot(
  # 数据、映射、颜色
  out, aes(x = logFC, y = -log10(`FDR.me1-me3`))) +
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
  geom_vline(xintercept=c(-0.7,0.7),lty=4,col="black",lwd=0.8) +
  # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (p-value)") +
  # xlim(-3,3)+
  # 图例
  theme(legend.position = "bottom",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme_bw()
outup <- out[which(out$`FDR.me1-me3`<0.05 & out$logFC>1.5),]
write.table(outup[,c("Chr","Start","End","Geneid")], file=paste0("data/samples/",tissue,"/H3K4me1_H3K4me3_logFC15_FDR005.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) 


############### plan B ###################
peak <- readPeakFile(paste0("data/samples/",tissue,"/combined_analysis2_enhancer/H3K4me1_H3K4me3_logFC15_FDR005.bed"))
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnno<-unique(as.data.frame(peakAnno))
peakAnno<-peakAnno[-which(str_detect(peakAnno$annotation,"Promoter")),]
write.table(peakAnno[,c(1:3,6)], file=paste0("data/samples/",tissue,"/combined_analysis2_enhancer/H3K4me1_H3K4me3_logFC15_FDR005_without_promoter.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) 

data <- read.table(paste0("data/samples/",tissue,"/combined_analysis2_enhancer/matrix/H3K4me1_H3K4me3_logFC15_FDR005_without_promoter.matrix"), header = F ,row.names = 1,) 
