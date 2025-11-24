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
tissue<-"testis"
state<-"E13"
antibody <-"H3K9me3"
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
peak <- readPeakFile(paste0("data/samples/",tissue,"/combined_analysis_enhancer/bed/15/merged_segments_",state,".bed"))
peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnno <- unique(as.data.frame(peakAnno))
colnames(peakAnno)[6]<-"Geneid"

tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_segments_",state,".counts"),skip=1)
counts = tab[,c(7:10)]
rownames(counts)= tab$Geneid
colnames(counts) = c("young_1","old_1","young_2","old_2")
group =c("young","old","young","old")
y= DGEList(counts=counts,group=group)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
y$samples$batch <- rep(c(rep("batch1", 2), rep("batch2", 2)), 1)
y$samples$year <- group
y$samples$year <- factor(y$samples$year,c("young","old"))
y <- calcNormFactors(y)
design <- model.matrix(~batch+year, y$samples)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)
lrt = glmLRT(fit_tag, coef = 3)
tab<-tab[keep,]
out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
              "PValue.old-young"=lrt$table$PValue,"FDR.old-young"= p.adjust(lrt$table$PValue,method="BH"),
              "LogFC.old-young"=lrt$table$logFC)
out <- merge(out, peakAnno[,c(6:7,16:17)], by = "Geneid") 
out$Significant <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$`LogFC.old-young`) >= 0, 
                            ifelse(out$`LogFC.old-young` > 0, "Up", "Down"), "Stable")
colour<- list(c("grey"),c("grey","red"),c("blue","grey","red"))
ggplot(
    # 数据、映射、颜色
    out, aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`))) +
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

H3K9me3_decrease <- out$Geneid[which(out$Significant=="Down")]
H3K9me3_increase <- out$Geneid[which(out$Significant=="Up")]
antibody <- "H3K27me3"
tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_segments_",state,".counts"),skip=1)
counts = tab[,c(7:10)]
rownames(counts)= tab$Geneid
colnames(counts) = c("young_1","old_1","young_2","old_2")
group =c("young","old","young","old")
y= DGEList(counts=counts,group=group)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
y$samples$batch <- rep(c(rep("batch1", 2), rep("batch2", 2)), 1)
y$samples$year <- group
y$samples$year <- factor(y$samples$year,c("young","old"))
y <- calcNormFactors(y)
design <- model.matrix(~batch+year, y$samples)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)
lrt = glmLRT(fit_tag, coef = 3)
tab<-tab[keep,]
out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
            "PValue.old-young"=lrt$table$PValue,"FDR.old-young"= p.adjust(lrt$table$PValue,method="BH"),
            "LogFC.old-young"=lrt$table$logFC)
out <- merge(out, peakAnno[,c(6:7,16:17)], by = "Geneid") 
out$Significant <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$`LogFC.old-young`) >= 0, 
                          ifelse(out$`LogFC.old-young` > 0, "Up", "Down"), "Stable")
out$H3K9me3_condition <-"stable"
out$H3K9me3_condition[which(out$Geneid %in% H3K9me3_increase)]<-"Up"
out$H3K9me3_condition[which(out$Geneid %in% H3K9me3_decrease)]<-"Down"

ggplot(
  # 数据、映射、颜色
  out, aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`))) +
  geom_point(aes(color = H3K9me3_condition), size=2) +
  scale_color_manual(values = colour[[nrow(as.data.frame(table(out$H3K9me3_condition)))]]) +
  geom_point(data=out[which(out$H3K9me3_condition=="Down"),], mapping=aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`)), size=2, color="blue")+
  geom_point(data=out[which(out$H3K9me3_condition=="Up"),], mapping=aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`)), size=2, color="red")+
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
ggplot(
  # 数据、映射、颜色
  out, aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`))) +
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
