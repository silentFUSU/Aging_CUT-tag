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
tissue <-"muscle"
bin_size="1kb"

tab = read.delim(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_",bin_size,"_bins.counts"),skip=1)
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
logCPMs <- cpm(y, log = TRUE)
rv <- apply(logCPMs, 1, var)
o <- order(rv, decreasing=TRUE)
# top1000 <- head(o, 1000)
# logCPM_top1000 <- logCPMs[top1000,]
pca <- prcomp(t(logCPMs))
to_plot <- data.frame(pca$x, y$samples)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
ggplot(to_plot, aes(x=PC1, y=PC2, color=batch, shape=year)) + 
  geom_point(size=5) +theme_bw()+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))
# ggsave(paste0("result/",tissue,"/pca/",antibody,"_",bin_size,"_bins.png"),width = 5,height =5)

batch <- factor(y$samples$batch)
logCPMs_corrected <- limma::removeBatchEffect(logCPMs, batch = batch)
# logCPM_corrected_top1000 <- logCPMs_corrected[top1000,]
pca <- prcomp(t(logCPMs_corrected))
to_plot <- data.frame(pca$x, y$samples)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))

ggplot(to_plot, aes(x=PC1, y=PC2, color=batch, shape=year)) + 
  geom_point(size=5) +theme_bw()+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))
# ggsave(paste0("result/",tissue,"/pca/",antibody,"_",bin_size,"_bins_pca_plot_after_remove_batch_effect.png"),width = 5,height =5)

design <- model.matrix(~batch+year, y$samples)
# y <- estimateDisp(y, design)

y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)
lrt = glmLRT(fit_tag, coef = 3)
tab<-tab[keep,]
# fit <- glmQLFit(y, design)
# fit  <- glmQLFTest(fit, coef = 3)
out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
            "PValue.old-young"=lrt$table$PValue,"FDR.old-young"= p.adjust(lrt$table$PValue,method="BH"),
            "LogFC.old-young"=lrt$table$logFC)

# txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
# # peak <- readPeakFile(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_macs2_mergecounts.bed"))
# peak <- readPeakFile(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_merge-W",window_size,"-G",gap_size,"-E100.bed"))
# peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
#                          TxDb=txdb, annoDb="org.Mm.eg.db")
# peakAnno<-unique(as.data.frame(peakAnno))
# 
# out$annotation <- peakAnno$annotation[which(peakAnno$V4 %in% out$Geneid)]
# out$ens_id <- peakAnno$geneId[which(peakAnno$V4 %in% out$Geneid)]
# out$symbol <- peakAnno$SYMBOL[which(peakAnno$V4 %in% out$Geneid)]
out$Significant <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$`LogFC.old-young`) >= 0, 
                          ifelse(out$`LogFC.old-young` > 0, "Up", "Down"), "Stable")
colour<- list(c("grey"),c("grey","red"),c("blue","grey","red"))
quadrants<-c("first","second","third","fourth")
quadrants_peaks <-list()
for (i in c(1:length(quadrants))){
  quadrant <- quadrants[i]
  peaks<- read.delim(paste0("data/samples/ATAC/",tissue,"/ATAC/bed/H3K27me3_H3K9me3_all_significant_",bin_size,"_",quadrant,"_quadrant_unique_sort.bed"),header = F)
  quadrants_peaks[[i]] <- peaks$V4
  names(quadrants_peaks)[i] <- quadrant
}
ggplot() +
  geom_point(out, mapping=aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`)), color = "grey",alpha=0.2) +
  geom_point(data=out[which(out$Geneid %in%quadrants_peaks[["fourth"]]),], mapping=aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`)),color ="#48466d",alpha=0.7)+
  geom_point(data=out[which(out$Geneid %in%quadrants_peaks[["first"]]),], mapping=aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`)),color = "#00b8a9",alpha=0.7)+
  geom_point(data=out[which(out$Geneid %in%quadrants_peaks[["third"]]),], mapping=aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`)),color ="#ffde7d",alpha=0.7)+
  geom_point(data=out[which(out$Geneid %in%quadrants_peaks[["second"]]),], mapping=aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`)),color ="#f6416c",alpha=0.7)+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="red",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (FDR)") +
  # xlim(-3,3)+
  # 图例
  theme_bw()+
  theme(text = element_text(size = 20))

ggplot() +
  geom_point(out, mapping=aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`)), color = "grey",alpha=0.2) +
  geom_point(data=out[which(out$Geneid %in%quadrants_peaks[["fourth"]] & out$`FDR.old-young`<0.05),], mapping=aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`)),color ="#48466d",alpha=0.7)+
  geom_point(data=out[which(out$Geneid %in%quadrants_peaks[["first"]] & out$`FDR.old-young`<0.05),], mapping=aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`)),color = "#00b8a9",alpha=0.7)+
  geom_point(data=out[which(out$Geneid %in%quadrants_peaks[["third"]] & out$`FDR.old-young`<0.05),], mapping=aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`)),color ="#ffde7d",alpha=0.7)+
  geom_point(data=out[which(out$Geneid %in%quadrants_peaks[["second"]] & out$`FDR.old-young`<0.05),], mapping=aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`)),color ="#f6416c",alpha=0.7)+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="red",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (FDR)") +
  # xlim(-3,3)+
  # 图例
  theme_bw()+
  theme(text = element_text(size = 20))
