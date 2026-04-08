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
library(dplyr)
library(reshape2)
args <- commandArgs(trailingOnly=TRUE)
arg_list <- list()
for (arg in args) {
  split_arg <- strsplit(arg, "=")[[1]]
  arg_name <- split_arg[1]
  arg_value <- split_arg[2]
  arg_list[[arg_name]] <- arg_value
}
tissue <- arg_list[["tissue"]]
window_size = "1000"
gap_size = "3000"
e_value = "100"
common_peak_diff <- function(tissue,antibody,window_size,gap_size,e_value){
tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,
                        "_merge-W",window_size,"-G",gap_size,"-E",e_value,".counts"),skip=1)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
commonpeak <- readPeakFile(paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/intersect_",antibody,"-W1000-G3000-E100_unique_sort.bed"))
peakAnno <- annotatePeak(commonpeak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnno <- unique(as.data.frame(peakAnno))
colnames(peakAnno)[6]<-"Geneid"
common_peak_chr <- as.data.frame(table(peakAnno$seqnames))
peak_chr <- as.data.frame(table(tab$Chr))
percent_chr <- data.frame(Chr=peak_chr$Var1,Freq=common_peak_chr$Freq/peak_chr$Freq)
percent_chr$Freq_not_overlap <- 1- percent_chr$Freq
percent_chr_melt <- melt(percent_chr)  
percent_chr_melt$Chr <- factor(percent_chr_melt$Chr,c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
                                                         "chr9","chr10","chr11","chr12","chr13","chr14","chr15",
                                                         "chr16","chr17","chr18","chr19","chrX","chrY"))
ggplot(percent_chr_melt, aes( x = Chr, weight = value, fill = variable))+
  geom_bar( position = "stack")+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0("result/",tissue,"/H3K27me3_H3K9me3_relationship/",antibody,"_commonpeak_diff_percent.png"),width=10,height = 10)
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
# fit <- glmQLFit(y, design)
# fit  <- glmQLFTest(fit, coef = 3)
out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
            "PValue.old-young"=lrt$table$PValue,"FDR.old-young"= p.adjust(lrt$table$PValue,method="BH"),
            "LogFC.old-young"=lrt$table$logFC)

commonout <- out[which(out$Geneid %in%peakAnno$Geneid),]
commonout$annotation <- peakAnno$annotation[which(peakAnno$Geneid %in% commonout$Geneid)]
commonout$ens_id <- peakAnno$geneId[which(peakAnno$Geneid %in% commonout$Geneid)]
commonout$symbol <- peakAnno$SYMBOL[which(peakAnno$Geneid %in% commonout$Geneid)]
commonout$Significant <- ifelse(commonout$`FDR.old-young` < 0.05 & abs(commonout$`LogFC.old-young`) >= 0, 
                          ifelse(commonout$`LogFC.old-young` > 0, "Up", "Down"), "Stable")

colour<- list(c("grey"),c("grey","red"),c("blue","grey","red"))
ggplot(
  # 数据、映射、颜色
  commonout, aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`))) +
  geom_point(aes(color = Significant), size=2) +
  scale_color_manual(values = colour[[nrow(as.data.frame(table(commonout$Significant)))]]) +
  # scale_color_manual(values = c("blue","grey")) +
  # 注释
  # geom_text_repel(
  #   data = subset(commonout,`FDR.FA-noFA` < 0.05 & abs(commonout$logFC) >= 1),
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
ggsave(paste0("result/",tissue,"/H3K27me3_H3K9me3_relationship/",antibody,"_commonpeak_diff_volcano.png"),width=10,height = 10)
outup <- commonout[which(commonout$Significant=="Up"),]
outdown <- commonout[which(commonout$Significant=="Down"),]
write.table(outdown[,c("Chr","Start","End","Geneid")], file=paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/intersect_",antibody,"-W1000-G3000-E100_down.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) 
write.table(outup[,c("Chr","Start","End","Geneid")], file=paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/intersect_",antibody,"-W1000-G3000-E100_up.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) 
}

common_peak_diff(tissue,"H3K9me3",window_size,gap_size,e_value)
common_peak_diff(tissue,"H3K27me3",window_size,gap_size,e_value)
