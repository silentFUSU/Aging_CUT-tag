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
antibody = "H3K9me3"
tissue = "colon"
window_size = "1000"
gap_size = "3000"
e_value = "100"
peak_preprocess_sicer <- function(tissue,antibody,window_size,gap_size,e_value){
  tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,
                          "_merge-W",window_size,"-G",gap_size,"-E",e_value,".counts"),skip=1)
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
  ggsave(paste0("result/",tissue,"/pca/",antibody,"_merge-W",window_size,"-G",gap_size,"-E",e_value,"_pca_plot_before_remove_batch_effect.png"),width = 5,height =5)
  
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
  ggsave(paste0("result/",tissue,"/pca/",antibody,"_merge-W",window_size,"-G",gap_size,"-E",e_value,"_pca_plot_after_remove_batch_effect.png"),width = 5,height =5)
  
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
  
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  # peak <- readPeakFile(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_macs2_mergecounts.bed"))
  peak <- readPeakFile(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_merge-W",window_size,"-G",gap_size,"-E100.bed"))
  peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                           TxDb=txdb, annoDb="org.Mm.eg.db")
  peakAnno<-unique(as.data.frame(peakAnno))
  
  out$annotation <- peakAnno$annotation[which(peakAnno$V4 %in% out$Geneid)]
  out$ens_id <- peakAnno$geneId[which(peakAnno$V4 %in% out$Geneid)]
  out$symbol <- peakAnno$SYMBOL[which(peakAnno$V4 %in% out$Geneid)]
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
  ggsave(paste0("result/",tissue,"/diffpeaks/",antibody,"_merge-W",window_size,"-G",gap_size,"-E",e_value,"_volcano_plot_after_remove_batch_effect.png"),width = 10,height = 10)
  write.csv(out,paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_merge-W",window_size,"-G",gap_size,"-E100_diff_after_remove_batch_effect.csv"),row.names = F)
  outup <- out[which(out$Significant=="Up"),]
  outdown <- out[which(out$Significant=="Down"),]
  write.table(outdown[,c("Chr","Start","End","Geneid")], file=paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_merge-W",window_size,"-G",gap_size,"-E",e_value,"_down.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) 
  write.table(outup[,c("Chr","Start","End","Geneid")], file=paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_merge-W",window_size,"-G",gap_size,"-E",e_value,"_up.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) 
}

peak_preprocess_macs <- function(tissue,antibody){
  tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_macs_narrowpeak.counts"),skip=1)
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
  ggsave(paste0("result/",tissue,"/pca/",antibody,"_pca_plot_before_remove_batch_effect.png"),width = 5,height =5)
  
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
  ggsave(paste0("result/",tissue,"/pca/",antibody,"_pca_plot_after_remove_batch_effect.png"),width = 5,height =5)
  
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
  
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  # peak <- readPeakFile(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_macs2_mergecounts.bed"))
  peak <- readPeakFile(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_macs_narrowpeak.bed"))
  peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                           TxDb=txdb, annoDb="org.Mm.eg.db")
  peakAnno<-unique(as.data.frame(peakAnno))
  
  out$annotation <- peakAnno$annotation[which(peakAnno$V4 %in% out$Geneid)]
  out$ens_id <- peakAnno$geneId[which(peakAnno$V4 %in% out$Geneid)]
  out$symbol <- peakAnno$SYMBOL[which(peakAnno$V4 %in% out$Geneid)]
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
  ggsave(paste0("result/",tissue,"/diffpeaks/",antibody,"_volcano_plot_after_remove_batch_effect.png"),width = 10,height = 10)
  write.csv(out,paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_macs_narrowpeak_diff_after_remove_batch_effect.csv"),row.names = F)
  outup <- out[which(out$Significant=="Up"),]
  outdown <- out[which(out$Significant=="Down"),]
  write.table(outdown[,c("Chr","Start","End","Geneid")], file=paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_macs_narrowpeak_down.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) 
  write.table(outup[,c("Chr","Start","End","Geneid")], file=paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_macs_narrowpeak_up.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) 
  }

# tissues <- c("brain","liver","testis","colon","kidney")
# tissues <- c("lung","spleen","muscle")
tissues <- c("pancreas")
antibodys <- c("H3K27me3","H3K9me3","H3K36me3")
for (j in c(1:length(tissues))){
  tissue <- tissues[j]
  for (i in c(1:length(antibodys))){
    antibody <- antibodys[i]
    peak_preprocess_sicer(tissue,antibody,window_size,gap_size,e_value)
  }
}

# tissues <- c("brain","liver","testis","colon","kidney")
antibodys <- c("H3K27ac","H3K4me3","H3K4me1")
# tissue <- "colon"
# antibodys <- c("H3K4me3")
for (j in c(1:length(tissues))){
  tissue <- tissues[j]
  for (i in c(1:length(antibodys))){
    antibody <- antibodys[i]
    peak_preprocess_macs(tissue,antibody)
  }
}
command <- paste("bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/bcv_report_after_remove_batch_effect.sh /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/", tissue, sep=" ") 
system(command, intern = TRUE)
bcv <- read.table(paste0("data/samples/",tissue,"/bcv_report_after_remove_batch_effect.txt"), header = FALSE, 
                  col.names = c("bcv", "antibody"),
                  colClasses = c("numeric", "character"))
ggplot(bcv,aes(x=antibody,y=bcv,fill = antibody))+geom_boxplot()+
  stat_summary(fun.y = mean, geom = "point", shape = 23, size=4)+
  # geom_jitter(shape=16, position = position_jitter(0.2))+
  scale_fill_brewer(palette="Set3")+guides(fill = FALSE) +
  ggtitle(paste0(tissue," BCV each antibody after remove batch effect"))+theme_bw()+theme(text = element_text(size = 18))+xlab("")+labs(fill = "", color = "") 
ggsave(paste0("result/",tissue,"/diffpeaks/bcv_report_after_remove_batch_effect.png"),width = 10,height = 8)
