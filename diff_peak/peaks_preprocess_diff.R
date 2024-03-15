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
antibody = "H3K27me3"
tissue = "testis"
window_size = "1000"
gap_size = "3000"
e_value = "100"
peak_preprocess_sicer <- function(tissue,antibody,window_size,gap_size,e_value){
  
  tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,
                          "_merge-W",window_size,"-G",gap_size,"-E",e_value,".counts"),skip=1)
  # tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_macs2_mergecounts.counts"),skip=1)
  # tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_PePr_old_increase_mergecounts.counts"),skip=1)
  # tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/W1000-G3000-W5000-G10000-intersect.counts"),skip=1)
  counts = tab[,c(7:10)]
  rownames(counts)= tab$Geneid
  colnames(counts) = c("young_1","old_1","young_2","old_2")
  counts <- counts[,c("young_1","young_2","old_1","old_2")]
  group =c("young","young","old","old")
  design = model.matrix(~0+group)
  lvls = levels(factor(group))
  colnames(design) = lvls
  len = length(lvls)
  myContrasts = c("old-young")
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
  out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion))
  # out[,c(7:10)] = sweep(out[,c(7:10)],1,out$Length,'/')*1000
  for (i in 1:length(myContrasts)){
    lrt = glmLRT(fit_tag, contrast = contrast.matrix[,i])
    out[,paste0("PValue.",myContrasts[i])] = lrt$table$PValue
    out[,paste0("FDR.",myContrasts[i])] = p.adjust(lrt$table$PValue,method="BH")
    out[,paste0("LogFC.",myContrasts[i])] = lrt$table$logFC
  }
  
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
  
  ggsave(paste0("result/",tissue,"/diffpeaks/",antibody,"_merge-W",window_size,"-G",gap_size,"-E",e_value,"_volcano_plot.png"),width = 10,height = 10)
  # out$Significant_pvalue <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$logFC) >= 0, 
  #                           ifelse(out$logFC > 0, "Up", "Down"), "Stable")
  out$Significant_pvalue <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$`LogFC.old-young`) >= 0, 
                                   ifelse(out$`LogFC.old-young` > 0, "Up", "Down"), "Stable")
  outup <- out[which(out$Significant_pvalue=="Up"),]
  outdown <- out[which(out$Significant_pvalue=="Down"),]
  
  # ggplot(outup, aes(x = Length, fill = Significant)) +  
  #   geom_density(alpha = 0.5) +  
  #   labs(title = "Distribution of Two Groups", x = "Value", y = "Density") +theme_bw()+theme(text = element_text(size = 18))
  # ggplot(outdown, aes(x = Length, fill = Significant)) +  
  #   geom_density(alpha = 0.5) +  
  #   labs(title = "Distribution of Two Groups", x = "Value", y = "Density") +theme_bw()+theme(text = element_text(size = 18))
  
  # ggplot() +
  #   geom_point(data=out, mapping=aes(logCPM, `LogFC.old-young`), size=2, color="grey") +
  #   geom_point(data=outup, mapping=aes(logCPM, `LogFC.old-young`), size=2, color="red") +
  #   geom_point(data=outdown, mapping=aes(logCPM, `LogFC.old-young`), size=2, color="blue") +
  # labs(x="log2(CPM)",
  #      y="log2(FC)") +
  #   # 图例
  #   theme(legend.position = "bottom",panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank())+
  #   theme_bw()
  write.table(outdown[,c("Chr","Start","End","Geneid")], file=paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_merge-W",window_size,"-G",gap_size,"-E",e_value,"_down.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) 
  write.table(outup[,c("Chr","Start","End","Geneid")], file=paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_merge-W",window_size,"-G",gap_size,"-E",e_value,"_up.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) 
  
  df <- as.data.frame(t(out[,c(7:10)]))
  df_pca <- prcomp(df) 
  df_pcs <-data.frame(df_pca$x,Species=rownames(df)) 
  percentage<-round(df_pca$sdev / sum(df_pca$sdev) * 100,2)
  
  percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))
  ggplot(df_pcs,aes(x=PC1,y=PC2,color=Species))+
    geom_point()+ 
    xlab(percentage[1]) +
    ylab(percentage[2])+    
    geom_text_repel(
      aes(label = rownames(df_pcs)),
      size = 5,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines"))+
    theme(legend.position = "bottom",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size = 20),
          axis.text = element_text(size = 20), 
          axis.title = element_text(size = 20))+
    guides(color = F)+
    theme_bw()
  ggsave(paste0("result/",tissue,"/pca/",antibody,"_merge-W",window_size,"-G",gap_size,"-E",e_value,"_pca_plot.png"),width = 5,height =5)
  write.csv(out,paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_merge-W",window_size,"-G",gap_size,"-E100_diff.csv"),row.names = F)
}
peak_preprocess_macs <- function(tissue,antibody){
  
  tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_macs_narrowpeak.counts"),skip=1)
  # tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_macs2_mergecounts.counts"),skip=1)
  # tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_PePr_old_increase_mergecounts.counts"),skip=1)
  # tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/W1000-G3000-W5000-G10000-intersect.counts"),skip=1)
  counts = tab[,c(7:10)]
  rownames(counts)= tab$Geneid
  colnames(counts) = c("young_1","old_1","young_2","old_2")
  counts <- counts[,c("young_1","young_2","old_1","old_2")]
  group =c("young","young","old","old")
  design = model.matrix(~0+group)
  lvls = levels(factor(group))
  colnames(design) = lvls
  len = length(lvls)
  myContrasts = c("old-young")
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
  out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion))
  # out[,c(7:10)] = sweep(out[,c(7:10)],1,out$Length,'/')*1000
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
  
  ggsave(paste0("result/",tissue,"/diffpeaks/",antibody,"_volcano_plot.png"),width = 10,height = 10)
  # out$Significant_pvalue <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$logFC) >= 0, 
  #                           ifelse(out$logFC > 0, "Up", "Down"), "Stable")
  out$Significant_pvalue <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$`LogFC.old-young`) >= 0, 
                                   ifelse(out$`LogFC.old-young` > 0, "Up", "Down"), "Stable")
  outup <- out[which(out$Significant_pvalue=="Up"),]
  outdown <- out[which(out$Significant_pvalue=="Down"),]
  
  # ggplot(outup, aes(x = Length, fill = Significant)) +  
  #   geom_density(alpha = 0.5) +  
  #   labs(title = "Distribution of Two Groups", x = "Value", y = "Density") +theme_bw()+theme(text = element_text(size = 18))
  # ggplot(outdown, aes(x = Length, fill = Significant)) +  
  #   geom_density(alpha = 0.5) +  
  #   labs(title = "Distribution of Two Groups", x = "Value", y = "Density") +theme_bw()+theme(text = element_text(size = 18))
  
  ggplot() +
    geom_point(data=out, mapping=aes(logCPM, `LogFC.old-young`), size=2, color="grey") +
    geom_point(data=outup, mapping=aes(logCPM, `LogFC.old-young`), size=2, color="red") +
    geom_point(data=outdown, mapping=aes(logCPM, `LogFC.old-young`), size=2, color="blue") +
    labs(x="log2(CPM)",
         y="log2(FC)") +
    # 图例
    theme(legend.position = "bottom",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    theme_bw()
  write.table(outdown[,c("Chr","Start","End","Geneid")], file=paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_macs_narrowpeak_down.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) 
  write.table(outup[,c("Chr","Start","End","Geneid")], file=paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_macs_narrowpeak_up.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) 
  
  df <- as.data.frame(t(out[,c(7:10)]))
  df_pca <- prcomp(df) 
  df_pcs <-data.frame(df_pca$x,Species=rownames(df)) 
  percentage<-round(df_pca$sdev / sum(df_pca$sdev) * 100,2)
  
  percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))
  ggplot(df_pcs,aes(x=PC1,y=PC2,color=Species))+
    geom_point()+ 
    xlab(percentage[1]) +
    ylab(percentage[2])+    
    geom_text_repel(
      aes(label = rownames(df_pcs)),
      size = 5,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines"))+
    theme(legend.position = "bottom",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size = 20),
          axis.text = element_text(size = 20), 
          axis.title = element_text(size = 20))+
    guides(color = F)+
    theme_bw()
  ggsave(paste0("result/",tissue,"/pca/",antibody,"_pca_plot.png"),width = 5,height =5)
  write.csv(out,paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_macs_narrowpeak_diff.csv"),row.names = F)
}

tissue <- "liver"
antibodys <- c("H3K27me3","H3K9me3","H3K36me3")
for (i in c(1:length(antibodys))){
  antibody <- antibodys[i]
  peak_preprocess_sicer(tissue,antibody,window_size,gap_size,e_value)
}

antibodys <- c("H3K27ac","H3K4me3","H3K4me1")
for (i in c(1:length(antibodys))){
  antibody <- antibodys[i]
  peak_preprocess_macs(tissue,antibody)
}
command <- paste("bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/bcv_report.sh /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/", tissue, sep=" ") 
system(command, intern = TRUE)
bcv <- read.table(paste0("data/samples/",tissue,"/bcv_report.txt"), header = FALSE, 
                           col.names = c("bcv", "antibody"),
                           colClasses = c("numeric", "character"))

ggplot(bcv,aes(x=antibody,y=bcv,fill = antibody))+geom_boxplot()+
  stat_summary(fun.y = mean, geom = "point", shape = 23, size=4)+
  # geom_jitter(shape=16, position = position_jitter(0.2))+
  scale_fill_brewer(palette="Set3")+guides(fill = FALSE) +
  ggtitle(paste0(tissue," BCV each antibody"))+theme_bw()+theme(text = element_text(size = 18))+xlab("")+labs(fill = "", color = "") 
ggsave(paste0("result/",tissue,"/diffpeaks/bcv_report.png"),width = 10,height = 8)
