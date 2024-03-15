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
antibody = "H3K27ac"
tissue = c("brain","liver","testis")
window_size = "1000"
gap_size = "3000"
e_value = "100"
group_label <- c("young_1","old_1","young_2","old_2")
group_label2 <- c("young","old","young","old")

peak_preprocess_sicer <- function(tissue,antibody,window_size,gap_size,e_value){
  tab = read.delim(paste0("data/samples/all/",antibody,"/merge-W",window_size,"-G",gap_size,"-E",e_value,".counts"),skip=1)
  # tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_macs2_mergecounts.counts"),skip=1)
  # tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_PePr_old_increase_mergecounts.counts"),skip=1)
  # tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/W1000-G3000-W5000-G10000-intersect.counts"),skip=1)
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  group <- c()
  for(i in c(1:length(tissue))){
    for(j in c(1:4)){
      colnames(counts)[(i-1)*4+j] <- paste0(tissue[i],"_",group_label[j])
      group <- c(group,paste0(tissue[i],"_",group_label2[j]))
    }
  }
  design = model.matrix(~0+group)
  lvls = levels(factor(group))
  colnames(design) = lvls
  len = length(lvls)
  
  myContrasts = c()
  for(i in c(1:length(tissue))){
    myContrasts <- c(myContrasts,paste0(tissue[i],"_old-",tissue[i],"_young"))
  }
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
  out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM)
  # out[,c(7:10)] = sweep(out[,c(7:10)],1,out$Length,'/')*1000
  for (i in 1:length(myContrasts)){
    lrt = glmLRT(fit_tag, contrast = contrast.matrix[,i])
    out[,paste0("PValue.",myContrasts[i])] = lrt$table$PValue
    out[,paste0("FDR.",myContrasts[i])] = p.adjust(lrt$table$PValue,method="BH")
    out[,paste0("LogFC.",myContrasts[i])] = lrt$table$logFC
  }
  
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  # peak <- readPeakFile(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_macs2_mergecounts.bed"))
  peak <- readPeakFile(paste0("data/samples/all/",antibody,"/bed/merge-W1000-G3000-E100.bed"))
  peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                           TxDb=txdb, annoDb="org.Mm.eg.db")
  peakAnno<-unique(as.data.frame(peakAnno))
  
  out$annotation <- peakAnno$annotation[which(peakAnno$V4 %in% out$Geneid)]
  out$ens_id <- peakAnno$geneId[which(peakAnno$V4 %in% out$Geneid)]
  out$symbol <- peakAnno$SYMBOL[which(peakAnno$V4 %in% out$Geneid)]
  
  df <- as.data.frame(t(out[,c(7:ncol(tab))]))
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
  ggsave(paste0("result/all/pca/",antibody,"_pca_plot.png"),width = 5,height =5)
  write.csv(out,paste0("data/samples/all/",antibody,"/merge-W",window_size,"-G",gap_size,"-E100_diff.csv"),row.names = F)
  write.csv(out[which(out$`FDR.brain_old-brain_young`<0.05 & out$`FDR.liver_old-liver_young`<0.05),],paste0("data/samples/all/",antibody,"/merge-W",window_size,"-G",gap_size,"-E100_diff_most_significant.csv"),row.names = F)
}

peak_preprocess_macs2 <- function(tissue,antibody){
  tab = read.delim(paste0("data/samples/all/",antibody,"/merge_macs_narrowpeak.counts"),skip=1)
  # tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_macs2_mergecounts.counts"),skip=1)
  # tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_PePr_old_increase_mergecounts.counts"),skip=1)
  # tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/W1000-G3000-W5000-G10000-intersect.counts"),skip=1)
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  group <- c()
  for(i in c(1:length(tissue))){
    for(j in c(1:4)){
      colnames(counts)[(i-1)*4+j] <- paste0(tissue[i],"_",group_label[j])
      group <- c(group,paste0(tissue[i],"_",group_label2[j]))
    }
  }
  design = model.matrix(~0+group)
  lvls = levels(factor(group))
  colnames(design) = lvls
  len = length(lvls)
  
  myContrasts = c()
  for(i in c(1:length(tissue))){
    myContrasts <- c(myContrasts,paste0(tissue[i],"_old-",tissue[i],"_young"))
  }
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
  out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM)
  # out[,c(7:10)] = sweep(out[,c(7:10)],1,out$Length,'/')*1000
  for (i in 1:length(myContrasts)){
    lrt = glmLRT(fit_tag, contrast = contrast.matrix[,i])
    out[,paste0("PValue.",myContrasts[i])] = lrt$table$PValue
    out[,paste0("FDR.",myContrasts[i])] = p.adjust(lrt$table$PValue,method="BH")
    out[,paste0("LogFC.",myContrasts[i])] = lrt$table$logFC
  }
  
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  # peak <- readPeakFile(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_macs2_mergecounts.bed"))
  peak <- readPeakFile(paste0("data/samples/all/",antibody,"/bed/merge_macs_narrowpeak.bed"))
  peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                           TxDb=txdb, annoDb="org.Mm.eg.db")
  peakAnno<-unique(as.data.frame(peakAnno))
  
  out$annotation <- peakAnno$annotation[which(peakAnno$V4 %in% out$Geneid)]
  out$ens_id <- peakAnno$geneId[which(peakAnno$V4 %in% out$Geneid)]
  out$symbol <- peakAnno$SYMBOL[which(peakAnno$V4 %in% out$Geneid)]
  
  df <- as.data.frame(t(out[,c(7:ncol(tab))]))
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
  ggsave(paste0("result/all/pca/",antibody,"_pca_plot.png"),width = 5,height =5)
  write.csv(out,paste0("data/samples/all/",antibody,"/merge_macs_narrowpeak_diff.csv"),row.names = F)
  write.csv(out[which(out$`FDR.brain_old-brain_young`<0.05 & out$`FDR.liver_old-liver_young`<0.05),],paste0("data/samples/all/",antibody,"/merge_macs_narrowpeak_diff_most_significant.csv"),row.names = F)
}

peak_preprocess_sicer_intersect <- function(tissue,antibody,window_size,gap_size,e_value){
  tab = read.delim(paste0("data/samples/all/",antibody,"/intersect-W",window_size,"-G",gap_size,"-E",e_value,".counts"),skip=1)
  # tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_macs2_mergecounts.counts"),skip=1)
  # tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_PePr_old_increase_mergecounts.counts"),skip=1)
  # tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/W1000-G3000-W5000-G10000-intersect.counts"),skip=1)
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  group <- c()
  for(i in c(1:length(tissue))){
    for(j in c(1:4)){
      colnames(counts)[(i-1)*4+j] <- paste0(tissue[i],"_",group_label[j])
      group <- c(group,paste0(tissue[i],"_",group_label2[j]))
    }
  }
  design = model.matrix(~0+group)
  lvls = levels(factor(group))
  colnames(design) = lvls
  len = length(lvls)
  
  myContrasts = c()
  for(i in c(1:length(tissue))){
    myContrasts <- c(myContrasts,paste0(tissue[i],"_old-",tissue[i],"_young"))
  }
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
  out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM)
  # out[,c(7:10)] = sweep(out[,c(7:10)],1,out$Length,'/')*1000
  for (i in 1:length(myContrasts)){
    lrt = glmLRT(fit_tag, contrast = contrast.matrix[,i])
    out[,paste0("PValue.",myContrasts[i])] = lrt$table$PValue
    out[,paste0("FDR.",myContrasts[i])] = p.adjust(lrt$table$PValue,method="BH")
    out[,paste0("LogFC.",myContrasts[i])] = lrt$table$logFC
  }
  
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  # peak <- readPeakFile(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_macs2_mergecounts.bed"))
  peak <- readPeakFile(paste0("data/samples/all/",antibody,"/bed/intersect-W1000-G3000-E100.bed"))
  peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                           TxDb=txdb, annoDb="org.Mm.eg.db")
  peakAnno<-unique(as.data.frame(peakAnno))
  
  out$annotation <- peakAnno$annotation[which(peakAnno$V4 %in% out$Geneid)]
  out$ens_id <- peakAnno$geneId[which(peakAnno$V4 %in% out$Geneid)]
  out$symbol <- peakAnno$SYMBOL[which(peakAnno$V4 %in% out$Geneid)]
  
  df <- as.data.frame(t(out[,c(7:ncol(tab))]))
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
  ggsave(paste0("result/all/pca/",antibody,"_intersect_pca_plot.png"),width = 5,height =5)
  write.csv(out,paste0("data/samples/all/",antibody,"/intersect-W",window_size,"-G",gap_size,"-E100_diff.csv"),row.names = F)
  write.csv(out[which(out$`FDR.brain_old-brain_young`<0.05 & out$`FDR.liver_old-liver_young`<0.05),],paste0("data/samples/all/",antibody,"/intersect-W",window_size,"-G",gap_size,"-E100_diff_most_significant.csv"),row.names = F)
}

peak_preprocess_macs2_intersect <- function(tissue,antibody){
  tab = read.delim(paste0("data/samples/all/",antibody,"/intersect_macs_narrowpeak.counts"),skip=1)
  # tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_macs2_mergecounts.counts"),skip=1)
  # tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_PePr_old_increase_mergecounts.counts"),skip=1)
  # tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/W1000-G3000-W5000-G10000-intersect.counts"),skip=1)
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  group <- c()
  for(i in c(1:length(tissue))){
    for(j in c(1:4)){
      colnames(counts)[(i-1)*4+j] <- paste0(tissue[i],"_",group_label[j])
      group <- c(group,paste0(tissue[i],"_",group_label2[j]))
    }
  }
  design = model.matrix(~0+group)
  lvls = levels(factor(group))
  colnames(design) = lvls
  len = length(lvls)
  
  myContrasts = c()
  for(i in c(1:length(tissue))){
    myContrasts <- c(myContrasts,paste0(tissue[i],"_old-",tissue[i],"_young"))
  }
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
  out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM)
  # out[,c(7:10)] = sweep(out[,c(7:10)],1,out$Length,'/')*1000
  for (i in 1:length(myContrasts)){
    lrt = glmLRT(fit_tag, contrast = contrast.matrix[,i])
    out[,paste0("PValue.",myContrasts[i])] = lrt$table$PValue
    out[,paste0("FDR.",myContrasts[i])] = p.adjust(lrt$table$PValue,method="BH")
    out[,paste0("LogFC.",myContrasts[i])] = lrt$table$logFC
  }
  
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  # peak <- readPeakFile(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_macs2_mergecounts.bed"))
  peak <- readPeakFile(paste0("data/samples/all/",antibody,"/bed/intersect_macs_narrowpeak.bed"))
  peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                           TxDb=txdb, annoDb="org.Mm.eg.db")
  peakAnno<-unique(as.data.frame(peakAnno))
  
  out$annotation <- peakAnno$annotation[which(peakAnno$V4 %in% out$Geneid)]
  out$ens_id <- peakAnno$geneId[which(peakAnno$V4 %in% out$Geneid)]
  out$symbol <- peakAnno$SYMBOL[which(peakAnno$V4 %in% out$Geneid)]
  
  df <- as.data.frame(t(out[,c(7:ncol(tab))]))
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
  ggsave(paste0("result/all/pca/",antibody,"_intersect_pca_plot.png"),width = 5,height =5)
  write.csv(out,paste0("data/samples/all/",antibody,"/intersect_macs_narrowpeak_diff.csv"),row.names = F)
  write.csv(out[which(out$`FDR.brain_old-brain_young`<0.05 & out$`FDR.liver_old-liver_young`<0.05),],paste0("data/samples/all/",antibody,"/intersect_macs_narrowpeak_diff_most_significant.csv"),row.names = F)
}

antibodys <- c("H3K27me3","H3K9me3","H3K36me3")
for(i in c(1:length(antibodys))){
  peak_preprocess_sicer(tissue,antibodys[i],window_size,gap_size,e_value)
}

antibodys <- c("H3K27ac","H3K4me3","H3K4me1")
for(i in c(1:length(antibodys))){
  peak_preprocess_macs2(tissue,antibodys[i])
}

antibodys <- c("H3K36me3")
for(i in c(1:length(antibodys))){
  peak_preprocess_sicer(tissue,antibodys[i],window_size,gap_size,e_value)
}

antibodys <- c("H3K27ac","H3K4me3","H3K4me1")
for(i in c(1:length(antibodys))){
  peak_preprocess_macs2_intersect(tissue,antibodys[i])
}
