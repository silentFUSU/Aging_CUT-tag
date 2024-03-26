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
library(tidyr)
library(MASS) 
antibody = "H3K9me3"
tissue = c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas","Hip","cecum","bonemarrow")
# brain liver testis colon kidney lung spleen muscle pancreas
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
  batch <- c()
  for(i in c(1:length(tissue))){
    for(j in c(1:4)){
      colnames(counts)[(i-1)*4+j] <- paste0(tissue[i],"_",group_label[j])
      group <- c(group,paste0(tissue[i],"_",group_label2[j]))
    }
    batch <- c(batch,rep(c(rep("batch1", 2), rep("batch2", 2)), 1))
  }
  y= DGEList(counts=counts,group=group)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  tab<-tab[keep,]
  y$samples$batch <- batch
  y$samples$year <- group
  y <- calcNormFactors(y)
  logCPMs <- cpm(y, log = TRUE)
  # rv <- apply(logCPMs, 1, var)
  # o <- order(rv, decreasing=TRUE)
  # top1000 <- head(o, 1000)
  # logCPM_top1000 <- logCPMs[top1000,]
  pca <- prcomp(t(logCPMs))
  to_plot <- data.frame(pca$x, y$samples)
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
  use.pcs <- c(1,2)
  labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
  ggplot(to_plot, aes(x=PC1, y=PC2, color=year, shape=batch)) + 
    geom_point(size=5) +theme_bw()+
    xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+    
    geom_text_repel(
      aes(label = rownames(to_plot)),
      size = 5,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines"))
  # ggsave(paste0("result/all/pca/",antibody,"_merge-W",window_size,"-G",gap_size,"-E",e_value,"_pca_plot_before_remove_batch_effect.png"),width = 10,height =8)
  
  for(i in c(1:length(tissue))){
    if(i == 1){
      t_logCPMs <- logCPMs[,c(((i-1)*4+1):(i*4))]
      t_batch <- batch[c(((i-1)*4+1):(i*4))]
      logCPMs_corrected <- limma::removeBatchEffect(t_logCPMs, batch = t_batch)
    }else{
      t_logCPMs <- logCPMs[,c(((i-1)*4+1):(i*4))]
      t_batch <- batch[c(((i-1)*4+1):(i*4))]
      t_logCPMs_corrected <- limma::removeBatchEffect(t_logCPMs, batch = t_batch)
      logCPMs_corrected <- cbind(logCPMs_corrected,t_logCPMs_corrected)
    }
  }

  # logCPM_corrected_top1000 <- logCPMs_corrected[top1000,]
  pca <- prcomp(t(logCPMs_corrected))
  to_plot <- data.frame(pca$x, y$samples)
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
  use.pcs <- c(1,2)
  labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
  
  ggplot(to_plot, aes(x=PC1, y=PC2, color=year, shape=batch)) + 
    geom_point(size=5) +theme_bw()+
    xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+    
    geom_text_repel(
      aes(label = rownames(to_plot)),
      size = 5,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines"))
  # ggsave(paste0("result/all/pca/",antibody,"_merge-W",window_size,"-G",gap_size,"-E",e_value,"_pca_plot_after_remove_batch_effect.png"),width = 10,height =8)
  
  to_plot <- separate(to_plot, year, into = c("tissue", "age"), sep = "_")  
  ggplot(to_plot, aes(x=PC1, y=PC2, color=tissue, shape=age)) + 
    geom_point(size=5) +theme_bw()+
    xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+    
    geom_text_repel(
      aes(label = group),
      size = 5,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines"),
      max.overlaps=Inf)
  ggsave(paste0("result/all/pca/",antibody,"_merge-W",window_size,"-G",gap_size,"-E",e_value,"_pca_plot_after_remove_batch_effect_tissue_age.png"),width = 10,height =8)
  # for(i in c(1:length(tissue))){
  #   t_counts <- counts[,c(((i-1)*4+1):(i*4))]
  #   t_batch <- batch[c(((i-1)*4+1):(i*4))]
  #   t_group <- group[c(((i-1)*4+1):(i*4))]
  #   t_y= DGEList(counts=t_counts,group=t_group)
  #   t_y=t_y[keep,]
  #   t_y$samples$batch <- t_batch
  #   t_y$samples$year <- t_group
  #   t_y <- calcNormFactors(t_y)
  #   t_y$samples$batch <- factor(t_y$samples$batch)
  #   t_y$samples$year <- factor(t_y$samples$year,levels = c(paste0(tissue[i],"_young"),paste0(tissue[i],"_old")))
  #   design <- model.matrix(~batch+year, t_y$samples)
  #   # y <- estimateDisp(y, design)
  #   y<-estimateCommonDisp(y)
  #   y<-estimateGLMTagwiseDisp(y,design)
  #   fit_tag = glmFit(y,design)
  #   lrt = glmLRT(fit_tag, coef = 3)
  #   
  #   if(i == 1){
  #     out = cbind(tab[,1:6],cpm(y),paste0(tissue[i],".logCPM")=lrt$table$logCPM, 
  #                 paste0(tissue[i],".bcv")=sqrt(fit_tag$dispersion),
  #                 paste0(tissue[i],".PValue.old-young")=lrt$table$PValue,
  #                 paste0(tissue[i],".FDR.old-young")= p.adjust(lrt$table$PValue,method="BH"),
  #                 paste0(tissue[i],".LogFC.old-young")=lrt$table$logFC)
  #   }else{
  #     t_logCPMs <- logCPMs[,c(((i-1)*4+1):(i*4))]
  #     t_batch <- batch[c(((i-1)*4+1):(i*4))]
  #     t_logCPMs_corrected <- limma::removeBatchEffect(t_logCPMs, batch = t_batch)
  #     logCPMs_corrected <- cbind(logCPMs_corrected,t_logCPMs_corrected)
  #   }
  # }
  # y$samples$batch <- factor(y$samples$batch)
  # y$samples$year <- factor(y$samples$year)
  # design <- model.matrix(~batch+year, y$samples, 
  #                        contrasts.arg = list(year=contrasts(y$samples$year, contrasts = FALSE)))
}

peak_preprocess_macs2 <- function(tissue,antibody){
  tab = read.delim(paste0("data/samples/all/",antibody,"/merge_macs_narrowpeak.counts"),skip=1)
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  group <- c()
  batch <- c()
  for(i in c(1:length(tissue))){
    for(j in c(1:4)){
      colnames(counts)[(i-1)*4+j] <- paste0(tissue[i],"_",group_label[j])
      group <- c(group,paste0(tissue[i],"_",group_label2[j]))
    }
    batch <- c(batch,rep(c(rep("batch1", 2), rep("batch2", 2)), 1))
  }
  y= DGEList(counts=counts,group=group)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  y$samples$batch <- batch
  y$samples$year <- group
  y <- calcNormFactors(y)
  
  logCPMs <- cpm(y, log = T)
  # rv <- apply(logCPMs, 1, var)
  # o <- order(rv, decreasing=TRUE)
  # top1000 <- head(o, 1000)
  # logCPM_top1000 <- logCPMs[top1000,]
  pca <- prcomp(t(logCPMs))
  to_plot <- data.frame(pca$x, y$samples)
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
  use.pcs <- c(1,2)
  labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
  ggplot(to_plot, aes(x=PC1, y=PC2, color=year, shape=batch)) + 
    geom_point(size=5) +theme_bw()+
    xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+    
    geom_text_repel(
      aes(label = rownames(to_plot)),
      size = 5,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines"))
  # ggsave(paste0("result/all/pca/",antibody,"_pca_plot_before_remove_batch_effect.png"),width = 10,height =8)
  
  batch <- factor(y$samples$batch)
  for(i in c(1:length(tissue))){
    if(i == 1){
      t_logCPMs <- logCPMs[,c(((i-1)*4+1):(i*4))]
      t_batch <- batch[c(((i-1)*4+1):(i*4))]
      logCPMs_corrected <- limma::removeBatchEffect(t_logCPMs, batch = t_batch)
    }else{
      t_logCPMs <- logCPMs[,c(((i-1)*4+1):(i*4))]
      t_batch <- batch[c(((i-1)*4+1):(i*4))]
      t_logCPMs_corrected <- limma::removeBatchEffect(t_logCPMs, batch = t_batch)
      logCPMs_corrected <- cbind(logCPMs_corrected,t_logCPMs_corrected)
    }
  }
  
  # logCPM_corrected_top1000 <- logCPMs_corrected[top1000,]
  pca <- prcomp(t(logCPMs_corrected))
  to_plot <- data.frame(pca$x, y$samples)
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
  use.pcs <- c(1,2)
  labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
  
  ggplot(to_plot, aes(x=PC1, y=PC2, color=year, shape=batch)) + 
    geom_point(size=5) +theme_bw()+
    xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+    
    geom_text_repel(
      aes(label = rownames(to_plot)),
      size = 5,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines"))
  # ggsave(paste0("result/all/pca/",antibody,"_pca_plot_after_remove_batch_effect.png"),width = 10,height =8)
  to_plot <- separate(to_plot, year, into = c("tissue", "age"), sep = "_")  
  ggplot(to_plot, aes(x=PC1, y=PC2, color=tissue, shape=age)) + 
    geom_point(size=5) +theme_bw()+
    xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))  
    geom_text_repel(
      aes(label = group),
      size = 5,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines"),
      max.overlaps=Inf)
  ggsave(paste0("result/all/pca/",antibody,"_pca_plot_after_remove_batch_effect_tissue_age.png"),width = 10,height =8)
}

bin_size <- "10kb"
peak_preprocess_bins <- function(tissue,antibody,bin_size){
  tab = read.delim(paste0("data/samples/all/",antibody,"/merge-",bin_size,"_bins.counts"),skip=1)
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  group <- c()
  batch <- c()
  tissue_group <- c()
  year <-c()
  for(i in c(1:length(tissue))){
    for(j in c(1:4)){
      colnames(counts)[(i-1)*4+j] <- paste0(tissue[i],"_",group_label[j])
      group <- c(group,paste0(tissue[i],"_",group_label2[j]))
      year <- c(year,paste0(group_label2[j]))
      tissue_group <-c(tissue_group,tissue[i])
    }
    batch <- c(batch,rep(c(rep("batch1", 2), rep("batch2", 2)), 1))
  }
  y= DGEList(counts=counts,group=group)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  y$samples$batch <- batch
  y$samples$year <- year
  y$samples$tissue <- tissue_group
  y <- calcNormFactors(y)
  
  logCPMs <- cpm(y, log = T)
  # rv <- apply(logCPMs, 1, var)
  # o <- order(rv, decreasing=TRUE)
  # top1000 <- head(o, 1000)
  # logCPM_top1000 <- logCPMs[top1000,]
  pca <- prcomp(t(logCPMs))
  to_plot <- data.frame(pca$x, y$samples)
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
  use.pcs <- c(1,2)
  labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
  ggplot(to_plot, aes(x=PC1, y=PC2, color=group, shape=batch)) + 
    geom_point(size=5) +theme_bw()+
    xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+    
    geom_text_repel(
      aes(label = rownames(to_plot)),
      size = 5,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines"))
  # ggsave(paste0("result/all/pca/",antibody,"_pca_plot_before_remove_batch_effect.png"),width = 10,height =8)
  
  batch <- factor(y$samples$batch)
  for(i in c(1:length(tissue))){
    if(i == 1){
      t_logCPMs <- logCPMs[,c(((i-1)*4+1):(i*4))]
      t_batch <- batch[c(((i-1)*4+1):(i*4))]
      logCPMs_corrected <- limma::removeBatchEffect(t_logCPMs, batch = t_batch)
    }else{
      t_logCPMs <- logCPMs[,c(((i-1)*4+1):(i*4))]
      t_batch <- batch[c(((i-1)*4+1):(i*4))]
      t_logCPMs_corrected <- limma::removeBatchEffect(t_logCPMs, batch = t_batch)
      logCPMs_corrected <- cbind(logCPMs_corrected,t_logCPMs_corrected)
    }
  }
  
  # logCPM_corrected_top1000 <- logCPMs_corrected[top1000,]
  pca <- prcomp(t(logCPMs_corrected))
  to_plot <- data.frame(pca$x, y$samples)
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
  use.pcs <- c(1,2)
  labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
  
  ggplot(to_plot, aes(x=PC1, y=PC2, color=group, shape=batch)) + 
    geom_point(size=5) +theme_bw()+
    xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+    
    geom_text_repel(
      aes(label = rownames(to_plot)),
      size = 5,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines"))
  # ggsave(paste0("result/all/pca/",antibody,"_pca_plot_after_remove_batch_effect.png"),width = 10,height =8)
  to_plot <- separate(to_plot, group, into = c("tissue", "age"), sep = "_")  
  to_plot$tissue[which(to_plot$tissue=="brain")] <- "brain_FC"
  to_plot$tissue[which(to_plot$tissue=="Hip")] <- "brain_Hip"
  to_plot$group <- paste0(to_plot$tissue,"_",to_plot$age)
  ggplot(to_plot, aes(x=PC1, y=PC2, color=tissue, shape=age)) + 
    geom_point(size=5) +theme_bw()+
    xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
  geom_text_repel(
    aes(label = group),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    max.overlaps=Inf)
  ggsave(paste0("result/all/pca/",antibody,"_bin_",bin_size,"_pca_plot_after_remove_batch_effect_tissue_age.png"),width = 10,height =8,type="cairo")
 
  y= DGEList(counts=counts,group=year)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  y$samples$batch <- batch
  y$samples$year <- year
  y$samples$tissue <- tissue_group
  y$samples$group <- factor(y$samples$group,c("young","old"))
  y <- calcNormFactors(y)
  design <- model.matrix(~tissue+group, y$samples)
  y<-estimateCommonDisp(y)
  y<-estimateGLMTagwiseDisp(y,design)
  fit_tag = glmFit(y,design)
  lrt = glmLRT(fit_tag, coef = ncol(design))
  tab<-tab[keep,]
  out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
              "PValue.old-young"=lrt$table$PValue,"FDR.old-young"= p.adjust(lrt$table$PValue,method="BH"),
              "LogFC.old-young"=lrt$table$logFC)
  out$Significant <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$`LogFC.old-young`) >= 0, 
                            ifelse(out$`LogFC.old-young` > 0, "Up", "Down"), "Stable")
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  peak <- GRanges(seqnames = out$Chr,   
                  ranges = IRanges(start = out$Start, end = out$End))
  peak_anno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                            TxDb=txdb, annoDb="org.Mm.eg.db")
  peak_anno<-unique(as.data.frame(peak_anno))
  peak_anno$V4 <- paste0(peak_anno$seqnames,"-",peak_anno$start,"-",peak_anno$end)
  out$V4 <- paste0(out$Chr,"-",out$Start,"-",out$End)
  out <- merge(out,peak_anno[,c(6:18)],by="V4")
  write.csv(out,paste0("data/samples/all/",antibody,"/",antibody,"_",bin_size,"_bins_diff_use_tissue_as_batch_effect.csv"),row.names = F)
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
  ggsave(paste0("result/all/diff/",antibody,"_all_tissue_as_batch_effect.png"),width = 10,height = 10,type="cairo")
  }

antibodys <- c("H3K27me3","H3K9me3","H3K36me3")
for (i in c(1:length(antibodys))){
  antibody <- antibodys[i]
  peak_preprocess_sicer(tissue,antibody,window_size,gap_size,e_value)
}
antibodys <- c("H3K27ac","H3K4me3","H3K4me1")
for (i in c(1:length(antibodys))){
  antibody <- antibodys[i]
  peak_preprocess_macs2(tissue,antibody)
}

dir.create("result/all/diff/")
antibodys <- c("H3K27me3","H3K9me3","H3K36me3")
bin_size <-"10kb"
for (i in c(1:length(antibodys))){
  antibody <- antibodys[i]
  peak_preprocess_bins(tissue,antibody,bin_size)
}
bin_size <-"1kb"
antibodys <- c("H3K27ac","H3K4me3","H3K4me1")
for (i in c(1:length(antibodys))){
  antibody <- antibodys[i]
  peak_preprocess_bins(tissue,antibody,bin_size)
}

chr_percent <- as.data.frame(table(out$Chr[which(out$Significant=="Up")]))
chr_percent$antibody <- "H3K9me3"
chr_percent$Var1 <- factor(chr_percent$Var1,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7",
                                                     "chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                                                     "chr17","chr18","chr19","chrX","chrY"))
colors <- c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#800000", "#008000", "#000080", "#808000", "#800080", "#008080", "#C0C0C0", "#808080", "#FFA500", "#800000", "#FFD700", "#FF69B4", "#00FA9A", "#9400D3", "#FF4500") 
ggplot(chr_percent, aes(fill=Var1, y=Freq, x=antibody)) + 
  geom_bar(position="stack", stat="identity")+scale_fill_manual(values =colors)+theme_bw()

chr_percent <- as.data.frame(table(out$Chr[which(out$Significant=="Down")]))
chr_percent$antibody <- "H3K27me3"
chr_percent$Var1 <- factor(chr_percent$Var1,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7",
                                                     "chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                                                     "chr17","chr18","chr19","chrX","chrY"))
colors <- c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#800000", "#008000", "#000080", "#808000", "#800080", "#008080", "#C0C0C0", "#808080", "#FFA500", "#800000", "#FFD700", "#FF69B4", "#00FA9A", "#9400D3", "#FF4500") 
ggplot(chr_percent, aes(fill=Var1, y=Freq, x=antibody)) + 
  geom_bar(position="stack", stat="identity")+scale_fill_manual(values =colors)+theme_bw()
