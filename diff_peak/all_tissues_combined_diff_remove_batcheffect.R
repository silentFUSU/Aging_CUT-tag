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
antibody = "H3K4me1"
tissue = c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas","Hip","cecum")
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
  ggsave(paste0("result/all/pca/",antibody,"_bin_",bin_size,"_pca_plot_after_remove_batch_effect_tissue_age.png"),width = 10,height =8)
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
