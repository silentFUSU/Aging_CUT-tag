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

peak_preprocess_bin_level <- function(tissue,antibody,bin_size){
  tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins.counts"),skip=1)
  counts = tab[,c(7:10)]
  rownames(counts)= tab$Geneid
  if (tissue=="ovary"){
    colnames(counts) = c("old_1","young_1","old_2","young_2")
    group =c("old","young","old","young")
  }else{
    colnames(counts) = c("young_1","old_1","young_2","old_2")
    group =c("young","old","young","old")
  }

  y= DGEList(counts=counts,group=group)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  y$samples$batch <- rep(c(rep("batch1", 2), rep("batch2", 2)), 1)
  y$samples$year <- group
  y$samples$year <- factor(y$samples$year,c("young","old"))
  y <- calcNormFactors(y)
  batch <- factor(y$samples$batch)
  design <- model.matrix(~batch+year, y$samples)
  y<-estimateCommonDisp(y)
  y<-estimateGLMTagwiseDisp(y,design)
  fit_tag = glmFit(y,design)
  lrt = glmLRT(fit_tag, coef = 3)
  tab<-tab[keep,]

  out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
              "PValue.old-young"=lrt$table$PValue,"FDR.old-young"= p.adjust(lrt$table$PValue,method="BH"),
              "LogFC.old-young"=lrt$table$logFC)
  
  out$Significant <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$`LogFC.old-young`) >= 0, 
                            ifelse(out$`LogFC.old-young` > 0, "Up", "Down"), "Stable")
  out$Significant_bar <- "Stable"
  out$Significant_bar[which(out$`FDR.old-young` < 0.05 & (out$old_1/out$young_1 > 1.2) & (out$old_2/out$young_2 > 1.2))] <- "Up"
  out$Significant_bar[which(out$`FDR.old-young` < 0.05 & (out$old_1/out$young_1 < 0.8) & (out$old_2/out$young_2 < 0.8))] <- "Down"
  colour<- list(c("grey"),c("grey","red"),c("blue","grey","red"))
  ggplot(
    out, aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`))) +
    geom_point(aes(color = Significant), size=2) +
    scale_color_manual(values = colour[[nrow(as.data.frame(table(out$Significant)))]]) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (p-value)") +
    theme_bw()+
    theme(text = element_text(size = 20))+
    annotate("text", x = min(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant_bar=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant_bar=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  # ggsave(paste0("result/",tissue,"/diffpeaks/",antibody,"_merge-W",window_size,"-G",gap_size,"-E",e_value,"_volcano_plot_after_remove_batch_effect.png"),width = 10,height = 10)
  write.csv(out,paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect.csv"),row.names = F)
  outup <- out[which(out$Significant_bar=="Up"),]
  outdown <- out[which(out$Significant_bar=="Down"),]
  write.table(outdown[,c("Chr","Start","End","Geneid")], file=paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect_down.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(outup[,c("Chr","Start","End","Geneid")], file=paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect_up.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

tissues <- c("muscle","brain","liver","testis","colon","kidney","lung","spleen","pancreas")
tissues <- c("Hip") 
tissues <- c("cecum") 
tissues <- c("bonemarrow")
tissues <- c("ileum")
tissues <- c("heart")
tissues <- c("thymus")
tissues <- c("stomach")
tissues <- c("skin")
tissues <- c("skin_H4K16ac")
tissues <- c("aorta","tongue")
tissues <- c("bladder")
tissues <- c("BAT")
tissues <- c("mammarygland")
antibodys <- c("H3K36me3","H3K27me3","H3K9me3","H3K27ac","H3K4me3","H3K4me1")
# antibodys <- c("ATAC")
# tissues<-c("Hip","testis", "colon", "kidney", "lung", "spleen", "muscle", "pancreas","cecum","bonemarrow","ileum","heart","thymus")
bin_size <-"10kb"
for (i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:lengåth(antibodys))){
    antibody <- antibodys[j]
    peak_preprocess_bin_level(tissue,antibody,bin_size)
  }
}
bin_size <-"1kb"
antibodys <- c("H3K27ac","H3K4me3","H3K4me1")
# antibodys <- c("H4K16ac")
for (i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(antibodys))){
    antibody <- antibodys[j]
    peak_preprocess_bin_level(tissue,antibody,bin_size)
  }
}
