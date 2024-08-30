rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(ChIPseeker)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
library(limma)
tissue <- "colon"
antibody <- "H3K36me3"
bin_size <- "10kb"
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
CUT_search_table <- read.csv("data/samples/all/CUTTag_search_table.csv")
ATAC_search_table <- read.csv("data/samples/all/ATAC_search_table.csv")
search_table <- rbind(CUT_search_table,ATAC_search_table)
tissue_label_change <- function(tissue){
  if(tissue=="brain"){
    tissue_label <- "Cortex"
  }else if(tissue == "Hip"){
    tissue_label <- "Hippocampus"
  }else if(tissue == "CB"){
    tissue_label <- "Cerebellum"
  }else{
    tissue_label <- str_to_title(tissue)
    if(tissue_label == "Bonemarrow"){
      tissue_label <- "Bone Marrow"
    }else if(tissue_label == "Bat"){
      tissue_label <- "BAT"
    }else if(tissue_label=="Mammarygland"){
      tissue_label <- "Mammary Gland"
    }
  }
  return(tissue_label)
}
peak_preprocess_bin_level <- function(tissue,antibody,bin_size){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    bin_size <- "10kb"
  }else{
    bin_size <- "1kb"
  }
  tab = read.delim(paste0("data/samples/intestine/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins.counts"),skip=1)
  counts = tab[,c(7:10)]
  rownames(counts)= tab$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(counts) <- gsub(pattern, "\\1", colnames(counts))
  t_search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]  
  t_search_table$sample_name <- factor(t_search_table$sample_name, levels = colnames(counts))
  t_search_table <- t_search_table[order(t_search_table$sample_name),]
  group <- t_search_table$age
  group[which(group=="3m")] <- "young"
  group[which(group=="24m")] <-"old"
  y= DGEList(counts=counts,group=group)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  if(tissue=="colon"){
    y$samples$batch <- c("batch1","batch1","batch2","batch2")
  }else if(tissue=="cecum"){
    y$samples$batch <- c("batch1","batch2","batch2","batch1")
  }

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
  if(tissue=="colon"){
    colnames(out)[7:10] <- c("old_1","young_1","young_2","old_2")
  }else if(tissue=="cecum"){
    colnames(out)[7:10] <- c("young_1","young_2","old_2","old_1")
  }
  out$Significant <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$`LogFC.old-young`) >= 0, 
                            ifelse(out$`LogFC.old-young` > 0, "Up", "Down"), "Stable")
  out$Significant_bar <- "Stable"
  out$Significant_bar[which(out$`FDR.old-young` < 0.05 & (out$old_1/out$young_1 > 1.2) & (out$old_2/out$young_2 > 1.2))] <- "Up"
  out$Significant_bar[which(out$`FDR.old-young` < 0.05 & (out$old_1/out$young_1 < 0.8) & (out$old_2/out$young_2 < 0.8))] <- "Down"
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  p <- ggplot(
    out, aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`))) +
    geom_point(aes(color = Significant_bar), size=2) +
    scale_color_manual(values = colour) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (p-value)") +
    theme_bw()+
    theme(text = element_text(size = 20))+
    ggtitle(paste0(tissue_label_change(tissue)," ",antibody))+
    annotate("text", x = min(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant_bar=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant_bar=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  # ggsave(paste0("result/",tissue,"/diffpeaks/",antibody,"_merge-W",window_size,"-G",gap_size,"-E",e_value,"_volcano_plot_after_remove_batch_effect.png"),width = 10,height = 10)
  write.csv(out,paste0("data/samples/intestine/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect.csv"),row.names = F)
  outup <- out[which(out$Significant_bar=="Up"),]
  outdown <- out[which(out$Significant_bar=="Down"),]
  write.table(outdown[,c("Chr","Start","End","Geneid")], file=paste0("data/samples/intestine/",tissue,"/",antibody,"/bed/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect_down.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(outup[,c("Chr","Start","End","Geneid")], file=paste0("data/samples/intestine/",tissue,"/",antibody,"/bed/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect_up.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  return(p)
}

p_list <- list()
antibodys <- c("H3K9me3","H3K27me3","H3K36me3","H3K4me3","H3K4me1","H3K27ac")
tissue <-"colon"
for(i in c(1:length(antibodys))){
  antibody <- antibodys[i]
  p_list[[i]] <- peak_preprocess_bin_level(tissue,antibody)
}
combined_plot <- plot_a_list(p_list,2,3)
ggsave(paste0("data/samples/intestine/",tissue,"/all_diff_volcano_plot.png"),combined_plot,width = 20,height = 10,type="cairo")
