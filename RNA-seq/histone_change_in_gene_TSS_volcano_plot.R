rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(edgeR)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
library(limma)
library(data.table)
library(ggrepel)
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
    }else if(tissue_label=="Iwat"){
      tissue_label <- "IWAT"
    }
  }
  return(tissue_label)
}
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}
tissue <- "MEF"
antibody <- "H3K27me3"
volcano_plot <- function(tissue,antibody){
  gene <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))  
  cuttag <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_gene_TSS_10kb_diff_after_remove_batch_effect.csv"))
  peaks <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_young_merge-W5000-G10000-E100.bed"))
  peaks <- as.data.table(peaks)
  setDT(peaks)
  setkey(peaks,V1,V2,V3)
  gene_tss <- cuttag[,c(1:4)]
  gene_tss <- as.data.table(gene_tss)
  gene_tss$Start <- as.numeric(gene_tss$Start)
  setDT(gene_tss)
  setkey(gene_tss,Chr,Start,End)
  overlaps <- foverlaps(gene_tss,peaks, type = "any", nomatch = 0L)
  
  cuttag <- cuttag[which(cuttag$Geneid %in% overlaps$Geneid),]
  colnames(cuttag)[ncol(cuttag)] <- "H3K27me3_Significant"
  gene$condition <- "Significant"
  gene$condition[which(gene$Significant=="Stable")] <-"Stable"
  gene$condition <- factor(gene$condition,levels = c("Stable","Significant"))
  to_plot <- merge(gene,cuttag[,c("Geneid","H3K27me3_Significant")],by.x="X",by.y="Geneid")
  to_plot <- to_plot[which((to_plot$Significant=="Up" & to_plot$H3K27me3_Significant == "Down")|(to_plot$Significant=="Down" & to_plot$H3K27me3_Significant == "Up")),]
  highlight_genes <- subset(to_plot, X %in% c("Cdkn2a", "H2-Q6", "H2-Q7","C4b"))
  
  p1 <- ggplot(
    gene, aes(x = logFC, y = -log10(fdr))) +
    geom_point(aes(color = condition),size=2) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(O/Y)",
         y="-log10 (fdr)") +
    theme_bw()+
    theme(text = element_text(size = 20))+
    ggtitle(paste0(tissue_label_change(tissue)))
  p2 <- p1 +
    geom_point(data = to_plot, aes(x = logFC, y = -log10(fdr), color = H3K27me3_Significant), size=2,alpha=0.5) +
    scale_color_manual(values = c("Significant"="#cecece","Stable"="#515151","Down" = "blue", "Up" = "red"))+
    geom_text_repel(data = highlight_genes, aes(label = X), size = 5, box.padding = 0.5, point.padding = 0.3, segment.color = "black")+xlim(-15,15)+
    geom_point(data = highlight_genes, aes(x = logFC, y = -log10(fdr)), color = "black", fill = "blue", shape = 21, size = 2.5, stroke = 1)

  p2
  # ggsave("result/figures/heart_H3K27me3_RNA_in_young_peaks_volcano_plot.pdf",p2,width = 5.5,height = 6)
  ggsave("result/Sup_figures/MEF_H3K27me3_RNA_in_young_peaks_volcano_plot.pdf",p2,width = 5.5,height = 6)
  return(p2)
}

tissues <- sort(c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
                  "thymus","skin","bladder","bonemarrow","Hip","heart",
                  "muscle","jejunum","uterus","ovary","liver","tongue",
                  "cecum","colon","testis","stomach","pancreas","iWAT","ileum"))
p_list <- list()
for(tissue in sort(tissues)){
  p_list[[tissue]] <- volcano_plot(tissue, antibody)
}
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}
volcano_plot_combined <- plot_a_list(p_list,no_of_rows = 4,no_of_cols = 7)
ggsave(paste0("tmp.png"),volcano_plot_combined,width = 42,height = 24,type="cairo")
