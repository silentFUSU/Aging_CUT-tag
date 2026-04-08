rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(ggrepel)
library(ChIPseeker)
senescence <- read.csv("data/samples/MEF/H3K27me3/H3K27me3_10kb_bins_diff_after_remove_batch_effect.csv")
senescence <- senescence[,c("Geneid","LogFC.old.young","Significant")]

antibody <- "H3K27me3"
conditions <- c("MEF_Bmi1", "MEF_Cbx2", "MEF_Cbx7")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
p_list <- list()
pvalue_summary <- data.frame()
for(condition in conditions){
  df <- read.csv(paste0("data/samples/MEF_OE/H3K27me3/",condition,"/H3K27me3_",condition,"_10kb_bins_diff_after_remove_batch_effect.csv"))
  # df <- df[,c("Geneid","LogFC.oe.vec","Significant")]
  to_plot <- merge(df,senescence,by="Geneid")
  to_plot$Significant.x[is.na(to_plot$Significant.x)] <- "Stable"
  to_plot$Significant.y[is.na(to_plot$Significant.y)] <- "Stable"
  to_plot$condition <- "Stable"
  to_plot$condition[which(to_plot$Significant.x=="Up" & to_plot$Significant.y=="Up")] <- "Up"
  to_plot$condition[which(to_plot$Significant.x=="Down" & to_plot$Significant.y=="Down")] <- "Down"
  to_plot$condition[which(to_plot$Significant.x=="Stable" & to_plot$Significant.y=="Stable")] <- "Stable"
  to_plot$condition[which(to_plot$Significant.x=="Up" & to_plot$Significant.y=="Down")] <- "Inconsistent"
  to_plot$condition[which(to_plot$Significant.x=="Down" & to_plot$Significant.y=="Up")] <- "Inconsistent"
  
  to_plot$quadrant <- "other"
  to_plot$quadrant[which(to_plot$condition=="Up")] <- "first"
  to_plot$quadrant[which(to_plot$condition=="Inconsistent" & to_plot$LogFC.oe.vec < 0)] <- "second"
  to_plot$quadrant[which(to_plot$condition=="Down")] <- "third"
  to_plot$quadrant[which(to_plot$condition=="Inconsistent" & to_plot$LogFC.oe.vec > 0)] <- "fourth"

  OE_RNA <- read.csv(paste0("data/samples/RNA/MEF_OE_RNA/MEF_Vector_",condition,"_diff_expression_gene_strict_filter_bar.csv"))
  SE_RNA <- read.csv(paste0("data/samples/RNA/MEF/diff_expression_gene_strict_filter_bar.csv"))
  colnames(SE_RNA)[1] <- "Geneid"
  logFC_summary <-data.frame()

  for(quadrant in c("first","second","third","fourth")){
    bins <- to_plot[which(to_plot$quadrant==quadrant),c(2:4)]
    gr <- GRanges(
      seqnames = bins$Chr,
      ranges   = IRanges(start = bins$Start + 1, end = bins$End)  # BED start是0-based，GRanges是1-based
    )
    peakAnno <- annotatePeak(
      gr,
      TxDb = txdb,
      tssRegion = c(-3000, 3000),
      annoDb = "org.Mm.eg.db"
    )
    peakAnno <- as.data.frame(peakAnno)
    peakAnno <- peakAnno[which(peakAnno$distanceToTSS==0),]
    # peakAnno <- peakAnno[grepl("Promoter", peakAnno$annotation), ]
    genes <- unique(peakAnno$SYMBOL)
    t_SE_RNA <- SE_RNA[which(SE_RNA$Geneid %in% genes),c("logCPM","logFC")]
    t_SE_RNA$condition <- "senescence"
    
    t_OE_RNA <- OE_RNA[which(OE_RNA$Geneid %in% genes),c("logCPM","LogFC.oe.vec")]
    colnames(t_OE_RNA) <- c("logCPM","logFC")
    t_OE_RNA$condition <- condition
    test <- t.test(t_SE_RNA$logFC,t_OE_RNA$logFC)
    t_pvalue_summary <- data.frame(condition=condition,quadrant=quadrant,pvalue=test$p.value)
    pvalue_summary <- rbind(pvalue_summary,t_pvalue_summary)
    t_logFC_summary <- rbind(t_OE_RNA,t_SE_RNA)
    t_logFC_summary$quadrant <- quadrant
    logFC_summary <- rbind(logFC_summary,t_logFC_summary)
  }
  logFC_summary$quadrant <- factor(logFC_summary$quadrant,levels=c("first","second","third","fourth"))
  p_list[[condition]] <- ggplot(logFC_summary, aes(x = quadrant, y = logFC, fill = condition)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
                alpha = 0.4, size = 0.8) +
    theme_bw() + 
    # ylim(-5,5) +
    ggtitle(condition) +
    labs(x = NULL, y = "logFC", color = "condition")
}

plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
combined_p <- plot_a_list(p_list,no_of_rows=1,no_of_cols=3)
ggsave(paste0("result/MEF_OE/H3K27me3_corr_with_senescence_gene_logFC.png"),combined_p,width=18,height=6)


