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
  CPM_summary <-data.frame()
  
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
    p2_cols <- grep("p2", colnames(SE_RNA), value = TRUE)
    p10_cols <- grep("p10", colnames(SE_RNA), value = TRUE)
    SE_RNA_p2 <- SE_RNA[which(SE_RNA$Geneid %in% genes), c("Geneid",p2_cols)]
    SE_RNA_p10 <- SE_RNA[which(SE_RNA$Geneid %in% genes), c("Geneid",p10_cols)]
    SE_RNA_p2$CPM <- rowMeans(SE_RNA_p2[, 2:3, drop = FALSE], na.rm = TRUE)
    SE_RNA_p2$condition <- "p2"
    
    SE_RNA_p10$CPM <- rowMeans(SE_RNA_p10[, 2:3, drop = FALSE], na.rm = TRUE)
    SE_RNA_p10$condition <- "p10"
    
    vec_cols <- grep(".vec.", colnames(OE_RNA), value = TRUE)
    oe_cols <- grep(".oe.u", colnames(OE_RNA), value = TRUE)
    OE_RNA_vec <- OE_RNA[which(OE_RNA$Geneid %in% genes), c("Geneid",vec_cols)]
    OE_RNA_oe <- OE_RNA[which(OE_RNA$Geneid %in% genes), c("Geneid",oe_cols)]
    OE_RNA_vec$CPM <- rowMeans(OE_RNA_vec[, 2:3, drop = FALSE], na.rm = TRUE)
    OE_RNA_vec$condition <- "vector"
    OE_RNA_oe$CPM <- rowMeans(OE_RNA_oe[, 2:3, drop = FALSE], na.rm = TRUE)
    OE_RNA_oe$condition <- condition
    
    t_CPM_summary <- rbind(SE_RNA_p2[,c("Geneid","CPM","condition")],SE_RNA_p10[,c("Geneid","CPM","condition")])
    t_CPM_summary <- rbind(t_CPM_summary,OE_RNA_vec[,c("Geneid","CPM","condition")])
    t_CPM_summary <- rbind(t_CPM_summary,OE_RNA_oe[,c("Geneid","CPM","condition")])
    t_CPM_summary$quadrant <- quadrant
    CPM_summary <- rbind(CPM_summary,t_CPM_summary)
    }
  CPM_summary$quadrant <- factor(CPM_summary$quadrant,levels=c("first","second","third","fourth"))
  CPM_summary$log2CPM <- log2(CPM_summary$CPM)
  CPM_summary$condition <- factor(CPM_summary$condition,levels=c("p2","p10","vector",condition))
  p_list[[condition]] <- ggplot(CPM_summary, aes(x = quadrant, y = log2CPM, fill = condition)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
    theme_bw() + 
    # ylim(-5,5) +
    ggtitle(condition) +
    labs(x = NULL, y = "log2(CPM)", color = "condition")
}

plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
combined_p <- plot_a_list(p_list,no_of_rows=1,no_of_cols=3)
ggsave(paste0("result/MEF_OE/H3K27me3_corr_with_senescence_gene_logCPM.png"),combined_p,width=18,height=6)


