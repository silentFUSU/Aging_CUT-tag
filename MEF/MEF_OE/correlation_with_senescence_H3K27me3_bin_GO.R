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
  p_list <- list()
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
    genes <- bitr(genes,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
    genes_GO <- enrichGO( genes$ENTREZID,#GO富集分析
                          OrgDb = GO_database,
                          keyType = "ENTREZID",#设定读取的gene ID类型
                          ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                          pvalueCutoff = 0.05,#设定p值阈值
                          qvalueCutoff = 0.05,#设定q值阈值
                          readable = T)
    # if(nrow(genes_GO) >0){
    #   p_list[[quadrant]] <- barplot(genes_GO,label_format = 50,showCategory = 10)+ggtitle(paste0(condition," ",quadrant))
    # }else{
    #   p_list[[quadrant]] <- ggplot()
    # }
    if(nrow(genes_GO) > 0){
      GO_result <- genes_GO@result
      GO_result <- GO_result[order(GO_result$p.adjust),]
      GO_result <- GO_result[c(1:3),]
      GO_result$p.adjust <- -log10(GO_result$p.adjust)
      GO_result$label <- paste0(GO_result$ID," ",GO_result$Description)
      GO_result <- GO_result[,c("label","p.adjust")]
      p_list[[quadrant]] <- ggplot(GO_result, aes(x = p.adjust, y = reorder(label, p.adjust))) +
        geom_bar(stat = "identity") +
        labs(x = "P Adjust", y = "Label") +
        theme_bw() +
        theme(axis.text.y = element_text(size = 10))
    }
    
  }
  plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
    
    patchwork::wrap_plots(master_list_with_plots, 
                          nrow = no_of_rows, ncol = no_of_cols)
  }
  combined_p <- plot_a_list(p_list,no_of_rows=length(p_list),no_of_cols=1)
  ggsave(paste0("result/MEF_OE/",condition,"_H3K27me3_corr_with_senescence_GO_TSS.pdf"),combined_p,width=8,height=15)
}




