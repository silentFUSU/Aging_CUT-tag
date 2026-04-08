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
quadrant_summary <- list(first=data.frame(),second=data.frame(),third=data.frame(),fourth=data.frame())
for(condition in conditions){
  df <- read.csv(paste0("data/samples/MEF_OE/H3K27me3/",condition,"/H3K27me3_",condition,"_10kb_bins_diff_after_remove_batch_effect.csv"))

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
    genes <- as.data.frame(unique(peakAnno$SYMBOL))
    genes$condition <- condition
    colnames(genes)[1] <- "gene"
    quadrant_summary[[quadrant]] <- rbind(quadrant_summary[[quadrant]],genes)
    }
}
quadrant_summary_count <- list()
quadrant_summary_tissue <- list()
for(quadrant in c("first","second","third","fourth")){
  quadrant_summary_count[[quadrant]] <- quadrant_summary[[quadrant]] %>%
    count(gene)
  quadrant_summary_tissue[[quadrant]] <- quadrant_summary[[quadrant]] %>%   
    group_by(gene) %>%   
    summarise(condition_content = paste(unique(condition), collapse = "/"))  
  quadrant_summary_count[[quadrant]] <- merge(quadrant_summary_count[[quadrant]],quadrant_summary_tissue[[quadrant]],by="gene")
  write.csv(quadrant_summary_count[[quadrant]],paste0("data/samples/MEF_OE/H3K27me3/correlation_with_senescence_",quadrant,"_quadrant_common_genes.csv"))
  }

