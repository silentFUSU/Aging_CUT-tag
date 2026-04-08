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
condition <- "MEF_Cbx7"
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'

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
quadrant <- "fourth"
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
genes <- unique(peakAnno$SYMBOL)

# common_gene <- read.csv("data/samples/MEF_OE/H3K27me3/correlation_with_senescence_fourth_quadrant_common_genes.csv")
# genes <- common_gene$gene[which(common_gene$n==3)]
OE_RNA <- read.csv(paste0("data/samples/RNA/MEF_OE_RNA/MEF_Vector_",condition,"_diff_expression_gene_strict_filter_bar.csv"))
SE_RNA <- read.csv(paste0("data/samples/RNA/MEF/diff_expression_gene_strict_filter_bar.csv"))
colnames(SE_RNA)[1] <- "Geneid"
logFC_summary <-data.frame()


t_SE_RNA <- SE_RNA[which(SE_RNA$Geneid %in% genes),c("Geneid","logCPM","logFC")]
t_SE_RNA$condition <- "senescence"

t_OE_RNA <- OE_RNA[which(OE_RNA$Geneid %in% genes),c("Geneid","logCPM","LogFC.oe.vec")]
colnames(t_OE_RNA) <- c("Geneid","logCPM","logFC")
t_OE_RNA$condition <- condition

t_logFC_summary <- rbind(t_OE_RNA,t_SE_RNA)
t_logFC_summary$quadrant <- quadrant
logFC_summary <- rbind(logFC_summary,t_logFC_summary)

genes_to_label <- c("Cdkn2a", "Msx2", "Hoxc13")

genes_to_label <- c("Cdkn2a", "Msx2", "Hoxc13")

df_lab <- logFC_summary %>%
  filter(Geneid %in% genes_to_label)

pd <- position_dodge(width = 0.8)
pjd <- position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8)

p <- ggplot(logFC_summary, aes(x = condition, y = logFC, fill = condition)) +
  geom_boxplot(outlier.shape = NA, position = pd) +
  geom_point(data = df_lab, position = pjd, size = 4, shape = 21,
             color = "black", stroke = 0.4) +
  geom_text_repel(
    data = df_lab,
    aes(label = Geneid),
    position = pjd,
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    min.segment.length = 0,
    max.overlaps = Inf
  ) +
  ylim(-5,5)+
  theme_bw() +
  labs(x = NULL, y = "logFC", color = "condition")
ggsave("result/figures/Cbx7_senescence_fourth_genes_logFC.pdf",p,width = 6,height = 8)
