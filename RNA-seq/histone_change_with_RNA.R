rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(reshape2)
library(ChIPseeker)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(dplyr)
tissue <- "CB"
antibody <- "H3K27me3"
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
histone <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_10kb_bins_diff_after_remove_batch_effect.csv"))
increase_histone <- histone[which(histone$Significant_bar=="Up"),]
peak_obj <- GRanges(seqnames = decrease_histone$Chr,   
                    ranges = IRanges(start = decrease_histone$Start, end = decrease_histone$End))
peak_anno <- annotatePeak(peak_obj, tssRegion=c(-3000, 3000),
                          TxDb=txdb, annoDb="org.Mm.eg.db")
peak_anno <- unique(as.data.frame(peak_anno))

rna <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_nodup.csv"))
rna <- rna[which(rna$X %in% peak_anno$SYMBOL),]

young_rna <- rna[,c(1,2,4)]
young_rna$year <- "young"
old_rna <- rna[,c(1,3,5)]
old_rna$year <- "old"
young_rna <- melt(young_rna)
old_rna <- melt(old_rna)
to_plot<- rbind(young_rna,old_rna)

ggplot(to_plot,aes(x=year,y=log10(value),fill = year))+
  geom_boxplot()
  scale_fill_brewer(palette="Set3")+
  guides(fill = FALSE) +
  ggtitle(paste0(antibody," TSSE"))+
  theme_bw()+theme(text = element_text(size = 18))+xlab("")+labs(fill = "", color = "") +ylab("tsse")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
