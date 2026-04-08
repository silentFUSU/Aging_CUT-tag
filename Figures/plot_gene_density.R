rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(karyoploteR)
library(GenomicRanges) 
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
all.genes <- genes(txdb)
head(all.genes)

pdf("result/Sup_figures/chr17_gene_density.pdf", width = 8, height = 3)  # 可按需要改尺寸
kp <- plotKaryotype(genome = "mm10", main="gene_density", chromosomes = "chr17",plot.type=6,cex=1.8)
kp <- kpDataBackground(kp, color = "#FFFFFFAA")
kp <- kpPlotDensity(kp, all.genes,window.size = 0.5e6, r0=0, r1=0.8,data.panel="ideogram", col="black", border="black")
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density)
dev.off()

