rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
library(ggrepel)
library(grid)
library(ggsignif)
library(data.table)
library(dplyr)
library(GenomeInfoDb)
library("GenomicRanges")
library(genomation)


peaks <- read.table("data/samples/MEF/H3K27me3/bed/H3K27me3_young_merge-W5000-G10000-E100.bed") 
peaks <- as.data.table(peaks)
setDT(peaks)
setkey(peaks,V1,V2,V3)
samples <- c("MEF_mid_age","MEF_EZH2_inhibit")
summary <- data.frame()
for(sample in samples){
  gene <- read.csv(paste0("data/samples/RNA/",sample,"/diff_expression_gene_change_filter_bar.csv"))
  gene_increased <- gene[which(gene$Significant=="Down"),]
  gene_increased <- gene_increased[order(gene_increased$fdr),]
  gene_increased <- gene_increased[1:1000,]
  TSS <- read.table("~/ref_data/for_normal_mapping/mm10/mm10_refseq_TSS.bed")
  TSS <- TSS[which(TSS$V6 %in% gene_increased$X),]
  
  TSS <- as.data.table(TSS[,c(1:3,6)])
  setDT(TSS)
  setkey(TSS,V1,V2,V3)
  
  overlaps <- foverlaps(TSS, peaks, type = "any", nomatch = 0L)  
  
  t_summary <- data.frame(sample=sample,
                          overlap=length(unique(overlaps$V6))/nrow(gene_increased)*100,
                          other=(100-length(unique(overlaps$V6))/nrow(gene_increased)*100))
  summary <- rbind(summary,t_summary)
}


to_plot <- reshape2::melt(summary)
ggplot(to_plot, aes(x = sample, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  labs(title = "increased gene overlap with H3K27me3 peak in p2 samples", x = "Sample", y = "Value") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")  # Optional: to set colors for fill


### increasd gene overlap with H3K27me3 peaks
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
gene <- read.csv(paste0("data/samples/RNA/",sample,"/diff_expression_gene_change_filter_bar.csv"))
gene_increased <- gene[which(gene$Significant=="Up"),]
TSS <- read.table("~/ref_data/for_normal_mapping/mm10/mm10_refseq_TSS.bed")
TSS <- TSS[which(TSS$V6 %in% gene_increased$X),]

TSS <- as.data.table(TSS[,c(1:3,6)])
setDT(TSS)
setkey(TSS,V1,V2,V3)

overlaps <- foverlaps(TSS, peaks, type = "any", nomatch = 0L)  
genelist <- unique(overlaps$V6)
genelist <- bitr(genelist,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_GO <- enrichGO( genelist$ENTREZID,#GO富集分析
                              OrgDb = GO_database,
                              keyType = "ENTREZID",#设定读取的gene ID类型
                              ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                              pvalueCutoff = 0.05,#设定p值阈值
                              qvalueCutoff = 0.05,#设定q值阈值
                              readable = T)
barplot(genelist_GO,label_format = 50,showCategory = 20)




