rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(stringr)
library(dplyr)
library(dbplyr)
library(clusterProfiler)
library(GSVA)
library(enrichplot)
options(scipen =0)  
genes <- read.table("~/ref_data/for_normal_mapping/mm10/mm10_refseq_TSS.bed")
genes <- genes[,c("V1","V2","V3","V6")]
genes <- as.data.table(genes)
setDT(genes)
setkey(genes,V1,V2,V3)
kmeans <- c(paste0("kmeans",1:4),"Stable")
for(kmean in kmeans){
  regions <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation_add_stable.csv",row.names = 1)
  split_chr <- strsplit(as.character(regions$label), ":")  
  chr_column <- sapply(split_chr, `[[`, 1)  
  split_start_end <- strsplit(sapply(split_chr, `[[`, 2), "-")  
  start_column <- sapply(split_start_end, `[[`, 1)  
  end_column <- sapply(split_start_end, `[[`, 2)  
  regions <- data.frame(chr = chr_column,start = start_column, end = end_column, cluster=regions$cluster)
  regions$start <- as.numeric(regions$start)
  regions$end <- as.numeric(regions$end)
  regions <-regions[which(regions$cluster==kmean),]
  regions <- as.data.table(regions)
  setDT(regions)
  setkey(regions,chr,start,end)
  
  overlaps <-  foverlaps(genes, regions, type = "any", nomatch = 0L)  
  write.csv(unique(overlaps$V6), paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/Gene_covered_by_",kmean,".csv"))
  
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  GO_database <- 'org.Mm.eg.db'
  genelist <- bitr(unique(overlaps$V6),fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  genelist_GO <- enrichGO( genelist$ENTREZID,#GO富集分析
                           OrgDb = GO_database,
                           keyType = "ENTREZID",#设定读取的gene ID类型
                           ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                           pvalueCutoff = 0.05,#设定p值阈值
                           qvalueCutoff = 0.05,#设定q值阈值
                           readable = T)
  result <- as.data.frame(genelist_GO@result)
  write.csv(result, paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/Gene_covered_by_",kmean,"_GO_pathway.csv"))
  genelist_GO <- pairwise_termsim(genelist_GO)  
  p <- emapplot(genelist_GO, showCategory = 30) 
  ggsave(paste0("result/Sup_figures/H3K9me3_",kmean,"_peak_cover_genes_GO_pathway.pdf"),p,width = 12,height = 10)
  # barplot(genelist_GO,title = paste0("Genes covered by kmeans3"),label_format = 50,showCategory = 30)
}




