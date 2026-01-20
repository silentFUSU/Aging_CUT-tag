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
to_plot <- to_plot[which(to_plot$Significant=="Up" & to_plot$H3K27me3_Significant == "Down"),]

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
genelist <- bitr(to_plot$X,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_GO <- enrichGO( genelist$ENTREZID,#GO富集分析
                              OrgDb = GO_database,
                              keyType = "ENTREZID",#设定读取的gene ID类型
                              ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                              pvalueCutoff = 0.05,#设定p值阈值
                              qvalueCutoff = 0.05,#设定q值阈值
                              readable = T)
barplot(genelist_GO,label_format = 50,showCategory = 20)
genelist_GO <- pairwise_termsim(genelist_GO)  
p <- emapplot(genelist_GO, showCategory = 30) 
ggsave("result/Sup_figures/MEF_H3K27me3_decreased_RNA_increased_GO_network.pdf",p,width = 7,height = 6)
