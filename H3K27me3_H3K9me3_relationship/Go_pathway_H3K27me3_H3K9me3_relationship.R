rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(edgeR)
library(ggplot2)
library(ChIPseeker)
library(EnsDb.Mmusculus.v79)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
library(limma)

GO_database <- 'org.Mm.eg.db'
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
tissues = c("brain","liver","testis","colon","kidney","lung","spleen","muscle","Hip","cecum","bonemarrow","heart","thymus")
peaks_list <-list()
GO_database <- 'org.Mm.eg.db'
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
quadrants <- c("second","fourth")
bin_size<-"10kb"
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  quadrant <-"second"
  peak <- read.table(paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/",bin_size,"_all_significant_",quadrant,"_quadrant.bed"))
  peak_obj <- GRanges(seqnames = peak$V1,   
                          ranges = IRanges(start = peak$V2, end = peak$V3))
  peak_anno <- annotatePeak(peak_obj, tssRegion=c(-3000, 3000),
                                TxDb=txdb, annoDb="org.Mm.eg.db")
  peak_anno <- unique(as.data.frame(peak_anno))
  genelist <- bitr(peak_anno$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  peaks_list[[i]]<-genelist
  names(peaks_list)[i] <- tissue
}

ck <- compareCluster(geneCluster = peaks_list, fun = enrichKEGG)
p<- dotplot(ck,show=3,label_format = 100)+theme(axis.text.x = element_text(angle = 45, hjust = 1)) +labs(x=NULL)
ck <- compareCluster(geneCluster = peaks_list[[2]], fun = enrichGO,OrgDb = GO_database, keyType = "ENTREZID")
ggsave(paste0("result/all/diff/",antibody,"/1kb_decrease_promoter_Go_dotplot.png"),p,width = 13,height = 10)
