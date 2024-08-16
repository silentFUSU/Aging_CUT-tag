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
antibody <- "H3K4me3"
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
histone <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_1kb_bins_diff_after_remove_batch_effect.csv"))
rna <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_nodup.csv"))

rna_gene_change <- rna$X[which(rna$Significant != "Stable")]
colnames(rna)[1]<-"SYMBOL"
peak_obj <- GRanges(seqnames = histone$Chr,   
                    ranges = IRanges(start = histone$Start, end = histone$End))
peak_anno <- annotatePeak(peak_obj, tssRegion=c(-3000, 3000),
                          TxDb=txdb, annoDb="org.Mm.eg.db")
peak_anno <- unique(as.data.frame(peak_anno))
peak_anno <- peak_anno[which(str_detect(peak_anno$annotation,"Promoter")),]
peak_anno <- peak_anno[which(peak_anno$SYMBOL %in% rna_gene_change),]

peak_anno$label <- paste0(peak_anno$seqnames,"-",peak_anno$start,"-",peak_anno$end)
histone$label <- paste0(histone$Chr,"-",histone$Start,"-",histone$End)
peak_anno <- merge(peak_anno,histone[,c("LogFC.old.young","label","Significant_bar")],by="label")
peak_anno <- merge(peak_anno,rna[,c("SYMBOL","logFC")],by="SYMBOL")
peak_anno <- peak_anno[which(peak_anno$Significant_bar != "Stable" & peak_anno$distanceToTSS==0),]

ggplot(peak_anno, aes(x = logFC, y = LogFC.old.young)) +
  # 坐标轴
  labs(x="RNA logFC",
       y=paste0(antibody," logFC")) +
  geom_point(color="grey")+
  theme(legend.position = "bottom",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme_bw()+theme(text = element_text(size = 18))+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  annotate("text",label = paste0(nrow(peak_anno[which(peak_anno$LogFC.old.young>0 & peak_anno$logFC>0),])),x=1, y=2,colour="#00b8a9",size=5)+
  annotate("text",label = paste0(nrow(peak_anno[which(peak_anno$LogFC.old.young<0 & peak_anno$logFC<0),])),x=-1, y=-2,colour="#ff9a00",size=5)+
  annotate("text",label = paste0(nrow(peak_anno[which(peak_anno$LogFC.old.young>0 & peak_anno$logFC<0),])),x=-1, y=2,colour="#f6416c",size=5)+
  annotate("text",label = paste0(nrow(peak_anno[which(peak_anno$LogFC.old.young<0 & peak_anno$logFC>0),])),x=1, y=-2,colour="#48466d",size=5)

