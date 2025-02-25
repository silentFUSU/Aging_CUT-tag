rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(edgeR)
library(ggrepel)
library(ChIPseeker)
library(GenomeInfoDb)
library(stringr)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(GenomicRanges)
library(reshape2)

antibody <- "H3K36me3"
bin_size <- ifelse(antibody %in% c("H3K36me3","H3K9me3","H3K27me3"), "10kb", "1kb")
tab <- read.table(paste0("/mnt/transposon1/zhangyanxiaoLab/suzhuojie/project/Aging_CUT_TAG/data/thymus_cut_tag_test/",antibody,"/",bin_size,"_bins.counts"),header=T)
tab2 <- read.table("data/raw_data_transposon2/20250218_LLX_Chip/H3K36me3_10kb_bins.counts",header=T)
tab <- merge(tab,tab2[,c(1,11,12)],by="Geneid")
rownames(tab) <- tab$Geneid
counts <- tab[,c(7:ncol(tab))]
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+|DYQ[0-9]+).*"
colnames(counts) <- gsub(pattern, "\\1", colnames(counts))
samples <- c("116","117","224","225","228","229","240","241","245","247","101","110","102","111","Chipseq","Chipseq")
age <- c("24m","24m","24m","24m","3m","3m","3m","3m","24m","24m","3m","24m","3m","24m","3m","24m")

group <- paste0(colnames(counts),"-",samples,"-",age)
y= DGEList(counts=counts,group = group)
y$samples$age <- age
keep = which(rowSums(cpm(y)>1)>=4)
y = y[keep,]
logCPMs <- cpm(y, log = TRUE)
logCPMs_corrected <- limma::removeBatchEffect(logCPMs,batch = c("batch1","batch1","batch1","batch1","batch1","batch1","batch4","batch4","batch4","batch4","batch2","batch2","batch3","batch3","batch5","batch5"))
colnames(logCPMs_corrected) <- paste0(colnames(logCPMs_corrected),"-",samples,"-",age)
pca <- prcomp(t(logCPMs_corrected))
to_plot <- data.frame(pca$x, tissue = paste0(y$samples$group))
to_plot$rownames <- rownames(to_plot)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
age <- factor(age,levels=c("3m","24m"))

ggplot(to_plot, aes(x=PC1, y=PC2, color=age)) + 
  geom_point(size=5) +theme_bw()+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
  ggtitle(antibody)+
  geom_text_repel(data = to_plot,
                  aes(x = PC1, y = PC2, label = group, color = age),  
                  size = 5,  
                  box.padding = unit(0.35, "lines"),  
                  point.padding = unit(0.3, "lines"),max.overlaps = Inf  
  )


# counts <- tab[,c(13:ncol(tab))]
# pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+|DYQ[0-9]+).*"
# colnames(counts) <- gsub(pattern, "\\1", colnames(counts))
# samples <- c("240","241","245","247","101","110","102","111","Chipseq","Chipseq")
# age <- c("3m","3m","24m","24m","3m","24m","3m","24m","3m","24m")
# 
# group <- paste0(colnames(counts),"-",samples,"-",age)
# y= DGEList(counts=counts,group = group)
# y$samples$age <- age
# keep = which(rowSums(cpm(y)>1)>=4)
# y = y[keep,]
# logCPMs <- cpm(y, log = TRUE)
# logCPMs_corrected <- limma::removeBatchEffect(logCPMs,batch = c("batch4","batch4","batch4","batch4","batch2","batch2","batch3","batch3","batch5","batch5"))
# colnames(logCPMs_corrected) <- paste0(colnames(logCPMs_corrected),"-",samples,"-",age)
# pca <- prcomp(t(logCPMs_corrected))
# to_plot <- data.frame(pca$x, tissue = paste0(y$samples$group))
# to_plot$rownames <- rownames(to_plot)
# percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
# use.pcs <- c(1,2)
# labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
# age <- factor(age,levels=c("3m","24m"))
# 
# ggplot(to_plot, aes(x=PC1, y=PC2, color=age)) + 
#   geom_point(size=5) +theme_bw()+
#   xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
#   ggtitle(antibody)+
#   geom_text_repel(data = to_plot,
#                   aes(x = PC1, y = PC2, label = group, color = age),  
#                   size = 5,  
#                   box.padding = unit(0.35, "lines"),  
#                   point.padding = unit(0.3, "lines")  
#   )

previous <- read.csv(paste0("data/samples/thymus/",antibody,"/",antibody,"_",bin_size,"_bins_diff.csv"))
increase <- previous$Geneid[which(previous$Significant=="Up")]
logCPMs_increase <- logCPMs_corrected[which(rownames(logCPMs_corrected) %in% increase),]
logCPMs_increase <- logCPMs_increase[,c(5:6,1:4,7:10,11,13,12,14,15,16)]
annotation <- data.frame(sample=colnames(logCPMs_corrected),age=age)
rownames(annotation) <- annotation$sample
annotation <- annotation[,"age",drop=F]
annotation$age <- factor(annotation$age, levels=c("3m","24m"))
pheatmap::pheatmap(logCPMs_increase,annotation = annotation,cluster_rows = T,cluster_cols = F,show_rownames = F,main = paste0(antibody," previous increase"),scale="row")

decrease <- previous$Geneid[which(previous$Significant=="Down")]
logCPMs_decrease <- logCPMs_corrected[which(rownames(logCPMs_corrected) %in% decrease),]
logCPMs_decrease <- logCPMs_decrease[,c(5:6,1:4,7:10,11,13,12,14,15,16)]
pheatmap::pheatmap(logCPMs_decrease,annotation = annotation,cluster_rows = T,cluster_cols = F,show_rownames = F,scale="row",main = paste0(antibody," previous decrease"))


to_plot <- logCPMs_decrease
to_plot <- reshape2::melt(to_plot)
to_plot$Var2 <- factor(to_plot$Var2, levels =unique(to_plot$Var2))
to_plot$age <- "young"
to_plot$age[which(to_plot$Var2 %in% unique(to_plot$Var2)[c(3,4,5,6,9,10,13,14,16)])] <- "old"
to_plot$age <- factor(to_plot$age,levels=c("young","old"))
to_plot$Var2 <- gsub("(^[A-Za-z0-9_]+)-.*", "\\1", to_plot$Var2) 
to_plot$Var2 <- factor(to_plot$Var2,levels = c("CKJ072","CKJ073","CKJ068","CKJ069","CKJ070","CKJ071","DYQ123","DYQ124","DYQ125","DYQ126","HJC_124","HJC_138","HJC_130","HJC_144","LLX923","LLX924"))
ggplot(to_plot, aes(x = Var2, y = value,fill=age)) +  
  geom_boxplot() +  
  labs(title = paste0("previous ",antibody," decrease site"),  
       x = NULL,  
       y = "log2(CPM)") +
  theme_bw() +
  theme(  
    axis.text.x = element_text(angle = 45, hjust = 1)  
  ) 
t.test(to_plot$value[which(to_plot$Var2 %in% unique(to_plot$Var2)[c(1,2)] & to_plot$age == "young")],
       to_plot$value[which(to_plot$Var2 %in% unique(to_plot$Var2)[c(3,4,5,6)] & to_plot$age == "old")])
t.test(to_plot$value[which(to_plot$Var2 %in% unique(to_plot$Var2)[c(7,9)] & to_plot$age == "young")],
       to_plot$value[which(to_plot$Var2 %in% unique(to_plot$Var2)[c(8,10)] & to_plot$age == "old")])