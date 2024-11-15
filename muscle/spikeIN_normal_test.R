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
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'

antibody <- "H3K36me3"
bin_size <- ifelse(antibody %in% c("H3K36me3","H3K9me3","H3K27me3"), "10kb", "1kb")
tab <- read.table(paste0("/mnt/transposon1/zhangyanxiaoLab/suzhuojie/project/Aging_CUT_TAG/data/muscle_cut_tag_test/",antibody,"/",antibody,"_",bin_size,"_bins.counts"), header = T)
counts <- tab[,c(7:ncol(tab))]
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+|DYQ[0-9]+).*"
colnames(counts) <- gsub(pattern, "\\1", colnames(counts))
samples <- c("228","229","116","117","225","237","226","230","234","235","101","110","102","111")
age <- c("3m","3m","24m","24m","24m","3m","24m","24m","3m","3m","3m","24m","3m","24m")

group <- paste0(colnames(counts),"-",samples,"-",age)
y= DGEList(counts=counts,group = group)
y$samples$age <- age
keep = which(rowSums(cpm(y)>1)>=4)
y = y[keep,]
logCPMs <- cpm(y, log = TRUE)
logCPMs_corrected <- limma::removeBatchEffect(logCPMs,
                                              batch = c("batch1","batch1","batch1","batch1","batch4","batch4","batch5","batch5","batch5","batch5","batch2","batch2","batch3","batch3"))

pca <- prcomp(t(logCPMs_corrected))
to_plot <- data.frame(pca$x, tissue = paste0(y$samples$group))
to_plot$rownames <- rownames(to_plot)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
to_plot$age <- factor(age, levels=c("3m","24m"))
ggplot(to_plot, aes(x=PC1, y=PC2, color=age)) + 
  geom_point(size=5) +theme_bw()+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
  ggtitle(antibody)+
  geom_text_repel(data = to_plot,
    aes(x = PC1, y = PC2, label = group, color = age),  
    size = 5,  
    box.padding = unit(0.35, "lines"),  
    point.padding = unit(0.3, "lines")  
  )

cor_matrix <- cor(logCPMs_corrected)
pheatmap::pheatmap(cor_matrix)

pca <- prcomp(t(logCPMs))
to_plot <- data.frame(pca$x, tissue = paste0(y$samples$group))
to_plot$rownames <- rownames(to_plot)
to_plot$age <- factor(age, levels=c("3m","24m"))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
ggplot(to_plot, aes(x=PC1, y=PC2, color=age)) + 
  geom_point(size=5) +theme_bw()+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
  ggtitle(antibody)+
  geom_text_repel(data = to_plot,
                  aes(x = PC1, y = PC2, label = group, color = age),  
                  size = 5,  
                  box.padding = unit(0.35, "lines"),  
                  point.padding = unit(0.3, "lines")  
  )
cor_matrix <- cor(logCPMs)
annotation <- data.frame(samples=colnames(counts),age=age)
rownames(annotation) <- annotation$samples
annotation <- annotation[,-1,drop=F]
annotation$age <- factor(annotation$age,levels=c("3m","24m"))
pheatmap::pheatmap(cor_matrix,annotation_row = annotation)

tab <- read.table(paste0("/mnt/transposon1/zhangyanxiaoLab/suzhuojie/project/Aging_CUT_TAG/data/muscle_cut_tag_test/",antibody,"/",antibody,"_",bin_size,"_bins.counts"), header = T)
counts <- tab[,c(7:ncol(tab))]
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+|DYQ[0-9]+).*"
colnames(counts) <- gsub(pattern, "\\1", colnames(counts))
counts <- counts[,-c(1:4)]
samples <- c("225","237","101","110","102","111")
age <- c("24m","3m","3m","24m","3m","24m")
group <- paste0(colnames(counts),"-",samples,"-",age)
y= DGEList(counts=counts,group = group)
y$samples$age <- age
keep = which(rowSums(cpm(y)>1)>=4)
y = y[keep,]
logCPMs <- cpm(y, log = TRUE)
logCPMs_corrected <- limma::removeBatchEffect(logCPMs,
                                              batch = c("batch4","batch4","batch2","batch2","batch3","batch3"))
pca <- prcomp(t(logCPMs))
to_plot <- data.frame(pca$x, tissue = paste0(y$samples$group))
to_plot$rownames <- rownames(to_plot)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
to_plot$age <- factor(age, levels=c("3m","24m"))
ggplot(to_plot, aes(x=PC1, y=PC2, color=age)) + 
  geom_point(size=5) +theme_bw()+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
  ggtitle(antibody)+
  geom_text_repel(data = to_plot,
                  aes(x = PC1, y = PC2, label = group, color = age),  
                  size = 5,  
                  box.padding = unit(0.35, "lines"),  
                  point.padding = unit(0.3, "lines")  
  )





annotation <- data.frame(samples=colnames(counts),age=age)
rownames(annotation) <- annotation$samples
annotation <- annotation[,-1,drop=F]
annotation$age <- factor(annotation$age,levels=c("3m","24m"))
previous <- read.csv(paste0("data/samples/muscle/",antibody,"/",antibody,"_",bin_size,"_bins_diff.csv"))
tab <- tab[keep,c(1:6)]
logCPMs_corrected <- cbind(tab,logCPMs_corrected)
increase <- previous$Geneid[which(previous$Significant=="Up")]
logCPMs_corrected_increase <- logCPMs_corrected[which(logCPMs_corrected$Geneid %in% increase),]
logCPMs_corrected_increase <- logCPMs_corrected_increase[,c(1:10,12,11,15,16,13,14,17,19,18,20)]
pheatmap::pheatmap(logCPMs_corrected_increase[,c(7:ncol(logCPMs_corrected_increase))],cluster_rows = T,annotation_col = annotation,cluster_cols = T,show_rownames = F,main = paste0(antibody," previous increase"))

decrease <- previous$Geneid[which(previous$Significant=="Down")]
logCPMs_corrected_decrease <- logCPMs_corrected[which(logCPMs_corrected$Geneid %in% decrease),]
logCPMs_corrected_decrease <- logCPMs_corrected_decrease[,c(1:10,12,11,15,16,13,14,17,19,18,20)]
pheatmap::pheatmap(logCPMs_corrected_decrease[,c(7:ncol(logCPMs_corrected_decrease))],cluster_rows = T,annotation_col = annotation,cluster_cols = F,show_rownames = F,main = paste0(antibody," previous decrease"))

to_plot <- logCPMs_corrected_increase[,c(1,7:ncol(logCPMs_corrected_increase))]
to_plot <- melt(to_plot)
to_plot$variable <- factor(to_plot$variable, levels =unique(to_plot$variable))
to_plot$age <- "young"
to_plot$age[which(to_plot$variable %in% unique(to_plot$variable)[c(3,4,6,9,10,13,14)])] <- "old"
to_plot$age <- factor(to_plot$age,levels=c("young","old"))
ggplot(to_plot, aes(x = variable, y = value,fill=age)) +  
  geom_boxplot() +  
  labs(title = "previous H3K36me3 increase site",  
       x = NULL,  
       y = "log2(CPM)") +
  theme_bw() +
  theme(  
    axis.text.x = element_text(angle = 45, hjust = 1)  
  ) 

to_plot <- logCPMs_corrected_decrease[,c(1,7:ncol(logCPMs_corrected_decrease))]
to_plot <- melt(to_plot)
to_plot$variable <- factor(to_plot$variable, levels =unique(to_plot$variable))
to_plot$age <- "young"
to_plot$age[which(to_plot$variable %in% unique(to_plot$variable)[c(3,4,6,9,10,13,14)])] <- "old"
to_plot$age <- factor(to_plot$age,levels=c("young","old"))
ggplot(to_plot, aes(x = variable, y = value,fill=age)) +  
  geom_boxplot() +  
  labs(title = "previous H3K36me3 decrease site",  
       x = NULL,  
       y = "log2(CPM)") +
  theme_bw() +
  theme(  
    axis.text.x = element_text(angle = 45, hjust = 1)  
  ) 
t.test(to_plot$value[which(to_plot$variable %in% unique(to_plot$variable)[c(2)] & to_plot$age == "young")],
       to_plot$value[which(to_plot$variable %in% unique(to_plot$variable)[c(3)] & to_plot$age == "old")])
t.test(to_plot$value[which(to_plot$variable %in% unique(to_plot$variable)[c(5:8)] & to_plot$age == "young")],
       to_plot$value[which(to_plot$variable %in% unique(to_plot$variable)[c(5:8)] & to_plot$age == "old")])

tab <- read.table(paste0("/mnt/transposon1/zhangyanxiaoLab/suzhuojie/project/Aging_CUT_TAG/data/muscle_cut_tag_test/",antibody,"/",antibody,"_",bin_size,"_bins.counts"), header = T)
counts <- tab[,c(7:ncol(tab))]
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+|DYQ[0-9]+).*"
colnames(counts) <- gsub(pattern, "\\1", colnames(counts))
counts <- counts[,c(7:10)]
age <- c("old","old","young","young")
y= DGEList(counts=counts,group=age)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
y$samples$group <- factor(y$samples$group,c("young","old"))
y <- calcNormFactors(y)

design <- model.matrix(~group, y$samples)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)
lrt = glmLRT(fit_tag, coef = 2)
tab<-tab[keep,]
out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
            "PValue.old-young"=lrt$table$PValue,"FDR.old-young"= p.adjust(lrt$table$PValue,method="BH"),
            "LogFC.old-young"=lrt$table$logFC)
out$Significant <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$`LogFC.old-young`) >= 0, 
                          ifelse(out$`LogFC.old-young` > 0, "Up", "Down"), "Stable")
colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
bin_in_peaks <- read.table("data/samples/muscle/H3K36me3/bed/H3K36me3_10kb_in_young_old_merge-W1000-G3000-E100.bed")
out <- out[which(out$Geneid %in% bin_in_peaks$V4),]
ggplot(
  out, aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`))) +
  geom_point(aes(color = Significant), size=2) +
  scale_color_manual(values = colour) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (p-value)") +
  theme_bw()+
  theme(text = element_text(size = 20),legend.position = "none")+
  ggtitle(paste0("Muscle H3K36me3"))+
  annotate("text", x = min(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
  annotate("text", x = max(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
# write.csv(out,paste0("data/raw_data/muscle_cut_tag_test/",antibody,"/",antibody,"_",bin_size,"_bins_diff_ConA.csv"),row.names = F)


table(out$Significant)
increase <- out$Geneid[which(out$Significant=="Up")]
logCPMs_corrected_increase <- logCPMs_corrected[which(logCPMs_corrected$Geneid %in% increase),]
logCPMs_corrected_increase <- logCPMs_corrected_increase[,c(1:10,11,13,12,14)]
pheatmap::pheatmap(logCPMs_corrected_increase[,c(7:ncol(logCPMs_corrected_increase))],cluster_rows = T,cluster_cols = F,show_rownames = F,main = paste0(antibody," increase"))

increase <- out[which(out$Significant=="Up"),]
peak_obj <- GRanges(seqnames = increase$Chr,   
                    ranges = IRanges(start = increase$Start, end = increase$End))
peak_anno <- annotatePeak(peak_obj, tssRegion=c(-3000, 3000),
                          TxDb=txdb, annoDb="org.Mm.eg.db")
peak_anno <- unique(as.data.frame(peak_anno))
genelist_up <- bitr(peak_anno$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_up_GO <- enrichGO( genelist_up$ENTREZID,#GO富集分析
                            OrgDb = GO_database,
                            keyType = "ENTREZID",#设定读取的gene ID类型
                            ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                            pvalueCutoff = 0.05,#设定p值阈值
                            qvalueCutoff = 0.05,#设定q值阈值
                            readable = T)
# genelist_up_GO_df <- genelist_up_GO@result
barplot(genelist_up_GO)



decrease <- out$Geneid[which(out$Significant=="Down")]
logCPMs_corrected_decrease <- logCPMs_corrected[which(logCPMs_corrected$Geneid %in% decrease),]
logCPMs_corrected_decrease <- logCPMs_corrected_decrease[,c(1:10,11,13,12,14)]
pheatmap::pheatmap(logCPMs_corrected_decrease[,c(7:ncol(logCPMs_corrected_decrease))],cluster_rows = T,cluster_cols = F,show_rownames = F,main = paste0(antibody," decrease"))


tab <- read.table(paste0("data/raw_data/muscle_cut_tag_test/",antibody,"/",bin_size,"_bins.counts"), header = T)
counts <- tab[,c(7:ncol(tab))]
counts <- counts[,c(3,4,6,8)]
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
colnames(counts) <- gsub(pattern, "\\1", colnames(counts))
group <- c("bing","bing","xiao","xiao")
y= DGEList(counts=counts,group=group)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
y$samples$group <- factor(y$samples$group,c("xiao","bing"))
y <- calcNormFactors(y)
design <- model.matrix(~group, y$samples)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)
lrt = glmLRT(fit_tag, coef = 2)
tab<-tab[keep,]
out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
            "PValue.old-young"=lrt$table$PValue,"FDR.old-young"= p.adjust(lrt$table$PValue,method="BH"),
            "LogFC.old-young"=lrt$table$logFC)
out$Significant <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$`LogFC.old-young`) >= 0, 
                          ifelse(out$`LogFC.old-young` > 0, "Up", "Down"), "Stable")

increase <- out[which(out$Significant=="Up"),]
peak_obj <- GRanges(seqnames = increase$Chr,   
                    ranges = IRanges(start = increase$Start, end = increase$End))
peak_anno <- annotatePeak(peak_obj, tssRegion=c(-3000, 3000),
                          TxDb=txdb, annoDb="org.Mm.eg.db")
peak_anno <- unique(as.data.frame(peak_anno))
genelist_up <- bitr(peak_anno$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_up_GO <- enrichGO( genelist_up$ENTREZID,#GO富集分析
                            OrgDb = GO_database,
                            keyType = "ENTREZID",#设定读取的gene ID类型
                            ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                            pvalueCutoff = 0.05,#设定p值阈值
                            qvalueCutoff = 0.05,#设定q值阈值
                            readable = T)
genelist_up_GO_df <- genelist_up_GO@result
barplot(genelist_up_GO)

decrease <- out[which(out$Significant=="Down"),]
peak_obj <- GRanges(seqnames = decrease$Chr,   
                    ranges = IRanges(start = decrease$Start, end = decrease$End))
peak_anno <- annotatePeak(peak_obj, tssRegion=c(-3000, 3000),
                          TxDb=txdb, annoDb="org.Mm.eg.db")
peak_anno <- unique(as.data.frame(peak_anno))
genelist_down <- bitr(peak_anno$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_down_GO <- enrichGO( genelist_down$ENTREZID,#GO富集分析
                            OrgDb = GO_database,
                            keyType = "ENTREZID",#设定读取的gene ID类型
                            ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                            pvalueCutoff = 0.05,#设定p值阈值
                            qvalueCutoff = 0.05,#设定q值阈值
                            readable = T)
genelist_down_GO_df <- genelist_down_GO@result
barplot(genelist_down_GO)
colour<- list(c("grey"),c("grey","red"),c("blue","grey","red"))
ggplot(
  # 数据、映射、颜色
  out, aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`))) +
  geom_point(aes(color = Significant), size=2,show.legend = FALSE) +
  scale_color_manual(values = colour[[nrow(as.data.frame(table(out$Significant)))]]) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (fdr)") +
  theme_bw()+
  ggtitle("ZhangBing old mouse vs YanXiao old mouse")+
  theme(text = element_text(size = 15))+
  annotate("text", x = min(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
  annotate("text", x = max(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)


tab <- read.table(paste0("data/raw_data/muscle_cut_tag_test/",antibody,"/",bin_size,"_bins.counts"), header = T)
counts <- tab[,c(7:ncol(tab))]
counts <- counts[,c(1,2,5,7)]
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
colnames(counts) <- gsub(pattern, "\\1", colnames(counts))
group <- c("new_data","new_data","previous_data","previous_data")
y= DGEList(counts=counts,group=group)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
y$samples$group <- factor(y$samples$group,c("previous_data","new_data"))
y <- calcNormFactors(y)
design <- model.matrix(~group, y$samples)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)
lrt = glmLRT(fit_tag, coef = 2)
tab<-tab[keep,]
out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
            "PValue.old-young"=lrt$table$PValue,"FDR.old-young"= p.adjust(lrt$table$PValue,method="BH"),
            "LogFC.old-young"=lrt$table$logFC)
out$Significant <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$`LogFC.old-young`) >= 0, 
                          ifelse(out$`LogFC.old-young` > 0, "Up", "Down"), "Stable")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
increase <- out[which(out$Significant=="Up"),]
peak_obj <- GRanges(seqnames = increase$Chr,   
                    ranges = IRanges(start = increase$Start, end = increase$End))
peak_anno <- annotatePeak(peak_obj, tssRegion=c(-3000, 3000),
                          TxDb=txdb, annoDb="org.Mm.eg.db")
peak_anno <- unique(as.data.frame(peak_anno))
genelist_up <- bitr(peak_anno$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_up_GO <- enrichGO( genelist_up$ENTREZID,#GO富集分析
                            OrgDb = GO_database,
                            keyType = "ENTREZID",#设定读取的gene ID类型
                            ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                            pvalueCutoff = 0.05,#设定p值阈值
                            qvalueCutoff = 0.05,#设定q值阈值
                            readable = T)
genelist_up_GO_df <- genelist_up_GO@result
barplot(genelist_up_GO)

decrease <- out[which(out$Significant=="Down"),]
peak_obj <- GRanges(seqnames = decrease$Chr,   
                    ranges = IRanges(start = decrease$Start, end = decrease$End))
peak_anno <- annotatePeak(peak_obj, tssRegion=c(-3000, 3000),
                          TxDb=txdb, annoDb="org.Mm.eg.db")
peak_anno <- unique(as.data.frame(peak_anno))
genelist_down <- bitr(peak_anno$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_down_GO <- enrichGO( genelist_down$ENTREZID,#GO富集分析
                              OrgDb = GO_database,
                              keyType = "ENTREZID",#设定读取的gene ID类型
                              ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                              pvalueCutoff = 0.05,#设定p值阈值
                              qvalueCutoff = 0.05,#设定q值阈值
                              readable = T)
genelist_down_GO_df <- genelist_down_GO@result
barplot(genelist_down_GO)
colour<- list(c("grey"),c("grey","red"),c("blue","grey","red"))
ggplot(
  # 数据、映射、颜色
  out, aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`))) +
  geom_point(aes(color = Significant), size=2,show.legend = FALSE) +
  scale_color_manual(values = colour[[nrow(as.data.frame(table(out$Significant)))]]) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (fdr)") +
  theme_bw()+
  ggtitle("new young mouse vs previous young mouse")+
  theme(text = element_text(size = 15))+
  annotate("text", x = min(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
  annotate("text", x = max(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)


FRiP<-data.frame(sample_name = character(),
                 antibody = character(),
                 FRiP = numeric(),
                 stringsAsFactors = FALSE)  
antibody <- "H3K4me3"
bin_size <- ifelse(antibody %in% c("H3K36me3","H3K9me3","H3K27me3"), "10kb", "1kb")
tab <- read.table(paste0("data/raw_data/muscle_cut_tag_test/",antibody,"/",bin_size,"_bins.counts"), header = T)
counts <- tab[,c(7:ncol(tab))]
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
colnames(counts) <- gsub(pattern, "\\1", colnames(counts))
samples <- colnames(counts)
for(sample in samples){
  df <- read.delim(paste0("data/raw_data/muscle_cut_tag_test/",antibody,"/featurecount/",sample,"_remove_blacklist.counts.summary"),row.names = 1)
  t_FRiP<-data.frame(sample_name = sample,
                     antibody = antibody,
                     FRiP = df[1,1]/sum(df[,1]),
                     stringsAsFactors = FALSE)  
  FRiP <- rbind(FRiP,t_FRiP)
}
FRiP$age <- c("3m","3m","24m","24m","3m","24m","3m","24m")
FRiP$FRiP <- FRiP$FRiP*100
FRiP$age <- factor(FRiP$age,levels=c("3m","24m"))
FRiP$condition <- c("ConA","ConA","ConA","ConA","Paired_CT","Paired_CT","Paired_CT","Paired_CT")
FRiP$mouse_ID <-  c("228","229","116","117","101","110","102","111")
ggplot(FRiP,aes(x=condition,y=FRiP,color = age))+
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.7)+
  geom_text(aes(label = mouse_ID), position = position_jitter(width = 0.2), vjust = -1, size = 3) +
  scale_fill_brewer(palette="Set3")+
  ggtitle(paste0(antibody," FRiP"))+ylim(0,100)+
  theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")+labs(fill = "", color = "") +ylab("FRiP(%)")


data_path <- "data/raw_data/muscle_cut_tag_test/"
tissue <- "muscle"
rna_name <- function(tissue){
  if(tissue=="brain"){
    return("FC")
  }
  return(tissue)
}
Histone_relationship_with_RNA <- function(antibody,tissue){
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  GO_database <- 'org.Mm.eg.db'
  bin_size <- ifelse(antibody %in% c("H3K36me3","H3K9me3","H3K27me3"), "10kb", "1kb")
  histone <- read.csv(paste0(data_path,"/",antibody,"/",antibody,"_",bin_size,"_bins_diff_ConA.csv"))
  rna <- read.csv(paste0("data/samples/RNA/",rna_name(tissue),"/diff_expression_gene_nodup.csv"))
  
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
  peak_anno <- merge(peak_anno,histone[,c("LogFC.old.young","label","Significant")],by="label")
  peak_anno <- merge(peak_anno,rna[,c("SYMBOL","logFC")],by="SYMBOL")
  peak_anno <- peak_anno[which(peak_anno$Significant != "Stable" & peak_anno$distanceToTSS==0),]
  p<- ggplot() +
    geom_point(data=peak_anno, mapping=aes(logFC, LogFC.old.young),color = "grey",alpha=0.5) +  
    # geom_point(data=peak_anno[which(peak_anno$LogFC.old.young<0 & peak_anno$logFC>0),], mapping=aes(logFC, LogFC.old.young),color = "#48466d") +
    geom_point(data=peak_anno[which(peak_anno$LogFC.old.young>0 & peak_anno$logFC>0),], mapping=aes(logFC, LogFC.old.young),color = "#00b8a9") +
    geom_point(data=peak_anno[which(peak_anno$LogFC.old.young<0 & peak_anno$logFC<0),], mapping=aes(logFC, LogFC.old.young),color = "#ff9a00") +
    # 坐标轴
    labs(x="RNA log2(Fold Change)",
         y=paste0(antibody," log2(Fold Change)")) +
    geom_point(color="grey")+
    theme(legend.position = "bottom",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    theme_bw()+theme(text = element_text(size = 18))+
    ggtitle(tissue)+
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    annotate("text",label = paste0(nrow(peak_anno[which(peak_anno$LogFC.old.young>0 & peak_anno$logFC>0),])),x=10, y=4,colour="#00b8a9",size=5)+
    annotate("text",label = paste0(nrow(peak_anno[which(peak_anno$LogFC.old.young<0 & peak_anno$logFC<0),])),x=-5, y=-4,colour="#ff9a00",size=5)+
    annotate("text",label = paste0(nrow(peak_anno[which(peak_anno$LogFC.old.young>0 & peak_anno$logFC<0),])),x=-5, y=4,colour="#f6416c",size=5)+
    annotate("text",label = paste0(nrow(peak_anno[which(peak_anno$LogFC.old.young<0 & peak_anno$logFC>0),])),x=10, y=-4,colour="#48466d",size=5)
  print(p)
  return(p)
}
