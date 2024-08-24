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
tab <- read.table(paste0("data/raw_data/muscle_cut_tag_test/",antibody,"/",bin_size,"_bins.counts"), header = T)
counts <- tab[,c(7:ncol(tab))]
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
colnames(counts) <- gsub(pattern, "\\1", colnames(counts))
samples <- c("228","229","116","117","101","110","102","111")
age <- c("3m","3m","24m","24m","3m","24m","3m","24m")

group <- paste0(colnames(counts),"-",samples,"-",age)
y= DGEList(counts=counts,group = group)
y$samples$age <- age
keep = which(rowSums(cpm(y)>1)>=4)
y = y[keep,]
logCPMs <- cpm(y, log = TRUE)
logCPMs_corrected <- limma::removeBatchEffect(logCPMs,
                                              batch = c("batch1","batch1","batch1","batch1","batch2","batch2","batch3","batch3"))

pca <- prcomp(t(logCPMs_corrected))
to_plot <- data.frame(pca$x, tissue = paste0(y$samples$group))
to_plot$rownames <- rownames(to_plot)
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
cor_matrix <- cor(logCPMs_corrected)
pheatmap::pheatmap(cor_matrix)

pca <- prcomp(t(logCPMs))
to_plot <- data.frame(pca$x, tissue = paste0(y$samples$group))
to_plot$rownames <- rownames(to_plot)
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
pheatmap::pheatmap(cor_matrix)



previous <- read.csv(paste0("data/samples/muscle/",antibody,"/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
tab <- tab[keep,c(1:6)]
logCPMs_corrected <- cbind(tab,logCPMs_corrected)
increase <- previous$Geneid[which(previous$Significant_bar=="Up")]
logCPMs_corrected_increase <- logCPMs_corrected[which(logCPMs_corrected$Geneid %in% increase),]
logCPMs_corrected_increase <- logCPMs_corrected_increase[,c(1:10,11,13,12,14)]
pheatmap::pheatmap(logCPMs_corrected_increase[,c(7:ncol(logCPMs_corrected_increase))],cluster_rows = T,cluster_cols = F,show_rownames = F,main = paste0(antibody," previous increase"))

decrease <- previous$Geneid[which(previous$Significant_bar=="Down")]
logCPMs_corrected_decrease <- logCPMs_corrected[which(logCPMs_corrected$Geneid %in% decrease),]
logCPMs_corrected_decrease <- logCPMs_corrected_decrease[,c(1:10,11,13,12,14)]
pheatmap::pheatmap(logCPMs_corrected_decrease[,c(7:ncol(logCPMs_corrected_decrease))],cluster_rows = T,cluster_cols = F,show_rownames = F,main = paste0(antibody," previous decrease"))
to_plot <- logCPMs_corrected_decrease[,c(1,7:ncol(logCPMs_corrected_decrease))]
to_plot <- melt(to_plot)
to_plot$variable <- factor(to_plot$variable, levels =unique(to_plot$variable))
to_plot$age <- "young"
to_plot$age[which(to_plot$variable %in% unique(to_plot$variable)[c(3,4,7,8)])] <- "old"
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

tab <- read.table(paste0("data/raw_data/muscle_cut_tag_test/",antibody,"/",bin_size,"_bins.counts"), header = T)
counts <- tab[,c(7:ncol(tab))]
counts <- counts[,c(1:4)]
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
colnames(counts) <- gsub(pattern, "\\1", colnames(counts))
age <- c("young","young","old","old")
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
