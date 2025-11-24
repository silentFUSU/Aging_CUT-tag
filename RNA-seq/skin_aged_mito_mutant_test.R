rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(ggrepel)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(clusterProfiler)
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
tab <-  read.delim(paste0("data/samples/RNA/skin/combined-chrM.counts"),skip=1)
tab2 <- read.delim(paste0("data/raw_data_transposon2/20250601_DYQ_RNA/combined-chrM.counts"),skip = 1)
tab <- merge(tab,tab2[,c(1,7:ncol(tab2))],by="Geneid")
rownames(tab) <- tab$Geneid
tab <- tab[,-1]
colnames <- colnames(tab)[6:length(tab)]
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|DYQ[0-9]+|HM[0-9]+).*"
new_colnames <- gsub(pattern, "\\1", colnames)
colnames(tab)[6:length(tab)] <- new_colnames
sorted_index <- order(new_colnames)
order_colnames <- new_colnames[sorted_index] 
counts <- tab[,order_colnames]   

counts <- counts[,c("DYQ188","HM028","HM029","LLX636","LLX637")]
group <- c("normal","normal","normal","mutant","mutant")
y= DGEList(counts=counts,group=group)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
logCPMs <- cpm(y, log = TRUE)
pca <- prcomp(t(logCPMs))
to_plot <- data.frame(pca$x)
to_plot$sample_name <- rownames(to_plot)

to_plot$group <- group
to_plot$rownames <- paste0(to_plot$sample_name,"-",to_plot$group)
to_plot$group <- factor(to_plot$group,levels=c("normal","mutant"))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
ggplot(to_plot, aes(x=PC1, y=PC2, color=group)) + 
  geom_point(size=5) +theme_bw()+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
  geom_text_repel(  
    data = to_plot,  
    aes(x = PC1, y = PC2, label = rownames, color = group),  
    size = 5,  
    box.padding = unit(0.35, "lines"),  
    point.padding = unit(0.3, "lines")  
  ) +
  ggtitle("Skin")

y$samples$group <- factor(y$samples$group, levels=c("normal","mutant"))
design <- model.matrix(~group, y$samples)
y <- calcNormFactors(y)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)
lrt = glmLRT(fit_tag, coef = 2)
tab<-tab[keep,]
out = cbind(cpm(y),lrt$table, "fdr"=p.adjust(lrt$table$PValue,method="BH"))
out$Significant <- ifelse(out$fdr< 0.05 & abs(out$logFC) >= 0, 
                          ifelse(out$logFC > 0, "Up", "Down"), "Stable")
colour=setNames(c("blue","grey","red"),c("Down","Stable","Up"))
ggplot(
  # 数据、映射、颜色
  out, aes(x = logFC, y = -log10(fdr))) +
  geom_point(aes(color = Significant), size=2) +
  scale_color_manual(values = colour) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (fdr)") +
  theme_bw()+
  theme(text = element_text(size = 20))+
  ggtitle("Skin Mut vs WT")+
  annotate("text", x = min(out$logFC), y = max(-log10(out$fdr)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
  annotate("text", x = max(out$logFC), y = max(-log10(out$fdr)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
genelist_up <- bitr(rownames(out)[which(out$Significant=="Up")],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_up_GO <- enrichGO( genelist_up$ENTREZID,#GO富集分析
                            OrgDb = GO_database,
                            keyType = "ENTREZID",#设定读取的gene ID类型
                            ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                            pvalueCutoff = 0.05,#设定p值阈值
                            qvalueCutoff = 0.05,#设定q值阈值
                            readable = T)
barplot(genelist_up_GO,title = paste0("Mut vs WT Increased gene GO pathway"),label_format = 50)

genelist_down <- bitr(rownames(out)[which(out$Significant=="Down")],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_down_GO <- enrichGO( genelist_down$ENTREZID,#GO富集分析
                              OrgDb = GO_database,
                              keyType = "ENTREZID",#设定读取的gene ID类型
                              ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                              pvalueCutoff = 0.05,#设定p值阈值
                              qvalueCutoff = 0.05,#设定q值阈值
                              readable = T)
barplot(genelist_down_GO,title = paste0("Mut vs WT Decreased gene GO pathway"),label_format = 50)

out$Geneid <- rownames(out)
signif_correlation_genes <- out[,c("Geneid","logFC")]
signif_correlation_genes <- signif_correlation_genes[order(signif_correlation_genes$logFC),]
genelist <- signif_correlation_genes[,c("Geneid","logFC")]
genelist <- setNames(genelist$logFC, genelist$Geneid)
genelist <- sort(genelist, decreasing = TRUE)
gse <- gseGO(geneList=genelist, 
             ont = "BP",
             keyType = "SYMBOL", 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = GO_database,
             pAdjustMethod = "none",eps = 1e-100)
results_df <- as.data.frame(gse@result)
positive_nes_results <- subset(results_df, NES > 0)
negative_nes_results <- subset(results_df, NES < 0)
positive_gsea <- new("gseaResult",
                     result = positive_nes_results,
                     geneSets = gse@geneSets,
                     geneList = gse@geneList,
                     params = gse@params,
                     setType = gse@setType,
                     organism = gse@organism)
negative_gsea <- new("gseaResult",
                     result = negative_nes_results,
                     geneSets = gse@geneSets,
                     geneList = gse@geneList,
                     params = gse@params,
                     setType = gse@setType,
                     organism = gse@organism)
positive_gsea <- pairwise_termsim(positive_gsea)  
emapplot(positive_gsea, showCategory = 50) 
negative_gsea <- pairwise_termsim(negative_gsea)  
emapplot(negative_gsea, showCategory = 50) 
