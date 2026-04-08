rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(edgeR)
library(ggplot2)
library(clusterProfiler)
library(stringr)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
conditions <- c("MEF_Bmi1", "MEF_Cbx2", "MEF_Cbx7")
p_list <- list()
for(condition in conditions){
  diff <- read.csv(paste0("data/samples/RNA/MEF_OE_RNA/MEF_Vector_",condition,"_diff_expression_gene.csv"))
  diff_Up <- diff[which(diff$Significant=="Up"),]
  genelist_up <- bitr(diff_Up$Geneid,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  genelist_up_GO <- enrichGO( genelist_up$ENTREZID,#GO富集分析
                              OrgDb = GO_database,
                              keyType = "ENTREZID",#设定读取的gene ID类型
                              ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                              pvalueCutoff = 0.05,#设定p值阈值
                              qvalueCutoff = 0.05,#设定q值阈值
                              readable = T)
  p_list[[paste0(condition," Up")]] <- barplot(genelist_up_GO,label_format = 50,showCategory = 10)+ggtitle(paste0(condition," Up"))
}
for(condition in conditions){
  diff <- read.csv(paste0("data/samples/RNA/MEF_OE_RNA/MEF_Vector_",condition,"_diff_expression_gene.csv"))
  diff_Down <- diff[which(diff$Significant=="Down"),]
  genelist_down <- bitr(diff_Down$Geneid,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  genelist_down_GO <- enrichGO( genelist_down$ENTREZID,#GO富集分析
                                OrgDb = GO_database,
                                keyType = "ENTREZID",#设定读取的gene ID类型
                                ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                                pvalueCutoff = 0.05,#设定p值阈值
                                qvalueCutoff = 0.05,#设定q值阈值
                                readable = T)
  p_list[[paste0(condition," Down")]] <- barplot(genelist_down_GO,label_format = 50,showCategory = 10)+ggtitle(paste0(condition," Down"))
}
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
combined_p <- plot_a_list(p_list,no_of_rows=2,no_of_cols=3)
ggsave("result/MEF_OE/DEG_GO.png",combined_p,width=30,height=15)
