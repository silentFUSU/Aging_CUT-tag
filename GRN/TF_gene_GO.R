rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
library(ggrepel)
library(tidyr)
library(dplyr)
library(tidyverse)
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
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
}
tissue <- "ovary"
TF <- "Nr5a2"
load("data/samples/GRN/grn_union_tissue.rdata")
load("data/samples/GRN/grn_union_skin.rdata")
grn_tissue[["skin"]] <- grn_union

df <- grn_tissue[[tissue]]
df <- df[which(df$TF == TF),]
df$gene_logFC <- log2(df$gene_old/df$gene_young)
df$condition <- "Up"
df$condition[which(df$gene_logFC < 0)] <- "Down"

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
genelist_up <- bitr(df$gene[which(df$condition=="Up")],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_up_GO <- enrichGO( genelist_up$ENTREZID,#GO富集分析
                            OrgDb = GO_database,
                            keyType = "ENTREZID",#设定读取的gene ID类型
                            ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                            pvalueCutoff = 0.05,#设定p值阈值
                            qvalueCutoff = 0.05,#设定q值阈值
                            readable = T)
barplot(genelist_up_GO,title = paste0(tissue_label_change(tissue)," ",TF," Increased gene GO pathway"),label_format = 50)


genelist_down <- bitr(df$gene[which(df$condition=="Down")],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_down_GO <- enrichGO( genelist_down$ENTREZID,#GO富集分析
                            OrgDb = GO_database,
                            keyType = "ENTREZID",#设定读取的gene ID类型
                            ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                            pvalueCutoff = 0.05,#设定p值阈值
                            qvalueCutoff = 0.05,#设定q值阈值
                            readable = T)
barplot(genelist_down_GO,title = paste0(tissue_label_change(tissue)," ",TF," decreased gene GO pathway"),label_format = 50)

RNA <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
df <- merge(df,RNA[,c("X","Significant")],by.x="gene",by.y="X")
genelist_up <- bitr(df$gene[which(df$Significant=="Up")],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_up_GO <- enrichGO( genelist_up$ENTREZID,#GO富集分析
                            OrgDb = GO_database,
                            keyType = "ENTREZID",#设定读取的gene ID类型
                            ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                            pvalueCutoff = 0.05,#设定p值阈值
                            qvalueCutoff = 0.05,#设定q值阈值
                            readable = T)
barplot(genelist_up_GO,title = paste0(tissue_label_change(tissue)," ",TF," Increased gene GO pathway"),label_format = 50)


genelist_down <- bitr(df$gene[which(df$Significant=="Down")],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_down_GO <- enrichGO( genelist_down$ENTREZID,#GO富集分析
                              OrgDb = GO_database,
                              keyType = "ENTREZID",#设定读取的gene ID类型
                              ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                              pvalueCutoff = 0.05,#设定p值阈值
                              qvalueCutoff = 0.05,#设定q值阈值
                              readable = T)
barplot(genelist_down_GO,title = paste0(tissue_label_change(tissue)," ",TF," decreased gene GO pathway"),label_format = 50)

to_plot <- genelist_down_GO@result
to_plot <- to_plot[order(to_plot$p.adjust,to_plot$pvalue),]
to_plot <- to_plot[c(1:5),]
to_plot$p.adjust <- -log10(to_plot$p.adjust)
to_plot$label <- paste0(to_plot$ID," ",to_plot$Description)
to_plot <- to_plot[,c("label","p.adjust")]
p <- ggplot(to_plot, aes(x = p.adjust, y = reorder(label, p.adjust))) +
  geom_bar(stat = "identity") +
  labs(x = "P Adjust", y = "Label") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10))
ggsave(paste0("result/Sup_figures/",TF,"_related_decreased_gene_GO.pdf"),p,width =15,height = 8)


