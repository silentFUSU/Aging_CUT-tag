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
library(biomaRt) 

tissues <- c("Hip","tongue","bladder","heart","stomach","CB","aorta")
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
chrM_issue_diff <- function(tissue){
  tab = read.delim(paste0("data/samples/RNA/",tissue,"/combined-chrM.counts"),skip=1)
  rownames(tab) <- tab$Geneid
  tab <- tab[,-1]
  colnames <- colnames(tab)[6:length(tab)]
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+).*"
  colnames(tab)[6:length(tab)] <- gsub(pattern, "\\1", colnames(tab)[6:length(tab)] )
  counts <- tab[6:length(tab)]
  search_table <- read.csv("data/samples/all/RNA_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue_label_change(tissue) & search_table$age=="3m"),]  
  counts <- counts[,search_table$sample_name]
  mouse <- as.character(search_table$mouse_ID)
  
  y= DGEList(counts=counts,group=mouse)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  y$samples$group <- factor(y$samples$group, levels=c(mouse[which(mouse != "100")],"100"))
  design <- model.matrix(~group, y$samples)
  y <- calcNormFactors(y)
  bcv <- 0.1
  fit_tag = glmFit(y,design,dispersion = bcv^2)
  lrt = glmLRT(fit_tag, coef = 2)
  t_tab<-tab[keep,]
  out = cbind(cpm(y),lrt$table, "fdr"=p.adjust(lrt$table$PValue,method="BH"))
  out$Significant <- ifelse(out$fdr< 0.05 & abs(out$logFC) >= 0, 
                            ifelse(out$logFC > 0, "Up", "Down"), "Stable")
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  
  p <- ggplot(
    out, aes(x = logFC, y = -log10(fdr))) +
    geom_point(aes(color = Significant), size=2) +
    scale_color_manual(values = colour) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (p-value)") +
    theme_bw()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ggtitle(paste0(tissue_label_change(tissue)," RNA 100 vs ",mouse[which(mouse != "100")]))+
    annotate("text", x = min(out$logFC), y = max(-log10(out$fdr)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$logFC), y = max(-log10(out$fdr)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  print(p)
  dir.create("data/chrM_inssue/RNA/")
  write.csv(out,paste0("data/chrM_inssue/RNA/",tissue,"_RNA_diff_bcv01.csv"))
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
  barplot(genelist_up_GO,title = paste0(tissue_label_change(tissue)," Increased gene GO pathway"),label_format = 50)
  GO_table <- genelist_up_GO@result
  write.csv(GO_table,paste0("data/chrM_inssue/RNA/",tissue,"_RNA_diff_bcv01_increase_GO.csv"))
  
  genelist_down <- bitr(rownames(out)[which(out$Significant=="Down")],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  genelist_down_GO <- enrichGO( genelist_down$ENTREZID,#GO富集分析
                                OrgDb = GO_database,
                                keyType = "ENTREZID",#设定读取的gene ID类型
                                ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                                pvalueCutoff = 0.05,#设定p值阈值
                                qvalueCutoff = 0.05,#设定q值阈值
                                readable = T)
  barplot(genelist_down_GO,title = paste0(tissue_label_change(tissue)," Decreased gene GO pathway"),label_format = 50)
  GO_table <- genelist_down_GO@result
  write.csv(GO_table,paste0("data/chrM_inssue/RNA/",tissue,"_RNA_diff_bcv01_decrease_GO.csv"))
}
for(tissue in tissues){
  chrM_issue_diff(tissue)
}

tissues <- c("Hip","tongue","bladder","heart","stomach","CB","aorta")
chrM_issue_diff_normal <- function(tissue){
  tab = read.delim(paste0("data/samples/RNA/",tissue,"/combined-chrM.counts"),skip=1)
  rownames(tab) <- tab$Geneid
  tab <- tab[,-1]
  colnames <- colnames(tab)[6:length(tab)]
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+).*"
  colnames(tab)[6:length(tab)] <- gsub(pattern, "\\1", colnames(tab)[6:length(tab)] )
  counts <- tab[6:length(tab)]
  search_table <- read.csv("data/samples/all/RNA_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue_label_change(tissue) & search_table$age=="24m"),]  
  counts <- counts[,search_table$sample_name]
  mouse <- search_table$mouse_ID
  
  y= DGEList(counts=counts,group=mouse)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  design <- model.matrix(~group, y$samples)
  y <- calcNormFactors(y)
  bcv <- 0.1
  fit_tag = glmFit(y,design,dispersion = bcv^2)
  lrt = glmLRT(fit_tag, coef = 2)
  t_tab<-tab[keep,]
  out = cbind(cpm(y),lrt$table, "fdr"=p.adjust(lrt$table$PValue,method="BH"))
  out$Significant <- ifelse(out$fdr< 0.05 & abs(out$logFC) >= 0, 
                            ifelse(out$logFC > 0, "Up", "Down"), "Stable")
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  
  p <- ggplot(
    out, aes(x = logFC, y = -log10(fdr))) +
    geom_point(aes(color = Significant), size=2) +
    scale_color_manual(values = colour) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (p-value)") +
    theme_bw()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ggtitle(paste0(tissue_label_change(tissue)," RNA ",levels(y$samples$group)[2]," vs ",levels(y$samples$group)[1]))+
    annotate("text", x = min(out$logFC), y = max(-log10(out$fdr)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$logFC), y = max(-log10(out$fdr)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  print(p)
  dir.create("data/chrM_inssue/RNA/")
  write.csv(out,paste0("data/chrM_inssue/RNA/",tissue,"_RNA_diff_bcv01_normal_",levels(y$samples$group)[2],"_",levels(y$samples$group)[1],".csv"))
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
  barplot(genelist_up_GO,title = paste0(tissue_label_change(tissue)," Increased gene GO pathway"),label_format = 50)
  GO_table <- genelist_up_GO@result
  write.csv(GO_table,paste0("data/chrM_inssue/RNA/",tissue,"_RNA_diff_bcv01_normal_",levels(y$samples$group)[2],"_",levels(y$samples$group)[1],"_increase_GO.csv"))
  
  genelist_down <- bitr(rownames(out)[which(out$Significant=="Down")],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  genelist_down_GO <- enrichGO( genelist_down$ENTREZID,#GO富集分析
                                OrgDb = GO_database,
                                keyType = "ENTREZID",#设定读取的gene ID类型
                                ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                                pvalueCutoff = 0.05,#设定p值阈值
                                qvalueCutoff = 0.05,#设定q值阈值
                                readable = T)
  barplot(genelist_down_GO,title = paste0(tissue_label_change(tissue)," Decreased gene GO pathway"),label_format = 50)
  GO_table <- genelist_down_GO@result
  write.csv(GO_table,paste0("data/chrM_inssue/RNA/",tissue,"_RNA_diff_bcv01_normal_",levels(y$samples$group)[2],"_",levels(y$samples$group)[1],"_decrease_GO.csv"))
}
for(tissue in tissues){
  chrM_issue_diff_normal(tissue)
}
