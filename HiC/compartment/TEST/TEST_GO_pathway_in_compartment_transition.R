rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
library(ggrepel)
library(grid)
library(data.table)
library(clusterProfiler)
library("AnnotationDbi")
library(org.Mm.eg.db)
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

TSS <- read.table("~/ref_data/for_normal_mapping/mm10/mm10_refseq_TSS.bed")
TSS <- TSS[,c("V1","V2","V3","V6")]
head(TSS)
TSS <- as.
data.table(TSS)
setDT(TSS)
setkey(TSS,V1,V2,V3)
tissues <- c("brain","CB", "kidney", "liver", "lung", "bonemarrow", "colon", "heart", "Hip", "mammarygland", "stomach", "thymus")
tissue <- "lung"
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
for(tissue in tissues){
  compartment <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/compartment_change_50000.csv"))
  compartment$start <- compartment$start + 1
  A_B <- compartment[which(compartment$condition=="A-B"),]
  B_A <- compartment[which(compartment$condition=="B-A"),]
  
  A_B <- as.data.table(A_B)
  setDT(A_B)
  setkey(A_B,chr,start,end)
  
  B_A <- as.data.table(B_A)
  setDT(B_A)
  setkey(B_A,chr,start,end)
  
  A_B_overlaps <- foverlaps(A_B, TSS, type = "any", nomatch = 0L)  
  B_A_overlaps <- foverlaps(B_A, TSS, type = "any", nomatch = 0L)  
  A_B_genes <- A_B_overlaps$V6
  B_A_genes <- B_A_overlaps$V6
  
  A_B_genes <- bitr(A_B_genes,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  B_A_genes <- bitr(B_A_genes,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  
  A_B_genes_GO <- enrichGO( A_B_genes$ENTREZID,
                              OrgDb = GO_database,
                              keyType = "ENTREZID",
                              ont = "BP",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.05,
                              readable = T)
  barplot(A_B_genes_GO,title = paste0(tissue_label_change(tissue)," Compartment A to B gene GO pathway"),label_format = 50)
  
  B_A_genes_GO <- enrichGO(B_A_genes$ENTREZID,
                            OrgDb = GO_database,
                            keyType = "ENTREZID",
                            ont = "BP",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05,
                            readable = T)
  barplot(B_A_genes_GO,title = paste0(tissue_label_change(tissue)," Compartment B to A gene GO pathway"),label_format = 50)
  }