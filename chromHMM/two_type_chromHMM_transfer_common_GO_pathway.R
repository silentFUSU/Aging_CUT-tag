rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(ggplot2)
library(ChIPseeker)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
type1 <- c("liver","bonemarrow","heart","skin","spleen","cecum","colon","lung","brain","Hip","aorta","muscle","stomach")
type2 <- c("ovary","mammarygland","tongue","uterus","thymus","jejunum","testis","iWAT","BAT","kidney","CB","bladder","pancreas")

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

tissues <- type1
state_num <- "14"
state <- "E14"
target_state <- "E11"

common_GO_pathway <-  function(tissues,state,target_state,state_num){
  GO_database <- 'org.Mm.eg.db'
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  gene_list <- list()
  for(tissue in tissues){
    young1 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_young1_",state_num,"_segments_1k.bed"),header = F)
    young2 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_young2_",state_num,"_segments_1k.bed"),header = F)
    young1$label <- paste(young1$V1,young1$V2,young1$V3,young1$V4,sep = "-")
    young2$label <- paste(young2$V1,young2$V2,young2$V3,young2$V4,sep = "-")
    young <- young1[which(young1$label %in% intersect(young1$label,young2$label)),]
    
    old1 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_old1_",state_num,"_segments_1k.bed"),header = F)
    old2 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_old2_",state_num,"_segments_1k.bed"),header = F)
    old1$label <- paste(old1$V1,old1$V2,old1$V3,old1$V4,sep = "-")
    old2$label <- paste(old2$V1,old2$V2,old2$V3,old2$V4,sep = "-")
    old <- old1[which(old1$label %in% intersect(old1$label,old2$label)),]
    
    young$label <- paste0(young$V1,":",young$V2,"-",young$V3)
    old$label <- paste0(old$V1,":",old$V2,"-",old$V3)
    young <- young[which(young$label %in% old$label),]  
    old <- old[which(old$label %in% young$label),]
    colnames(young)[4] <- "young_state"
    colnames(old)[4] <- "old_state"
    state_change <- merge(young[,c(1:5)],old[,c(4:5)],by="label")
    state_change <- state_change[which(state_change$V1 %in% c(paste0("chr",c(1:19,"X","Y")))),]
    state_change <- state_change[,c(-1)]
    state_change$condition <- paste0(state_change$young_state,"-",state_change$old_state)  
    state_change <- state_change[which(state_change$young_state==state & state_change$old_state==target_state),] 
    state_change_region <- GRanges(seqnames = state_change$V1,   
                              ranges = IRanges(start =state_change$V2, end = state_change$V2))
    state_change_region_anno  <- annotatePeak(state_change_region, tssRegion=c(-3000, 3000),
                                   TxDb=txdb, annoDb="org.Mm.eg.db")
    state_change_region_anno <- unique(as.data.frame(state_change_region_anno))
    gene_list[[tissue]] <- bitr(state_change_region_anno$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  }
  gene_list <- readRDS("data/t.rds")
  ck <- compareCluster(gene_list, fun = enrichGO,OrgDb = GO_database, keyType = "ENTREZID",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
}
