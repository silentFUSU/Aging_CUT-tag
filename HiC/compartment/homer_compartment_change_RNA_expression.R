rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
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
tissue <- "lung"
resolution <- "50000"
RNA_in_compartment <- function(tissue,resolution){
  RNA <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene.csv"),row.names = 1)
  ref <- read.table("/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/TSS/refBed/mm10_refGene.bed")
  # ref <- read.table("/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/mm10/mm10_refseq_TSS.bed")
  ref$length <- ref$V3-ref$V2+1
  ref <- ref %>%
    group_by(V5) %>%
    filter(length == max(length)) %>%
    ungroup()
  colnames(ref)[5] <- "Geneid"
  RNA$Geneid <- rownames(RNA)
  RNA <- merge(RNA,ref[,c("Geneid","V1","V2","V3")])
  RNA <- RNA[,c("Geneid","V1","V2","V3","logFC","fdr")]
  RNA <- as.data.table(RNA)
  setDT(RNA)  
  setkey(RNA,V1,V2,V3)
  compartment <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/compartment_change_",resolution,".csv"))
  compartment <- as.data.table(compartment)
  setDT(compartment)  
  setkey(compartment, chr, start, end) 
  overlaps <- foverlaps(RNA, compartment, type = "any", nomatch = 0L)  
  overlaps <- as.data.frame(overlaps)
  overlaps$condition <- factor(overlaps$condition,levels = c("A-A","B-B","A-B","B-A"))
  p <- ggplot(overlaps, aes(x = condition, y = logFC, fill=condition)) +  
    geom_boxplot() +
    # scale_fill_manual(values = colours) +
    coord_flip() +
    theme_minimal()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ylab("log2(Fold change)") +
    xlab(NULL)+
    ggtitle(tissue_label_change(tissue))
  print(p)
}
for(tissue in c("lung","liver","CB","brain")){
    RNA_in_compartment(tissue,resolution)
}
