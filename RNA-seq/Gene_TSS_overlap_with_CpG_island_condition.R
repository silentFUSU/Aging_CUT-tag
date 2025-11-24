rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(stringr)
library(ggrepel)
library(grid)
library(clusterProfiler)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(tidyverse)
library(data.table)
library(GO.db)
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
cpg <- read.table("~/ref_data/for_normal_mapping/mm10/cpgi.mm10.bed.txt")
cpg <- cpg[which(cpg$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
cpg <- as.data.table(cpg[,c(1:3)])
setDT(cpg)
setkey(cpg,V1,V2,V3)

genes <- read.table("~/ref_data/for_normal_mapping/mm10/mm10_refseq_TSS.bed")
genes <- as.data.table(genes[,c(1:3,6)])
setDT(genes)
setkey(genes,V1,V2,V3)
overlaps <- foverlaps(genes, cpg, type = "any", nomatch = 0L)  
summary <- data.frame(tissue="all",
                      all_genes=nrow(genes),
                      overlap_genes=length(unique(overlaps$V6)),
                      overlap=length(unique(overlaps$V6))/nrow(genes) *100)
summary$out <- 100 - summary$overlap
to_plot <- reshape2::melt(summary[,c(1,4,5)])
to_plot$variable <- factor(to_plot$variable,levels=c("out","overlap"))
to_plot$position <- 100
to_plot$position[which(to_plot$variable=="overlap")] <- to_plot[which(to_plot$variable=="overlap"),"value"]

color <- setNames(c("red","gray"),c("overlap","out"))
ggplot(to_plot, aes(x = tissue, y = value, fill = variable)) +  
  geom_bar(stat = 'identity') +   
  theme_minimal() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank())+
  geom_text(data = to_plot,
            aes(label = sprintf("%.1f", value), y = position),
            color = "black", size = 5, vjust = 0.5)

### CPM filtered genes
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
tissue_summary <- data.frame()
for(tissue in tissues){
  genes <- read.table("~/ref_data/for_normal_mapping/mm10/mm10_refseq_TSS.bed")
  df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
  genes <- genes[which(genes$V6 %in% df$X),]
  genes <- as.data.table(genes[,c(1:3,6)])
  setDT(genes)
  setkey(genes,V1,V2,V3)
  
  cpg <- read.table("~/ref_data/for_normal_mapping/mm10/cpgi.mm10.bed.txt")
  cpg <- cpg[which(cpg$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  cpg <- as.data.table(cpg[,c(1:3)])
  setDT(cpg)
  setkey(cpg,V1,V2,V3)
  overlaps <- foverlaps(genes, cpg, type = "any", nomatch = 0L)  
  t_tissue_summary <- data.frame(tissue=tissue_label_change(tissue),
                                 all_genes=nrow(genes),
                                 overlap_genes=length(unique(overlaps$V6)),
                                 overlap=length(unique(overlaps$V6))/nrow(genes) *100)
  tissue_summary <- rbind(tissue_summary,t_tissue_summary)
}
tissue_summary$out <- 100 - tissue_summary$overlap
to_plot <- reshape2::melt(tissue_summary[,c(1,4,5)])
to_plot$variable <- factor(to_plot$variable,levels=c("out","overlap"))
to_plot$position <- 100
to_plot$position[which(to_plot$variable=="overlap")] <- to_plot[which(to_plot$variable=="overlap"),"value"]

color <- setNames(c("red","gray"),c("overlap","out"))
ggplot(to_plot, aes(x = tissue, y = value, fill = variable)) +  
  geom_bar(stat = 'identity') +   
  theme_minimal() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank())+
  geom_text(data = to_plot,
            aes(label = sprintf("%.1f", value), y = position),
            color = "black", size = 5, vjust = 0.5)

### CPM filtered genes diff condition
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
tissue_summary <- data.frame()
condition <- "Down"
for(tissue in tissues){
  genes <- read.table("~/ref_data/for_normal_mapping/mm10/mm10_refseq_TSS.bed")
  df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
  df <- df[which(df$Significant==condition),]
  genes <- genes[which(genes$V6 %in% df$X),]
  genes <- as.data.table(genes[,c(1:3,6)])
  setDT(genes)
  setkey(genes,V1,V2,V3)
  
  cpg <- read.table("~/ref_data/for_normal_mapping/mm10/cpgi.mm10.bed.txt")
  cpg <- cpg[which(cpg$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  cpg <- as.data.table(cpg[,c(1:3)])
  setDT(cpg)
  setkey(cpg,V1,V2,V3)
  overlaps <- foverlaps(genes, cpg, type = "any", nomatch = 0L)  
  t_tissue_summary <- data.frame(tissue=tissue_label_change(tissue),
                                 all_genes=nrow(genes),
                                 overlap_genes=length(unique(overlaps$V6)),
                                 overlap=length(unique(overlaps$V6))/nrow(genes) *100)
  tissue_summary <- rbind(tissue_summary,t_tissue_summary)
}
tissue_summary$out <- 100 - tissue_summary$overlap
to_plot <- reshape2::melt(tissue_summary[,c(1,4,5)])
to_plot$variable <- factor(to_plot$variable,levels=c("out","overlap"))
to_plot$position <- 100
to_plot$position[which(to_plot$variable=="overlap")] <- to_plot[which(to_plot$variable=="overlap"),"value"]

color <- setNames(c("red","gray"),c("overlap","out"))
ggplot(to_plot, aes(x = tissue, y = value, fill = variable)) +  
  geom_bar(stat = 'identity') +   
  theme_minimal() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank())+
  geom_text(data = to_plot,
            aes(label = sprintf("%.1f", value), y = position),
            color = "black", size = 5, vjust = 0.5)
