rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(ggrepel)
library(ggsci)
library(gridExtra)
library(data.table)
library(karyoploteR)
library(GenomicRanges) 
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
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
                  "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT","ileum")) 
increased_summary <- data.frame()
decreased_summary <- data.frame()
for(tissue in tissues){
  DML <- read.table(paste0("data/samples/WGBS/",tissue,"/DSS_table/",tissue,"_DML_delta01.txt"),header = T)
  DML <- DML[which(DML$fdr < 0.05),c("chr","pos","fdr","diff")]
  DML$label <- paste0(DML$chr,"-",DML$pos)
  DML$tissue <- tissue_label_change(tissue)
  increased_summary <- rbind(increased_summary,DML[which(DML$diff >0),c("label","tissue")])
  decreased_summary <- rbind(decreased_summary,DML[which(DML$diff <0),c("label","tissue")])
}

increased_summary_count <- as.data.frame(table(increased_summary$label))
increased_summary_count <- increased_summary_count[order(increased_summary_count$Freq,decreasing = T),]
ggplot(increased_summary_count, aes(x = Freq)) +
  geom_histogram(binwidth = 1, boundary = 0, closed = "left",
                 fill = "steelblue", color = "white") +
  scale_y_log10() +
  labs(x = "Freq", y = "Count (log10)") +
  theme_bw()

decreased_summary_count <- as.data.frame(table(decreased_summary$label))
decreased_summary_count <- decreased_summary_count[order(decreased_summary_count$Freq,decreasing = T),]
ggplot(decreased_summary_count, aes(x = Freq)) +
  geom_histogram(binwidth = 1, boundary = 0, closed = "left",
                 fill = "steelblue", color = "white") +
  scale_y_log10() +
  labs(x = "Freq", y = "Count (log10)") +
  theme_bw()

increased_common_DML <- increased_summary_count[which(increased_summary_count$Freq >= 10),]
increased_common_DML <- increased_common_DML %>%
  separate(Var1, into = c("chr", "start"), sep = "-", convert = TRUE)
increased_common_DML$end <- increased_common_DML$start+1
increased_gr <- GRanges(
  seqnames = increased_common_DML$chr,
  ranges   = IRanges(start = increased_common_DML$start, end = increased_common_DML$end)
)
mouse.chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", 
                       "chr6", "chr7", "chr8", "chr9", "chr10",
                       "chr11", "chr12", "chr13", "chr14", "chr15",
                       "chr16", "chr17", "chr18", "chr19", "chrX", "chrY")
kp <- plotKaryotype(genome = "mm10", main="common hyper-DML", chromosomes = mouse.chromosomes,plot.type=6,cex=1.8)
kpDataBackground(kp, color = "#FFFFFFAA")
kpPlotDensity(kp, increased_gr,window.size = 0.5e6, data.panel="ideogram", col="#3c5488", border="#3c5488")


decreased_common_DML <- decreased_summary_count[which(decreased_summary_count$Freq >= 10),]
decreased_common_DML <- decreased_common_DML %>%
  separate(Var1, into = c("chr", "start"), sep = "-", convert = TRUE)
decreased_common_DML$end <- decreased_common_DML$start+1
decreased_gr <- GRanges(
  seqnames = decreased_common_DML$chr,
  ranges   = IRanges(start = decreased_common_DML$start, end = decreased_common_DML$end)
)
mouse.chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", 
                       "chr6", "chr7", "chr8", "chr9", "chr10",
                       "chr11", "chr12", "chr13", "chr14", "chr15",
                       "chr16", "chr17", "chr18", "chr19", "chrX", "chrY")
kp <- plotKaryotype(genome = "mm10", main="common hyper-DML", chromosomes = mouse.chromosomes,plot.type=6,cex=1.8)
kpDataBackground(kp, color = "#FFFFFFAA")
kpPlotDensity(kp, decreased_gr,window.size = 0.5e6, data.panel="ideogram", col="#3c5488", border="#3c5488")


increase_count <- increased_summary %>%   
  count(label)
increase_tissue <- increased_summary %>%   
  group_by(label) %>%   
  summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
