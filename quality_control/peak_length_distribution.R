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
library(corrplot)  
antibody = "H3K27me3"
tissue = "brain"
window_size = "1000"
gap_size = "3000"
e_value = "100"
# peak_short <- as.data.frame(readPeakFile(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_merge-W",window_size,"-G",gap_size,"-E",e_value,".bed")))
peak_short <- as.data.frame(readPeakFile(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_macs2_mergecounts.bed")))
# peak_short$group <- paste0("W",window_size,"-G",gap_size,"-E",e_value)
peak_short$group <- paste0("Macs2")
# window_size = "5000"
# gap_size = "10000"
# e_value = "100"
peak_long <- as.data.frame(readPeakFile(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_merge-W",window_size,"-G",gap_size,"-E",e_value,".bed")))
peak_long$group <- paste0("W",window_size,"-G",gap_size,"-E",e_value)

distribution <- rbind(peak_short[,c(4,7)],peak_long[,c(4,8)])
# distribution <- rbind(peak_short[,c(4,8)],peak_long[,c(4,8)])
group <- as.vector(as.data.frame(table(distribution$group))[,1])
ggplot(distribution, aes(x = width, fill = group)) +  
   geom_density(alpha = 0.5) +  
   labs(title = "Distribution of Two Groups", x = "Value", y = "Density") +  xlim(0,10000) +theme_bw()+theme(text = element_text(size = 18))
nrow(peak_short[which(peak_short$width>250000),])

nrow(peak_long[which(peak_long$width>250000),])

ggplot(distribution[which(distribution$width<250000),], aes(x = width, fill = group)) +  
  geom_density(alpha = 0.5) +  
  labs(title = "Distribution of Two Groups", x = "Value", y = "Density") +  
  scale_fill_manual(values = c("W1000-G3000-E100" = "blue", "W5000-G10000-E100" = "red"))+theme_bw()+theme(text = element_text(size = 18))
