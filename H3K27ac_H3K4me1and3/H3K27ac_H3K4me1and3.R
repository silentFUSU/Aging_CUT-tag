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
tissue = "brain"
tab = read.delim(paste0("data/samples/brain/H3K27ac_H3K4me1and3/H3K27ac_H3K4me1and3.counts"),skip=1)
counts = tab[,c(7:18)]
rownames(counts)= tab$Geneid
colnames(counts) = c("H3K27ac_young_1","H3K27ac_old_1","H3K27ac_young_2","H3K27ac_old_2",
                     "H3K4me3_young_1","H3K4me3_old_1","H3K4me3_young_2","H3K4me3_old_2",
                     "H3K4me1_young_1","H3K4me1_old_1","H3K4me1_young_2","H3K4me1_old_2")
out <- cbind(tab[,1:6],cpm(counts))

df <- as.data.frame(t(out[,c(7:18)]))
df_pca <- prcomp(df) 
df_pcs <-data.frame(df_pca$x,Species=rownames(df)) 
percentage<-round(df_pca$sdev / sum(df_pca$sdev) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))
ggplot(df_pcs,aes(x=PC1,y=PC2,color=Species))+
  geom_point()+ 
  xlab(percentage[1]) +
  ylab(percentage[2])+    
  geom_text_repel(
    aes(label = rownames(df_pcs)),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))+
  theme(legend.position = "bottom",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 20),
        axis.text = element_text(size = 20), 
        axis.title = element_text(size = 20))+
  guides(color = F)+
  theme_bw()

