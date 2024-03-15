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
library(limma)
library(gg.gap)
library(scales)
depth <- read.csv("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/Haplox_sequencing_depth.csv",sep = "\t")
ATAC <- depth[str_detect(depth$Info,"ATAC"),]
CUT <- depth[!str_detect(depth$Info,"ATAC"),]
CUT$antibody<-regmatches(depth$Info, regexpr("H3K\\w+", depth$Info))
ATAC $ antibody <-"ATAC"
depth <- rbind(CUT,ATAC)
depth$antibody[which(depth$antibody=="H3K4me")]<-"H3K4me1"
depth$antibody[which(depth$antibody=="H3K27AC")]<-"H3K27ac"
depth <- depth[!str_detect(depth$Info,"Test"),]
ggplot(depth,aes(x=antibody,y=Depth,fill=antibody))+
  geom_violin()+
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
  geom_jitter(width = 0.1, alpha = 0.7)+
  # geom_jitter(shape=16,size=1,position = position_jitter(0.2))+
  scale_fill_brewer(palette="Set3")+
  theme_bw()+theme(text = element_text(size = 18))+ylab("Depth (G)")+
  xlab("")+labs(fill = "", color = "")+ggtitle(paste0("Haplox sequencing depth"))+ylim(0,35)

depth <- read.csv("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/Jiangbei_sequencing_depth.csv",sep = "\t")
depth <- depth[-which(depth$Depth==2),]
ATAC <- depth[str_detect(depth$Info,"ATAC"),]
CUT <- depth[!str_detect(depth$Info,"ATAC"),]
CUT$antibody<-regmatches(depth$Info, regexpr("H3K\\w+", depth$Info))
ATAC $ antibody <-"ATAC"
depth <- rbind(CUT,ATAC)
depth$antibody[which(depth$antibody=="H3K4me")]<-"H3K4me1"

ggplot(depth,aes(x=antibody,y=Real_depth,fill=antibody))+
  geom_violin()+
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
  geom_jitter(width = 0.1, alpha = 0.7)+
  # geom_jitter(shape=16,size=1,position = position_jitter(0.2))+
  scale_fill_brewer(palette="Set3")+
  theme_bw()+theme(text = element_text(size = 18))+ylab("Depth (G)")+
  xlab("")+labs(fill = "", color = "")+ggtitle(paste0("Jiangbei sequencing depth"))+ylim(0,35)
