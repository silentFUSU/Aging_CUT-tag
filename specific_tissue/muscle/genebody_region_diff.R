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
tissue <-"muscle"
antibody <- "H3K36me3"
tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_gene_body.counts"),skip=1)
counts = tab[,c(7:10)]
rownames(counts)= tab$Geneid
colnames(counts) = c("young_1","old_1","young_2","old_2")
group =c("young","old","young","old")
y= DGEList(counts=counts,group=group)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
y$samples$batch <- rep(c(rep("batch1", 2), rep("batch2", 2)), 1)
y$samples$year <- group
y$samples$year <- factor(y$samples$year,c("young","old"))
y <- calcNormFactors(y)
logCPMs <- cpm(y, log = TRUE)
batch <- factor(y$samples$batch)
design <- model.matrix(~batch+year, y$samples)
# y <- estimateDisp(y, design)

y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)
lrt = glmLRT(fit_tag, coef = 3)
tab<-tab[keep,]
# fit <- glmQLFit(y, design)
# fit  <- glmQLFTest(fit, coef = 3)
out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
            "PValue.old-young"=lrt$table$PValue,"FDR.old-young"= p.adjust(lrt$table$PValue,method="BH"),
            "LogFC.old-young"=lrt$table$logFC)
out$Significant <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$`LogFC.old-young`) >= 0, 
                          ifelse(out$`LogFC.old-young` > 0, "Up", "Down"), "Stable")
colour<- list(c("grey"),c("grey","red"),c("blue","grey","red"))
ggplot(
  # 数据、映射、颜色
  out, aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`))) +
  geom_point(aes(color = Significant), size=2) +
  scale_color_manual(values = colour[[nrow(as.data.frame(table(out$Significant)))]]) +
  # scale_color_manual(values = c("blue","grey")) +
  # 注释
  # geom_text_repel(
  #   data = subset(out,`FDR.FA-noFA` < 0.05 & abs(out$logFC) >= 1),
  #   aes(label = Geneid),
  #   size = 5,max.overlaps = 100,
  #   box.padding = unit(0.35, "lines"),
  #   point.padding = unit(0.3, "lines")) +
  # 辅助线
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (FDR)") +
  # xlim(-3,3)+
  # 图例
  theme_bw()+
  theme(text = element_text(size = 20))

ggplot() +
  geom_point(data=out[which(out$Significant=="Stable"),], mapping=aes(logCPM,`LogFC.old-young`), size=2, color="grey") +
  geom_point(data=out[which(out$Significant=="Up"),], mapping=aes(logCPM, `LogFC.old-young`), size=2, color="red") +
  geom_point(data=out[which(out$Significant=="Down"),], mapping=aes(logCPM,`LogFC.old-young`), size=2, color="blue") +
  labs(x="log2(CPM)",
       y="log2(FC)") +
  theme(legend.position = "bottom",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme_bw()

ggplot(out[which(out$Significant !="Stable"),],aes(x=Significant,y=Length,fill = Significant))+
  geom_boxplot()+
  # geom_jitter(shape=16,size=1,position = position_jitter(0.2))+
  scale_fill_brewer(palette="Set3")+
  theme_bw()+theme(text = element_text(size = 18))+ylab("gene body length")+
  xlab("")+labs(fill = "", color = "")
t.test(out$Length[which(out$Significant=="Up")],out$Length[which(out$Significant=="Down")])

window_size = "1000"
gap_size = "3000"
e_value = "100"
tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,
                        "_merge-W",window_size,"-G",gap_size,"-E",e_value,".counts"),skip=1)
counts = tab[,c(7:10)]
rownames(counts)= tab$Geneid
colnames(counts) = c("young_1","old_1","young_2","old_2")
group =c("young","old","young","old")
y= DGEList(counts=counts,group=group)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
y$samples$batch <- rep(c(rep("batch1", 2), rep("batch2", 2)), 1)
y$samples$year <- group
y$samples$year <- factor(y$samples$year,c("young","old"))
y <- calcNormFactors(y)
logCPMs <- cpm(y, log = TRUE)
batch <- factor(y$samples$batch)
design <- model.matrix(~batch+year, y$samples)
# y <- estimateDisp(y, design)

y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)
lrt = glmLRT(fit_tag, coef = 3)
tab<-tab[keep,]
# fit <- glmQLFit(y, design)
# fit  <- glmQLFTest(fit, coef = 3)
out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
            "PValue.old-young"=lrt$table$PValue,"FDR.old-young"= p.adjust(lrt$table$PValue,method="BH"),
            "LogFC.old-young"=lrt$table$logFC)
out$Significant <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$`LogFC.old-young`) >= 0, 
                          ifelse(out$`LogFC.old-young` > 0, "Up", "Down"), "Stable")
colour<- list(c("grey"),c("grey","red"),c("blue","grey","red"))
out$genebody<-"no"
genebody_peak <- read.table("data/samples/muscle/H3K36me3/bed/H3K36me3_merge-W1000-G3000-E100_intersect_with_genebody_unique_sort.bed")
out$genebody[which(out$Geneid %in% genebody_peak$V4)] <- "yes"

ggplot() +
  geom_point(out, mapping = aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`)), size=2, color="grey") +
  geom_point(out[which(out$genebody=="yes" & out$Significant=="Up"),], mapping = aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`)), size=2, color="red") +
  geom_point(out[which(out$genebody=="yes" & out$Significant=="Down"),], mapping = aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`)), size=2, color="blue")
  