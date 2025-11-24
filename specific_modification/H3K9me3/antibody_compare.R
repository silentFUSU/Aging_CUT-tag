rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(edgeR)
library(MASS) 
library(gridExtra)
library(dplyr)
library(stringr)
tab = read.delim(paste0("data/raw_data_transposon2/20250110_DYQ_CUTTag/H3K9me3_10kb_bins.counts"),skip=1)  
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+|DYQ[0-9]+).*"
colnames(tab)[7:length(tab)] <- gsub(pattern, "\\1", colnames(tab)[7:length(tab)])
counts <- tab[7:length(tab)]
age <- c("young","young","young","young","old","old","old","old","young","old","young","old")
y= DGEList(counts=counts,group = age)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
logCPMs <- cpm(y, log = TRUE)
pca <- prcomp(t(logCPMs))
to_plot <- data.frame(pca$x, age = paste0(y$samples$group))
to_plot$rownames <- rownames(to_plot)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
to_plot$age <- factor(to_plot$age, levels = c("young","old"))
ggplot(to_plot, aes(x=PC1, y=PC2, color=age)) + 
  geom_point(size=5) +theme_bw()+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
  geom_text_repel(  
    data = to_plot,  
    aes(x = PC1, y = PC2, label = rownames, color = age),  
    size = 5,  
    box.padding = unit(0.35, "lines"),  
    point.padding = unit(0.3, "lines")  
  )

batch=c("batch1","batch1","batch2","batch2","batch1","batch1","batch2","batch2","batch3","batch3","batch3","batch3")
logCPMs_corrected <- limma::removeBatchEffect(logCPMs, batch = batch)
pca <- prcomp(t(logCPMs_corrected))
to_plot <- data.frame(pca$x, age = paste0(y$samples$group))
to_plot$rownames <- rownames(to_plot)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
to_plot$age <- factor(to_plot$age, levels = c("young","old"))
ggplot(to_plot, aes(x=PC1, y=PC2, color=age)) + 
  geom_point(size=5) +theme_bw()+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
  geom_text_repel(  
    data = to_plot,  
    aes(x = PC1, y = PC2, label = rownames, color = age),  
    size = 5,  
    box.padding = unit(0.35, "lines"),  
    point.padding = unit(0.3, "lines")  
  )

## 8898
tab = read.delim(paste0("data/raw_data_transposon2/20250110_DYQ_CUTTag/H3K9me3_10kb_bins.counts"),skip=1)  
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+|DYQ[0-9]+).*"
colnames(tab)[7:length(tab)] <- gsub(pattern, "\\1", colnames(tab)[7:length(tab)])
counts <- tab[7:length(tab)]
# counts <- counts[,c(1,3,5,7)]
counts <- counts[,c(2,4,6,8)]
age <- c("young","young","old","old")
y= DGEList(counts=counts,group=age)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
y$samples$year <- age
y$samples$year <- factor(y$samples$year,c("young","old"))
y <- calcNormFactors(y)
design <- model.matrix(~year, y$samples)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)
lrt = glmLRT(fit_tag, coef = 2)
tab<-tab[keep,]

out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
            "PValue.old-young"=lrt$table$PValue,"FDR.old-young"= p.adjust(lrt$table$PValue,method="BH"),
            "LogFC.old-young"=lrt$table$logFC)
out$Significant <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$`LogFC.old-young`) >= log2(1.2), 
                          ifelse(out$`LogFC.old-young` > log2(1.2), "Up", "Down"), "Stable")
colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))

ggplot(
  out, aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`))) +
  geom_point(aes(color = Significant), size=2) +
  scale_color_manual(values = colour) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (p-value)") +
  theme_bw()+
  theme(text = element_text(size = 20),legend.position = "none")+
  ggtitle("Lung H3K9me3")+
  annotate("text", x = min(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
  annotate("text", x = max(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
antibody8898 <- out
antibody176916 <- out
to_plot <- merge(antibody8898[,c("Geneid","LogFC.old-young")],antibody176916[,c("Geneid","LogFC.old-young")],by="Geneid")
colnames(to_plot)[c(2,3)]<-c("a8898","a176916")
ggplot(to_plot, aes(x = a8898, y = a176916)) +
  # 坐标轴
  labs(x="log2(Fold change) #8898",
       y="log2(Fold change) #176916") +
  geom_point(color="grey")+
  geom_smooth(method=lm,color = "#c9d6df")+
  xlim(-3,3)+
  ylim(-3,3)+
  # 图例
  theme(legend.position = "bottom",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme_bw()+theme(text = element_text(size = 18))+
  geom_abline(intercept = 0, slope = 1,color="red") 

peaks <- read.table("data/samples/lung/H3K9me3/bed/H3K9me3_10kb_in_young_old_merge-W1000-G3000-E100.bed")

ggplot(to_plot[which(to_plot$Geneid %in% peaks$V4),], aes(x = a8898, y = a176916)) +
  # 坐标轴
  labs(x="log2(Fold change) #8898",
       y="log2(Fold change) #176916") +
  geom_point(color="grey")+
  geom_smooth(method=lm,color = "#c9d6df")+
  xlim(-3,3)+
  ylim(-3,3)+
  # 图例
  theme(legend.position = "bottom",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme_bw()+theme(text = element_text(size = 18))+
  geom_abline(intercept = 0, slope = 1,color="red") 



