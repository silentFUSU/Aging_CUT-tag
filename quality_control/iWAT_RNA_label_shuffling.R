rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(stringr)
tissue <- "iWAT"
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
    }else if(tissue_label == "Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
}
tab = read.delim(paste0("data/samples/RNA/iWAT/combined-chrM.counts"),skip=1)
counts = tab[,c(7:ncol(tab))]
rownames(counts)= tab$Geneid
colnames(counts) = c("young_1","old_1","young_2","old_2")
group =c("label1","label2","label2","label1")
y= DGEList(counts=counts,group=group)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
# y$samples$batch <- rep(c(rep("batch1", 2), rep("batch2", 2)), 1)
y <- calcNormFactors(y)
design <- model.matrix(~group, y$samples)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)
lrt = glmLRT(fit_tag, coef = 2)
tab<-tab[keep,]

out = as.data.frame(cbind(cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
            "PValue"=lrt$table$PValue,"FDR"= p.adjust(lrt$table$PValue,method="BH"),
            "LogFC"=lrt$table$logFC))

out$Significant <- ifelse(out$FDR < 0.05 & abs(out$LogFC) >= 0, 
                          ifelse(out$LogFC > 0, "Up", "Down"), "Stable")

color <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
ggplot(
  out, aes(x = LogFC, y = -log10(FDR))) +
  geom_point(aes(color = Significant), size=2) +
  scale_color_manual(values = color) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (FDR)") +
  theme_bw()+
  ggtitle(tissue_label_change(tissue))+
  theme(text = element_text(size = 20))+
  annotate("text", x = min(out$LogFC), y = max(-log10(out$FDR)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
  annotate("text", x = max(out$LogFC), y = max(-log10(out$FDR)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
tab = read.delim(paste0("data/samples/RNA/iWAT/combined-chrM.counts"),skip=1)
counts = tab[,c(7:ncol(tab))]
rownames(counts)= tab$Geneid
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+).*"
colnames(counts) <- gsub(pattern, "\\1", colnames(counts))
y= DGEList(counts=counts)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
logCPMs <- cpm(y, log = TRUE)
pheatmap::pheatmap(logCPMs,show_rownames = F)
