rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(edgeR)
df1 <- read.delim("data/samples/ATAC/brain/ATAC/revoked_data/ATAC_1kb_bins.counts",skip=1)
df2 <- read.delim("data/samples/ATAC/brain/ATAC/ATAC_1kb_bins.counts",skip = 1)
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY_[0-9]+).*"
colnames(df1)[7:ncol(df1)] <-  gsub(pattern, "\\1",colnames(df1)[7:ncol(df1)])
colnames(df2)[7:ncol(df1)] <-  gsub(pattern, "\\1",colnames(df2)[7:ncol(df2)])

df <- merge(df1[,c(1,7,8)],df2,by="Geneid")
rownames(df) <- df$Geneid
counts <- df[,c(2,3,9:12)]
y= DGEList(counts=counts)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
logCPMs <- cpm(y, log = TRUE)
pca <- prcomp(t(logCPMs))
to_plot <- data.frame(pca$x)
to_plot$age <- c("3m","3m","24m","24m","3m","3m")
to_plot$label <- paste0(rownames(to_plot),"-",to_plot$age)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))

ggplot(to_plot, aes(x=PC1, y=PC2, color=age)) + 
  geom_point(size=5) +theme_bw()+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
  geom_text_repel(  
    data = to_plot,  
    aes(x = PC1, y = PC2, label = label, color = age),  
    size = 5,  
    box.padding = unit(0.35, "lines"),  
    point.padding = unit(0.3, "lines")  
  )

logCPMs_corrected <- limma::removeBatchEffect(logCPMs,batch = c("batch1","batch1","batch1","batch1","batch2","batch2"))
pca <- prcomp(t(logCPMs_corrected))
to_plot <- data.frame(pca$x)
to_plot$age <- c("3m","3m","24m","24m","3m","3m")
to_plot$label <- paste0(rownames(to_plot),"-",to_plot$age)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))

ggplot(to_plot, aes(x=PC1, y=PC2, color=age)) + 
  geom_point(size=5) +theme_bw()+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
  geom_text_repel(  
    data = to_plot,  
    aes(x = PC1, y = PC2, label = label, color = age),  
    size = 5,  
    box.padding = unit(0.35, "lines"),  
    point.padding = unit(0.3, "lines")  
  )
age <- c("young","young","old","old","young","young")
y= DGEList(counts=counts,group=age)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
y$samples$group <- factor(y$samples$group,levels=c("young","old"))
y <- calcNormFactors(y)
design <- model.matrix(~group, y$samples)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)
lrt = glmLRT(fit_tag, coef = 2)
df<-df[keep,]
out = cbind(df[,c(1,4:6)],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
            "PValue.old-young"=lrt$table$PValue,"FDR.old-young"= p.adjust(lrt$table$PValue,method="BH"),
            "LogFC.old-young"=lrt$table$logFC)
out$Significant <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$`LogFC.old-young`) >= 0, 
                          ifelse(out$`LogFC.old-young` > 0, "Up", "Down"), "Stable")
colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
p <- ggplot(
  out, aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`))) +
  geom_point(aes(color = Significant), size=2) +
  scale_color_manual(values = colour) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (p-value)") +
  theme_bw()+
  theme(text = element_text(size = 20),legend.position = "none")+
  ggtitle(paste0("Cortex ATAC"))+
  annotate("text", x = min(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
  annotate("text", x = max(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)


# out <- read.csv("data/samples/ATAC/brain/ATAC/ATAC_1kb_bins_diff.csv")
out <- read.csv("data/samples/ATAC/brain/ATAC/revoked_data/ATAC_1kb_bins_diff_after_remove_batch_effect.csv")
colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
ggplot(
  out, aes(x = `LogFC.old.young`, y = -log10(`FDR.old.young`))) +
  geom_point(aes(color = Significant), size=2) +
  scale_color_manual(values = colour) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (p-value)") +
  theme_bw()+
  theme(text = element_text(size = 20),legend.position = "none")+
  ggtitle(paste0("Cortex ATAC"))+
  annotate("text", x = min(out$`LogFC.old.young`), y = max(-log10(out$`FDR.old.young`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
  annotate("text", x = max(out$`LogFC.old.young`), y = max(-log10(out$`FDR.old.young`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
peak <- read.table("data/samples/ATAC/brain/ATAC/bed/ATAC_1kb_in_young_old_merge_macs_narrowpeak.bed")
out <- out[which(out$Geneid %in% peak$V4),]
ggplot(
  out, aes(x = `LogFC.old.young`, y = -log10(`FDR.old.young`))) +
  geom_point(aes(color = Significant), size=2) +
  scale_color_manual(values = colour) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (p-value)") +
  theme_bw()+
  theme(text = element_text(size = 20),legend.position = "none")+
  ggtitle(paste0("Cortex ATAC"))+
  annotate("text", x = min(out$`LogFC.old.young`), y = max(-log10(out$`FDR.old.young`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
  annotate("text", x = max(out$`LogFC.old.young`), y = max(-log10(out$`FDR.old.young`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)


ATAC <- read.csv("data/samples/ATAC/Hip/ATAC/ATAC_1kb_bins_diff.csv")
peak <- read.table("data/samples/ATAC/Hip/ATAC/bed/ATAC_1kb_in_young_old_merge_macs_narrowpeak.bed")
ATAC <- ATAC[which(ATAC$Geneid %in% peak$V4),]
colnames(ATAC)[15] <- "ATAC"
antibodys <- c("H3K27ac","H3K4me3","H3K4me1")
for(antibody in antibodys){
  df <- read.csv(paste0("data/samples/Hip/",antibody,"/",antibody,"_1kb_bins_diff.csv"))
  peak <- read.table(paste0("data/samples/Hip/",antibody,"/bed/",antibody,"_1kb_in_young_old_merge_macs_narrowpeak.bed"))
  df <- df[which(df$Geneid %in% peak$V4),]
  colnames(df)[15] <- "other"
  to_plot <- merge(df[,c(1,15)],ATAC[,c(1,15)],by="Geneid")
  p <- ggplot()+
    geom_point(data=to_plot, mapping=aes(other,ATAC),color = "grey",alpha=0.5) +  
    geom_point(data=to_plot[which(to_plot[,2]>0 & to_plot[,3]>0),], mapping=aes(other,ATAC),color = "#00b8a9") +
    geom_point(data=to_plot[which(to_plot[,2]<0 & to_plot[,3]<0),], mapping=aes(other,ATAC),color = "#ff9a00") +
    labs(x=antibody,
         y="ATAC") +
    theme_bw()+theme(text = element_text(size = 18))+
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ggtitle(paste0("Hippocampus ATAC and ",antibody))+
    annotate("text",label = paste0(nrow(to_plot[which(to_plot[,2]>0 & to_plot[,3]>0),])),x = Inf, y = Inf,hjust = 1.1, vjust = 1.2,colour="#00b8a9",size=5)+
    annotate("text",label = paste0(nrow(to_plot[which(to_plot[,2]<0 & to_plot[,3]>0),])),x = -Inf, y = Inf,hjust = -0.1, vjust = 1.2,colour="#ff9a00",size=5)+
    annotate("text",label = paste0(nrow(to_plot[which(to_plot[,2]<0 & to_plot[,3]<0),])),x = -Inf, y = -Inf,hjust = -0.1, vjust = -1.2,colour="#f6416c",size=5)+
    annotate("text",label = paste0(nrow(to_plot[which(to_plot[,2]>0 & to_plot[,3]<0),])),x = Inf, y = -Inf,hjust = 1.1, vjust = -1.2,colour="#48466d",size=5)
  print(p)
  }
