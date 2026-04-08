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
antibody = "H3K36me3"
tissue = "liver"
# bamcoverage <- list()
# sample_names <- list.files(paste0("data/samples/",tissue,"/",antibody,"/bam/"))
# sample_names_without_extension <- sub("\\..*", "", sample_names)  
# sample_names_without_extension <- unique(sample_names_without_extension)
# mapped_reads_samples <- c(31810943,27562610,29639278,25575270)
# for(i in c(1:length(sample_names_without_extension))){
#   bamcoverage[[i]] <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",sample_names_without_extension[i],"_macs_broad_peak_5k_merge_coverage_output.txt"), sep = "\t")
#   bamcoverage[[i]]$RPKM <- (10^9*bamcoverage[[i]]$V5)/(bamcoverage[[i]]$V7*mapped_reads_samples[i])
#   if(i==1){
#     QC<-as.data.frame(bamcoverage[[i]]$RPKM)
#     colnames(QC)[i] <- sample_names_without_extension[i]
#     rownames(QC) <- bamcoverage[[i]]$V4
#   }
#   else{
#     QC <- cbind(QC,as.data.frame(bamcoverage[[i]]$RPKM))  
#     colnames(QC)[i] <- sample_names_without_extension[i]
#   }
# }
# 
# QC$average <- rowMeans(QC[,1:4])
tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_macs_broad_peak_5k_merge.counts"),skip=1)
tab <- tab[which(tab$Length >500),]
# tab <- tab[which(tab$Geneid %in% rownames(QC)[which(QC$average > 1)]),]

counts = tab[,c(7:10)]
rownames(counts)= tab$Geneid
colnames(counts) = c("young_1","old_1","young_2","old_2")
counts <- counts[,c("young_1","young_2","old_1","old_2")]
group =c("young","young","old","old")
design = model.matrix(~0+group)
lvls = levels(factor(group))
colnames(design) = lvls
len = length(lvls)
myContrasts = c("old-young")
contrast.matrix = eval(as.call(c(as.symbol("makeContrasts"),as.list(myContrasts),levels=list(design))))
y= DGEList(counts=counts,group=group)
keep = which(rowSums(cpm(y)>1)>=4)
y = y[keep,]
y =  calcNormFactors(y)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)
lrt = glmLRT(fit_tag, contrast = contrast.matrix)
qBH = p.adjust(lrt$table$PValue,method="BH")
tab<-tab[keep,]
out = cbind(tab[,1:6],cpm(y),lrt$table,qBH)
out[,c(7:10)] = sweep(out[,c(7:10)],1,out$Length,'/')*1000
for (i in 1:length(myContrasts)){
  lrt = glmLRT(fit_tag, contrast = contrast.matrix[,i])
  out[,paste0("PValue.",myContrasts[i])] = lrt$table$PValue
  out[,paste0("FDR.",myContrasts[i])] = p.adjust(lrt$table$PValue,method="BH")
}
write.csv(out,paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_macs_broad_peak_5k_diff.csv"),row.names = F)

out$Significant <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$logFC) >= 0, 
                          ifelse(out$logFC > 0, "Up", "Down"), "Stable")
colour<- list(c("grey"),c("grey","red"),c("blue","grey","red"))
ggplot(
  # 数据、映射、颜色
  out, aes(x = logFC, y = -log10(`FDR.old-young`))) +
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
  geom_vline(xintercept=c(-0.7,0.7),lty=4,col="black",lwd=0.8) +
  # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (p-value)") +
  # xlim(-3,3)+
  # 图例
  theme(legend.position = "bottom",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme_bw()

df <- as.data.frame(t(out[,c(7:10)]))
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
ggsave(paste0("result/",tissue,"/pca/",antibody,"_pca_plot.png"),width = 5,height =5)
