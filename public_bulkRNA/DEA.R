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
library(biomaRt) 
tab = read.csv("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/GSE132040/GSE132040_190214_A00111_0269_AHH3J3DSXX_190214_A00111_0270_BHHMFWDSXX.csv")
tab <- tab[-c(54353:54357),]
table = read.delim("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/GSE132040/GSE132040_MACA_Bulk_metadata.csv",sep = ',')
table <- table[which(str_detect(table$source.name,"Muscle")),]
table <- table[which(table$characteristics..sex=="m"),]
table$characteristics..age <- as.numeric(table$characteristics..age)
table <- table[order(table$characteristics..age), ]
new_column_name <- sub("\\..*", "", colnames(tab[,-1])) 
colnames(tab)[2:ncol(tab)] <- new_column_name
counts <- tab[,which(colnames(tab) %in% table$Sample.name)]
table <- table[which(table$Sample.name %in% colnames(counts)),]

# counts <- cbind(tab[,which(colnames(tab) %in% table$Run[which(table$Age %in% c(3,21,24))])])
# table <- table[which(table$Age %in% c(3,21,24)),]
counts <- cbind(tab[,which(colnames(tab) %in% table$Sample.name[which(table$characteristics..age %in% c(3,24))])])
table <- table[which(table$characteristics..age %in% c(3,24)),]
rownames(counts) <- tab$gene
counts <- counts[,match(table$Sample.name,colnames(counts),)]
# colnames(counts) <- c("m3_1","m3_2","m3_3","m3_4","m21_1","m21_2","m21_3","m21_4","m24_1","m24_2","m24_3")
# group =c("m3","m3","m3","m3","m21","m21","m21","m21","m24","m24","m24")
colnames(counts) <- c("m3_1","m3_2","m3_3","m3_4","m24_1","m24_2","m24_3")
group =c("m3","m3","m3","m3","m24","m24","m24")


design = model.matrix(~0+group)
lvls = levels(factor(group))
colnames(design) = lvls
len = length(lvls)
# myContrasts = c("m24-m3","m21-m3")
myContrasts = c("m24-m3")

contrast.matrix = eval(as.call(c(as.symbol("makeContrasts"),as.list(myContrasts),levels=list(design))))
y= DGEList(counts=counts,group=group)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
y =  calcNormFactors(y)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)
lrt = glmLRT(fit_tag, contrast = contrast.matrix)
# qBH = p.adjust(lrt$table$PValue,method="BH")
tab<-tab[keep,]
out = cbind(tab[,1:2],cpm(y))
# out[,c(3:13)] = sweep(out[,c(3:13)],1,out$Length,'/')*1000
for (i in 1:length(myContrasts)){
  lrt = glmLRT(fit_tag, contrast = contrast.matrix[,i])
  out[,paste0("PValue.",myContrasts[i])] = lrt$table$PValue
  out[,paste0("FDR.",myContrasts[i])] = p.adjust(lrt$table$PValue,method="BH")
  out[,paste0("logFC.",myContrasts[i])] = lrt$table$logFC
}

out$Significant <- ifelse(out$`FDR.m24-m3` < 0.05 & abs(out$`logFC.m24-m3`) >= 0, 
                          ifelse(out$`logFC.m24-m3` > 0, "Up", "Down"), "Stable")
colour<- list(c("grey"),c("grey","red"),c("blue","grey","red"))
ggplot(
  # 数据、映射、颜色
  out, aes(x =`logFC.m24-m3` , y = -log10(`FDR.m24-m3`))) +
  geom_point(aes(color = Significant), size=2) +
  scale_color_manual(values = colour[[nrow(as.data.frame(table(out$Significant)))]]) +
  # scale_color_manual(values = c("blue","grey")) +
  # 注释
  geom_text_repel(
    data = subset(out,-log10(`FDR.m24-m3`) > 15 |`logFC.m24-m3`>5 |`logFC.m24-m3` < -4 ),
    aes(label = gene),
    size = 5,max.overlaps = 100,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")) +
  # 辅助线
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
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
write.csv(out,"data/public_data/GSE132040/Limb_muscle.csv",row.names = F)

df <- as.data.frame(t(out[,c(3:(3+length(group)-1))]))
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
