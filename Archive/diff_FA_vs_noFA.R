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
tissues <- c("pancreas")
antibodys <- c("H3K27ac","H3K4me1")
for (i in c(1:length(tissues))){
  tissue <- tissues[i]
  for (j in c(1:length(antibodys))){
    antibody <- antibodys[j]
    tab = read.delim(paste0("data/FA_vs_noFA/",tissue,"/",antibody,"/",antibody,"_macs_narrowpeak.counts"),skip=1)
    counts = tab[,c(7:8)]
    rownames(counts)= tab$Geneid
    # counts <- cbind(counts,counts)
    colnames(counts) = c("FA_1","noFA_1")
    # colnames(counts) = c("noFA_1","FA_1")
    colnames(tab)[7:8] <- c("FA_1","noFA_1")
    # colnames(tab)[7:8] <- c("noFA_1","FA_1")
    # counts <- counts[,c("noFA_1","noFA_2","FA_1","FA_2")]
    group =c("noFA","FA")
    design = model.matrix(~0+group)
    lvls = levels(factor(group))
    colnames(design) = lvls
    len = length(lvls)
    myContrasts = c("FA-noFA")
    contrast.matrix = eval(as.call(c(as.symbol("makeContrasts"),as.list(myContrasts),levels=list(design))))
    y= DGEList(counts=counts,group=group)
    keep = which(rowSums(cpm(y)>1)>=2)
    y = y[keep,]
    y =  calcNormFactors(y, method ="TMM")
    bcv <- 0.1
    # bcv <- 0.2
    y <- exactTest(y,dispersion = bcv^2, pair = c("noFA","FA"))
    # y<-estimateCommonDisp(y)
    # y<-estimateGLMTagwiseDisp(y,design)
    # fit_tag = glmFit(y,design)
    # lrt = glmLRT(fit_tag, contrast = contrast.matrix)
    # qBH = p.adjust(lrt$table$PValue,method="BH")
    tab<-tab[keep,]
    out = cbind(tab[,1:6],cpm(tab[,7:8]),logFC=y$table$logFC)
    out[,c(7:8)] = sweep(out[,c(7:8)],1,out$Length,'/')*1000
    for (i in 1:length(myContrasts)){
      out[,paste0("PValue.",myContrasts[i])] = y$table$PValue
      out[,paste0("FDR.",myContrasts[i])] = p.adjust(y$table$PValue,method="BH")
    }
    out$Significant <- ifelse(out$`FDR.FA-noFA` < 0.05 & abs(out$logFC) >= 1, 
                              ifelse(out$logFC > 1, "Up", "Down"), "Stable")
    
    
    colour<- list(c("grey"),c("grey","red"),c("blue","grey","red"))
    ggplot(
      # 数据、映射、颜色
      out, aes(x = logFC, y = -log10(`FDR.FA-noFA`))) +
      geom_point(aes(color = Significant), size=2) +
      scale_color_manual(values = colour[[nrow(as.data.frame(table(out$Significant)))]]) +
      # scale_color_manual(values = c("blue","grey")) +
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
      theme_bw()+theme(text = element_text(size = 18))
    
    ggsave(paste0("result/FA_vs_noFA/",tissue,"/",tissue,"_",antibody,"_volcano_plot_bcv01.png"),width = 10,height = 10)
    write.csv(out, paste0("result/FA_vs_noFA/",tissue,"/",tissue,"_",antibody,"_diff_bcv01.csv"))
    cordata <- out[,c(7,8)]
    cor_matrix <- cor(cordata)
    png(height=1800, width=1800, file=paste0("result/FA_vs_noFA/",tissue,"/",tissue,"_",antibody,"_corrplot_bcv01.png"), type = "cairo")
    corrplot(cor_matrix, method = "number",number.cex = 10, tl.cex=5)
    dev.off()
    
    scatter_data <- tab[,c(7:8)]
    ggplot(out, aes(x = log10(FA_1), y = log10(noFA_1))) +
      # 坐标轴
      labs(x="log10(reads(FA))",
           y="log10(reads(noFA))") +
      geom_point(color="grey")+
      geom_smooth(method=lm,color = "#c9d6df")+
      # xlim(-3,3)+
      # 图例
      theme(legend.position = "bottom",panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
      theme_bw()+theme(text = element_text(size = 18))+
      geom_abline(intercept = 0, slope = 1,color="red") 
    ggsave(paste0("result/FA_vs_noFA/",tissue,"/",tissue,"_",antibody,"_scatter_plot_bcv01.png"),width = 10,height = 10)
    
    scatter_data <- out[,c(7,8)]
    ggplot(out, aes(x = log10(FA_1), y = log10(noFA_1))) +
      # 坐标轴
      labs(x="log10(cpm(FA))",
           y="log10(cpm(noFA))") +
      geom_point(color="grey")+
      geom_smooth(method=lm,color = "#c9d6df")+
      # xlim(-3,3)+
      # 图例
      theme(legend.position = "bottom",panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
      theme_bw()+theme(text = element_text(size = 18))+
      geom_abline(intercept = 0, slope = 1,color="red") 
    ggsave(paste0("result/FA_vs_noFA/",tissue,"/",tissue,"_",antibody,"_cpm_scatter_plot_bcv01.png"),width = 10,height = 10)
  }
}

#H3K27ac outliers
tissue <- "liver"
antibody <- "H3K27ac"
tab = read.delim(paste0("data/FA_vs_noFA/",tissue,"/",antibody,"/",antibody,"_macs2callsplit_mergecounts.counts"),skip=1)
counts = tab[,c(7:8)]
rownames(counts)= tab$Geneid
# counts <- cbind(counts,counts)
colnames(counts) = c("FA_1","noFA_1")
colnames(tab)[7:8] <- c("FA_1","noFA_1")
# counts <- counts[,c("noFA_1","noFA_2","FA_1","FA_2")]
group =c("FA","noFA")
design = model.matrix(~0+group)
lvls = levels(factor(group))
colnames(design) = lvls
len = length(lvls)
myContrasts = c("FA-noFA")
contrast.matrix = eval(as.call(c(as.symbol("makeContrasts"),as.list(myContrasts),levels=list(design))))
y= DGEList(counts=counts,group=group)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
y =  calcNormFactors(y, method ="TMM")
bcv <- 0.2
y <- exactTest(y,dispersion = bcv^2, pair = c("noFA","FA"))
# y<-estimateCommonDisp(y)
# y<-estimateGLMTagwiseDisp(y,design)
# fit_tag = glmFit(y,design)
# lrt = glmLRT(fit_tag, contrast = contrast.matrix)
# qBH = p.adjust(lrt$table$PValue,method="BH")
tab<-tab[keep,]
out = cbind(tab[,1:6],cpm(tab[,7:8]),logFC=y$table$logFC)
out[,c(7:8)] = sweep(out[,c(7:8)],1,out$Length,'/')*1000
for (i in 1:length(myContrasts)){
  out[,paste0("PValue.",myContrasts[i])] = y$table$PValue
  out[,paste0("FDR.",myContrasts[i])] = p.adjust(y$table$PValue,method="BH")
}
out$Significant <- ifelse(out$`FDR.FA-noFA` < 0.05 & abs(out$logFC) >= 1, 
                          ifelse(out$logFC > 1, "Up", "Down"), "Stable")

# scatter_data <- out[,c(1:7,9)]
# scatter_data[,c(7:8)] <- log10(scatter_data[,c(7:8)])
# # scatter_data_outliers <- scatter_data[which(scatter_data$noFA_1>1.5 & scatter_data$FA_1 <1),]
# scatter_data_outliers <- scatter_data[which(scatter_data$noFA_1>1.5 & scatter_data$FA_1 > 2.4),]

GO_database <- 'org.Mm.eg.db'
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
peak <- readPeakFile(paste0("data/FA_vs_noFA/",tissue,"/",antibody,"/bed/",antibody,"_macs2callsplit_mergecounts.bed"))
peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnno<-as.data.frame(peakAnno)
out$annotation <- peakAnno$annotation[which(peakAnno$V4 %in% out$Geneid)]
out$ens_id <- peakAnno$geneId[which(peakAnno$V4 %in% out$Geneid)]
out$symbol <- peakAnno$SYMBOL[which(peakAnno$V4 %in% out$Geneid)]
out_gene<-out[-which(is.na(out$symbol)),]
out_gene<-out_gene[which(str_detect(out_gene$annotation,"Promoter")),]
table(out_gene$annotation)
out_gene$Significant <- ifelse(out_gene$`FDR.FA-noFA` < 0.05 & abs(out_gene$logFC) >= 1, 
                          ifelse(out_gene$logFC > 1, "Up", "Down"), "Stable")

genelist_up <-bitr(out_gene$symbol[which(out_gene$`FDR.FA-noFA` < 0.05 & out_gene$logFC >= 0)],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_up_GO <-enrichGO( genelist_up$ENTREZID,#GO富集分析
                           OrgDb = GO_database,
                           keyType = "ENTREZID",#设定读取的gene ID类型
                           ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                           pvalueCutoff = 0.05,#设定p值阈值
                           qvalueCutoff = 0.05,#设定q值阈值
                           readable = T)
dotplot(genelist_up_GO,font.size=15,label_format = 50)
genelist_up_GO_result <- genelist_up_GO@result
write.csv(genelist_up_GO_result,paste0("result/FA_vs_noFA/",tissue,"/",tissue,"_",antibody,"_increase_logFC02_GO.csv"))

genelist_down <-bitr(out_gene$symbol[which(out_gene$`FDR.FA-noFA` < 0.05 & out_gene$logFC <= -0.2)],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database) 
genelist_down_GO <-enrichGO( genelist_down$ENTREZID,#GO富集分析
                             OrgDb = GO_database,
                             keyType = "ENTREZID",#设定读取的gene ID类型
                             ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                             pvalueCutoff = 0.05,#设定p值阈值
                             qvalueCutoff = 0.05,#设定q值阈值
                             readable = T)
dotplot(genelist_down_GO,font.size=15,label_format = 50)

genelist_down_GO_result <- genelist_down_GO@result
write.csv(genelist_down_GO_result,paste0("result/FA_vs_noFA/",tissue,"/",tissue,"_",antibody,"_decrease_logFC02_GO.csv"))

# ggplot(
#   # 数据、映射、颜色
#   out_gene, aes(x = logFC, y = -log10(`FDR.FA-noFA`))) +
#   geom_point(aes(color = Significant), size=2) +
#   scale_color_manual(values = colour[[nrow(as.data.frame(table(out_gene$Significant)))]]) +
#   geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
#   # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
#   geom_hline(yintercept = 1,lty=4,col="black",lwd=0.8) +
#   # 坐标轴
#   labs(x="log2(fold change)",
#        y="-log10 (p-value)") +
#   # xlim(-3,3)+
#   # 图例
#   theme(legend.position = "bottom",panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())+
#   theme_bw()+theme(text = element_text(size = 18))

# ggplot(
#   # 数据、映射、颜色
#   out_gene[which(out_gene$annotation=="Promoter (1-2kb)"),], aes(x = logFC, y = -log10(`FDR.FA-noFA`))) +
#   geom_point(aes(color = Significant), size=2) +
#   scale_color_manual(values = colour[[nrow(as.data.frame(table(out_gene$Significant)))]]) +
#   geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
#   # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
#   geom_hline(yintercept = 1,lty=4,col="black",lwd=0.8) +
#   # 坐标轴
#   labs(x="log2(fold change)",
#        y="-log10 (p-value)") +
#   # xlim(-3,3)+
#   # 图例
#   theme(legend.position = "bottom",panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())+
#   theme_bw()+theme(text = element_text(size = 18))
# 
# ggplot(
#   # 数据、映射、颜色
#   out_gene[which(out_gene$annotation=="Promoter (1-2kb)"),], aes(x = logFC, y = -log10(`FDR.FA-noFA`))) +
#   geom_point(aes(color = Significant), size=2) +
#   scale_color_manual(values = colour[[nrow(as.data.frame(table(out_gene$Significant)))]]) +
#   geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
#   # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
#   geom_hline(yintercept = 1,lty=4,col="black",lwd=0.8) +
#   # 坐标轴
#   labs(x="log2(fold change)",
#        y="-log10 (p-value)") +
#   # xlim(-3,3)+
#   # 图例
#   theme(legend.position = "bottom",panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())+
#   theme_bw()+theme(text = element_text(size = 18))
# 
# ggplot(
#   # 数据、映射、颜色
#   out_gene[which(out_gene$annotation=="Promoter (2-3kb)"),], aes(x = logFC, y = -log10(`FDR.FA-noFA`))) +
#   geom_point(aes(color = Significant), size=2) +
#   scale_color_manual(values = colour[[nrow(as.data.frame(table(out_gene$Significant)))]]) +
#   geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
#   # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
#   geom_hline(yintercept = 1,lty=4,col="black",lwd=0.8) +
#   # 坐标轴
#   labs(x="log2(fold change)",
#        y="-log10 (p-value)") +
#   # xlim(-3,3)+
#   # 图例
#   theme(legend.position = "bottom",panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())+
#   theme_bw()+theme(text = element_text(size = 18))

# genelist_up <-bitr(out_gene$symbol[which(out_gene$`FDR.FA-noFA` < 0.05 & out_gene$logFC >= 1 & out_gene$annotation=="Promoter (<=1kb)")],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
# genelist_up_GO <-enrichGO( genelist_up$ENTREZID,#GO富集分析
#                            OrgDb = GO_database,
#                            keyType = "ENTREZID",#设定读取的gene ID类型
#                            ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
#                            pvalueCutoff = 0.05,#设定p值阈值
#                            qvalueCutoff = 0.05,#设定q值阈值
#                            readable = T)
# dotplot(genelist_up_GO,font.size=15,label_format = 50)
# genelist_up_GO_result <- genelist_up_GO@result



# genelist_down <-bitr(out_gene$symbol[which(out_gene$`FDR.FA-noFA` < 0.05 & out_gene$logFC <= -1 & out_gene$annotation=="Promoter (<=1kb)")],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database) 
# genelist_down_GO <-enrichGO( genelist_down$ENTREZID,#GO富集分析
#                              OrgDb = GO_database,
#                              keyType = "ENTREZID",#设定读取的gene ID类型
#                              ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
#                              pvalueCutoff = 0.05,#设定p值阈值
#                              qvalueCutoff = 0.05,#设定q值阈值
#                              readable = T)
# dotplot(genelist_down_GO,font.size=15,label_format = 50)


