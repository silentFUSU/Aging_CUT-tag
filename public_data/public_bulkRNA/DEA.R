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
table <- table[which(str_detect(table$source.name,"Brain")),]
table <- table[which(table$characteristics..sex=="m"),]
table$characteristics..age <- as.numeric(table$characteristics..age)
table <- table[order(table$characteristics..age), ]
new_column_name <- sub("\\..*", "", colnames(tab[,-1])) 
colnames(tab)[2:ncol(tab)] <- new_column_name
counts <- tab[,which(colnames(tab) %in% table$Sample.name)]
table <- table[which(table$Sample.name %in% colnames(counts)),]

counts <- cbind(tab[,which(colnames(tab) %in% table$Sample.name[which(table$characteristics..age %in% c(3,27))])])
table <- table[which(table$characteristics..age %in% c(3,27)),]
rownames(counts) <- tab$gene
counts <- counts[,match(table$Sample.name,colnames(counts),)]
colnames(counts) <- c(paste0("m",table$characteristics..age,"_",c(1:nrow(table))))
counts <- counts[,c(1,2,5,6)]
colnames(counts) <- c("Young1","Young2","Old1","Old2")
counts$Geneid <- rownames(counts)
counts <- counts[,c("Geneid","Young1","Young2","Old1","Old2")]
write.table(counts,"~/projects/bioinfo_class/homework/lab2/lab2_counts.txt",row.names = F,append = F,sep = "\t")
group =c(paste0("m",table$characteristics..age))



y= DGEList(counts=counts,group=group)
y$samples$group <- factor(y$samples$group,levels = c("m3","m24"))
design <- model.matrix(~group, y$samples)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
y =  calcNormFactors(y)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)
lrt = glmLRT(fit_tag, coef = ncol(design))
tab<-tab[keep,]

out = cbind(cpm(y),lrt$table, "fdr"=p.adjust(lrt$table$PValue,method="BH"))
out$Significant <- ifelse(out$fdr< 0.05 & abs(out$logFC) >= 0, 
                          ifelse(out$logFC > 0, "Up", "Down"), "Stable")
colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
top_genes <- out[which(out$Significant=="Up"),] %>%  
  arrange(desc(logFC)) %>%  
  head(10)  
bottom_genes <- out[which(out$Significant=="Down"),] %>%  
  arrange(logFC) %>%  
  head(10)  
highlight_genes <- rbind(top_genes, bottom_genes) 
highlight_genes$Gene <- rownames(highlight_genes)
ggplot(
  out, aes(x = logFC, y = -log10(fdr))) +
  geom_point(aes(color = Significant), size=2) +
  scale_color_manual(values = colour) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (fdr)") +
  theme_bw()+
  theme(text = element_text(size = 20))+
  ggtitle("27m-6m")+
  annotate("text", x = min(out$logFC), y = max(-log10(out$fdr)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
  annotate("text", x = max(out$logFC), y = max(-log10(out$fdr)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)+
  geom_text_repel(data = highlight_genes, aes(logFC,-log10(fdr), label = Gene), max.overlaps=100,
                  size = 5, # 字体大小  
                  nudge_y = 0.2)

write.csv(out,"data/public_data/GSE132040/diff/Skin_27-6.csv",)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
genelist_up <- bitr(rownames(out)[which(out$Significant=="Up")],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_up_GO <- enrichGO( genelist_up$ENTREZID,#GO富集分析
                            OrgDb = GO_database,
                            keyType = "ENTREZID",#设定读取的gene ID类型
                            ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                            pvalueCutoff = 0.05,#设定p值阈值
                            qvalueCutoff = 0.05,#设定q值阈值
                            readable = T)
barplot(genelist_up_GO,title = paste0("Skin Increased gene GO pathway"),label_format = 50)
GO_table <- genelist_up_GO@result
write.csv(GO_table,paste0("data/public_data/GSE132040/diff/Skin_27-6_increase_GO.csv"))

genelist_down <- bitr(rownames(out)[which(out$Significant=="Down")],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_down_GO <- enrichGO( genelist_down$ENTREZID,#GO富集分析
                              OrgDb = GO_database,
                              keyType = "ENTREZID",#设定读取的gene ID类型
                              ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                              pvalueCutoff = 0.05,#设定p值阈值
                              qvalueCutoff = 0.05,#设定q值阈值
                              readable = T)
barplot(genelist_down_GO,title = paste0("Skin Decreased gene GO pathway"),label_format = 50)
GO_table <- genelist_down_GO@result
write.csv(GO_table,paste0("data/public_data/GSE132040/diff/Skin_27-6_decrease_GO.csv"))

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
