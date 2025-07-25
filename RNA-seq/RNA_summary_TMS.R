rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(stringr)
library(edgeR)
gene <- "Cdkn2a"
tab = read.csv("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/GSE132040/GSE132040_190214_A00111_0269_AHH3J3DSXX_190214_A00111_0270_BHHMFWDSXX.csv")
tab <- tab[-c(54353:54357),]
tab <- as.data.frame(t(tab))
colnames(tab) <- tab[1,]
tab <- tab[-1,]
new_row_name <- sub("\\..*", "", rownames(tab)) 
tab$Sample.name <- new_row_name
table = read.delim("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/GSE132040/GSE132040_MACA_Bulk_metadata.csv",sep = ',')
table$tissue <- sub("_\\d.*$", "", table$source.name)
tab <- merge(tab,table[,c(1,3,5,7)],by="Sample.name")
tab <- tab[!grepl("^NA", tab$source.name), ] 
rownames(tab) <- tab$Sample.name
tab <- tab[,-1]
tab <- tab[which(tab$characteristics..sex=="m"),]
tab <- tab[which(tab$characteristics..age %in% c("3","24","27")),]
count <- as.data.frame(t(tab[,-c(54353:54355)])) 
count <- as.data.frame(apply(count, 2, as.numeric))
rownames(count) <- colnames(tab[,-c(54353:54355)])
y= DGEList(counts=count)
cpm <- as.data.frame(cpm(y))
to_plot <- cpm[gene,]
to_plot <- as.data.frame(t(to_plot))
to_plot$Sample.name <- rownames(to_plot)

to_plot <- merge(to_plot,table[,c("Sample.name","characteristics..age","tissue")],by ="Sample.name")
colnames(to_plot)[2:3]<-c("CPM","Age")
to_plot$Age <- factor(to_plot$Age,levels=c("1","3","6","9","12","15","18","21","24","27"))
ggplot(to_plot,aes(x=tissue,y=CPM,fill=Age))+    
  geom_boxplot()+
  ggtitle(paste0(gene," CPM"))+
  theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")+labs(fill = "", color = "") +ylab("CPM")
