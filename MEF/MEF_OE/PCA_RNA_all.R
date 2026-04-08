rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(patchwork)
library(edgeR)
library(MASS) 
library(Seurat)
library(gridExtra)
library(ggrepel)
library(stringr)

tab <- read.table(paste0("data/samples/RNA/MEF_OE_RNA/all_merge.counts"),header = T)
new_names <- sub("^.*\\.bam\\.([A-Za-z]+\\d+\\.\\d+)(?:[_.].*)?$", "\\1", colnames(tab)[7:ncol(tab)])
new_names <- gsub("\\.", "-", new_names)
colnames(tab)[7:ncol(tab)] <- new_names
tab <- tab[!grepl("chrY", tab$Chr), ]
rownames(tab) <- tab$Geneid
tab <- tab[,-1]
tab <- tab[,c(1:13,20,21)]
tab_MEF_senescence <- read.table(paste0("data/samples/RNA/MEF/combined-chrM.counts"),header = T)
new_names <- sub("^.*\\bTK(\\d+).*$", "TK\\1",colnames(tab_MEF_senescence)[7:ncol(tab_MEF_senescence)])
new_names <- gsub("\\.", "-", new_names)
colnames(tab_MEF_senescence)[7:ncol(tab_MEF_senescence)] <- new_names
tab_MEF_senescence <- tab_MEF_senescence[!grepl("chrY", tab_MEF_senescence$Chr), ]
rownames(tab_MEF_senescence) <- tab_MEF_senescence$Geneid
tab_MEF_senescence <- tab_MEF_senescence[,-1]

# tab_MEF_inhibit <- read.table(paste0("data/samples/RNA/MEF_EZH2_inhibit/combined-chrM.counts"),header = T)
# new_names <- sub("^.*\\.bam\\.(XM\\d+).*$", "\\1",colnames(tab_MEF_inhibit)[7:ncol(tab_MEF_inhibit)])
# new_names <- gsub("\\.", "-", new_names)
# colnames(tab_MEF_inhibit)[7:ncol(tab_MEF_inhibit)] <- new_names
# tab_MEF_inhibit <- tab_MEF_inhibit[!grepl("chrY", tab_MEF_inhibit$Chr), ]
# rownames(tab_MEF_inhibit) <- tab_MEF_inhibit$Geneid
# tab_MEF_inhibit <- tab_MEF_inhibit[,-1]

tab <- merge(tab[,c(6:ncol(tab))],tab_MEF_senescence[,c(6:ncol(tab_MEF_senescence))],by="row.names")
# tab <- merge(tab,tab_MEF_inhibit[,c(6:ncol(tab_MEF_inhibit))],by.x="Row.names",by.y="row.names")

rownames(tab) <- tab$Row.names
counts <- tab[,-1]

search_table <- read.csv("data/samples/all/RNA_search_table.csv")
search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
search_table <- search_table[!duplicated(search_table$sample_name), ]
search_table$mouse_ID[7:16] <- c("vec","vec","Bmi1","Bmi1","Cbx2","Cbx2","Cbx7","Cbx7","Cbx8","Cbx8")
# search_table$mouse_ID[7:22] <- c("vec","vec","Bmi1","Bmi1","Cbx2","Cbx2","Cbx7","Cbx7","vec","vec","mEzh2","mEzh2","hEzh2","hEzh2","mCbx8","mCbx8")
counts <- counts[,search_table$sample_name]
y= DGEList(counts=counts,group=search_table$mouse_ID)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]

logCPMs <- cpm(y, log = TRUE)
pca <- prcomp(t(logCPMs))
to_plot <- data.frame(pca$x)
to_plot$sample_name <- rownames(to_plot)

to_plot <- merge(to_plot,search_table,by="sample_name")
# to_plot$age <- factor(to_plot$age,levels=c("vec","oe"))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
to_plot$label <- paste0(to_plot$sample_name,"-",to_plot$mouse_ID)
p <-ggplot(to_plot, aes(x=PC1, y=PC2, color=mouse_ID)) + 
  geom_point(size=5) +theme_bw()+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
  geom_text_repel(  
    data = to_plot,  
    aes(x = PC1, y = PC2, label = label),  
    size = 5,  
    box.padding = unit(0.35, "lines"),  
    point.padding = unit(0.3, "lines")  
  )
ggsave("result/figures/MEF_OE_RNA_PCA.pdf",p,height = 6,width = 7.5)

batch <- c(1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2)
logCPMs_corrected <- limma::removeBatchEffect(logCPMs, batch = batch)
pca <- prcomp(t(logCPMs_corrected))
to_plot <- data.frame(pca$x)
to_plot$sample_name <- rownames(to_plot)

to_plot <- merge(to_plot,search_table,by="sample_name")
# to_plot$age <- factor(to_plot$age,levels=c("vec","oe"))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
to_plot$label <- paste0(to_plot$sample_name,"-",to_plot$mouse_ID)
ggplot(to_plot, aes(x=PC1, y=PC2, color=mouse_ID)) + 
  geom_point(size=5) +theme_bw()+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
  geom_text_repel(  
    data = to_plot,  
    aes(x = PC1, y = PC2, label = label),  
    size = 5,  
    box.padding = unit(0.35, "lines"),  
    point.padding = unit(0.3, "lines")  
  )
