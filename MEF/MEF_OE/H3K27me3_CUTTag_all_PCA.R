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
tab <- read.table(paste0("data/samples/MEF_OE/H3K27me3/H3K27me3_all_10kb.counts"),header = T)
new_names <- sub("^.*\\.bam\\.([A-Za-z]+\\d+\\.\\d+)\\.bam.*$", "\\1", colnames(tab)[7:ncol(tab)])
new_names <- gsub("\\.", "-", new_names)
colnames(tab)[7:ncol(tab)] <- new_names
tab <- tab[!grepl("chrY", tab$Chr), ]
rownames(tab) <- tab$Geneid
tab <- tab[,-1]

tab_MEF_senescence <- read.table(paste0("data/samples/MEF/H3K27me3/H3K27me3_10kb_bins.counts"),header = T)
new_names <- sub("^.*\\.bam\\.(NTY\\d+).*$", "\\1",colnames(tab_MEF_senescence)[7:ncol(tab_MEF_senescence)])
new_names <- gsub("\\.", "-", new_names)
colnames(tab_MEF_senescence)[7:ncol(tab_MEF_senescence)] <- new_names
tab_MEF_senescence <- tab_MEF_senescence[!grepl("chrY", tab_MEF_senescence$Chr), ]
rownames(tab_MEF_senescence) <- tab_MEF_senescence$Geneid
tab_MEF_senescence <- tab_MEF_senescence[,-1]

# tab_MEF_inhibit <- read.table(paste0("data/samples/MEF_EZH2_inhibit/H3K27me3/H3K27me3_10kb_bins.counts"),header = T)
# new_names <- sub("^.*\\.bam\\.(DYQ\\d+).*$", "\\1",colnames(tab_MEF_inhibit)[7:ncol(tab_MEF_inhibit)])
# new_names <- gsub("\\.", "-", new_names)
# colnames(tab_MEF_inhibit)[7:ncol(tab_MEF_inhibit)] <- new_names
# tab_MEF_inhibit <- tab_MEF_inhibit[!grepl("chrY", tab_MEF_inhibit$Chr), ]
# rownames(tab_MEF_inhibit) <- tab_MEF_inhibit$Geneid
# tab_MEF_inhibit <- tab_MEF_inhibit[,-1]

tab <- merge(tab[,c(6:ncol(tab))],tab_MEF_senescence[,c(6:ncol(tab_MEF_senescence))],by="row.names")
# tab <- merge(tab,tab_MEF_inhibit[,c(6:ncol(tab_MEF_inhibit))],by.x="Row.names",by.y="row.names")
rownames(tab) <- tab$Row.names
counts <- tab[,-1]

search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
search_table <- search_table[!duplicated(search_table$sample_name), ]
# search_table$mouse_ID <- c("p6","p6","p10","p10","DMSO","DMSO","GSK126","GSK126","vec","vec","Bmi1","Bmi1","Cbx2","Cbx2","Cbx7","Cbx7")
search_table$mouse_ID <- c("p2","p2","p10","p10","vec","vec","Bmi1","Bmi1","Cbx2","Cbx2","Cbx7","Cbx7")

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

# batch <- c(1,1,1,1,2,2,2,2,3,3,3,3,3,3,3,3)
batch <- c(1,1,1,1,2,2,2,2,2,2,2,2)
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
