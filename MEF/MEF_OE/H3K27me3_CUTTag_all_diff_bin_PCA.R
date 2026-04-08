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


tab <- merge(tab[,c(6:ncol(tab))],tab_MEF_senescence[,c(6:ncol(tab_MEF_senescence))],by="row.names")
rownames(tab) <- tab$Row.names
counts <- tab[,-1]

# mm10 <- read.table("~/ref_data/mm10_10kb_bins.bed")
# mm10 <- as.data.table(mm10)
# setDT(mm10)
# setkey(mm10,V1,V2,V3)
# H3K9me3_peaks <- read.table("data/samples/MEF/H3K9me3/bed/H3K9me3_young_merge-W5000-G10000-E100.bed")
# H3K9me3_peaks <- as.data.table(H3K9me3_peaks[,c(1:3)])
# setDT(H3K9me3_peaks)
# setkey(H3K9me3_peaks,V1,V2,V3)
# overlaps <- as.data.frame(foverlaps(mm10,H3K9me3_peaks, type = "any", nomatch = 0L))
# counts <- counts[which(!rownames(counts) %in% overlaps$V4),]

search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
search_table <- search_table[!duplicated(search_table$sample_name), ]
search_table$mouse_ID <- c("p2","p2","p10","p10","vec","vec","Bmi1","Bmi1","Cbx2","Cbx2","Cbx7","Cbx7")

diff_bin <- c()
for(condition in c("Bmi1","Cbx2","Cbx7")){
  diff <- read.csv(paste0("data/samples/MEF_OE/H3K27me3/MEF_",condition,"/H3K27me3_MEF_",condition,"_10kb_bins_diff_after_remove_batch_effect.csv"))
  diff_bin <- c(diff_bin,diff$Geneid[which(diff$Significant !="Stable")])
}
diff <- read.csv(paste0("data/samples/MEF/H3K27me3/H3K27me3_10kb_bins_diff_after_remove_batch_effect.csv"))
diff_bin <- c(diff_bin,diff$Geneid[which(diff$Significant !="Stable")])

# diff_bin <- c()
# for(condition in c("Bmi1","Cbx2","Cbx7")){
#   diff <- read.csv(paste0("data/samples/MEF_OE/H3K27me3/MEF_",condition,"/H3K27me3_MEF_",condition,"_10kb_bins_diff_after_remove_batch_effect.csv"))
#   if(condition == "Bmi1"){
#     diff_bin <- diff$Geneid[which(diff$Significant !="Stable")]
#   }else{
#     diff_bin <- intersect(diff_bin,diff$Geneid[which(diff$Significant !="Stable")])
#   }
# }
# 
# diff <- read.csv(paste0("data/samples/MEF/H3K27me3/H3K27me3_10kb_bins_diff_after_remove_batch_effect.csv"))
# diff_bin <- intersect(diff_bin,diff$Geneid[which(diff$Significant !="Stable")])

diff_bin <- unique(diff_bin)

counts <- counts[,search_table$sample_name]
counts <- counts[which(rownames(counts) %in% diff_bin),]
y= DGEList(counts=counts,group=search_table$mouse_ID)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]

logCPMs <- cpm(y, log = TRUE)
pca <- prcomp(t(logCPMs))
to_plot <- data.frame(pca$x)
to_plot$sample_name <- rownames(to_plot)

to_plot <- merge(to_plot,search_table,by="sample_name")
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
to_plot$label <- paste0(to_plot$sample_name,"-",to_plot$mouse_ID)
p <- ggplot(to_plot, aes(x=PC1, y=PC2, color=mouse_ID)) + 
  geom_point(size=5) +theme_bw()+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
  geom_text_repel(  
    data = to_plot,  
    aes(x = PC1, y = PC2, label = label),  
    size = 5,  
    box.padding = unit(0.35, "lines"),  
    point.padding = unit(0.3, "lines")  
  )
ggsave("result/figures/MEF_OE_H3K27me3_PCA.pdf",p,height = 6,width = 7.5)

batch <- c(1,1,1,1,2,2,2,2,2,2,2,2)
logCPMs_corrected <- limma::removeBatchEffect(logCPMs, batch = batch)
pca <- prcomp(t(logCPMs_corrected))
to_plot <- data.frame(pca$x)
to_plot$sample_name <- rownames(to_plot)

to_plot <- merge(to_plot,search_table,by="sample_name")
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
