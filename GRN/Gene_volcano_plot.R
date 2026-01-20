rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
library(ggrepel)
library(tidyr)
library(dplyr)
library(tidyverse)
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
    }else if(tissue_label == "Bat"){
      tissue_label <- "BAT"
    }else if(tissue_label=="Mammarygland"){
      tissue_label <- "Mammary Gland"
    }else if(tissue_label=="Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
}
tissue <- "ovary"
TF_gene <- "Esr2"
load("data/samples/GRN/grn_union_tissue.rdata")
load("data/samples/GRN/grn_union_skin.rdata")
grn_tissue[["skin"]] <- grn_union
TF <- grn_tissue[[tissue]]
cre <- TF[which(TF$TF==TF_gene | TF$gene %in% c("Htra1")),]
cre <- as.data.frame(do.call(rbind, strsplit(as.character(TF$peak), "_")))
cre$V2 <- as.numeric(cre$V2)
cre$V3 <- cre$V2+501
write.table(cre,paste0("data/samples/ATAC/",tissue,"/ATAC/bed/",TF_gene,"_cre.bed"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)
TF <- TF[which(TF$TF==TF_gene),]

RNA <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
RNA_TF <- RNA[which(RNA$X %in% TF$gene),]
RNA_TF$Significant[which(RNA_TF$Significant=="Down")] <- "Significant decrease"
RNA_TF$Significant[which(RNA_TF$Significant=="Up")] <- "Significant increase"
RNA_TF$Significant[which(RNA_TF$Significant=="Stable" & RNA_TF$logFC > 0 )] <- "increase"
RNA_TF$Significant[which(RNA_TF$Significant=="Stable" & RNA_TF$logFC < 0 )] <- "decrease"
highlight_genes <- subset(RNA_TF, X %in% c("Lama1","Hao2","Cst8","Idh1"))


p1 <- ggplot(
  RNA, aes(x = logFC, y = -log10(fdr))) +
  geom_point(color="gray",size=2) +
  # geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=2,col="black",lwd=0.8) +
  labs(x="log2(O/Y)",
       y="-log10 (fdr)") +
  theme_bw()+
  theme(text = element_text(size = 20))+
  ggtitle(paste0(tissue_label_change(tissue),"-",TF_gene))
p1
p2 <- p1 +
  geom_point(data=RNA_TF, aes(x = logFC, y = -log10(fdr), color = Significant), size=2,alpha=0.5) +xlim(-15,15)+
  scale_color_manual(values = c("Significant decrease" = "blue", "Significant increase" = "red","increase"="dark red","decrease"="dark blue"))+
  geom_text_repel(data = highlight_genes, aes(label = X), size = 5, box.padding = 0.5, point.padding = 0.3, segment.color = "black")+
  geom_point(data = highlight_genes, aes(x = logFC, y = -log10(fdr)), color = "black", fill = "blue", shape = 21, size = 2.5, stroke = 1)+
  annotate("text", x = min(RNA$logFC), y = max(-log10(RNA$fdr)), label = paste0(nrow(RNA_TF[which(RNA_TF$Significant=="Significant decrease"),]),"/",nrow(RNA[which(RNA$Significant=="Down"),])), vjust = 5, hjust = 0,colour="dark blue",size=5)+
  annotate("text", x = max(RNA$logFC), y = max(-log10(RNA$fdr)), label = paste0(nrow(RNA_TF[which(RNA_TF$Significant=="Significant increase"),]),"/",nrow(RNA[which(RNA$Significant=="Up"),])), vjust = 5, hjust = 1.5,colour="dark red",size=5)
ggsave(paste0("result/figures/ovary_",TF_gene,"_volcano_plot.pdf"),p2,width = 7,height = 6)

proportion <- as.data.frame(table(RNA_TF$Significant))
proportion$percentage <- proportion$Freq / sum(proportion$Freq) *100
