rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
tissue <- "MEF_mid_age"
out <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"),row.names = 1)
top_genes <- out[which(out$Significant=="Up"),] %>%  
  arrange(desc(logFC)) %>%  
  head(10)  
bottom_genes <- out[which(out$Significant=="Down"),] %>%  
  arrange(logFC) %>%  
  head(10)  
highlight_genes <- rbind(top_genes, bottom_genes) 
highlight_genes$Gene <- rownames(highlight_genes)
colour=setNames(c("blue","grey","red"),c("Down","Stable","Up"))
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
  annotate("text", x = min(out$logFC), y = max(-log10(out$fdr)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
  annotate("text", x = max(out$logFC), y = max(-log10(out$fdr)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)+
  geom_text_repel(data = highlight_genes, aes(logFC,-log10(fdr), label = Gene), max.overlaps=100,
                  size = 5, 
                  nudge_y = 0.2)
