rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GO.db)
library(DOSE)
library(data.table)
data(geneList)
options(bitmapType="cairo")  
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
      tissue_label <- "IWAT"
    }
  }
  return(tissue_label)
}

tissues <- sort(c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
                  "thymus","skin","bladder","bonemarrow","Hip","heart",
                  "muscle","jejunum","uterus","ovary","liver","tongue",
                  "cecum","colon","testis","stomach","pancreas","iWAT","ileum"))

second_quadrant_summary <- data.frame()
fourth_quadrant_summary <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_gene_TSS_10kb_diff_after_remove_batch_effect.csv"))
  
  peaks <- read.table(paste0("data/samples/",tissue,"/H3K27me3/bed/H3K27me3_young_merge-W5000-G10000-E100.bed"))
  peaks <- as.data.table(peaks)
  setDT(peaks)
  setkey(peaks,V1,V2,V3)
  gene_tss <- df[,c(1:4)]
  gene_tss <- as.data.table(gene_tss)
  gene_tss$Start <- as.numeric(gene_tss$Start)
  setDT(gene_tss)
  setkey(gene_tss,Chr,Start,End)
  overlaps <- foverlaps(gene_tss,peaks, type = "any", nomatch = 0L)
  df <- df[which(df$Geneid %in% overlaps$Geneid),]

  df <- df[which(df$Significant != "Stable"),]
  RNA <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
  RNA <- RNA[which(RNA$Significant != "Stable"),]
  colnames(RNA)[1] <- "Geneid"
  df <- merge(df[,c("Geneid","LogFC.old.young")],RNA[,c("Geneid","logFC")],by="Geneid")
  second_quadrant <- df[which(df$LogFC.old.young > 0 & df$logFC < 0),]
  fourth_quadrant <- df[which(df$LogFC.old.young < 0 & df$logFC > 0),]
  
  if(nrow(second_quadrant) >0){
    second_quadrant$tissue <- tissue_label_change(tissue)
    second_quadrant_summary <- rbind(second_quadrant_summary,second_quadrant[,c(1,4)])
  }
  if(nrow(fourth_quadrant >0)){
    fourth_quadrant$tissue <- tissue_label_change(tissue)
    fourth_quadrant_summary <- rbind(fourth_quadrant_summary,fourth_quadrant[,c(1,4)])
  }
}
second_quadrant_summary_count <- second_quadrant_summary %>%   
  count(Geneid)
second_quadrant_summary_tissue <- second_quadrant_summary %>%   
  group_by(Geneid) %>%   
  summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
second_quadrant_summary_count <- merge(second_quadrant_summary_count,second_quadrant_summary_tissue,by="Geneid")
# write.csv(second_quadrant_summary_count[order(second_quadrant_summary_count$n,decreasing = T),],"data/samples/all/H3K27me3/with_RNA/H3K27me3_increased_gene_expression_decreased_common_genes.csv")

fourth_quadrant_summary_count <- fourth_quadrant_summary %>%   
  count(Geneid)
fourth_quadrant_summary_tissue <- fourth_quadrant_summary %>%   
  group_by(Geneid) %>%   
  summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
fourth_quadrant_summary_count <- merge(fourth_quadrant_summary_count,fourth_quadrant_summary_tissue,by="Geneid")
# write.csv(fourth_quadrant_summary_count[order(fourth_quadrant_summary_count$n,decreasing = T),],"data/samples/all/H3K27me3/with_RNA/H3K27me3_decreased_gene_expression_increased_common_genes.csv")

to_plot <- as.data.frame(table(second_quadrant_summary_count$n))
ggplot(data = to_plot, aes(x = Var1, y = Freq)) +  
  geom_bar(stat = "identity") +  
  labs(  
    title = "H3K27me3 increased and Gene expression decreased",  
    x = NULL,  
    y = "Count"  
  ) +  
  theme_minimal() +  
  xlab(NULL) +  
  geom_text(  
    aes(label = Freq),   
    vjust = 0,   
    size = 3.5  # Adjust text size as needed  
  ) +
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
    axis.text.y = element_text(size = 12),  
    plot.title = element_text(size = 16, face = "bold")  
  ) 

to_plot <- as.data.frame(table(fourth_quadrant_summary_count$n))
p <- ggplot(data = to_plot, aes(x = Var1, y = Freq)) +  
  geom_bar(stat = "identity") +  
  labs(  
    title = "H3K27me3 decreased and Gene expression increased",  
    x = NULL,  
    y = "Count"  
  ) +  
  theme_bw() +  
  xlab(NULL) +  
  geom_text(  
    aes(label = Freq),   
    vjust = 0,   
    size = 3.5  # Adjust text size as needed  
  ) +
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
    axis.text.y = element_text(size = 12),  
    plot.title = element_text(size = 16, face = "bold")  
  ) 
ggsave("result/figures/H3K27me3_decreased_RNA_increased_in_young_peaks_commom_genes.pdf",p,width = 6,height = 4)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
genelist <- fourth_quadrant_summary_count$Geneid[which(fourth_quadrant_summary_count$n >=4)]
genelist <- bitr(genelist,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelistGO <- enrichGO( genelist$ENTREZID,#GO富集分析
                        OrgDb = GO_database,
                        keyType = "ENTREZID",#设定读取的gene ID类型
                        ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                        pvalueCutoff = 0.05,#设定p值阈值
                        qvalueCutoff = 0.05,#设定q值阈值
                        readable = T)
p <- barplot(genelistGO,label_format = 50)
ggsave("result/figures/H3K27me3_decreasd_gene_increased_in_H3K27me3_young_peaks_larger4_GO.pdf",p,width = 6,height = 8)
