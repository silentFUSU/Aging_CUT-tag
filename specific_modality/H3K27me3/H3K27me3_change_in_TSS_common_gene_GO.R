rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GO.db)
library(DOSE)
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

increase_summary <- data.frame()
decrease_summary <- data.frame()
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
  df <- df[which(df$Geneid %in% overlaps$Geneid & df$Significant != "Stable"),]

  
  increase <- df[which(df$Significant=="Up"),c("Geneid","Significant")]
  decrease <- df[which(df$Significant=="Down"),c("Geneid","Significant")]
  
  if(nrow(increase) >0){
    increase$tissue <- tissue_label_change(tissue)
    increase_summary <- rbind(increase_summary,increase)
  }
  if(nrow(decrease >0)){
    decrease$tissue <- tissue_label_change(tissue)
    decrease_summary <- rbind(decrease_summary,decrease)
  }
}

increase_summary_count <- increase_summary %>%   
  count(Geneid)
increase_summary_tissue <- increase_summary %>%   
  group_by(Geneid) %>%   
  summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
increase_summary_count <- merge(increase_summary_count,increase_summary_tissue,by="Geneid")


decrease_summary_count <- decrease_summary %>%   
  count(Geneid)
decrease_summary_tissue <- decrease_summary %>%   
  group_by(Geneid) %>%   
  summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
decrease_summary_count <- merge(decrease_summary_count,decrease_summary_tissue,by="Geneid")

to_plot <- as.data.frame(table(increase_summary_count$n))
p <- ggplot(data = to_plot, aes(x = Var1, y = Freq)) +  
  geom_bar(stat = "identity") +  
  labs(  
    title = "H3K27me3 increased",  
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
    axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.1,size = 12),  
    axis.text.y = element_text(size = 12),  
    plot.title = element_text(size = 16, face = "bold")  
  ) +ylim(0,6500)
p
ggsave("result/Sup_figures/H3K27me3_increased_in_young_peaks_common_gene.pdf",p,width = 8,height = 6)
to_plot <- as.data.frame(table(decrease_summary_count$n))
p <- ggplot(data = to_plot, aes(x = Var1, y = Freq)) +  
  geom_bar(stat = "identity") +  
  labs(  
    title = "H3K27me3 decreased",  
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
    axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.1,size = 12),  
    axis.text.y = element_text(size = 12),  
    plot.title = element_text(size = 16, face = "bold")  
  ) +ylim(0,6500)
ggsave("result/Sup_figures/H3K27me3_decreased_in_young_peaks_common_gene.pdf",p,width = 8,height = 6)


txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
genelist <- increase_summary_count$Geneid[which(increase_summary_count$n >=5)]
# genelist <- decrease_summary_count$Geneid[which(decrease_summary_count$n >=20)]
genelist <- bitr(genelist,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelistGO <- enrichGO( genelist$ENTREZID,#GO富集分析
                        OrgDb = GO_database,
                        keyType = "ENTREZID",#设定读取的gene ID类型
                        ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                        pvalueCutoff = 0.05,#设定p值阈值
                        qvalueCutoff = 0.05,#设定q值阈值
                        readable = T)
p <- barplot(genelistGO,label_format = 50)
p
to_plot <- genelistGO@result
to_plot <- to_plot[order(to_plot$p.adjust),]
to_plot <- to_plot[c(1:10),]
to_plot$p.adjust <- -log10(to_plot$p.adjust)
to_plot$label <- paste0(to_plot$ID," ",to_plot$Description)
to_plot <- to_plot[,c("label","p.adjust")]
p <- ggplot(to_plot, aes(x = p.adjust, y = reorder(label, p.adjust))) +
  geom_bar(stat = "identity") +
  labs(x = "P Adjust", y = "Label") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10))
p

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'

genelist <- decrease_summary_count$Geneid[which(decrease_summary_count$n >=20)]
genelist <- bitr(genelist,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelistGO <- enrichGO( genelist$ENTREZID,#GO富集分析
                        OrgDb = GO_database,
                        keyType = "ENTREZID",#设定读取的gene ID类型
                        ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                        pvalueCutoff = 0.05,#设定p值阈值
                        qvalueCutoff = 0.05,#设定q值阈值
                        readable = T)
p <- barplot(genelistGO,label_format = 50)
p
to_plot <- genelistGO@result
to_plot <- to_plot[order(to_plot$p.adjust),]
to_plot <- to_plot[c(1:10),]
to_plot$p.adjust <- -log10(to_plot$p.adjust)
to_plot$label <- paste0(to_plot$ID," ",to_plot$Description)
to_plot <- to_plot[,c("label","p.adjust")]
p <- ggplot(to_plot, aes(x = p.adjust, y = reorder(label, p.adjust))) +
  geom_bar(stat = "identity") +
  labs(x = "P Adjust", y = "Label") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10))
p
ggsave("result/Sup_figures/H3K27me3_decreasd_in_H3K27me3_young_peaks_larger20_GO.pdf",p,width =15,height = 8)

