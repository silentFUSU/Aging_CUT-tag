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
genes <- vector()
second_quadrant_summary <- data.frame()
fourth_quadrant_summary <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_gene_TSS_diff_after_remove_batch_effect.csv"))
  df <- df[which(df$Significant != "Stable"),]
  RNA <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
  RNA <- RNA[which(RNA$Significant != "Stable"),]
  colnames(RNA)[1] <- "Geneid"
  df <- merge(df[,c("Geneid","LogFC.old.young")],RNA[,c("Geneid","logFC")],by="Geneid")
  second_quadrant <- df[which(df$LogFC.old.young > 0 & df$logFC < 0),]
  if(nrow(second_quadrant)!=0){
    second_quadrant$tissue <- tissue
    second_quadrant_summary <- rbind(second_quadrant_summary,second_quadrant)
  }
  
  fourth_quadrant <- df[which(df$LogFC.old.young < 0 & df$logFC > 0),]
  if(nrow(fourth_quadrant) != 0){
    fourth_quadrant$tissue <- tissue
    fourth_quadrant_summary <- rbind(fourth_quadrant_summary,fourth_quadrant)
  }
  genes <- c(genes,fourth_quadrant$Geneid)
}
genes <- unique(genes)

fourth_quadrant_summary_count <- fourth_quadrant_summary %>%   
  count(Geneid)
fourth_quadrant_summary_tissue <- fourth_quadrant_summary %>%   
  group_by(Geneid) %>%   
  summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
fourth_quadrant_summary_count <- merge(fourth_quadrant_summary_count,fourth_quadrant_summary_tissue,by="Geneid")
# genes <- fourth_quadrant_summary_count$Geneid[which(fourth_quadrant_summary_count$n >=4)]

summary <- data.frame()
for(tissue in tissues){
  RNA <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
  colnames(RNA)[1] <- "Geneid"
  RNA <- RNA[which(RNA$Geneid %in% genes),c("Geneid","logFC")]
  colnames(RNA)[2] <- tissue_label_change(tissue)
  if(nrow(summary)==0){
    summary <- RNA
  }else{
    summary <- merge(summary,RNA,by="Geneid",all=T)
  }
}
rownames(summary) <- summary$Geneid
summary[is.na(summary)] <- 0
summary <- summary[,-1]
set.seed(1)
k <- 4
kmeans_result <- kmeans(summary, centers=k)
summary$cluster <- kmeans_result$cluster
annotation <- summary[,"cluster",drop=F]
# write.csv(annotation,"data/samples/RNA/H3K27me3_decreased_RNA_increased_union_kmeans_change_filter_bar.csv")
# annotation <- read.csv("data/samples/RNA/H3K27me3_decreased_RNA_increased_union_kmeans_change_filter_bar.csv",row.names = 1)
colnames(annotation)[1] <- "cluster"
summary$cluster <- annotation$cluster  
summary_sorted <- summary[order(summary$cluster),]  
data_for_heatmap <- summary_sorted[, -ncol(summary_sorted)]  
 
annotation$cluster <- as.character(annotation$cluster)
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)  
breaks <- c(seq(-1, -0.29, length.out = 40), seq(-0.3, 0.3, length.out = 20), seq(0.31, 1, length.out = 40))  
# data_for_heatmap <- data_for_heatmap[,tissues_order]
data_for_heatmap <- summary_sorted[,-ncol(summary_sorted)]
pheatmap::pheatmap(data_for_heatmap,cluster_rows = F,cluster_cols = T,show_rownames = F,breaks = breaks, annotation_row = annotation, color = color_palette, clustering_distance_cols="manhattan")

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
genelist <- bitr(rownames(annotation)[which(annotation$cluster==4)],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelistGO <- enrichGO( genelist$ENTREZID,#GO富集分析
                            OrgDb = GO_database,
                            keyType = "ENTREZID",#设定读取的gene ID类型
                            ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                            pvalueCutoff = 0.05,#设定p值阈值
                            qvalueCutoff = 0.05,#设定q值阈值
                            readable = T)
barplot(genelistGO,label_format = 70)

# annotation <- read.csv("data/samples/RNA/H3K27me3_decreased_RNA_increased_union_kmeans.csv")
summary <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_gene_TSS_diff_after_remove_batch_effect.csv"))
  df <- df[which(df$Geneid %in% genes),c("Geneid","LogFC.old.young")]
  colnames(df)[2] <- tissue_label_change(tissue)
  if(nrow(summary)==0){
    summary <- df
  }else{
    summary <- merge(summary,df,by="Geneid",all=TRUE)
  }
}
rownames(summary) <- summary$Geneid
summary[is.na(summary)] <- 0
summary <- summary[,-1]
set.seed(1)
k <- 4
kmeans_result <- kmeans(summary, centers=k)
summary$cluster <- kmeans_result$cluster
annotation <- summary[,"cluster",drop=F]
# write.csv(annotation,"data/samples/RNA/H3K27me3_decreased_RNA_increased_H3K27me3_union_kmeans_tissue_lager4.csv")
annotation <- read.csv("data/samples/RNA/H3K27me3_decreased_RNA_increased_H3K27me3_union_kmeans.csv",row.names = 1)
colnames(annotation)[1] <- "cluster"
summary$cluster <- annotation$cluster  
summary_sorted <- summary[order(summary$cluster),]  
data_for_heatmap <- summary_sorted[, -ncol(summary_sorted)]  

annotation$cluster <- as.character(annotation$cluster)
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)  
breaks <- c(seq(-1, -0.29, length.out = 40), seq(-0.3, 0.3, length.out = 20), seq(0.31, 1, length.out = 40))  
# data_for_heatmap <- data_for_heatmap[,tissues_order]
data_for_heatmap <- summary_sorted[,-ncol(summary_sorted)]

color_palette <- colorRampPalette(c("blue", "white", "red"))(100)  
breaks <- c(seq(-1, -0.29, length.out = 40), seq(-0.3, 0.3, length.out = 20), seq(0.31, 1, length.out = 40))  

annotation$cluster <- as.character(annotation$cluster)
data_for_heatmap <- summary_sorted[, -ncol(summary_sorted)]  

tissues_order <- c("Ovary","Uterus","Thymus","Mammary Gland","Aorta","Pancreas","Skin","Ileum",
                   "Jejunum","Kidney","Liver","Muscle","BAT","IWAT","Lung","Cecum","Colon","Bone Marrow",
                   "Spleen","Heart","Testis","Cerebellum","Cortex","Hippocampus","Stomach","Bladder","Tongue")
data_for_heatmap <- data_for_heatmap[,tissues_order]
pheatmap::pheatmap(data_for_heatmap,cluster_rows = F,cluster_cols = T,show_rownames = F,breaks = breaks, annotation_row = annotation, color = color_palette, clustering_distance_cols="manhattan")

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
genelist <- bitr(rownames(annotation)[which(annotation$cluster==1)],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelistGO <- enrichGO( genelist$ENTREZID,#GO富集分析
                        OrgDb = GO_database,
                        keyType = "ENTREZID",#设定读取的gene ID类型
                        ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                        pvalueCutoff = 0.05,#设定p值阈值
                        qvalueCutoff = 0.05,#设定q值阈值
                        readable = T)
barplot(genelistGO,label_format = 70)

