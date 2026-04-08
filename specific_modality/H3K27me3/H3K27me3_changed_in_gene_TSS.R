rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(DOSE)
library(clusterProfiler)
options(scipen = 0)
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
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'

for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_gene_TSS_diff_after_remove_batch_effect.csv"))
  
  peaks <- read.table(paste0("data/samples/",tissue,"/H3K27me3/bed/H3K27me3_young_old_merge-W5000-G10000-E100.bed"))
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
  
  df <- df[which(df$Significant != "Stable"),c("Geneid","LogFC.old.young")]
  fourth_quadrant <- df[which(df$LogFC.old.young < 0),]
  second_quadrant <- df[which(df$LogFC.old.young > 0),] 
  if(nrow(fourth_quadrant) >= 100){
    genelist <- bitr(fourth_quadrant$Geneid,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
    genelistGO <- enrichGO( genelist$ENTREZID,#GO富集分析
                            OrgDb = GO_database,
                            keyType = "ENTREZID",#设定读取的gene ID类型
                            ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                            pvalueCutoff = 0.05,#设定p值阈值
                            qvalueCutoff = 0.05,#设定q值阈值
                            readable = T)
    result <- genelistGO@result
    write.csv(result,paste0("data/samples/RNA/",tissue,"/H3K27me3_decreased_in_gene_TSS_GO.csv"))
    p <- barplot(genelistGO,label_format = 70)+ ggtitle(tissue_label_change(tissue),"Genes with decreased H3K27me3 in the TSS region")
    ggsave(paste0("result/RNA/",tissue,"/H3K27me3_decreased_in_gene_TSS_GO.png"),p,width=12,height=5,type="cairo")
  }

  if(nrow(second_quadrant) >= 100){
    genelist <- bitr(second_quadrant$Geneid,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
    genelistGO <- enrichGO( genelist$ENTREZID,#GO富集分析
                            OrgDb = GO_database,
                            keyType = "ENTREZID",#设定读取的gene ID类型
                            ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                            pvalueCutoff = 0.05,#设定p值阈值
                            qvalueCutoff = 0.05,#设定q值阈值
                            readable = T)
    result <- genelistGO@result
    write.csv(result,paste0("data/samples/RNA/",tissue,"/H3K27me3_increased_in_gene_TSS_GO.csv"))
    p <- barplot(genelistGO,label_format = 70)+ ggtitle(tissue_label_change(tissue),"Genes with increased H3K27me3 in the TSS region")
    ggsave(paste0("result/RNA/",tissue,"/H3K27me3_increased_in_gene_TSS_GO.png"),p,width=12,height=5,type="cairo")
  }
}

increase_summary <- data.frame()
decrease_summary <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_gene_TSS_10kb_diff_after_remove_batch_effect.csv"))
  
  peaks <- read.table(paste0("data/samples/",tissue,"/H3K27me3/bed/H3K27me3_young_old_merge-W5000-G10000-E100.bed"))
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
  
  df <- df[which(df$Significant != "Stable"),c("Geneid","LogFC.old.young")]
  decrease <- df[which(df$LogFC.old.young < 0),]
  decrease$tissue <- tissue_label_change(tissue)
  increase <- df[which(df$LogFC.old.young > 0),] 
  increase$tissue <- tissue_label_change(tissue)
  
  increase_summary <- rbind(increase_summary,increase)
  decrease_summary <- rbind(decrease_summary,decrease)
  }
increase_summary_count <- increase_summary %>%   
  count(Geneid)
increase_summary_tissue <- increase_summary %>%   
  group_by(Geneid) %>%   
  summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
increase_summary_count <- merge(increase_summary_count,increase_summary_tissue,by="Geneid")
to_plot <- as.data.frame(table(increase_summary_count$n))
ggplot(data = to_plot, aes(x = Var1, y = Freq)) +  
  geom_bar(stat = "identity") +  
  labs(  
    title = "Distribution of tissues number in H3K27me3 increased gene",  
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

decrease_summary_count <- decrease_summary %>%   
  count(Geneid)
decrease_summary_tissue <- decrease_summary %>%   
  group_by(Geneid) %>%   
  summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
decrease_summary_count <- merge(decrease_summary_count,decrease_summary_tissue,by="Geneid")
to_plot <- as.data.frame(table(decrease_summary_count$n))
ggplot(data = to_plot, aes(x = Var1, y = Freq)) +  
  geom_bar(stat = "identity") +  
  labs(  
    title = "Distribution of tissues number in H3K27me3 decreased gene",  
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
genelist <- increase_summary_count$Geneid[which(increase_summary_count$n >=15)]
genelist <- bitr(genelist,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelistGO <- enrichGO( genelist$ENTREZID,#GO富集分析
                        OrgDb = GO_database,
                        keyType = "ENTREZID",#设定读取的gene ID类型
                        ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                        pvalueCutoff = 0.05,#设定p值阈值
                        qvalueCutoff = 0.05,#设定q值阈值
                        readable = T)
barplot(genelistGO,label_format = 70,font.size = 14)

genelist <- decrease_summary_count$Geneid[which(decrease_summary_count$n >=20)]
genelist <- bitr(genelist,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelistGO <- enrichGO( genelist$ENTREZID,#GO富集分析
                        OrgDb = GO_database,
                        keyType = "ENTREZID",#设定读取的gene ID类型
                        ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                        pvalueCutoff = 0.05,#设定p值阈值
                        qvalueCutoff = 0.05,#设定q值阈值
                        readable = T)
barplot(genelistGO,label_format = 70,font.size = 14)

summary <- data.frame()
gene <- unique(c(increase_summary_count$Geneid,decrease_summary_count$Geneid))
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_gene_TSS_diff_after_remove_batch_effect.csv"))
  df <- df[which(df$Geneid %in% gene),c("Geneid","LogFC.old.young")]
  colnames(df)[2] <- tissue_label_change(tissue)
  if(nrow(summary) == 0){
    summary <- df
  }else{
    summary <- merge(summary,df,by="Geneid",all=T)
  }
}
summary[is.na(summary)] <- 0
rownames(summary) <- summary$Geneid
summary <- summary[,-1]
k <- 10
set.seed(1)
kmeans_result <- kmeans(summary, centers=k)
summary$cluster <- kmeans_result$cluster
annotation <- summary[,"cluster",drop=F]
annotation$cluster <- factor(annotation$cluster,levels=c(1:k))
write.csv(annotation,"data/samples/all/H3K27me3/H3K27me3_changed_in_gene_TSS_all_tissues_kmeans.csv")
summary$cluster <- factor(summary$cluster,levels=c(1:k))
summary_sorted <- summary[order(summary$cluster),]  
data_for_heatmap <- summary_sorted[, -ncol(summary_sorted)]  
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)  
breaks <- c(seq(-1, -0.29, length.out = 40), seq(-0.3, 0.3, length.out = 20), seq(0.31, 1, length.out = 40))  
pheatmap::pheatmap(data_for_heatmap,cluster_rows = F,cluster_cols = T,show_rownames = F,breaks = breaks, annotation_row = annotation, color = color_palette)

genelist <- list()
for(kmean in c(1:k)){
  df <- bitr(rownames(annotation)[which(annotation$cluster==kmean)],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  genelist[[paste0("kmeans",kmean)]] <- df$ENTREZID
  genelistGO <- enrichGO( df$ENTREZID,#GO富集分析
                          OrgDb = GO_database,
                          keyType = "ENTREZID",#设定读取的gene ID类型
                          ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                          pvalueCutoff = 0.05,#设定p值阈值
                          qvalueCutoff = 0.05,#设定q值阈值
                          readable = T)
  p <- barplot(genelistGO,label_format = 70,font.size = 7)
  ggsave(paste0("result/RNA/histone_relationship_with_RNA/H3K27me3/H3K27me3_changed_in_Gene_TSS_kmeans",kmean,"_GO.png"),p,width = 8,height = 4,type= "cairo")
}
genelistGO <- compareCluster(genelist,
                             fun="enrichGO",
                             OrgDb = GO_database,
                             ont = "BP",
                             keyType = "ENTREZID")
p <- dotplot(genelistGO,label_format = 100,showCategory = 10)+theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
ggsave("result/RNA/histone_relationship_with_RNA/H3K27me3/H3K27me3_changed_in_Gene_TSS_kmeans_GO.png",p,width = 12,height = 15,type= "cairo")
