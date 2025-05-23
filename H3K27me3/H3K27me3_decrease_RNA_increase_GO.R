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
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_gene_TSS_diff_after_remove_batch_effect.csv"))
  df <- df[which(df$Significant != "Stable"),]
  RNA <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
  RNA <- RNA[which(RNA$Significant != "Stable"),]
  colnames(RNA)[1] <- "Geneid"
  df <- merge(df[,c("Geneid","LogFC.old.young")],RNA[,c("Geneid","logFC")],by="Geneid")
  fourth_quadrant <- df[which(df$LogFC.old.young < 0 & df$logFC > 0),]
  second_quadrant <- df[which(df$LogFC.old.young > 0 & df$logFC < 0),] 
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
    write.csv(result,paste0("data/samples/RNA/",tissue,"/H3K27me3_decreased_gene_expression_increased_GO.csv"))
    p <- barplot(genelistGO,label_format = 70)+ ggtitle(tissue_label_change(tissue),"Genes with decreased H3K27me3 in the TSS region and increased self-expression")
    ggsave(paste0("result/RNA/",tissue,"/H3K27me3_decreased_in_TSS_gene_expression_increased_GO.png"),p,width=12,height=5,type="cairo")
  }

  # if(nrow(second_quadrant) >= 100){
  #   genelist <- bitr(second_quadrant$Geneid,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  #   genelistGO <- enrichGO( genelist$ENTREZID,#GO富集分析
  #                           OrgDb = GO_database,
  #                           keyType = "ENTREZID",#设定读取的gene ID类型
  #                           ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
  #                           pvalueCutoff = 0.05,#设定p值阈值
  #                           qvalueCutoff = 0.05,#设定q值阈值
  #                           readable = T)
  #   result <- genelistGO@result
  #   write.csv(result,paste0("data/samples/RNA/",tissue,"/H3K27me3_increased_gene_expression_decreased_GO.csv"))
  #   p <- barplot(genelistGO,label_format = 70)+ ggtitle(tissue_label_change(tissue),"Genes with increased H3K27me3 in the TSS region and decreased self-expression")
  #   ggsave(paste0("result/RNA/",tissue,"/H3K27me3_increased_in_TSS_gene_expression_decreased_GO.png"),p,width=12,height=5,type="cairo")
  # }
}

genelist_summary <- list()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_gene_TSS_diff_after_remove_batch_effect.csv"))
  df <- df[which(df$Significant != "Stable"),]
  RNA <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
  RNA <- RNA[which(RNA$Significant != "Stable"),]
  colnames(RNA)[1] <- "Geneid"
  df <- merge(df[,c("Geneid","LogFC.old.young")],RNA[,c("Geneid","logFC")],by="Geneid")
  fourth_quadrant <- df[which(df$LogFC.old.young < 0 & df$logFC > 0),]
  if(nrow(fourth_quadrant) >= 100){
    genelist <- bitr(fourth_quadrant$Geneid,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
    genelist_summary[[tissue_label_change(tissue)]] <- genelist$ENTREZID
  }
}
genelistGO  <- compareCluster(geneCluster = genelist_summary, fun = enrichGO, OrgDb = GO_database, keyType = "ENTREZID",ont = "BP")
result <- genelistGO@compareClusterResult
cluster_result <- reshape2::dcast(result, Cluster ~ID, value.var = "p.adjust")  
cluster_result[is.na(cluster_result)] <- 1
rownames(cluster_result) <- cluster_result$Cluster
cluster_result <- cluster_result[,-1]
distance_matrix <- dist(cluster_result, method = "euclidean")  
hc <- hclust(distance_matrix, method = "ward.D2")  
plot(hc)
plot(hc,labels = FALSE)  

plot.phylo(phylo_hc, type = "cladogram", edge.color = "black")  
plot(hc, hang = -1, labels = rownames(cluster_result))  
order_from_left_to_right <- hc$order  
sorted_data <- cluster_result[order_from_left_to_right, ]  
genelistGO@compareClusterResult$Cluster<- factor(genelistGO@compareClusterResult$Cluster, levels = rownames(sorted_data))  

p <- dotplot(genelistGO,showCategory = 8,label_format = 100)+theme(axis.text.x = element_text(angle = 90, hjust = 1,size=15),axis.text.y = element_text(size=20)) + coord_flip() + xlab(NULL)  
ggsave("result/RNA/GO/plot/all_tissues_H3K27me3_decreased_RNA_increased_GO_pathway.png",p,width = 40,height = 20)


GO_summary <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_gene_TSS_diff_after_remove_batch_effect.csv"))
  df <- df[which(df$Significant != "Stable"),]
  RNA <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
  RNA <- RNA[which(RNA$Significant != "Stable"),]
  colnames(RNA)[1] <- "Geneid"
  df <- merge(df[,c("Geneid","LogFC.old.young")],RNA[,c("Geneid","logFC")],by="Geneid")
  fourth_quadrant <- df[which(df$LogFC.old.young < 0 & df$logFC > 0),]
  if(nrow(fourth_quadrant) >=100){
    df <- read.csv(paste0("data/samples/RNA/",tissue,"/H3K27me3_decreased_gene_expression_increased_GO.csv"))
    df <- df[which(df$p.adjust < 0.05),c(2,3)]
    df$tissue <- tissue_label_change(tissue)
    GO_summary <- rbind(GO_summary,df)
  }
}
GO_summary_count <- GO_summary %>%   
  count(ID)
GO_summary_tissue <- GO_summary %>%   
  group_by(ID) %>%   
  summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
GO_summary_count <- merge(GO_summary_count,GO_summary_tissue,by="ID")
GO_summary_count$Description <- sapply(GO_summary_count$ID, function(go_id) {  
  Term(GOTERM[[go_id]])  
})  
to_plot <- as.data.frame(table(GO_summary_count$n))
ggplot(data = to_plot, aes(x = Var1, y = Freq)) +  
  geom_bar(stat = "identity") +  
  labs(  
    title = "Distribution of tissues number in common changed GO terms",  
    x = NULL,  
    y = "Count"  
  ) +  
  geom_text(  
    aes(label = Freq),   
    vjust = 0,   
    size = 3.5  # Adjust text size as needed  
  ) +
  theme_minimal() +  
  xlab(NULL) +  
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
    axis.text.y = element_text(size = 12),  
    plot.title = element_text(size = 16, face = "bold")  
  ) 



GO_term <- GO_summary_count$ID[which(GO_summary_count$n>=6)]

GO_summary <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_gene_TSS_diff_after_remove_batch_effect.csv"))
  df <- df[which(df$Significant != "Stable"),]
  RNA <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
  RNA <- RNA[which(RNA$Significant != "Stable"),]
  colnames(RNA)[1] <- "Geneid"
  df <- merge(df[,c("Geneid","LogFC.old.young")],RNA[,c("Geneid","logFC")],by="Geneid")
  fourth_quadrant <- df[which(df$LogFC.old.young < 0 & df$logFC > 0),]
  if(nrow(fourth_quadrant) >=100){
    df <- read.csv(paste0("data/samples/RNA/",tissue,"/H3K27me3_decreased_gene_expression_increased_GO.csv"))
    df <- df[which(df$ID %in% GO_term),c(2,3,7)]
    df$p.adjust <- -log10(df$p.adjust)
    colnames(df)[3] <- tissue_label_change(tissue)
    if(nrow(GO_summary)==0){
      GO_summary <- df
    }else{
      GO_summary <- merge(GO_summary,df[,c(1,3)],by="ID",all=T)
    }
  }
}
GO_summary$Description <- sapply(GO_summary$ID, function(go_id) {  
  Term(GOTERM[[go_id]])  
})  
rownames(GO_summary) <- GO_summary$Description
GO_summary[is.na(GO_summary)] <- 0  
GO_summary <- GO_summary[,-c(1:2)]
set.seed(1)
k <- 4
kmeans_result <- kmeans(GO_summary, centers=k)
GO_summary$cluster<- kmeans_result$cluster
GO_summary_sorted <- GO_summary[order(GO_summary$cluster),]  
annotation <- GO_summary[,"cluster",drop=F]
annotation$cluster <- as.character(annotation$cluster)
color_palette <- colorRampPalette(c("white", "red"))(100)  
breaks <- c(seq(-log10(0.05),10, length.out = 100))  
data_for_heatmap <- GO_summary_sorted[, -ncol(GO_summary_sorted)]  
pheatmap::pheatmap(data_for_heatmap,cluster_rows = F,show_rownames = T,breaks = breaks, annotation_row = annotation, color = color_palette, clustering_distance_cols="manhattan",filename = "result/RNA/histone_relationship_with_RNA/H3K27me3/H3K27me3_decreased_RNA_inreased_GO_heatmap.png",width=15,height=40)
data_for_heatmap2 <- data_for_heatmap[which(rownames(data_for_heatmap) %in% rownames(annotation)[which(annotation$cluster %in% c(2,3))]),]
pheatmap::pheatmap(data_for_heatmap2,cluster_rows = F,show_rownames = T,breaks = breaks, annotation_row = annotation, color = color_palette, clustering_distance_cols="manhattan",filename = "result/RNA/histone_relationship_with_RNA/H3K27me3/H3K27me3_decreased_RNA_inreased_GO_heatmap_subclusters.png",width=15,height=20)

GO_summary <- vector()
top <- 8
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_gene_TSS_diff_after_remove_batch_effect.csv"))
  df <- df[which(df$Significant != "Stable"),]
  RNA <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
  RNA <- RNA[which(RNA$Significant != "Stable"),]
  colnames(RNA)[1] <- "Geneid"
  df <- merge(df[,c("Geneid","LogFC.old.young")],RNA[,c("Geneid","logFC")],by="Geneid")
  fourth_quadrant <- df[which(df$LogFC.old.young < 0 & df$logFC > 0),]
  if(nrow(fourth_quadrant) >=100){
    df <- read.csv(paste0("data/samples/RNA/",tissue,"/H3K27me3_decreased_gene_expression_increased_GO.csv"))
    df <- df[order(df$p.adjust),]
    GO_summary <- c(GO_summary,df$ID[1:top])
  }
}
GO_terms <- unique(GO_summary)
GO_summary <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_gene_TSS_diff_after_remove_batch_effect.csv"))
  df <- df[which(df$Significant != "Stable"),]
  RNA <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
  RNA <- RNA[which(RNA$Significant != "Stable"),]
  colnames(RNA)[1] <- "Geneid"
  df <- merge(df[,c("Geneid","LogFC.old.young")],RNA[,c("Geneid","logFC")],by="Geneid")
  fourth_quadrant <- df[which(df$LogFC.old.young < 0 & df$logFC > 0),]
  if(nrow(fourth_quadrant) >=100){
    df <- read.csv(paste0("data/samples/RNA/",tissue,"/H3K27me3_decreased_gene_expression_increased_GO.csv"))
    df <- df[which(df$ID %in% GO_terms),c("ID","p.adjust")]
    colnames(df)[2] <- tissue_label_change(tissue)
    if(nrow(GO_summary)==0){
      GO_summary <- df
    }else{
      GO_summary <- merge(GO_summary,df,by="ID",all=T)
    }
  }
}
GO_summary[is.na(GO_summary)] <- 1
GO_summary$Description <- sapply(GO_summary$ID, function(go_id) {  
  Term(GOTERM[[go_id]])  
})  
rownames(GO_summary) <- GO_summary$Description

GO_summary <- GO_summary[,-c(1,ncol(GO_summary))]
GO_summary[] <- lapply(GO_summary, function(x) -log10(x))  
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)  
breaks <- c(seq(0, -log10(0.05), length.out = 40), seq(-log10(0.05)+0.1, 10, length.out = 60)) 
pheatmap::pheatmap(t(GO_summary),breaks = breaks,color = color_palette)
