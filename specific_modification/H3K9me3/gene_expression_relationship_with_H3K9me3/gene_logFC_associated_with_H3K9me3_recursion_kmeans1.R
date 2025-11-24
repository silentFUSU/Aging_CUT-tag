rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(stringr)
library(dplyr)
library(dbplyr)
library(clusterProfiler)
library(GSVA)
library(enrichplot)
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
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
}
antibody <- "H3K9me3"
regions <- read.table("data/samples/all/H3K9me3/recursion_peaks_diff_table/bed/kmeans1_uinon_recursion_peaks.bed")
annotation <- data.frame(label = paste0(regions$V1,":",regions$V2,"-",regions$V3),chr=regions$V1)
rownames(annotation) <- annotation$label
annotation <- annotation[,-1,drop=F]
annotation$chr <- factor(annotation$chr,levels=paste0("chr",c(1:19,"X","Y")))

regions <- paste0(regions$V1,":",regions$V2,"-",regions$V3)
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver","lung","muscle","pancreas","skin","spleen","stomach","testis","thymus","tongue","iWAT")
diff_summary <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_young_old_merge-W5000-G10000-E100_recursion_diff_after_remove_batch_effect.csv"))
  df <- df[which(df$Geneid %in% regions),c("Geneid","LogFC.old.young","Significant")]
  df <- df[,-3]
  colnames(df)[2] <- tissue_label_change(tissue)
  if(nrow(diff_summary) == 0){
    diff_summary <- df
  }else{
    diff_summary <- merge(diff_summary,df,by="Geneid",all=T)    
  }
}
rownames(diff_summary) <- diff_summary$Geneid
diff_summary <- diff_summary[,-1]
diff_summary$row_mean <- rowMeans(diff_summary, na.rm = TRUE)
col_means <- colMeans(as.matrix(diff_summary[, -ncol(diff_summary)]), na.rm = TRUE)
diff_summary_sorted_rows <- diff_summary %>%
  arrange(row_mean)
diff_summary_sorted_rows <- diff_summary_sorted_rows[ , -ncol(diff_summary_sorted_rows)]
diff_summary_sorted <- diff_summary_sorted_rows[, order(col_means)]
diff_summary <- diff_summary_sorted
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- c(seq(-1, -0.29, length.out = 40), seq(-0.3, 0.3, length.out = 20), seq(0.31, 1, length.out = 40))
color <- read.table("data/samples/30_distinct_color.txt")
annotation_color <- list(chr=setNames(color$V1[1:21],paste0("chr",c(1:19,"X","Y"))))
p <- pheatmap::pheatmap(diff_summary,show_rownames = F,cluster_rows = F,cluster_cols = F,breaks = breaks, color = color_palette, clustering_distance_cols="manhattan",clustering_distance_rows="manhattan",annotation_row = annotation,annotation_colors = annotation_color)

medians <- apply(diff_summary, 2, median, na.rm = TRUE)
medians <- data.frame(tissues=colnames(diff_summary),median=medians)
medians <- medians[order(medians$median),]
medians$rank <- 1:nrow(medians)
colnames(medians) <- c("tissue","logFC","rank")

logFC <- data.frame() 
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"),header = T)
  df <- df[,c("X","logFC")]
  colnames(df) <- c("Geneid","logFC")
  colnames(df)[2] <- tissue_label_change(tissue)
  if(nrow(logFC)==0){
    logFC <- df
  }else{
    logFC <- merge(logFC,df,by="Geneid")
  }
}
rownames(logFC) <- logFC$Geneid
logFC <- logFC[,-1]

Indicator <- "logFC"
summary_antibody <- medians[,c("tissue",Indicator)]
correlation_summary <- data.frame()
method <- "spearman"

for(i in c(1:nrow(logFC))){
  gene <- rownames(logFC)[i]
  t_logFC <- as.data.frame(t(logFC[i,]))
  colnames(t_logFC)[1] <- "logFC"
  t_logFC <- as.data.frame(t_logFC)
  t_logFC <- merge(t_logFC,summary_antibody,by.x="row.names",by.y="tissue")
  colnames(t_logFC)[3] <- "histone"
  cortest <- cor.test(t_logFC$logFC,t_logFC$histone,method = method)
  t_correlation_summary <- data.frame(p_value=cortest$p.value[1], cor = as.numeric(cortest$estimate),gene=gene,histone=antibody)
  correlation_summary <- rbind(correlation_summary,t_correlation_summary)
}
correlation_summary_output <- correlation_summary 
correlation_summary_output$cor <- -correlation_summary_output$cor
correlation_summary_output <- correlation_summary_output[order(correlation_summary_output$cor,decreasing = T),]
# write.csv(correlation_summary_output,"data/samples/all/H3K9me3/recursion_peaks_diff_table/gene_expression_logFC_associated_with_kmeans1.csv",row.names = F)

positive_correlation <- correlation_summary[which(correlation_summary$p_value < 0.05 & correlation_summary$cor < 0),]
negative_correlation <- correlation_summary[which(correlation_summary$p_value < 0.05 & correlation_summary$cor > 0),]

positive_correlation <- positive_correlation[order(positive_correlation$cor),]
positive_correlation_genes <- positive_correlation$gene[1:min(20,nrow(positive_correlation))]
negative_correlation <- negative_correlation[order(negative_correlation$cor,decreasing = T),]
negative_correlation_genes <- negative_correlation$gene[1:min(20,nrow(negative_correlation))]

to_plot <- data.frame()
for(i in c(1:length(positive_correlation_genes))){
  gene <- positive_correlation_genes[i]
  t_logFC<- as.data.frame(t(logFC[gene,]))
  colnames(t_logFC)[1] <- "logFC"
  colnames(t_logFC)[1] <- gene
  t_logFC$tissue <- rownames(t_logFC)
  if(nrow(to_plot)==0){
    to_plot <- t_logFC
  } else{
    to_plot <- merge(to_plot,t_logFC,by="tissue")
  }
}

to_plot <- reshape2::melt(to_plot)
to_plot <- merge(to_plot,summary_antibody,by="tissue")
colnames(to_plot)[4] <- "histone" 
if(Indicator=="rank"){
  to_plot$histone <- factor(to_plot$histone,c(27:1))
  ggplot(to_plot,aes(x=histone,y=value,color =variable))+    
    geom_jitter(size = 3, alpha = 0.7)+
    ggtitle(paste0("Gene logFC positively correlated with the degree of downregulation of ",antibody))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab(Indicator)+labs(fill = "", color = "") +ylab(paste0("Gene expression log2(Fold Change)"))
}else{
  ggplot(to_plot,aes(x=histone,y=value,color =variable))+    
    geom_jitter(size = 3, alpha = 0.7)+
    ggtitle(paste0("Gene logFC positively correlated with the degree of downregulation of ",antibody))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("H3K9me3 log2(Fold change)")+labs(fill = "", color = "") +ylab(paste0("Gene expression log2(Fold Change)"))+ scale_x_reverse()
}


to_plot <- data.frame()
for(i in c(1:length(negative_correlation_genes))){
  gene <- negative_correlation_genes[i]
  t_logFC<- as.data.frame(t(logFC[gene,]))
  colnames(t_logFC)[1] <- "logFC"
  colnames(t_logFC)[1] <- gene
  t_logFC$tissue <- rownames(t_logFC)
  if(nrow(to_plot)==0){
    to_plot <- t_logFC
  } else{
    to_plot <- merge(to_plot,t_logFC,by="tissue")
  }
}

to_plot <- reshape2::melt(to_plot)
to_plot <- merge(to_plot,summary_antibody,by="tissue")
colnames(to_plot)[4] <- "histone" 
if(Indicator=="rank"){
  to_plot$histone <- factor(to_plot$histone,c(27:1))
  ggplot(to_plot,aes(x=histone,y=value,color =variable))+    
    geom_jitter(size = 3, alpha = 0.7)+
    ggtitle(paste0("Gene logFC negatively correlated with the degree of downregulation of ",antibody))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab(Indicator)+labs(fill = "", color = "") +ylab(paste0("log2(Fold Change)"))
}else{
  ggplot(to_plot,aes(x=histone,y=value,color =variable))+    
    geom_jitter(size = 3, alpha = 0.7)+
    ggtitle(paste0("Gene logFC negatively correlated with the degree of downregulation of ",antibody))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("H3K9me3 log2(Fold change)")+labs(fill = "", color = "") +ylab(paste0("Gene expression log2(Fold Change)"))+ scale_x_reverse()
}


GO_database <- 'org.Mm.eg.db'
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
genelist_up <- bitr(positive_correlation$gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_up_GO <- enrichGO( genelist_up$ENTREZID,#GO富集分析
                            OrgDb = GO_database,
                            keyType = "ENTREZID",#设定读取的gene ID类型
                            ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                            pvalueCutoff = 0.05,#设定p值阈值
                            qvalueCutoff = 0.05,#设定q值阈值
                            readable = T)
result <- as.data.frame(genelist_up_GO@result)
barplot(genelist_up_GO,title = paste0("Gene logFC positive correlaiton with ",antibody," change"),label_format = 50,showCategory = 30)

genelist_down <- bitr(negative_correlation$gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_down_GO <- enrichGO( genelist_down$ENTREZID,#GO富集分析
                              OrgDb = GO_database,
                              keyType = "ENTREZID",#设定读取的gene ID类型
                              ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                              pvalueCutoff = 0.05,#设定p值阈值
                              qvalueCutoff = 0.05,#设定q值阈值
                              readable = T)
result <- as.data.frame(genelist_down_GO@result)
barplot(genelist_down_GO,title = paste0("Gene expression negative correlaiton with ",antibody," change"),label_format = 70,showCategory = 30)


signif_correlation_genes <- correlation_summary
signif_correlation_genes <- signif_correlation_genes[order(signif_correlation_genes$cor),]
signif_correlation_genes_rank <- signif_correlation_genes$gene
# genelist <- bitr(signif_correlation_genes_rank,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
# genelist <- merge(genelist,signif_correlation_genes,by.x="SYMBOL",by.y="gene")
genelist <- signif_correlation_genes[,c("gene","cor")]
genelist$cor <- -genelist$cor
genelist <- setNames(genelist$cor, genelist$gene)
genelist <- sort(genelist, decreasing = TRUE)

gse <- gseGO(geneList=genelist, 
             ont = "BP",
             keyType = "SYMBOL", 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = GO_database,
             pAdjustMethod = "none",eps = 1e-100)
results_df <- as.data.frame(gse@result)
write.csv(results_df,"data/samples/all/H3K9me3/recursion_peaks_diff_table/gene_expression_logFC_associated_with_kmeans1_GO_terms.csv",row.names = F)

positive_nes_results <- subset(results_df, NES > 0)
negative_nes_results <- subset(results_df, NES < 0)
positive_gsea <- new("gseaResult",
                     result = positive_nes_results,
                     geneSets = gse@geneSets,
                     geneList = gse@geneList,
                     params = gse@params,
                     setType = gse@setType,
                     organism = gse@organism)
negative_gsea <- new("gseaResult",
                     result = negative_nes_results,
                     geneSets = gse@geneSets,
                     geneList = gse@geneList,
                     params = gse@params,
                     setType = gse@setType,
                     organism = gse@organism)

positive_gsea <- pairwise_termsim(positive_gsea)  
p <- emapplot(positive_gsea, showCategory = 50) 
ggsave(paste0("result/Sup_figures/H3K9me3_kmeans1_gene_expression_logFC_correlation_GO_pathway.pdf"),p,width = 12,height = 10)
negative_gsea <- pairwise_termsim(negative_gsea)  
emapplot(negative_gsea, showCategory = 50) 

GO_id <- "GO:0050778"
p <-gseaplot2( 
  gse,
  geneSetID = GO_id,
  title = paste0("Enrichment Plot for ",gse@result$Description[which(gse@result$ID==GO_id)]),
  pvalue_table = TRUE
)
ggsave("result/figures/H3K9me3_Immune_logFC_kmeans1.pdf",p,width = 8,height = 5)
