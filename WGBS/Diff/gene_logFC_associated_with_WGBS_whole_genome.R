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

tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
                  "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT","ileum")) 

tissue_summary <- read.csv("data/samples/WGBS/all_tissues_delta_in_200kb_bins_cross_comparison.csv",row.names = 1)
H3K9me3_tissue_order_label <- c() 
annotation_col <- data.frame()
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  t_search_table <- search_table[which(search_table$tissue==tissue),]
  young_samples <- t_search_table$sample_name[which(t_search_table$age=="3M")]
  old_samples <- t_search_table$sample_name[which(t_search_table$age=="24M")]
  combinations <- as.data.frame(expand.grid(young = young_samples, old = old_samples))
  if(tissue=="bonemarrow"){
    H3K9me3_tissue_order_label <- c(H3K9me3_tissue_order_label, paste0("Bone.Marrow",".",paste0(combinations$old,".",combinations$young)))  
    t_annotation_col <- data.frame(tissue=tissue_label_change(tissue),sample=paste0("Bone.Marrow",".",paste0(combinations$old,".",combinations$young)))  
  }else if(tissue=="mammarygland"){
    H3K9me3_tissue_order_label <- c(H3K9me3_tissue_order_label, paste0("Mammary.Gland",".",paste0(combinations$old,".",combinations$young)))  
    t_annotation_col <- data.frame(tissue=tissue_label_change(tissue),sample=paste0("Mammary.Gland",".",paste0(combinations$old,".",combinations$young)))  
  }
  else{
    H3K9me3_tissue_order_label <- c(H3K9me3_tissue_order_label, paste0(tissue_label_change(tissue),".",paste0(combinations$old,".",combinations$young)))  
    t_annotation_col <- data.frame(tissue=tissue_label_change(tissue),sample=paste0(tissue_label_change(tissue),".",paste0(combinations$old,".",combinations$young)))  
  }
  annotation_col <- rbind(annotation_col,t_annotation_col)
}
rownames(tissue_summary) <- tissue_summary$label

tissue_mean_summary <- data.frame() 
annotation <- annotation_col
for(tissue in tissues){
  t_annotation <- annotation[which(annotation$tissue==tissue_label_change(tissue)),]
  t_tissue_mean_summary <- tissue_summary[,t_annotation$sample]
  t_tissue_mean_summary$mean_delta <- rowMeans(t_tissue_mean_summary) 
  t_tissue_mean_summary$label <- rownames(t_tissue_mean_summary)
  t_tissue_mean_summary <- t_tissue_mean_summary[,c("label","mean_delta")]
  colnames(t_tissue_mean_summary)[2] <- tissue_label_change(tissue)
  if(nrow(tissue_mean_summary)==0){
    tissue_mean_summary <- t_tissue_mean_summary  
  }else{
    tissue_mean_summary <- merge(tissue_mean_summary,t_tissue_mean_summary,by="label")
  }
}

to_plot_WGBS_order_long <- reshape2::melt(tissue_mean_summary)

medians <- to_plot_WGBS_order_long %>%
  group_by(variable) %>%
  summarise(median_value = median(value, na.rm = TRUE))
medians <- medians[-which(medians$variable %in% c("Mammary Gland","Uterus","Ovary")),]
medians <- medians[order(medians$median_value),]
medians$rank <- c(1:nrow(medians))
medians <- as.data.frame(medians)

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

colnames(medians)[1:2] <- c("tissue","delta")
Indicator <- "delta"
summary_antibody <- medians[,c("tissue",Indicator)]
correlation_summary <- data.frame()
method <- "spearman"

colnames(medians)[1:2] <- c("tissue","delta")
Indicator <- "delta"
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
  t_correlation_summary <- data.frame(p_value=cortest$p.value[1], cor = as.numeric(cortest$estimate),gene=gene,histone="WGBS")
  correlation_summary <- rbind(correlation_summary,t_correlation_summary)
}

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
    ggtitle(paste0("Gene logFC positively correlated with the degree of downregulation of WGBS"))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab(Indicator)+labs(fill = "", color = "") +ylab(paste0("log2(Fold Change)"))
}else{
  ggplot(to_plot,aes(x=histone,y=value,color =variable))+    
    geom_jitter(size = 3, alpha = 0.7)+
    ggtitle(paste0("Gene logFC positively correlated with the degree of downregulation of WGBS"))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("Delta")+labs(fill = "", color = "") +ylab(paste0("log2(Fold Change)"))+
    scale_x_reverse()
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
    ggtitle(paste0("Gene logFC negatively correlated with the degree of downregulation of WGBS"))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab(Indicator)+labs(fill = "", color = "") +ylab(paste0("log2(Fold Change)"))
}else{
  ggplot(to_plot,aes(x=histone,y=value,color =variable))+    
    geom_jitter(size = 3, alpha = 0.7)+
    ggtitle(paste0("Gene logFC negatively correlated with the degree of downregulation of WGBS"))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("Delta")+labs(fill = "", color = "") +ylab(paste0("log2(Fold Change)"))+
    scale_x_reverse()
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
barplot(genelist_up_GO,title = paste0("Gene logFC positive correlaiton with DNA methylation change"),label_format = 50,showCategory = 30)

genelist_down <- bitr(negative_correlation$gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_down_GO <- enrichGO( genelist_down$ENTREZID,#GO富集分析
                              OrgDb = GO_database,
                              keyType = "ENTREZID",#设定读取的gene ID类型
                              ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                              pvalueCutoff = 0.05,#设定p值阈值
                              qvalueCutoff = 0.05,#设定q值阈值
                              readable = T)
result <- as.data.frame(genelist_down_GO@result)
barplot(genelist_down_GO,title = paste0("Gene expression negative correlaiton with DNA methylation change"),label_format = 70,showCategory = 30)

signif_correlation_genes <- correlation_summary
signif_correlation_genes <- signif_correlation_genes[order(signif_correlation_genes$cor),]
signif_correlation_genes_rank <- signif_correlation_genes$gene
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
emapplot(positive_gsea, showCategory = 50) 
negative_gsea <- pairwise_termsim(negative_gsea)  
emapplot(negative_gsea, showCategory = 50) 
