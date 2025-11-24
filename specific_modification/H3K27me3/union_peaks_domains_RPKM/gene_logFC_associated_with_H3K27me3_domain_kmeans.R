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

antibody <- "H3K27me3"
annotation <- read.csv("data/samples/all/H3K27me3/edd_domain_merged/kmeans_annotation_RPKM.csv")
regions <- annotation[which(annotation$cluster=="3"),]
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver","lung","pancreas","skin","spleen","stomach","testis","thymus","tongue","iWAT","muscle")
diff_summary <- data.frame()
rpkm_log2FC_summary <- data.frame()
condition <- "domain"
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  if(condition=="domain"){
    tab <- read.table(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_edd_domain_merged.counts"),header = T)
    summary <- read.table(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_edd_domain_merged.counts.summary"),header = T,row.names = 1)
  }else{
    tab <- read.table(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_young_merge-W5000-G10000-E100_peak_merged_exist_in_larger_20_tissues.counts"),header = T)
    summary <- read.table(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_young_merge-W5000-G10000-E100_peak_merged_exist_in_larger_20_tissues.counts.summary"),header = T,row.names = 1)
  }
  
  if(tissue %in% c("mammarygland","uterus","ovary")){
    tab <- tab[which(tab$Chr %in% paste0("chr",c(1:19,"X"))),]
  }
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+|DYQ[0-9]+).*"
  colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
  colnames(summary) <- gsub(pattern,"\\1",colnames(summary))
  search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
  counts <- counts[,search_table$sample_name]
  summary <- summary[-2,search_table$sample_name]
  total_reads <- colSums(summary)
  length <- as.numeric(tab$Length)
  
  rpkm <- sweep(counts,2,total_reads,"/")
  rpkm <- sweep(rpkm,1,length,"/") * 1000000000
  
  rpkm_young <- rpkm[,search_table$sample_name[which(search_table$age=="3m")]]
  rpkm_old <- rpkm[,search_table$sample_name[which(search_table$age=="24m")]]
  
  rpkm_young$mean_young <- rowMeans(rpkm_young)
  rpkm_old$mean_old <- rowMeans(rpkm_old)
  
  rpkm_mean_summary <- merge(rpkm_young[,"mean_young",drop=F],rpkm_old[,"mean_old",drop=F],by="row.names")
  rpkm_mean_summary$log2FC <- log2(rpkm_mean_summary$mean_old/rpkm_mean_summary$mean_young)
  colnames(rpkm_mean_summary)[which(colnames(rpkm_mean_summary)=="log2FC")] <- tissue_label_change(tissue)
  colnames(rpkm_mean_summary)[1] <- "Geneid"
  rpkm_mean_summary <- rpkm_mean_summary[which(rpkm_mean_summary$Geneid %in% regions$X),]
  if(nrow(rpkm_log2FC_summary )==0){
    rpkm_log2FC_summary <- rpkm_mean_summary[,c(1,4)]
  }else{
    rpkm_log2FC_summary <- merge(rpkm_log2FC_summary,rpkm_mean_summary[,c(1,4)],by="Geneid",all=T)
  }
}
diff_summary <- rpkm_log2FC_summary
rownames(diff_summary) <- diff_summary$Geneid
diff_summary <- diff_summary[,-1]
to_plot <- diff_summary
to_plot[is.na(to_plot)] <- 0
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- c(seq(-5, -0.51, length.out = 40), seq(-0.5, 0.5, length.out = 20), seq(0.51, 5, length.out = 40))
color <- read.table("data/samples/30_distinct_color.txt")
annotation_color <- list(chr=setNames(color$V1[1:21],paste0("chr",c(1:19,"X","Y"))))
p <- pheatmap::pheatmap(to_plot,show_rownames = F,breaks = breaks, color = color_palette, clustering_distance_cols="manhattan",clustering_distance_rows="manhattan")

medians <- apply(diff_summary, 2, median, na.rm = TRUE)
medians <- data.frame(tissues=colnames(diff_summary),median=medians)
medians <- medians[order(medians$median,decreasing = T),]
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
positive_correlation <- correlation_summary[which(correlation_summary$p_value < 0.05 & correlation_summary$cor > 0),]
negative_correlation <- correlation_summary[which(correlation_summary$p_value < 0.05 & correlation_summary$cor < 0),]
positive_correlation <- positive_correlation[order(positive_correlation$cor,decreasing = T),]
positive_correlation_genes <- positive_correlation$gene[1:min(20,nrow(positive_correlation))]
negative_correlation <- negative_correlation[order(negative_correlation$cor),]
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
    ggtitle(paste0("Gene expression relationship with ",antibody," change"))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab(Indicator)+labs(fill = "", color = "") +ylab(paste0("log2(CPM)"))
}else{
  ggplot(to_plot,aes(x=histone,y=value,color =variable))+    
    geom_jitter(size = 3, alpha = 0.7)+
    ggtitle(paste0("Gene expression changed positively correlated with the degree of upregulation of H3K27me3"))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab(Indicator)+labs(fill = "", color = "") +ylab(paste0("log2(Fold change)"))
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
    ggtitle(paste0("Gene expression relationship with ",antibody," change"))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab(Indicator)+labs(fill = "", color = "") +ylab(paste0("log2(CPM)"))
}else{
  ggplot(to_plot,aes(x=histone,y=value,color =variable))+    
    geom_jitter(size = 3, alpha = 0.7)+
    ggtitle(paste0("Gene expression changed positively correlated with the degree of upregulation of H3K27me3"))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab(Indicator)+labs(fill = "", color = "") +ylab(paste0("log2(Fold change)"))
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
barplot(genelist_up_GO,title = paste0("Gene expression change positive correlaiton with ",antibody," change"),label_format = 50,showCategory = 30)

genelist_down <- bitr(negative_correlation$gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_down_GO <- enrichGO( genelist_down$ENTREZID,#GO富集分析
                              OrgDb = GO_database,
                              keyType = "ENTREZID",#设定读取的gene ID类型
                              ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                              pvalueCutoff = 0.05,#设定p值阈值
                              qvalueCutoff = 0.05,#设定q值阈值
                              readable = T)
result <- as.data.frame(genelist_down_GO@result)
barplot(genelist_down_GO,title = paste0("Gene expression change negative correlaiton with ",antibody," change"),label_format = 70,showCategory = 30)

###GSEA
signif_correlation_genes <- correlation_summary
signif_correlation_genes <- signif_correlation_genes[order(signif_correlation_genes$cor),]
signif_correlation_genes_rank <- signif_correlation_genes$gene
genelist <- signif_correlation_genes[,c("gene","cor")]
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
