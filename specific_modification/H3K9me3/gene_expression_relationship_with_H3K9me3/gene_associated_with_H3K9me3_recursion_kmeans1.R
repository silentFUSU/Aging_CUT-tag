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
# regions <- regions[-which(regions$V1=="chrY"),]
annotation <- data.frame(label = paste0(regions$V1,":",regions$V2,"-",regions$V3),chr=regions$V1)
rownames(annotation) <- annotation$label
annotation <- annotation[,-1,drop=F]
annotation$chr <- factor(annotation$chr,levels=paste0("chr",c(1:19,"X","Y")))
# regions <- regions[-which(regions$V1 == "chrY"),]
regions <- paste0(regions$V1,":",regions$V2,"-",regions$V3)
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver","lung","muscle","pancreas","skin","spleen","stomach","testis","thymus","tongue","iWAT")
# tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver","lung","muscle","pancreas","skin","spleen","stomach","testis","thymus","tongue","iWAT")
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
counts <- data.frame() 
search_table <- read.csv("data/samples/all/RNA_search_table.csv")
for(tissue in tissues){
  df <- read.table(paste0("data/samples/RNA/",tissue,"/combined-chrM.counts"),header = T)
  df <- df[,c(1,7:ncol(df))]
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|HM[0-9]+).*"
  colnames(df)[-1] <- gsub(pattern, "\\1",colnames(df)[-1])
  df <- df[,c("Geneid",search_table$sample_name[which(search_table$tissue==tissue_label_change(tissue))])]
  if(nrow(counts)==0){
    counts <- df
  }else{
    counts <- merge(counts,df,by="Geneid")
  }
}
rownames(counts) <- counts$Geneid
counts <- counts[,-1]
CPM <- as.data.frame(edgeR::cpm(counts))

Indicator <- "logFC"
summary_antibody <- medians[,c("tissue",Indicator)]

correlation_summary <- data.frame()
method <- "spearman"
# method <- "pearson"
for(i in c(1:nrow(CPM))){
  gene <- rownames(CPM)[i]
  t_CPM <- as.data.frame(t(CPM[i,]))
  search_table <- read.csv("data/samples/all/RNA_search_table.csv")
  colnames(t_CPM)[1] <- "CPM"
  t_CPM <- merge(t_CPM,search_table,by.x="row.names",by.y="sample_name")
  average_CPM <- t_CPM %>%
    group_by(tissue) %>%
    summarize(mean_CPM = mean(CPM))
  average_CPM <- as.data.frame(average_CPM)
  average_CPM <- merge(average_CPM,summary_antibody,by="tissue")
  colnames(average_CPM)[3] <- "histone"
  cortest <- cor.test(average_CPM$mean_CPM,average_CPM$histone,method = method)
  t_correlation_summary <- data.frame(p_value=cortest$p.value[1], cor = as.numeric(cortest$estimate),gene=gene,histone=antibody)
  correlation_summary <- rbind(correlation_summary,t_correlation_summary)
}
correlation_summary_output <- correlation_summary 
correlation_summary_output$cor <- -correlation_summary_output$cor
correlation_summary_output <- correlation_summary_output[order(correlation_summary_output$cor,decreasing = T),]
write.csv(correlation_summary_output,"data/samples/all/H3K9me3/recursion_peaks_diff_table/gene_expression_associated_with_kmeans1.csv",row.names = F)
positive_correlation <- correlation_summary[which(correlation_summary$p_value < 0.05 & correlation_summary$cor < 0),]
negative_correlation <- correlation_summary[which(correlation_summary$p_value < 0.05 & correlation_summary$cor > 0),]

positive_correlation <- positive_correlation[order(positive_correlation$cor),]
positive_correlation_genes <- positive_correlation$gene[1:min(20,nrow(positive_correlation))]
negative_correlation <- negative_correlation[order(negative_correlation$cor,decreasing = T),]
negative_correlation_genes <- negative_correlation$gene[1:min(20,nrow(negative_correlation))]

to_plot <- data.frame()
for(i in c(1:length(positive_correlation_genes))){
  gene <- positive_correlation_genes[i]
  t_CPM <- as.data.frame(t(CPM[gene,]))
  search_table <- read.csv("data/samples/all/RNA_search_table.csv")
  colnames(t_CPM)[1] <- "CPM"
  t_CPM <- merge(t_CPM,search_table,by.x="row.names",by.y="sample_name")
  average_CPM <- t_CPM %>%
    group_by(tissue) %>%
    summarize(mean_CPM = mean(CPM))
  average_CPM <- as.data.frame(average_CPM)
  colnames(average_CPM)[2] <- gene
  if(nrow(to_plot)==0){
    to_plot <- average_CPM
  } else{
    to_plot <- merge(to_plot,average_CPM,by="tissue")
  }
}
to_plot <- reshape2::melt(to_plot)
to_plot <- merge(to_plot,summary_antibody,by="tissue")
colnames(to_plot)[4] <- "histone" 
if(Indicator=="rank"){
  to_plot$histone <- factor(to_plot$histone,c(27:1))
  ggplot(to_plot,aes(x=histone,y=log2(value),color =variable))+    
    geom_jitter(size = 3, alpha = 0.7)+
    ggtitle(paste0("Gene expression relationship with ",antibody," change"))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab(Indicator)+labs(fill = "", color = "") +ylab(paste0("log2(CPM)"))
}else{
  ggplot(to_plot,aes(x=histone,y=log2(value),color =variable))+    
    geom_jitter(size = 3, alpha = 0.7)+
    ggtitle(paste0("Gene expression relationship with ",antibody," change"))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab(Indicator)+labs(fill = "", color = "") +ylab(paste0("log2(CPM)"))+
    scale_x_reverse()
  
}

to_plot <- data.frame() 
for(i in c(1:length(negative_correlation_genes))){
  gene <- negative_correlation_genes[i]
  t_CPM <- as.data.frame(t(CPM[gene,]))
  search_table <- read.csv("data/samples/all/RNA_search_table.csv")
  colnames(t_CPM)[1] <- "CPM"
  t_CPM <- merge(t_CPM,search_table,by.x="row.names",by.y="sample_name")
  average_CPM <- t_CPM %>%
    group_by(tissue) %>%
    summarize(mean_CPM = mean(CPM))
  average_CPM <- as.data.frame(average_CPM)
  colnames(average_CPM)[2] <- gene
  if(nrow(to_plot)==0){
    to_plot <- average_CPM
  } else{
    to_plot <- merge(to_plot,average_CPM,by="tissue")
  }
}
to_plot <- reshape2::melt(to_plot)
to_plot <- merge(to_plot,summary_antibody,by="tissue")
colnames(to_plot)[4] <- "histone" 
if(Indicator=="rank"){
  to_plot$histone <- factor(to_plot$histone,c(27:1))
  ggplot(to_plot,aes(x=histone,y=log2(value),color =variable))+    
    geom_jitter(position = position_jitter(width = 0.2), size = 3, alpha = 0.7)+
    ggtitle(paste0("Gene expression relationship with ",antibody," change"))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab(Indicator)+labs(fill = "", color = "") +ylab(paste0("log2(CPM)"))
}else{
  ggplot(to_plot,aes(x=histone,y=log2(value),color =variable))+    
    geom_jitter( size = 3, alpha = 0.7)+
    ggtitle(paste0("Gene expression relationship with ",antibody," change"))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab(Indicator)+labs(fill = "", color = "") +ylab(paste0("log2(CPM)")) +
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
barplot(genelist_up_GO,title = paste0("Gene expression positive correlaiton with ",antibody," change"),label_format = 50,showCategory = 30)


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

immunoglobin_production <- read.csv("data/public_data/GO_term_summary_0002377.csv")
immunoglobin_production <- unique(immunoglobin_production$Symbol)

####### ssGSEA
rpkm <- data.frame()
for(tissue in tissues){
  df <- read.table(paste0("data/samples/RNA/",tissue,"/combined-chrM.counts"),header = T)
  df <- df[,c(1,6,7:ncol(df))]
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|HM[0-9]+).*"
  colnames(df)[-c(1,2)] <- gsub(pattern, "\\1",colnames(df)[-c(1,2)])
  gene_lengths <- df$Length
  total_mapped_reads <- colSums(df[, 3:ncol(df)])
  rpkm_df <- data.frame(Geneid = df$Geneid)
  for (i in 3:ncol(df)) {
    counts <- df[[i]]
    t_rpkm <- (counts / (gene_lengths / 1000)) / (total_mapped_reads[i - 2] / 1e6)
    rpkm_df[[colnames(df)[i]]] <- t_rpkm
  }
  
  if(nrow(rpkm)==0){
    rpkm <- rpkm_df
  }else{
    rpkm <- merge(rpkm,rpkm_df,by="Geneid")
  }
}
rownames(rpkm) <- rpkm$Geneid
rpkm <- rpkm[,-1]

mitotic_nuclear_division <- read.csv("data/public_data/GO_term_summary_0140014.csv")
mitotic_nuclear_division <- unique(mitotic_nuclear_division$Symbol)
target_genes <- mitotic_nuclear_division
GO_id <- "GO:0140014"
description <- "mitotic nuclear division"
# immunoglobin_production <- read.csv("data/public_data/GO_term_summary_0002377.csv")
# immunoglobin_production <- unique(immunoglobin_production$Symbol)
# target_genes <- immunoglobin_production
# GO_id <- "GO:0002377"
# description <- "immunoglobin production"
# description <- result$Description[which(result$ID==GO_id)]
search_table <- read.csv("data/samples/all/RNA_search_table.csv")
genelist <- list(score=target_genes)
rpkm_matrix <- as.matrix(rpkm)
re <- gsva(rpkm_matrix,genelist , method="ssgsea",ssgsea.norm=TRUE) 
re <- as.data.frame(t(re))
to_plot <- merge(re,search_table,by.x="row.names",by.y="sample_name")
color <- read.table("data/samples/30_distinct_color.txt")
tissues_label <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
tissues_label <- sapply(tissues_label, tissue_label_change)
color <- setNames(color$V1,sort(tissues_label))
color <- color[which(names(color) %in% sapply(tissues, tissue_label_change))]

average_scores <-aggregate(score ~ tissue, data = to_plot, FUN = mean)
average_scores <- average_scores[order(average_scores$score),]
to_plot$tissue <- factor(to_plot$tissue,levels = average_scores$tissue)
to_plot$age[which(to_plot$age=="3m")] <- "Young"
to_plot$age[which(to_plot$age=="24m")] <- "Old"
to_plot$age <- factor(to_plot$age,levels=c("Young","Old"))
ggplot(to_plot,aes(x=tissue,y=score,color = tissue,shape=age))+    
  geom_jitter( size = 3, alpha = 0.7)+
  scale_color_manual(values = color)+
  ggtitle(paste0(description))+
  theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("Tissues")+labs(fill = "", color = "") 

to_plot <- merge(to_plot,summary_antibody,by="tissue")
colnames(to_plot)[ncol(to_plot)] <- "histone"
cor_test <- cor.test(to_plot$score,to_plot$histone,method="spearman")
average_scores <- merge(average_scores,summary_antibody,by="tissue")
colnames(average_scores)[ncol(average_scores)] <- "histone"
cor_test <- cor.test(average_scores$score,average_scores$histone,method="spearman")
lm_model <- lm(score ~ histone, data = average_scores)
r_squared <- summary(lm_model)$r.squared

if(Indicator=="rank"){
  to_plot$histone <- factor(to_plot$histone,levels = c(27:1))
  ggplot(to_plot,aes(x=histone,y=score,color = tissue,shape=age))+    
    geom_jitter( size = 3, alpha = 0.7)+
    scale_color_manual(values = color)+
    ggtitle(paste0(description," with H3K9me3"))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("Rank")+labs(fill = "", color = "")
}else{
 p <-  ggplot(to_plot,aes(x=histone,y=score,color = tissue,shape=age))+    
    geom_jitter( size = 3, alpha = 0.7)+
    scale_color_manual(values = color,breaks = sort(sapply(tissues, tissue_label_change)))+
    ggtitle(paste0(description," with H3K9me3"))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("H3K9m3 log2(Fold change)")+labs(fill = "", color = "")+
    scale_x_reverse()
 
 p <-  ggplot(average_scores,aes(x=histone,y=score))+    
   geom_jitter(size = 3, alpha = 0.7,color="#4DBBD5")+
   geom_smooth(method = "lm", color = "#3C5488", se = TRUE, level = 0.95) +
   # scale_color_manual(values = color,breaks = sort(sapply(tissues, tissue_label_change)))+
   ggtitle(paste0(description," with H3K9me3"))+
   theme_bw()+theme(text = element_text(size = 18))+
   xlab("H3K9m3 log2(Fold change)")+labs(fill = "", color = "")+
   scale_x_reverse()
 ggsave("result/figures/H3K9me3_mitotic_nuclear_division_kmeans1.pdf",p,width = 6,height = 4)
}

GO_result <- result[which(result$ID==GO_id),]
GO_result_gene <- strsplit(GO_result$geneID, split = "/")
GO_result_gene <- GO_result_gene[[1]]
genelist <- list(score=GO_result_gene)
re <- gsva(counts_matrix,genelist , method="ssgsea",ssgsea.norm=TRUE) 
re <- as.data.frame(t(re))
to_plot <- merge(re,search_table,by.x="row.names",by.y="sample_name")
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(unique(as.character(to_plot$tissue))))

average_scores <-aggregate(score ~ tissue, data = to_plot, FUN = mean)
average_scores <- average_scores[order(average_scores$score),]
to_plot$tissue <- factor(to_plot$tissue,levels = average_scores$tissue)
to_plot$age[which(to_plot$age=="3m")] <- "Young"
to_plot$age[which(to_plot$age=="24m")] <- "Old"
to_plot$age <- factor(to_plot$age,levels=c("Young","Old"))
ggplot(to_plot,aes(x=tissue,y=score,color = tissue,shape=age))+    
  geom_jitter(position = position_jitter(width = 0.2), size = 3, alpha = 0.7)+
  scale_color_manual(values = color)+
  ggtitle(paste0(description))+
  theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("Tissues")+labs(fill = "", color = "") 

to_plot <- merge(to_plot,summary_antibody,by="tissue")
colnames(to_plot)[ncol(to_plot)] <- "histone"
cor_test <- cor.test(to_plot$score,to_plot$histone)
average_scores <-aggregate(score ~ histone, data = to_plot, FUN = mean)
cor_test <- cor.test(average_scores$score,average_scores$histone)
colnames(to_plot)[ncol(to_plot)] <- "histone"
ggplot(to_plot,aes(x=histone,y=score,color = tissue,shape=age))+    
  geom_jitter(position = position_jitter(width = 0.2), size = 3, alpha = 0.7)+
  scale_color_manual(values = color)+
  ggtitle(paste0(description," Score with H3K9me3"))+
  theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("Rank")+labs(fill = "", color = "") +
  scale_x_continuous(breaks = 1:27) 

#### GSEA 
signif_correlation_genes <- correlation_summary
signif_correlation_genes <- signif_correlation_genes[order(signif_correlation_genes$cor),]
signif_correlation_genes_rank <- signif_correlation_genes$gene
# genelist <- bitr(signif_correlation_genes_rank,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
# genelist <- merge(genelist,signif_correlation_genes,by.x="SYMBOL",by.y="gene")
genelist <- signif_correlation_genes[,c("gene","cor")]
genelist$cor <- -genelist$cor
# genelist <- genelist[,c("ENTREZID","cor")]
mitotic_nuclear_division <- read.csv("data/public_data/GO_term_summary_0140014.csv")
mitotic_nuclear_division <- unique(mitotic_nuclear_division$Symbol)
target_genes <- mitotic_nuclear_division
# target_genes <- bitr(target_genes,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
description <- "mitotic nuclear division"
target_genes <- data.frame(TERM = description, ENTREZID = target_genes)
genelist <- setNames(genelist$cor, genelist$gene)
genelist <- sort(genelist, decreasing = TRUE)
GSEA_result <- GSEA(genelist, TERM2GENE = target_genes)
gseaplot2(
  GSEA_result,
  geneSetID = "mitotic nuclear division", # 使用一个实际存在的 TERM
  title = "Enrichment Plot for Mitotic Nuclear Division",
  pvalue_table = TRUE
)

#### GSEA GO term

gse <- gseGO(geneList=genelist, 
             ont = "BP",
             keyType = "SYMBOL", 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = GO_database,
             pAdjustMethod = "none",eps = 1e-100)
results_df <- as.data.frame(gse@result)
write.csv(results_df,"data/samples/all/H3K9me3/recursion_peaks_diff_table/gene_expression_associated_with_kmeans1_GO_terms.csv",row.names = F)
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
p <- emapplot(negative_gsea, showCategory = 50) 
ggsave(paste0("result/Sup_figures/H3K9me3_kmeans1_gene_expression_correlation_GO_pathway.pdf"),p,width = 12,height = 10)
