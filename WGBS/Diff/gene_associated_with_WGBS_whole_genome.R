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
library(MASS)  
library(RANSAC)
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

# tissue_summary <- read.csv("data/samples/WGBS/all_tissues_delta_in_200kb_bins_cross_comparison.csv",row.names = 1)
# H3K9me3_tissue_order_label <- c() 
# annotation_col <- data.frame()
# for(tissue in tissues){
#   search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
#   t_search_table <- search_table[which(search_table$tissue==tissue),]
#   young_samples <- t_search_table$sample_name[which(t_search_table$age=="3M")]
#   old_samples <- t_search_table$sample_name[which(t_search_table$age=="24M")]
#   combinations <- as.data.frame(expand.grid(young = young_samples, old = old_samples))
#   if(tissue=="bonemarrow"){
#     H3K9me3_tissue_order_label <- c(H3K9me3_tissue_order_label, paste0("Bone.Marrow",".",paste0(combinations$old,".",combinations$young)))  
#     t_annotation_col <- data.frame(tissue=tissue_label_change(tissue),sample=paste0("Bone.Marrow",".",paste0(combinations$old,".",combinations$young)))  
#   }else if(tissue=="mammarygland"){
#     H3K9me3_tissue_order_label <- c(H3K9me3_tissue_order_label, paste0("Mammary.Gland",".",paste0(combinations$old,".",combinations$young)))  
#     t_annotation_col <- data.frame(tissue=tissue_label_change(tissue),sample=paste0("Mammary.Gland",".",paste0(combinations$old,".",combinations$young)))  
#   }
#   else{
#     H3K9me3_tissue_order_label <- c(H3K9me3_tissue_order_label, paste0(tissue_label_change(tissue),".",paste0(combinations$old,".",combinations$young)))  
#     t_annotation_col <- data.frame(tissue=tissue_label_change(tissue),sample=paste0(tissue_label_change(tissue),".",paste0(combinations$old,".",combinations$young)))  
#   }
#   annotation_col <- rbind(annotation_col,t_annotation_col)
# }
# rownames(tissue_summary) <- tissue_summary$label
# 
# tissue_mean_summary <- data.frame() 
# annotation <- annotation_col
# for(tissue in tissues){
#   t_annotation <- annotation[which(annotation$tissue==tissue_label_change(tissue)),]
#   t_tissue_mean_summary <- tissue_summary[,t_annotation$sample]
#   t_tissue_mean_summary$mean_delta <- rowMeans(t_tissue_mean_summary) 
#   t_tissue_mean_summary$label <- rownames(t_tissue_mean_summary)
#   t_tissue_mean_summary <- t_tissue_mean_summary[,c("label","mean_delta")]
#   colnames(t_tissue_mean_summary)[2] <- tissue_label_change(tissue)
#   if(nrow(tissue_mean_summary)==0){
#     tissue_mean_summary <- t_tissue_mean_summary  
#   }else{
#     tissue_mean_summary <- merge(tissue_mean_summary,t_tissue_mean_summary,by="label")
#   }
# }
# 
# to_plot_WGBS_order_long <- reshape2::melt(tissue_mean_summary)

summary <- read.csv("data/samples/WGBS/CG_manual.csv")
search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
summary$tissue <- sapply(summary$tissue,tissue_label_change)
summary <- merge(summary,search_table[,c("sample_name","age")],by.x="sample",by.y="sample_name")
result <- summary %>%
  group_by(tissue, age) %>%
  summarise(mean_CG = mean(CG, na.rm = TRUE), .groups = 'drop')
result$mean_CG <- result$mean_CG *100

result_young <- result[which(result$age=="3M"),]
result_old <- result[which(result$age=="24M"),]
to_plot <- merge(result_young,result_old,by="tissue")
to_plot$delta <- to_plot$mean_CG.y - to_plot$mean_CG.x
to_plot <- to_plot[order(to_plot$delta),]

# medians <- to_plot_WGBS_order_long %>%
#   group_by(variable) %>%
#   summarise(median_value = median(value, na.rm = TRUE))
medians <- to_plot[,c("tissue","delta")]
colnames(medians) <- c("variable","median_value")
medians <- medians[-which(medians$variable %in% c("Mammary Gland","Uterus","Ovary")),]
medians <- medians[order(medians$median_value),]
medians$rank <- c(1:nrow(medians))
medians <- as.data.frame(medians)

counts <- data.frame() 
search_table <- read.csv("data/samples/all/RNA_search_table.csv")
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
                  "lung","muscle","pancreas","skin","spleen","stomach","testis","thymus","tongue","iWAT","ileum")) 
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

colnames(medians)[1:2] <- c("tissue","delta")
Indicator <- "delta"
summary_antibody <- medians[,c("tissue",Indicator)]
correlation_summary <- data.frame()
method <- "spearman"
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
    ggtitle(paste0("Gene expression positively correlated with the degree of downregulation of WGBS"))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab(Indicator)+labs(fill = "", color = "") +ylab(paste0("log2(CPM)"))
}else{
  ggplot(to_plot,aes(x=histone,y=log2(value),color =variable))+    
    geom_jitter(size = 3, alpha = 0.7)+
    ggtitle(paste0("Gene expression positively correlated with the degree of downregulation of WGBS"))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("Delta")+labs(fill = "", color = "") +ylab(paste0("log2(CPM)"))+
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
    ggtitle(paste0("Gene expression negatively correlated with the degree of downregulation of WGBS"))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab(Indicator)+labs(fill = "", color = "") +ylab(paste0("log2(CPM)"))
}else{
  ggplot(to_plot,aes(x=histone,y=log2(value),color =variable))+    
    geom_jitter( size = 3, alpha = 0.7)+
    ggtitle(paste0("Gene expression negatively correlated with the degree of downregulation of WGBS"))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("Delta")+labs(fill = "", color = "") +ylab(paste0("log2(CPM)")) +
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
barplot(genelist_up_GO,title = paste0("Gene expression positive correlaiton with WGBS change"),label_format = 50,showCategory = 30)

genelist_down <- bitr(negative_correlation$gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_down_GO <- enrichGO( genelist_down$ENTREZID,#GO富集分析
                              OrgDb = GO_database,
                              keyType = "ENTREZID",#设定读取的gene ID类型
                              ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                              pvalueCutoff = 0.05,#设定p值阈值
                              qvalueCutoff = 0.05,#设定q值阈值
                              readable = T)
result <- as.data.frame(genelist_down_GO@result)
barplot(genelist_down_GO,title = paste0("Gene expression negative correlaiton with WGBS change"),label_format = 70,showCategory = 30)

####### ssGSEA
mitotic_nuclear_division <- read.csv("data/public_data/GO_term_summary_0140014.csv")
mitotic_nuclear_division <- unique(mitotic_nuclear_division$Symbol)
target_genes <- mitotic_nuclear_division
GO_id <- "GO:0140014"
# target_genes <- immunoglobin_production
# GO_id <- "GO:0002377"
# description <- result$Description[which(result$ID==GO_id)]
description <- "mitotic nuclear division"
search_table <- read.csv("data/samples/all/RNA_search_table.csv")
genelist <- list(score=target_genes)

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
rpkm_matrix <- as.matrix(rpkm)

re <- gsva(rpkm_matrix,genelist , method="ssgsea",ssgsea.norm=TRUE) 
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
cor_test <- cor.test(to_plot$score,to_plot$histone,method="spearman")
average_scores <- merge(average_scores,summary_antibody,by="tissue")
colnames(average_scores)[ncol(average_scores)] <- "histone"
average_scores_test  <- average_scores[which(average_scores$tissue !="Pancreas"),]
cor_test <- cor.test(average_scores_test$score,average_scores_test$histone,method="spearman")
if(Indicator=="rank"){
  to_plot$histone <- factor(to_plot$histone,levels = c(27:1))
}
p <- ggplot(to_plot,aes(x=histone,y=score,color = tissue,shape=age))+    
  geom_jitter(size = 3, alpha = 0.7)+
  scale_color_manual(values = color)+
  ggtitle(paste0(description," with DNA methylation"))+
  theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("Delta")+labs(fill = "", color = "")+
  scale_x_reverse() 

#### Ransac regression
# set.seed(1)
# model <- ransac_reg(score ~ histone, data = average_scores, n_min = 2, tol = 0.001,n_iter = 10000,verbose = T)
# coefficients <- coef(model)
# intercept <- coefficients[1]
# slope <- -coefficients[2]
# predictions <-average_scores %>%
#   mutate(fitted = intercept + slope * histone)
# 
# 
# # Calculate the confidence interval (we'll use 95% CI here)
# alpha <- 0.05
# ci_multiplier <- qt(1 - alpha/2, df = length(average_scores$score) - 2)
# average_scores <- average_scores %>%
#   mutate(
#     fitted = predictions$fit,
#     se = predictions$se.fit,
#     lower = fitted - ci_multiplier * se,
#     upper = fitted + ci_multiplier * se
#   )

p <- ggplot(average_scores,aes(x=histone,y=score))+    
  geom_jitter(size = 3, alpha = 0.7,color="#f39b7f")+
  geom_smooth(data = average_scores[-which(average_scores$tissue %in% c("Pancreas")),], aes(x = histone, y = score),
              method = "lm", color = "#e64b35", se = TRUE, level = 0.95) +
  # geom_abline(intercept = intercept, slope = slope, color = "#e64b35", size = 1) +
  # geom_ribbon(aes(ymin = lower, ymax = upper))+
  ggtitle(paste0(description))+
  theme_bw()+theme(text = element_text(size = 18))+
  xlab("Delta")+labs(fill = "", color = "")+
  scale_x_reverse()
p
ggsave("result/figures/WGBS_mitotic_nuclear_division_ransac.pdf",p,width = 6,height = 4)

#### GSEA 
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
p <- emapplot(positive_gsea, showCategory = 50) 
ggsave(paste0("result/Sup_figures/WGBS_gene_expression_correlation_GO_pathway_positive.pdf"),p,width = 12,height = 10)

negative_gsea <- pairwise_termsim(negative_gsea)  
p <-emapplot(negative_gsea, showCategory = 50) 
ggsave(paste0("result/Sup_figures/WGBS_gene_expression_correlation_GO_pathway_negative.pdf"),p,width = 12,height = 10)
