rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
library(ggplot2)
library(stringr)
library(dplyr)
library(dbplyr)
library(clusterProfiler)
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
      tissue_label <- "Mammary gland"
    }else if(tissue_label=="Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
}
tissues <- c("skin","CB","spleen","heart","bladder","tongue","uterus","aorta","thymus","stomach","Hip","brain","BAT","iWAT","muscle","bonemarrow","lung","kidney","liver","testis","colon","cecum","ileum","jejunum","ovary","mammarygland","pancreas")
# antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me3","H3K4me1")
antibodys <- c("H3K9me3")
summary <- data.frame()
for(antibody in antibodys){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    bin_size <- "10kb"
  }else{
    bin_size <- "1kb"
  }
  summary_per_antibody <- data.frame()
  for(tissue in tissues){
    df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
    df <- df[which(df$Significant != "Stable"),]
    if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
      peaks <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_",bin_size,"_in_young_old_merge-W1000-G3000-E100.bed"))    
    }else{
      peaks <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_",bin_size,"_in_young_old_merge_macs_narrowpeak.bed"))
    }
    df <- df[which(df$Geneid %in% peaks$V4),]
    t_summary_per_antibody <- data.frame(tissue=tissue_label_change(tissue),count=nrow(df))
    summary_per_antibody <- rbind(summary_per_antibody,t_summary_per_antibody)
  }
  summary_per_antibody <- summary_per_antibody[order(summary_per_antibody$count,decreasing = TRUE),]
  summary_per_antibody$rank <- c(1:length(tissues))
  colnames(summary_per_antibody)[which(colnames(summary_per_antibody)=="rank")] <- antibody
  summary_per_antibody <- summary_per_antibody[,c(1,3)]
  if(nrow(summary)==0){
    summary <- summary_per_antibody
  }else{
    summary <- merge(summary,summary_per_antibody,by="tissue")
  }
}
antibody <- "H3K9me3"
summary_antibody <- summary[,c("tissue",antibody)]

logFC<- data.frame() 
search_table <- read.csv("data/samples/all/RNA_search_table.csv")
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"),header = T)
  df <- df[,c("X","logFC")]
  colnames(df) <- c("Geneid",tissue_label_change(tissue))
  if(nrow(logFC)==0){
    logFC <- df
  }else{
    logFC <- merge(logFC,df,by="Geneid")
  }
}
rownames(logFC) <- logFC$Geneid
logFC <- logFC[,-1]

correlation_summary <- data.frame()
for(i in c(1:nrow(logFC))){
  gene <- rownames(logFC)[i]
  t_logFC <- as.data.frame(t(logFC[i,]))
  search_table <- read.csv("data/samples/all/RNA_search_table.csv")
  colnames(t_logFC)[1] <- "logFC"
  t_logFC <- merge(t_logFC,summary_antibody,by.x="row.names",by.y="tissue")
  colnames(t_logFC)[3] <- "histone"
  cortest <- cor.test(t_logFC$logFC,t_logFC$histone)
  t_correlation_summary <- data.frame(p_value=cortest$p.value[1], cor = as.numeric(cortest$estimate),gene=gene,histone=antibody)
  correlation_summary <- rbind(correlation_summary,t_correlation_summary)
}
positive_correlation <- correlation_summary[which(correlation_summary$p_value < 0.05 & correlation_summary$cor < 0),]
negative_correlation <- correlation_summary[which(correlation_summary$p_value < 0.05 & correlation_summary$cor > 0),]

positive_correlation <- positive_correlation[order(positive_correlation$cor),]
positive_correlation_genes <- positive_correlation$gene[1:20]
to_plot <- data.frame()

for(i in c(1:length(positive_correlation_genes))){
  gene <- positive_correlation_genes[i]
  t_logFC <- as.data.frame(t(logFC[gene,]))
  search_table <- read.csv("data/samples/all/RNA_search_table.csv")
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
to_plot$histone <- factor(to_plot$histone, levels = 27:1)  
ggplot(to_plot,aes(x=histone,y=value,color =variable))+    
  geom_jitter(position = position_jitter(width = 0.2), size = 3, alpha = 0.7)+
  ggtitle(paste0("Gene expression relationship with ",antibody," change"))+
  theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("Rank")+labs(fill = "", color = "") +ylab(paste0("log2(Fold Change)"))

negative_correlation <- negative_correlation[order(negative_correlation$cor,decreasing = T),]
negative_correlation_genes <- negative_correlation$gene[1:20]
to_plot <- data.frame()
for(i in c(1:length(negative_correlation_genes))){
  gene <- negative_correlation_genes[i]
  t_logFC <- as.data.frame(t(logFC[gene,]))
  search_table <- read.csv("data/samples/all/RNA_search_table.csv")
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
to_plot$histone <- factor(to_plot$histone, levels = 27:1)  
ggplot(to_plot,aes(x=histone,y=value,color =variable))+    
  geom_jitter(position = position_jitter(width = 0.2), size = 3, alpha = 0.7)+
  ggtitle(paste0("Gene expression relationship with ",antibody," change"))+
  theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("Rank")+labs(fill = "", color = "") +ylab(paste0("log2(CPM)"))

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
genelist_up <- bitr(positive_correlation$gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_up_GO <- enrichGO( genelist_up$ENTREZID,#GO富集分析
                            OrgDb = GO_database,
                            keyType = "ENTREZID",#设定读取的gene ID类型
                            ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                            pvalueCutoff = 0.05,#设定p值阈值
                            qvalueCutoff = 0.05,#设定q值阈值
                            readable = T)
barplot(genelist_up_GO,title = paste0("Gene expression positive correlaiton with ",antibody," change"),label_format = 50,showCategory = 30)

genelist_down <- bitr(negative_correlation$gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_down_GO <- enrichGO( genelist_down$ENTREZID,#GO富集分析
                              OrgDb = GO_database,
                              keyType = "ENTREZID",#设定读取的gene ID类型
                              ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                              pvalueCutoff = 0.05,#设定p值阈值
                              qvalueCutoff = 0.05,#设定q值阈值
                              readable = T)
barplot(genelist_down_GO,title = paste0("Gene expression negative correlaiton with ",antibody," change"),label_format = 50,showCategory = 30)
