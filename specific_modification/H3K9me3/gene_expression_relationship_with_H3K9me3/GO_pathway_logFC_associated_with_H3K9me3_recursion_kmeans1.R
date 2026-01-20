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

medians <- apply(diff_summary, 2, median, na.rm = TRUE)
medians <- data.frame(tissues=colnames(diff_summary),median=medians)
medians <- medians[order(medians$median),]
medians$rank <- 1:nrow(medians)
colnames(medians) <- c("tissue","logFC","rank")
search_table <- read.csv("data/samples/all/RNA_search_table.csv")
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

Indicator <- "logFC"
summary_antibody <- medians[,c("tissue",Indicator)]

#### ssGSEA score logFC
# immunoglobin_production <- read.csv("data/public_data/GO_term_summary_0002377.csv")
# immunoglobin_production <- unique(immunoglobin_production$Symbol)
# target_genes <- immunoglobin_production
# description <- "immunoglobin production"
# activation_of_immune_response <- read.csv("data/public_data/GO_term_summary_0002253.csv")
# activation_of_immune_response <- unique(activation_of_immune_response$Symbol)
# target_genes <- activation_of_immune_response
# description <- "activation of immune response"
positive_regulation_of_innate_immune_response <- read.csv("data/public_data/GO_term_summary_0045089.csv")
positive_regulation_of_innate_immune_response <- unique(positive_regulation_of_innate_immune_response$Symbol)
target_genes <- positive_regulation_of_innate_immune_response
description <- "positive regulation of innate immune response"

search_table <- read.csv("data/samples/all/RNA_search_table.csv")
genelist <- list(score=target_genes)
rpkm_matrix <- as.matrix(rpkm)
re <- gsva(rpkm_matrix,genelist , method="ssgsea",ssgsea.norm=TRUE) 
re <- as.data.frame(t(re))
to_plot <- merge(re,search_table,by.x="row.names",by.y="sample_name")

to_plot_logFC <- data.frame()
for(tissue in tissues){
  t_to_plot <- to_plot[which(to_plot$tissue==tissue_label_change(tissue)),]
  mean_3m <- mean(t_to_plot$score[t_to_plot$age == "3m"])
  mean_24m <- mean(t_to_plot$score[t_to_plot$age == "24m"])
  logFC <- log2(mean_24m/mean_3m)
  t_to_plot <- data.frame(tissue=tissue_label_change(tissue),logFC_score=logFC)
  to_plot_logFC <- rbind(to_plot_logFC,t_to_plot)
}
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(unique(as.character(to_plot$tissue))))

to_plot_logFC <- merge(to_plot_logFC,medians,by="tissue")
cor_test <- cor.test(to_plot_logFC$logFC_score,to_plot_logFC$logFC,method = "spearman")
p <- ggplot(to_plot_logFC,aes(x=logFC,y=logFC_score,color = tissue))+    
  geom_jitter(size = 3, alpha = 0.7)+
  scale_color_manual(values = color)+
  ggtitle(paste0(description))+
  theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("H3K9me3 log2(Fold change)")+ylab("GO pathway score log2(Fold change)")+labs(fill = "", color = "")+scale_x_reverse()

linear_model <- lm(logFC_score ~ logFC, data = to_plot_logFC)
model_summary <- summary(linear_model)
p <-  ggplot(to_plot_logFC,aes(x=logFC,y=logFC_score))+    
  geom_jitter(size = 3, alpha = 0.7,color="#F39B7F")+
  geom_smooth(method = "lm", color = "#E64B35", se = TRUE, level = 0.95) +
  ggtitle(paste0(description))+
  theme_bw()+theme(text = element_text(size = 18))+
  xlab("H3K9m3 log2(Fold change)")+ylab("GO pathway score log2(Fold change)")+labs(fill = "", color = "")+
  scale_x_reverse()

ggsave("result/figures/H3K9me3_positive_regulation_of_innate_immune_response_logFC_kmeans1.pdf",p,width = 6,height = 4)
