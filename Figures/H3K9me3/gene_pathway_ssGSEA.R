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
Indicator <- "logFC"
summary_antibody <- medians[,c("tissue",Indicator)]

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

# mitotic_nuclear_division <- read.csv("data/public_data/GO_term_summary_0140014.csv")
# mitotic_nuclear_division <- unique(mitotic_nuclear_division$Symbol)
# target_genes <- mitotic_nuclear_division
# GO_id <- "GO:0140014"
# description <- "mitotic nuclear division"
immunoglobin_production <- read.csv("data/public_data/GO_term_summary_0002377.csv")
immunoglobin_production <- unique(immunoglobin_production$Symbol)
target_genes <- immunoglobin_production
GO_id <- "GO:0002377"
description <- "immunoglobin production"
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
  # ggsave("result/figures/H3K9me3_mitotic_nuclear_division_kmeans1.pdf",p,width = 10,height = 5)
}
