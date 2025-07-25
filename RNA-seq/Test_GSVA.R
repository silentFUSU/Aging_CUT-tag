rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
library(ggplot2)
library(stringr)
library(dplyr)
library(dbplyr)
library(clusterProfiler)
library(GSVA)
library(org.Mm.eg.db)
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
#GO:0044772 mitotic cell cycle phase transition
#GO:0019827 stem cell population maintenance

mitotic_cell_cycle_phase_transition <- read.csv("data/public_data/GO_term_summary_0044772.csv")
mitotic_cell_cycle_phase_transition <- unique(mitotic_cell_cycle_phase_transition$Symbol)
stem_cell_population_maintenance<- read.csv("data/public_data/GO_term_summary_0019827.csv")
stem_cell_population_maintenance <- unique(stem_cell_population_maintenance$Symbol)

target_genes <- mitotic_cell_cycle_phase_transition
GO_id <- "GO:0044772"
description <- "mitotic cell cycle phase transition"
# target_genes <- stem_cell_population_maintenance
# GO_id <- "GO:0019827"
# description <- "stem cell population maintenance"
counts <- data.frame() 
search_table <- read.csv("data/samples/all/RNA_search_table.csv")
tissues <- c("skin","CB","spleen","heart","bladder","tongue","uterus","aorta","thymus","stomach","Hip","brain","BAT","iWAT","muscle","bonemarrow","lung","kidney","liver","testis","colon","cecum","ileum","jejunum","ovary","mammarygland","pancreas")
for(tissue in tissues){
  df <- read.table(paste0("data/samples/RNA/",tissue,"/combined-chrM.counts"),header = T)
  df <- df[,c(1,7:ncol(df))]
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+).*"
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
counts <- as.matrix(counts)
genelist <- list(score=target_genes)
re <- gsva(counts,genelist , method="ssgsea",ssgsea.norm=TRUE) 
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
to_plot <- merge(to_plot,summary,by="tissue")
colnames(to_plot)[ncol(to_plot)] <- "histone"
cor_test <- cor.test(to_plot$score,to_plot$histone)
average_scores <-aggregate(score ~ histone, data = to_plot, FUN = mean)
cor_test <- cor.test(average_scores$score,average_scores$histone)
ggplot(to_plot,aes(x=histone,y=score,color = tissue,shape=age))+    
  geom_jitter(position = position_jitter(width = 0.2), size = 3, alpha = 0.7)+
  scale_color_manual(values = color)+
  ggtitle(paste0(description," with H3K9me3"))+
  theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("Rank")+labs(fill = "", color = "") +ylab("stem cell population maintenance")+
  scale_x_continuous(breaks = 1:27) 

GO_result <- read.csv("result/RNA/histone_relationship_with_RNA/H3K9me3/gene_expression_relationship_with_H3K9me3_rank_positive_GO.csv")

GO_result <- GO_result[which(GO_result$ID==GO_id),]
GO_result_gene <- strsplit(GO_result$geneID, split = "/")
GO_result_gene <- GO_result_gene[[1]]
genelist <- list(score=GO_result_gene)
re <- gsva(counts,genelist , method="ssgsea",ssgsea.norm=TRUE) 
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

to_plot <- merge(to_plot,summary,by="tissue")
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

correlation_summary <- read.csv("result/RNA/histone_relationship_with_RNA/H3K9me3/gene_expression_relationship_with_H3K9me3_rank.csv")
correlation_summary <- correlation_summary[which(correlation_summary$gene %in% GO_result_gene),]



