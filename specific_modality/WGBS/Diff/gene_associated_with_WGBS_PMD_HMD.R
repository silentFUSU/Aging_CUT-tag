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
region <-"HMD"
tissue_overall_summary <- read.csv("data/samples/WGBS/all/DNA_methylation_change_in_PMD_HMD_overall_summary.csv",row.names = 1)
tissue_overall_summary <- tissue_overall_summary[which(tissue_overall_summary$V5==region),]
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
                  "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT","ileum"))
summary_antibody <- data.frame()
for(tissue in tissues){
  df <- tissue_overall_summary[which(tissue_overall_summary$tissue==tissue_label_change(tissue)),]
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  young <- mean(df$methylation[which(df$sample %in% c(search_table$sample_name[which(search_table$age=="3M")]))])
  old <- mean(df$methylation[which(df$sample %in% c(search_table$sample_name[which(search_table$age=="24M")]))])
  t_summary_antibody <- data.frame(tissue=tissue_label_change(tissue),delta=old-young)  
  summary_antibody <- rbind(summary_antibody,t_summary_antibody)
  }

mitotic_nuclear_division <- read.csv("data/public_data/GO_term_summary_0140014.csv")
mitotic_nuclear_division <- unique(mitotic_nuclear_division$Symbol)
target_genes <- mitotic_nuclear_division
GO_id <- "GO:0140014"

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

to_plot$tissue <- as.character(to_plot$tissue)
to_plot$tissue[which(to_plot$tissue=="Mammary gland")] <- "Mammary Gland"
to_plot <- merge(to_plot,summary_antibody,by="tissue")
colnames(to_plot)[ncol(to_plot)] <- "histone"
cor_test <- cor.test(to_plot$score,to_plot$histone,method="spearman")

average_scores$tissue[which(average_scores$tissue=="Mammary gland")] <- "Mammary Gland"
average_scores <- merge(average_scores,summary_antibody,by="tissue")
colnames(average_scores)[ncol(average_scores)] <- "histone"

cor_test <- cor.test(average_scores$score,average_scores$histone,method="spearman")
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(unique(as.character(to_plot$tissue))))
ggplot(to_plot,aes(x=histone,y=score,color = tissue,shape=age))+    
  geom_jitter(size = 3, alpha = 0.7)+
  scale_color_manual(values = color)+
  ggtitle(paste0(description," with DNA methylation"))+
  theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("Delta")+labs(fill = "", color = "")+
  scale_x_reverse() 

to_plot <- average_scores[which(!average_scores$tissue%in% c("Mammary Gland","Uterus","Ovary")),]
p <- ggplot(to_plot,aes(x=histone,y=score))+    
  geom_jitter(size = 3, alpha = 0.7,color="#f39b7f")+
  geom_smooth(data = to_plot[which(to_plot$tissue!="Pancreas"),], aes(x = histone, y = score),
              method = "lm", color = "#e64b35", se = TRUE, level = 0.95) +
  # geom_smooth(data = to_plot, aes(x = histone, y = score),
  #             method = "lm", color = "#e64b35", se = TRUE, level = 0.95) +
  ggtitle(paste0(description," ",region))+
  theme_bw()+theme(text = element_text(size = 18))+
  xlab("Delta")+labs(fill = "", color = "")+
  scale_y_continuous(limits = c(2.5, 3.8), breaks = seq(2.5, 3.8, by = 0.3))+
  scale_x_reverse(limits = c(10,-6))
p
ggsave(paste0("result/Sup_figures/",region,"_mitosis.pdf"),p,width = 6,height = 4)
cor_test <- cor.test(to_plot$score,to_plot$histone,method="spearman")
cor_test <- cor.test(to_plot$score[which(to_plot$tissue!="Pancreas")],to_plot$histone[which(to_plot$tissue!="Pancreas")],method="spearman")

