rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
library(CytoTRACE2)
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
tissue <- "lung"
search_table <- read.csv("data/samples/all/RNA_search_table.csv")
tissues <- c("skin","CB","spleen","heart","bladder","tongue","uterus","aorta","thymus","stomach","Hip","brain","BAT","iWAT","muscle","bonemarrow","lung","kidney","liver","testis","colon","cecum","ileum","jejunum","ovary","mammarygland","pancreas")
summary <- data.frame()
for(tissue in tissues){
  df <- read.table(paste0("data/samples/RNA/",tissue,"/combined-chrM.counts"),header = T)
  df <- df[,c(1,7:ncol(df))]
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+).*"
  colnames(df)[-1] <- gsub(pattern, "\\1",colnames(df)[-1])
  df <- df[,c("Geneid",search_table$sample_name[which(search_table$tissue==tissue_label_change(tissue))])]
  if(nrow(summary)==0){
    summary <- df
  }else{
    summary <- merge(summary,df,by="Geneid")
  }
}

rownames(summary) <- summary$Geneid
summary <- summary[,-1]
cytotrace2_result <- cytotrace2(summary,species = "mouse",slot_type="counts",seed=1,ncores=10)
cytotrace2_result$sample_name <- rownames(cytotrace2_result)
cytotrace2_result <- merge(cytotrace2_result,search_table,by="sample_name")

to_plot <- cytotrace2_result[,c("sample_name","CytoTRACE2_Score","tissue","age")]
average_scores <-aggregate(CytoTRACE2_Score ~ tissue, data = to_plot, FUN = mean)
average_scores <- average_scores[order(average_scores$CytoTRACE2_Score),]
to_plot$tissue <- factor(to_plot$tissue,levels = average_scores$tissue)
to_plot$age <- factor(to_plot$age,levels = c("3m","24m"))
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(unique(as.character(to_plot$tissue))))
ggplot(to_plot,aes(x=tissue,y=CytoTRACE2_Score,color = tissue,shape=age))+    
  geom_jitter(position = position_jitter(width = 0.2), size = 3, alpha = 0.7)+
  scale_color_manual(values = color)+
  ggtitle(paste0("Stemness Score"))+ylim(0,1)+
  theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("")+labs(fill = "", color = "") +ylab("CytoTrace2 stemness score")


antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me3","H3K4me1")
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
to_plot_rank <- merge(to_plot,summary, by="tissue")
antibodys <- c("H3K27me3","H3K36me3","H3K9me3","H3K4me3","H3K4me1","H3K27ac")
p_value_summary <- data.frame()
for(antibody in antibodys){
  t_to_plot_rank <- to_plot_rank[,c("tissue","sample_name","CytoTRACE2_Score","age",antibody)]
  colnames(t_to_plot_rank)[5] <- "histone"
  t_to_plot_rank$tissue <- as.character(t_to_plot_rank$tissue)
  correlation_test <- cor.test(t_to_plot_rank$histone, t_to_plot_rank$CytoTRACE2_Score)
  t_p_value_summary <- data.frame(antibody=antibody,p_value=correlation_test$p.value)
  p_value_summary <- rbind(p_value_summary,t_p_value_summary)
  p <- ggplot(t_to_plot_rank,aes(x=histone,y=CytoTRACE2_Score,color = tissue,shape=age))+    
    geom_jitter(position = position_jitter(width = 0.2), size = 3, alpha = 0.7)+
    scale_color_manual(values = color)+
    ggtitle(paste0("Stemness Score with ",antibody))+ylim(0,1)+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("Rank")+labs(fill = "", color = "") +ylab("CytoTrace2 stemness score")+
    scale_x_continuous(breaks = 1:27) 
  print(p)
  }
