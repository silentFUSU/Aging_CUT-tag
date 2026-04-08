rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
library(scran)

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
gene <- read.table("~/ref_data/for_normal_mapping/TSS/refBed/mm10_refGene.bed")
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
summary <- summary[,c("Geneid",sort(colnames(summary)[-1]))]
summary <- merge(summary,gene[,c("V4","V5")],by.x="Geneid",by.y="V5")
summary$V4 <- sub("\\..*", "", summary$V4)
rownames(summary) <- summary$V4
summary <- summary[,-1]
summary <- summary[,-ncol(summary)]
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
assignments <- cyclone(summary, mm.pairs, gene.names=rownames(summary))
assignments_df <- as.data.frame(assignments$normalized.scores)
rownames(assignments_df) <- colnames(summary)

cyclone_result <- merge(assignments_df,search_table,by.x="row.names",by.y="sample_name")
conditions <- c("G1","S","G2M")
for(condition in conditions){
  to_plot <- cyclone_result[,c("Row.names","tissue","age",condition)]
  colnames(to_plot)[ncol(to_plot)] <- "score"
  color <- read.table("data/samples/30_distinct_color.txt")
  color <- setNames(color$V1,sort(unique(to_plot$tissue)))
  average_scores <-aggregate(score ~ tissue, data = to_plot, FUN = mean)
  average_scores <- average_scores[order(average_scores$score),]
  to_plot$tissue <- factor(to_plot$tissue,levels = average_scores$tissue)
  to_plot$age <- factor(to_plot$age,levels = c("3m","24m"))

  p <- ggplot(to_plot,aes(x=tissue,y=score,color = tissue,shape=age))+    
    geom_jitter(position = position_jitter(width = 0.2), size = 3, alpha = 0.7)+
    scale_color_manual(values = color)+
    ggtitle(paste0(condition," Score"))+ylim(0,1)+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("")+labs(fill = "", color = "") +ylab(paste0(condition," score"))
  print(p)
  }

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

conditions <- c("G1","S","G2M")
p_value_summary <- data.frame()
for(condition in conditions){
  to_plot <- cyclone_result[,c("Row.names","tissue","age",condition)]
  colnames(to_plot)[ncol(to_plot)] <- "score"
  to_plot_rank <- merge(to_plot,summary, by="tissue")
  antibodys <- c("H3K27me3","H3K36me3","H3K9me3","H3K4me3","H3K4me1","H3K27ac")
  for(antibody in antibodys){
    t_to_plot_rank <- to_plot_rank[,c("tissue","Row.names","score","age",antibody)]
    colnames(t_to_plot_rank)[5] <- "histone"
    t_to_plot_rank$tissue <- as.character(t_to_plot_rank$tissue)
    correlation_test <- cor.test(t_to_plot_rank$histone, t_to_plot_rank$score)
    t_p_value_summary <- data.frame(condition=condition,antibody=antibody,p_value=correlation_test$p.value)
    p_value_summary <- rbind(p_value_summary,t_p_value_summary)
    t_to_plot_rank$age <- factor(t_to_plot_rank$age, levels=c("3m","24m"))
    p <- ggplot(t_to_plot_rank,aes(x=histone,y=score,color = tissue,shape=age))+    
      geom_jitter(position = position_jitter(width = 0.2), size = 3, alpha = 0.7)+
      scale_color_manual(values = color)+
      ggtitle(paste0(condition," Score with ",antibody))+ylim(0,1)+
      theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
      xlab("Rank")+labs(fill = "", color = "") +ylab(paste0(condition," score"))+
      scale_x_continuous(breaks = 1:27) 
    print(p)
  }
}


