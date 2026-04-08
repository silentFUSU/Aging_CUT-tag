rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(factoextra)
library(cluster)
library(umap)
library(data.table)
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

regions <- read.table("data/samples/all/H3K9me3/recursion_peaks_diff_table/bed/kmeans1_uinon_recursion_peaks.bed")
annotation <- data.frame(label = paste0(regions$V1,":",regions$V2,"-",regions$V3),chr=regions$V1)
rownames(annotation) <- annotation$label
annotation <- annotation[,-1,drop=F]
annotation$chr <- factor(annotation$chr,levels=paste0("chr",c(1:19,"X","Y")))
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
antibody <- "H3K36me3"
diff_summary <- data.frame()
# antibodys <- c("H3K9me3","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K27ac","ATAC")
antibodys <- c("RNA")
for(antibody in antibodys){
  for(tissue in tissues){
    search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
    if(antibody == "H3K9me3"){
      tab <- read.table(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_young_old_merge-W5000-G10000-E100_recursion.counts"),header = T)
      summary <- read.table(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_young_old_merge-W5000-G10000-E100_recursion.counts.summary"),header = T,row.names = 1)
    }else if(antibody == "RNA"){
      search_table <- read.csv("data/samples/all/RNA_search_table.csv")
      tab <- read.table(paste0("data/samples/RNA/",tissue,"/counts/",tissue,"_H3K9me3_peaks.counts"),header = T) 
      summary <- read.table(paste0("data/samples/RNA/",tissue,"/counts/",tissue,"_H3K9me3_peaks.counts.summary"),header = T)  
    }else if(antibody == "ATAC"){
      search_table <- read.csv("data/samples/all/ATAC_search_table_batch.csv")
      tab <- read.table(paste0("data/samples/ATAC/",tissue,"/ATAC/",tissue,"_H3K9me3_peaks.counts"),header = T) 
      summary <- read.table(paste0("data/samples/ATAC/",tissue,"/ATAC/",tissue,"_H3K9me3_peaks.counts.summary"),header=T)
    }else{
      tab <- read.table(paste0("data/samples/",tissue,"/",antibody,"/",tissue,"_H3K9me3_peaks.counts"),header = T)
      summary <- read.table(paste0("data/samples/",tissue,"/",antibody,"/",tissue,"_H3K9me3_peaks.counts.summary"),header = T)
    }
    if(tissue %in% c("mammarygland","uterus","ovary")){
      tab <- tab[which(tab$Chr %in% paste0("chr",c(1:19,"X"))),]
    }
    counts = tab[,c(7:ncol(tab))]
    rownames(counts)= tab$Geneid
    pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+|DYQ[0-9]+|HM[0-9]+).*"
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
    colnames(rpkm_mean_summary)[1] <- "Geneid"
    H3K9me3_diff <- read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_young_old_merge-W5000-G10000-E100_recursion_diff_after_remove_batch_effect.csv"))
    H3K9me3_diff <- H3K9me3_diff[,c("Geneid","Length","LogFC.old.young","Significant")]
    colnames(H3K9me3_diff)[4] <- "histone_significant"
    rpkm_mean_summary <- merge(H3K9me3_diff,rpkm_mean_summary,by="Geneid",all=T)
    colnames(rpkm_mean_summary)[c(3,7)] <- c("LogFC.old.young","logFC")
    rpkm_mean_summary$tissue <- tissue_label_change(tissue)
    rpkm_mean_summary$antibody <- antibody
    diff_summary <- rbind(diff_summary,rpkm_mean_summary)
  }
}
diff_summary <- diff_summary[which(diff_summary$Length > 200000),]



if(antibody == "RNA"){
  median_results <- diff_summary[,c("histone_significant", "tissue","antibody","mean_young","mean_old")]
  median_results <- diff_summary %>%
    group_by(histone_significant, tissue,antibody) %>%
    summarise(
      mean_young = mean(mean_young, na.rm = TRUE),
      mean_old = mean(mean_old, na.rm = TRUE)
    )
  p_value_summary <- data.frame()
  for(antibody in antibodys){
    for(condition in c("Up","Stable","Down")){
      t_df <- median_results[which(median_results$antibody=="RNA" & median_results$histone_significant==condition),]    
      test <- wilcox.test(t_df$mean_young,t_df$mean_old,paired = T)
      t_p_value_summary <- data.frame(antibody=antibody,condition=condition,p_value=test$p.value)  
      p_value_summary <- rbind(p_value_summary,t_p_value_summary)
    }
  }
  to_plot <- median_results
  colnames(to_plot)[4:5]<- c("young","old")
  to_plot <- reshape2::melt(to_plot)
  to_plot$variable <- factor(to_plot$variable, levels = c("young","old"))
  color <- setNames(c("#f39b7f","#4dbbd5"),c("young","old"))
  p <- ggplot(to_plot[-which(to_plot$histone_significant=="Stable"),], aes(x = histone_significant, y =value, fill = variable)) +  
    geom_boxplot(alpha = 0.7,outliers = F) +  
    scale_fill_manual(values = color) +
    labs(  
      title = paste0("Gene expression in H3K9me3 peaks"),  
      x = "Cluster",  
      y="RPKM"
    ) + 
    theme_bw() +  
    theme(  
      axis.title.x = element_text(size = 14),  
      axis.title.y = element_text(size = 14),  
      axis.text.x = element_text(size = 14),  
      axis.text.y = element_text(size = 14),  
      plot.title = element_text(size = 16, face = "bold")
    ) 
  ggsave(paste0("result/Sup_figures/H3K9me3_broad_recursion_peaks_RNA_change.pdf"),p,height = 8,width = 6)
}else{
  median_results <- diff_summary[,c("histone_significant", "tissue","antibody","mean_young","mean_old")]
  p_value_summary <- data.frame()
  for(antibody in antibodys){
    for(condition in c("Up","Stable","Down")){
      t_df <- median_results[which(median_results$antibody==antibody & median_results$histone_significant==condition),]    
      test <- t.test(t_df$mean_young,t_df$mean_old,paired = T)
      t_p_value_summary <- data.frame(antibody=antibody,condition=condition,p_value=test$p.value)  
      p_value_summary <- rbind(p_value_summary,t_p_value_summary)
    }
  }
  for(condition in c("Up","Stable","Down")){
    to_plot <- median_results[which(median_results$histone_significant==condition),]
    colnames(to_plot)[4:5]<- c("young","old")
    to_plot <- reshape2::melt(to_plot)
    to_plot$variable <- factor(to_plot$variable, levels = c("young","old"))
    color <- setNames(c("#f39b7f","#4dbbd5"),c("young","old"))
    to_plot$antibody <- factor(to_plot$antibody,levels = c("H3K9me3","H3K27me3","H3K36me3","H3K27ac","H3K4me1","H3K4me3","ATAC"))
    p <- ggplot(to_plot, aes(x = antibody, y =value, fill = variable)) +  
      geom_boxplot(alpha = 0.7,outliers = F) +  
      scale_fill_manual(values = color) +
      labs(  
        title = paste0(condition," H3K9me3 peaks"),  
        x = "Cluster",  
        y="RPKM"
      ) + 
      theme_bw() +  
      theme(  
        axis.title.x = element_text(size = 14),  
        axis.title.y = element_text(size = 14),  
        axis.text.x = element_text(size = 14),  
        axis.text.y = element_text(size = 14),  
        plot.title = element_text(size = 16, face = "bold")
      ) 
    ggsave(paste0("result/Sup_figures/H3K9me3_broad_recursion_peaks_",condition,"_other_modification_change.pdf"),p,height = 8,width = 6)
  }
}

### WGBS
tissue_summary <- data.frame()
ages <- c("3M","24M")
for(age in ages){
  for(tissue in tissues){
    annotation <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation_add_stable.csv",row.names = 1)
    split_names <- strsplit(annotation$label, "[:-]")
    annotation_df <- data.frame(
      chr = sapply(split_names, "[", 1),
      start = sapply(split_names, "[", 2),
      end = sapply(split_names, "[", 3),
      cluster = annotation$cluster
    )
    if(tissue %in% c("mammarygland","ovary","uterus")){
      annotation_df <- annotation_df[which(annotation_df$chr %in% paste0("chr",c(1:19,"X"))),]
    }
    annotation_df$start <- as.numeric(annotation_df$start)
    annotation_df$end <- as.numeric(annotation_df$end)
    annotation_df$label <- paste0(annotation_df$chr,":",annotation_df$start,"-",annotation_df$end)
    annotation_df <- as.data.table(annotation_df)
    setDT(annotation_df)
    setkey(annotation_df,chr,start,end)
    search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
    t_search_table <- search_table[which(search_table$tissue==tissue & search_table$age==age),]
    summary <- data.frame()
    for(sample in t_search_table$sample_name){
      df <- fread(paste0("data/samples/WGBS/",tissue,"/bdg/",sample,"_CpG.bdg"),sep = "\t")
      setDT(df)
      setkey(df,V1,V2,V3)  
      overlaps <- foverlaps(df, annotation_df, type = "any", nomatch = 0L)  
      
      result <- overlaps[, .(V4_sum = sum(V4), V5_sum = sum(V5)), by = label]
      result <- as.data.frame(result)
      result$methylation <- result$V4_sum/result$V5_sum * 100
      result <- result[,c("label","methylation")]
      colnames(result) <- c("label",sample)
      if(nrow(summary)==0){
        summary <- result
      }else{
        summary <- merge(summary,result,by="label")
      }
    }
    summary$mean <- rowMeans(summary[,-1])
    summary <- summary[,c("label","mean")]
    summary$tissue <- tissue_label_change(tissue)
    summary$age <- age
    H3K9me3_diff <- read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_young_old_merge-W5000-G10000-E100_recursion_diff_after_remove_batch_effect.csv"))
    H3K9me3_diff <- H3K9me3_diff[,c("Geneid","Length","LogFC.old.young","Significant")]
    colnames(H3K9me3_diff)[4] <- "histone_significant"
    summary <- merge(summary,H3K9me3_diff,by.x="label",by.y="Geneid")
    tissue_summary <- rbind(tissue_summary,summary)
  }
}
# median_results <- tissue_summary[,c("histone_significant", "tissue","mean","age")]
median_results <-  tissue_summary %>%
  group_by(tissue,age,histone_significant) %>%
  summarise(
    mean = mean(mean, na.rm = TRUE)
  )

p_value_summary <- data.frame()

for(condition in c("Up","Stable","Down")){
  t_df <- median_results[which(median_results$histone_significant==condition),]    
  # test <- wilcox.test(t_df$median_mean_young,t_df$median_mean_old)
  test <- t.test(t_df$mean[which(t_df$age=="3M")],t_df$mean[which(t_df$age=="24M")])
  t_p_value_summary <- data.frame(antibody="WGBS",condition=condition,p_value=test$p.value)  
  p_value_summary <- rbind(p_value_summary,t_p_value_summary)
}
median_results$age[which(median_results$age == "3M")] <- "young"
median_results$age[which(median_results$age == "24M")] <- "old"
ggplot(median_results, aes(x = histone_significant, y =mean, fill = age)) +  
  geom_boxplot(alpha = 0.7,outliers = F) +  
  scale_fill_manual(values = color) +
  labs(  
    title = paste0(condition," H3K9me3 peaks"),  
    x = "Cluster",  
    y="RPKM"
  ) + 
  theme_bw() +  
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 14),  
    axis.text.y = element_text(size = 14),  
    plot.title = element_text(size = 16, face = "bold")
  ) 
