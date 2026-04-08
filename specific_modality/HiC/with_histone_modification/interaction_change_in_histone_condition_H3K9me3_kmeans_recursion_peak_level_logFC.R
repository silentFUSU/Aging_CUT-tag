rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(tidyverse)  
library(dplyr) 
library(data.table)
library(ggsignif)
## kmeans1-kmeans1 kmeans2-kmeans2 kmeans3-kmeans3 kmeans4-kmeans4
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}

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
    }else if(tissue_label == "Mammarygland"){
      tissue_label <- "Mammary Gland"
    }else if(tissue_label == "Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
} 
tissue <- "lung"
antibody <- "H3K9me3"
resolution <- "200000"
kmeans <- "kmeans1"
interaction_change_in_histone_condition_peak_level_logFC <- function(tissue,resolution,kmeans){
  window_size <- "5000"
  gap_size <- "10000"
  histone <- read.table(paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/bed/",kmeans,"_uinon_recursion_peaks.bed"))
  # histone <- histone[which((histone$V3 - histone$V2 + 1)>200000),]
  histone$peaks <- paste0("peaks",c(1:nrow(histone)))
  histone$V2 <- histone$V2+1
  histone <- as.data.table(histone)
  setDT(histone)
  setkey(histone,V1,V2,V3)
  HiC_search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  HiC_search_table <- HiC_search_table[which(HiC_search_table$tissue==tissue),]
  
  bed <- read.table(paste0("data/samples/HiC/",tissue,"/raw_matrix/",HiC_search_table$sample_name[1],"_",resolution,"_abs.bed"))
  bed$V2 <- bed$V2+1
  bed <- as.data.table(bed)
  setDT(bed)
  setkey(bed,V1,V2,V3)
  overlaps <- foverlaps(histone, bed, type = "any", nomatch = 0L)  
  overlaps$condition <- "within_peaks"
  
  bed$condition <- "out_of_peaks"
  bed$condition[which(bed$V4 %in% overlaps$V4)] <- "within_peaks"
  summary <- data.frame()
  for(sample in HiC_search_table$sample_name){
    df <- fread(paste0("data/samples/HiC/",tissue,"/raw_matrix/",sample,"_",resolution,".matrix"))
    if(resolution == "200000"){
      df <- df %>%
        # dplyr::filter(abs(V1 - V2) >= 5 & abs(V1 - V2) <= 200)
        dplyr::filter(abs(V1 - V2) >= 5 & abs(V1 - V2) <= 300)
    }else if(resolution=="10000"){
      df <- df %>%
        # dplyr::filter(abs(V1 - V2) >= 20 & abs(V1 - V2) <= 600)
        dplyr::filter(abs(V1 - V2) >= 100 & abs(V1 - V2) <= 6000)
    }else{
      print("need change resolution")
    }
    df <- df %>%
      mutate(label = paste(V1, V2, sep = "-"))
    if(nrow(summary)==0){
      colnames(df)[3] <- sample
      summary <- df
    }else{
      df_subset <- df[, .(V3, label)]
      colnames(df_subset)[1]<-sample
      summary <- left_join(summary,df_subset,by="label")
    }
  }
  
  HiC_search_table$age <- factor(HiC_search_table$age,levels=c("3M","24M"))
  HiC_search_table <- HiC_search_table[order(HiC_search_table$age),]
  setcolorder(summary,c("V1","V2","label",HiC_search_table$sample_name))
  summary[is.na(summary)] <- 0
  sample_columns <- names(summary)[4:ncol(summary)]
  if(resolution == "10000"){
    summary[, (sample_columns) := lapply(.SD, function(x) x + 1), .SDcols = sample_columns]
  }
  CPM <- copy(summary)
  CPM[,(sample_columns):= lapply(.SD, function(column) {
    total_counts <- sum(column, na.rm = TRUE)
    cpm <- (column / total_counts) * 1e6
    return(cpm)
  }), .SDcols = sample_columns]
  
  young_samples <- HiC_search_table$sample_name[which(HiC_search_table$age=="3M")]
  old_samples <- HiC_search_table$sample_name[which(HiC_search_table$age=="24M")]
  young_CPM <- CPM[,c("V1","V2","label",young_samples), with = FALSE]
  old_CPM <- CPM[,c("V1","V2","label",old_samples),with=FALSE]
  young_CPM[, young_average := rowMeans(.SD, na.rm = TRUE), .SDcols = young_samples]
  old_CPM[, old_average := rowMeans(.SD, na.rm = TRUE), .SDcols = old_samples]
  young_CPM <- young_CPM[,.(V1,V2,label,young_average)]
  old_CPM <- old_CPM[,.(label,old_average)]
  CPM_average <- left_join(young_CPM,old_CPM,by="label")
  CPM_average <- copy(CPM_average)
  # CPM_average[, log2_ratio := log2((old_average + epsilon) / (young_average + epsilon))]
  if(resolution == "200000"){
    epsilon <- .Machine$double.eps
    CPM_average[, log2_ratio := log2((old_average + epsilon) / (young_average + epsilon))]
  }else if (resolution=="10000"){
    CPM_average[, log2_ratio := log2((old_average) / (young_average))]
  }else{
    print("need change resolution")
  }
  CPM_average[, `:=`(bin1_condition = "out", bin2_condition = "out")]
  CPM_average[, bin1_condition := ifelse(as.character(V1) %in% overlaps$V4, "within", bin1_condition)]
  CPM_average[, bin2_condition := ifelse(as.character(V2) %in% overlaps$V4, "within", bin2_condition)]
  CPM_average[, combined_condition := paste(bin1_condition, bin2_condition, sep = "-")]
  CPM_average[, combined_condition := ifelse(combined_condition == "out-within", "within-out", combined_condition)]
  to_plot <- as.data.frame(CPM_average[,.(combined_condition,log2_ratio)])
  if(resolution == "200000"){
    median_values <- to_plot %>%
      filter(is.finite(log2_ratio)) %>%
      group_by(combined_condition) %>%
      summarise(median_log2_ratio = median(log2_ratio, na.rm = TRUE))
    write.csv(median_values,paste0("result/HiC/",tissue,"/differential_analysis/with_histone/",tissue,"_",resolution,"_interaction_change_in_",antibody,"_histone_condition_",kmeans,"_recursion_peak_level_logFC.csv"),row.names = F)
  }else{
    mean_values <- to_plot %>%
      filter(is.finite(log2_ratio)) %>%
      group_by(combined_condition) %>%
      summarise(mean_log2_ratio = mean(log2_ratio, na.rm = TRUE))
    write.csv(mean_values,paste0("result/HiC/",tissue,"/differential_analysis/with_histone/",tissue,"_",resolution,"_interaction_change_in_",antibody,"_histone_condition_",kmeans,"_recursion_peak_level_logFC.csv"),row.names = F)
  }
  p <- ggplot(to_plot,aes(x=combined_condition,y=log2_ratio,color = combined_condition))+
    geom_boxplot(outliers = F) +
    ggtitle(tissue_label_change(tissue),paste0("Interaction change relationship with ",antibody," recursion peaks"))+
    theme_bw()+xlab("")+ylab("log2(Fold change)") + 
    theme(
      plot.title = element_text(size = 16, hjust = 0.5),
      axis.title.y = element_text(size = 14),
      axis.text = element_text(size = 12),
      strip.text = element_text(size = 14),
      legend.title = element_blank(),
      legend.text = element_text(size = 10)
    )
  return(p)
}
tissues <- c("brain","CB","kidney","liver","lung","bonemarrow","colon","heart","Hip","mammarygland","stomach","thymus")
for(kmeans in paste0("kmeans",1:4)){
  p_list <- list()
  for(tissue in tissues){
    p_list[[tissue]] <- interaction_change_in_histone_condition_peak_level_logFC(tissue,resolution,kmeans)
  }
  combined_plot <- plot_a_list(p_list,no_of_rows = 3, no_of_cols = 4)
  ggsave(paste0("result/HiC/all_tissues_",resolution,"_interaction_change_in_",antibody,"_histone_condition_",kmeans,"_recursion_peak_level_logFC.png"),combined_plot,width = 18,height = 20,type="cairo")
  
}

to_plot <- data.frame(
  kmeans1 = rep(NA, length(tissues)),
  kmeans2 = rep(NA, length(tissues)),
  kmeans3 = rep(NA, length(tissues)),
  kmeans4 = rep(NA, length(tissues))
)
rownames(to_plot) <- tissues
for(kmeans in paste0("kmeans",1:4)){
  for(tissue in tissues){
    df <- read.csv(paste0("result/HiC/",tissue,"/differential_analysis/with_histone/",tissue,"_",resolution,"_interaction_change_in_",antibody,"_histone_condition_",kmeans,"_recursion_peak_level_logFC.csv"))
    to_plot[tissue,kmeans] <- df$median_log2_ratio[which(df$combined_condition=="within-within")]   
  }
}
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- c(seq(-0.2, -0.06, length.out = 40), seq(-0.05, 0.05, length.out = 20), seq(0.06, 0.2, length.out = 40))
pheatmap::pheatmap(to_plot,breaks = breaks,color = color_palette,cluster_cols = F)
