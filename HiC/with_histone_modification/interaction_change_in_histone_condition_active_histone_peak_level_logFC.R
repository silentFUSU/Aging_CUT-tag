rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(tidyverse)  
library(dplyr) 
library(data.table)
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
antibody <- "H3K27ac"
resolution <- "10000"

interaction_change_in_histone_condition_peak_level <- function(tissue,resolution,antibody){
  histone <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_macs_young_old_narrowpeak.bed"))

  histone$peaks <- paste0("peaks",c(1:nrow(histone)))
  histone$length <- histone$V3 - histone$V2 +1
  ggplot(histone, aes(x = length)) +
    geom_histogram(binwidth = 100, fill = "blue", color = "black", alpha = 0.7) +
    labs(title = paste0("Histogram of ",tissue_label_change(tissue)," ",antibody," Length"), x = "Length", y = "Frequency") +
    theme_minimal()
  
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
    df <- df %>%
      dplyr::filter(abs(V1 - V2) >= 100 & abs(V1 - V2) <= 6000)
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
  summary[, (sample_columns) := lapply(.SD, function(x) x + 1), .SDcols = sample_columns]
  
  CPM <- copy(summary)
  CPM[,(sample_columns):= lapply(.SD, function(column) {
    total_counts <- sum(column, na.rm = TRUE)
    cpm <- (column / total_counts) * 1e6
    return(cpm)
  }), .SDcols = sample_columns]
  
  CPM <- CPM[rowSums(CPM[, ..sample_columns] >= 1) >= 1]
  
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
  CPM_average[, log2_ratio := log2((old_average) / (young_average))]
  CPM_average[, `:=`(bin1_condition = "out", bin2_condition = "out")]
  CPM_average[, bin1_condition := ifelse(as.character(V1) %in% overlaps$V4, "within", bin1_condition)]
  CPM_average[, bin2_condition := ifelse(as.character(V2) %in% overlaps$V4, "within", bin2_condition)]
  CPM_average[, combined_condition := paste(bin1_condition, bin2_condition, sep = "-")]
  CPM_average[, combined_condition := ifelse(combined_condition == "out-within", "within-out", combined_condition)]
  to_plot <- as.data.frame(CPM_average[,.(combined_condition,log2_ratio)])
  mean_values <- to_plot %>%
    filter(is.finite(log2_ratio)) %>%
    group_by(combined_condition) %>%
    summarise(mean_log2_ratio = mean(log2_ratio, na.rm = TRUE))
  write.csv(mean_values,paste0("result/HiC/",tissue,"/differential_analysis/with_histone/",tissue,"_",resolution,"_interaction_change_in_",antibody,"_histone_condition_peak_level_logFC.csv"),row.names = F)
  median_values <- to_plot %>%
    filter(is.finite(log2_ratio)) %>%
    group_by(combined_condition) %>%
    summarise(median_log2_ratio = median(log2_ratio, na.rm = TRUE))
  write.csv(median_values,paste0("result/HiC/",tissue,"/differential_analysis/with_histone/",tissue,"_",resolution,"_interaction_change_in_",antibody,"_histone_condition_peak_level_logFC.csv"),row.names = F)
  p <- ggplot(to_plot,aes(x=combined_condition,y=log2_ratio,color = combined_condition))+
    # geom_boxplot(outliers = F) +
    geom_boxplot() +
    ggtitle(tissue_label_change(tissue),paste0("Interaction change relationship with ",antibody))+
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
p_list <- list()
for(tissue in tissues){
  p_list[[tissue]] <- interaction_change_in_histone_condition_peak_level(tissue,resolution,antibody)
}
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}
combined_plot <- plot_a_list(p_list,no_of_rows = 3,no_of_cols = 4)
ggsave(paste0("result/HiC/all_tissues_",resolution,"_interaction_change_in_",antibody,"_histone_condition_peak_level.png"),combined_plot,width = 18,height = 20,type="cairo")
