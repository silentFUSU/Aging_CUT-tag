rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)  
tissue <- "ovary"

compress_to_peak <- function(tissue, antibody){
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  samples <- search_table$sample_name[which(search_table$tissue==tissue)]
  peak_summary <- data.frame()
  for(sample in samples){
    df <- fread(paste0("data/samples/WGBS/",tissue,"/bdg/",sample,"_CpG.bdg"),sep = "\t")
    df <- df[which(df$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
    df$V1 <- factor(df$V1,paste0("chr",c(1:19,"X","Y")))
    df$V3 <- df$V2
    df <- df[order(df$V1,df$V2),]
    if(antibody %in% c("H3K9me3","H3K27me3","H3K36me3")){
      ref <- fread(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_young_old_merge-W1000-G3000-E100.bed")) 
    }
    ref$V1 <- factor(ref$V1,paste0("chr",c(1:19,"X","Y")))
    ref <- ref[order(ref$V1,ref$V2),]
    ref$label <- paste0("peak",c(1:nrow(ref)))
    ref[, `:=`(total_V4 = NA, total_V5 = NA)]  
    
    setDT(df)  
    setDT(ref)  
    setkey(df, V1, V2, V3) 
    setkey(ref, V1, V2, V3) 
    overlaps <- foverlaps(df, ref, type = "any", nomatch = 0L)  
    results <- overlaps[, .(total_V4 = sum(V4), total_V5 = sum(V5)), by = .(label)]  
    results$percent <- results$total_V4/results$total_V5
    colnames(results)[1] <- "label"
    results$tissue <- tissue
    results$sample <- sample
    results <- merge(results,ref[,c(1:4)],by="label")
    peak_summary <- rbind(peak_summary,results)
    # t_df <- df %>% filter(V2 > "27570001" & V2 <= "27580000" & V1 == "chr1") #check
  }
  write.csv(peak_summary,paste0("data/samples/WGBS/",tissue,"/compress2bin/",antibody,"_peaks_all_depth.csv"),row.names = F)
}
compress_to_peak(tissue,"H3K27me3")
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
    }
  }
  return(tissue_label)
}
WGBS_change_in_histone_peak <- funcion(tissue,antibody){
  df <- read.csv(paste0("data/samples/WGBS/",tissue,"/compress2bin/",antibody,"_peaks_all_depth.csv"))
  colnames(df)[6] <- "sample_name"
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  df <- merge(df,search_table[,c(3,4,5)],by="sample_name")  
  df$age[which(df$age=="3M")] <- "young"
  df$age[which(df$age=="24M")] <- "old"
  df$age <- factor(df$age,levels=c("young","old"))
  df$percent <- 100*(df$total_V4/df$total_V5)
  t <- t.test(df$percent[which(df$age=="old")],df$percent[which(df$age=="young")])
  ggplot(df, aes(x = age, y = percent,fill=sample_name)) +  
    scale_fill_brewer(palette = "Pastel1") +
    geom_boxplot(outlier.shape = NA) +  
    theme_minimal() +
    theme(text = element_text(size = 20)) +
    labs(title = paste0(tissue_label_change(tissue),"\n",antibody," peak region"), x = NULL, y = "CG%")+
    annotate("text", x = Inf, y = -Inf, label = paste("p-value =",  format(t$p.value, scientific = TRUE, digits = 3)  ),   
             hjust = 1.1, vjust = -1.1, size = 5, colour = "red") +
    annotate("text", x = Inf, y = Inf, label = paste("old mean =",  round(t$estimate[[1]],2)  ),   
             hjust = 1.1, vjust = 1.1, size = 5, colour = "red") +
    annotate("text", x = -Inf, y = Inf, label = paste("young mean =",  round(t$estimate[[2]],2)  ),   
             hjust = 0, vjust = 1.1, size = 5, colour = "red")
}

WGBS_change_in_histone_bin <- function(tissue,antibody){
  if(antibody %in% c("H3K9me3","H3K27me3","H3K36me3")){
    bin_size <- "10kb"
  }else{
    bin_size <- "1kb"
  }
  df <- read.csv(paste0("data/samples/WGBS/",tissue,"/compress2bin/",bin_size,"_bins_all_depth.csv"))
  df <- df[which(df$total_V5>15),]
  colnames(df)[6] <- "sample_name"
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  df <- merge(df,search_table[,c(3,4,5)],by="sample_name")  
  df$age[which(df$age=="3M")] <- "young"
  df$age[which(df$age=="24M")] <- "old"
  df$age <- factor(df$age,levels=c("young","old"))
  df$percent <- 100*(df$total_V4/df$total_V5)
  peak <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_10kb_in_young_old_merge-W1000-G3000-E100.bed"))
  
  t <- t.test(df$percent[which(df$age=="old" & df$label %in% peak$V4)],df$percent[which(df$age=="young" & df$label %in% peak$V4)])
  ggplot(df[which(df$label %in% peak$V4),], aes(x = age, y = percent,fill=sample_name)) +  
    scale_fill_brewer(palette = "Pastel1") +
    geom_boxplot(outlier.shape = NA) +  
    theme_minimal() +
    theme(text = element_text(size = 20)) +
    labs(title = paste0(tissue_label_change(tissue),"\n","bins in ",antibody," peak region"), x = NULL, y = "CG%")+
    annotate("text", x = Inf, y = -Inf, label = paste("p-value =",  format(t$p.value, scientific = TRUE, digits = 3)  ),   
             hjust = 1.1, vjust = -1.1, size = 5, colour = "red") +
    annotate("text", x = Inf, y = Inf, label = paste("old mean =",  round(t$estimate[[1]],2)  ),   
             hjust = 1.1, vjust = 1.1, size = 5, colour = "red") +
    annotate("text", x = -Inf, y = Inf, label = paste("young mean =",  round(t$estimate[[2]],2)  ),   
             hjust = 0, vjust = 1.1, size = 5, colour = "red")
  
  t <- t.test(df$percent[which(df$age=="old" & !(df$label %in% peak$V4))],df$percent[which(df$age=="young" & !(df$label %in% peak$V4))])
  ggplot(df[which(!df$label %in% peak$V4),], aes(x = age, y = percent,fill=sample_name)) +  
    scale_fill_brewer(palette = "Pastel1") +
    geom_boxplot(outlier.shape = NA) +  
    theme_minimal() +
    theme(text = element_text(size = 20)) +
    labs(title = paste0(tissue_label_change(tissue),"\n","bins outside ",antibody," peak region"), x = NULL, y = "CG%")+
    annotate("text", x = Inf, y = -Inf, label = paste("p-value =",  format(t$p.value, scientific = TRUE, digits = 3)  ),   
             hjust = 1.1, vjust = -1.1, size = 5, colour = "red") +
    annotate("text", x = Inf, y = Inf, label = paste("old mean =",  round(t$estimate[[1]],2)  ),   
             hjust = 1.1, vjust = 1.1, size = 5, colour = "red") +
    annotate("text", x = -Inf, y = Inf, label = paste("young mean =",  round(t$estimate[[2]],2)  ),   
             hjust = 0, vjust = 1.1, size = 5, colour = "red")

  result <- df %>% 
    group_by(age,label) %>%
    summarise(sum_total_V4 = sum(total_V4),
              sum_total_V5 = sum(total_V5))
  result$percent <- result$sum_total_V4/result$sum_total_V5*100
  result <- result[,c("age","label","percent")]
  result_to_plot <- result %>%   
    pivot_wider(names_from = age, values_from = percent)
  result_to_plot <- na.omit(result_to_plot)
  result_to_plot$difference <- result_to_plot$old - result_to_plot$young
  result_to_plot$condition <- "peak region"
  result_to_plot$condition[which(! result_to_plot$label %in% peak$V4)] <- "outside region"
  result_to_plot$condition <- factor(result_to_plot$condition, levels = c("peak region","outside region"))
  t <- t.test(result_to_plot$difference[which(result_to_plot$condition=="peak region" & result_to_plot$difference <0)],result_to_plot$difference[which(result_to_plot$condition=="outside region" & result_to_plot$difference <0)])
  ggplot(result_to_plot[which(result_to_plot$difference < 0),], aes(x = condition, y = abs(difference),fill=condition)) +  
    scale_fill_brewer(palette = "Pastel1") +
    geom_boxplot() +  
    theme_minimal() +
    theme(text = element_text(size = 20)) +
    guides(fill = FALSE)+
    labs(title = paste0(tissue_label_change(tissue),"\n","Hypo bins difference"), x = NULL, y = "CG% difference")+
    annotate("text", x = Inf, y = -Inf, label = paste("p-value =",  format(t$p.value, scientific = TRUE, digits = 3)  ),   
             hjust = 1.1, vjust = -1.1, size = 5, colour = "red")
  }
