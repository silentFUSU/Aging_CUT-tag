rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(tidyverse)  
library(dplyr) 
library(limma)
library(data.table)
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
antibody <- "H3K27me3"
resolution <- "200000"
peak_overlap_summary <- data.frame()
interaction_change_in_histone_condition_peak_level <- function(tissue,resolution,antibody){
  window_size <- "5000"
  gap_size <- "10000"
  histone <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_young_old_merge-W",window_size,"-G",gap_size,"-E100.bed"))
  # histone <- read.table(paste0("data/samples/all/H3K9me3/bed/H3K9me3_young_old_merge-W5000-G10000-E100_recursion.bed"))
  # histone <- read.table(paste0("data/samples/",tissue,"/",antibody,"/peaks/edd/edd_peaks_fdr05.bed"))
  histone <- histone[which((histone$V3-histone$V2 +1)>200000),]
  histone$peaks <- paste0("peaks",c(1:nrow(histone)))
  histone <- as.data.table(histone)
  setDT(histone)
  setkey(histone,V1,V2,V3)
  re <-  read.table(paste0("data/samples/HiC/",tissue,"/differential_analysis/",tissue,"_",resolution,".FDR"))
  re$Significant <- "Stable"
  re$Significant[which((re$V5 < 0.05 & re$V6 < 0.05) & re$V4 < 0)] <- "Down"
  re$Significant[which((re$V5 < 0.05 & re$V6 < 0.05) & re$V4 > 0)] <- "Up"
  
  all_dynamic <- read.table(paste0("data/samples/HiC/",tissue,"/dynamics/",tissue,"_",resolution,".dynamics"))
  all_dynamic <- all_dynamic[which(abs(all_dynamic$V2 - all_dynamic$V3) > 4 & abs(all_dynamic$V2 - all_dynamic$V3) <= 200),c(1:3)]
  all_dynamic$Significant <- "Stable"
  all_dynamic$label <- paste(all_dynamic$V1,all_dynamic$V2,all_dynamic$V3,sep = "-")
  
  re_sig <- re[which(re$Significant!="Stable"),]
  re_sig <- re_sig[which(abs(re_sig$V2 - re_sig$V3)>4 & abs(re_sig$V2 - re_sig$V3) <= 200),]
  re_sig$label <- paste(re_sig$V1,re_sig$V2,re_sig$V3,sep = "-")
  all_dynamic$Significant[which(all_dynamic$label %in% re_sig$label[which(re_sig$Significant=="Up")])] <- "Up"
  all_dynamic$Significant[which(all_dynamic$label %in% re_sig$label[which(re_sig$Significant=="Down")])] <- "Down"
  
  re_sig <- all_dynamic[,c(1:4)]
  re_sig$V1 <- paste0("chr",re_sig$V1)
  re_sig$V1[which(re_sig$V1=="chr20")] <- "chrX"
   
  re_sig_df <- re_sig
  re_sig_df$label1 <- paste(re_sig_df$V1,re_sig_df$V2,sep = "-")
  re_sig_df$label2 <- paste(re_sig_df$V1,re_sig_df$V3,sep = "-")
  
  HiC_search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  HiC_search_table <- HiC_search_table[which(HiC_search_table$tissue==tissue),]
  bed <- read.table(paste0("data/samples/HiC/",tissue,"/raw_matrix/",HiC_search_table$sample_name[1],"_",resolution,"_abs.bed"))
  bed$V2 <- bed$V2+1
  count <- 1
  bed[1,"V5"] <- count
  for(i in c(2:nrow(bed))){
    if(bed[i,"V1"] != bed[i-1,"V1"]){
      count <- 1
      bed[i,"V5"] <- count
    }else{
      count <- count+1
      bed[i,"V5"] <- count
    }
  }
  bed <- as.data.table(bed)
  setDT(bed)
  setkey(bed,V1,V2,V3)
  overlaps <- foverlaps(histone, bed, type = "any", nomatch = 0L)  
  overlaps$condition <- "within_peaks"
  overlaps$label <- paste0(overlaps$V1,"-",overlaps$V5)
  bed$condition <- "out_of_peaks"
  bed$condition[which(bed$V4 %in% overlaps$V4)] <- "within_peaks"
  t_peak_overlap_summary <- as.data.frame(table(bed$condition))
  t_peak_overlap_summary$percent <- t_peak_overlap_summary$Freq/sum(t_peak_overlap_summary$Freq)*100
  t_peak_overlap_summary$tissue <- tissue_label_change(tissue)
  peak_overlap_summary <<- rbind(peak_overlap_summary,t_peak_overlap_summary)
  re_sig_df$condition_bin1 <- "out"
  re_sig_df$condition_bin1[which(re_sig_df$label1 %in% overlaps$label)] <- "within"
  re_sig_df$condition_bin2 <- "out"
  re_sig_df$condition_bin2[which(re_sig_df$label2 %in% overlaps$label)] <- "within"
  re_sig_df$condition <- paste0(re_sig_df$condition_bin1,"-",re_sig_df$condition_bin2)
  re_sig_df$condition[which(re_sig_df$condition=="out-within")] <- "within-out"
  result <- re_sig_df %>%
    group_by(condition, Significant) %>%
    summarize(count = n(), .groups = 'drop') %>%
    group_by(condition) %>%
    mutate(total = sum(count), proportion = count / total) %>%
    ungroup()
  result$proportion <- result$proportion*100
  result <- as.data.frame(result)
  result$Significant <- factor(result$Significant, levels = c("Up","Down","Stable"))
  color <- setNames(c("#fc5185","#00adb5","grey"),c("Up","Down","Stable"))
  p <- ggplot(result, aes(x = condition, y = proportion, fill = Significant)) +  
    geom_bar(stat = 'identity',colour = "white") +   
    scale_fill_manual(values = color) +
    theme_minimal() +   
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 20),legend.title = element_blank()) +
    ylab("Proportion")+
    ggtitle(paste0(tissue_label_change(tissue)," ",antibody))
  ggsave(paste0("result/HiC/",tissue,"/differential_analysis/with_histone/",tissue,"_",resolution,"_interaction_change_in_histone_condition_",antibody,"_peak_level.png"),p,width = 5,height = 7,type="cairo")
  write.csv(result,paste0("result/HiC/",tissue,"/differential_analysis/with_histone/",tissue,"_",resolution,"_interaction_change_in_histone_condition_",antibody,"_peak_level.csv"),row.names = F)
  # ggsave(paste0("result/HiC/",tissue,"/differential_analysis/with_histone/",tissue,"_",resolution,"_interaction_change_in_histone_condition_",antibody,"_domain_level.png"),p,width = 5,height = 7,type="cairo")
  # write.csv(result,paste0("result/HiC/",tissue,"/differential_analysis/with_histone/",tissue,"_",resolution,"_interaction_change_in_histone_condition_",antibody,"_domain_level.csv"),row.names = F)
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
ggsave(paste0("result/HiC/all_tissues_",resolution,"_interaction_change_in_",antibody,"_histone_condition_domain_level.png"),combined_plot,width = 18,height = 20,type="cairo")

ggplot(peak_overlap_summary, aes(x = tissue, y = percent, fill = Var1)) +  
  geom_bar(stat = 'identity',colour = "white") +   
  theme_minimal() +   
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Proportion")+
  ggtitle(paste0("HiC ",resolution," bins overlap with ",antibody," domains"))

summary <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("result/HiC/",tissue,"/differential_analysis/with_histone/",tissue,"_",resolution,"_interaction_change_in_histone_condition_peak_level.csv"))
  df$tissue <- tissue_label_change(tissue)
  summary <- rbind(summary,df)
  }

change_percent <- summary
change_percent$Significant[which(change_percent$Significant %in% c("Up","Down"))] <- "Changed"
result <- change_percent %>%
  filter(Significant == "Changed") %>%  
  group_by(tissue, condition) %>%      
  summarize(total_proportion = sum(proportion), .groups = 'drop')  
result <- as.data.frame(result)
to_plot <- reshape2::dcast(result, tissue ~ condition, value.var = "total_proportion")
rownames(to_plot) <- to_plot$tissue
to_plot <- to_plot[,-1]
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)  
pheatmap::pheatmap(to_plot,scale="row",cluster_cols = F,fontsize = 15, color = color_palette)

change_condition_foldchange <- summary
change_condition_foldchange <- change_condition_foldchange[which(change_condition_foldchange$Significant %in% c("Up","Down")),]
result <- change_condition_foldchange %>%
  filter(Significant %in% c("Down", "Up")) %>% 
  group_by(tissue, condition, Significant) %>%
  summarize(count = sum(count), .groups = 'drop') %>%
  spread(Significant, count) %>% 
  mutate(ratio = ifelse(is.na(Down) | is.na(Up), NA, Down / Up)) %>% 
  ungroup()
result <- as.data.frame(result)
to_plot <- reshape2::dcast(result, tissue ~ condition, value.var = "ratio")
rownames(to_plot) <- to_plot$tissue
to_plot <- to_plot[,-1]
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)  
breaks <- c(seq(0, 0.8, length.out = 40), seq(0.81, 1.2, length.out = 20), seq(1.21, 2, length.out = 40))  
pheatmap::pheatmap(to_plot,cluster_cols = F,fontsize = 15,breaks = breaks, color = color_palette)
