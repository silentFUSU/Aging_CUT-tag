rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
tissue <- "Hip"
antibody <- "H3K27me3"
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols, guides = "collect")
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
capitalize <- function(string) {  
  paste0(toupper(substr(string, 1, 1)), substring(string, 2))  
}  
age <- "3m"
interaction_histone_RPKM_three_type <- function(tissue,antibody,age){
  p_list <- list()
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    bin_size <- "10kb"
  }else{
    bin_size <- "1kb"
  }
  histone <- read.table(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins.counts"),header = T)
  counts <- histone[,c(7:ncol(histone))]
  rownames(counts)= histone$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff.csv")
  search_table <- search_table[which(search_table$tissue==tissue & search_table$age==age & search_table$antibody==antibody),]
  counts <- counts[,c(search_table$sample_name)]
  lengths_kb <- histone$Length/1000
  rpkm_matrix <- matrix(0, nrow = nrow(counts), ncol = ncol(counts))  
  rownames(rpkm_matrix) <- rownames(counts)  
  colnames(rpkm_matrix) <- colnames(counts) 
  total_mapped_reads <- colSums(counts)
  for (i in 1:ncol(counts)) {  
    rpkm_matrix[, i] <- (counts[, i] / lengths_kb) / (total_mapped_reads[i] / 1e6)  
  } 
  rpkm_matrix <- as.data.frame(rpkm_matrix)
  row_avg <- rowMeans(rpkm_matrix)  
  histone <- cbind(histone[,c(1:6)],row_avg)
  if(age=="3m"){
    age_label <- "young"
  }else{
    age_label <- "old"
  }
  if(antibody %in% c("H3K27me3","H3K36me3","H3K9me3")){
    peaks <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_10kb_in_",age_label,"_merge-W1000-G3000-E100.bed"))
  }else{
    peaks <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_1kb_in_",age_label,"_merge_macs_narrowpeak.bed"))
  }
  
  histone$peak_condition <- "Out of peaks"
  histone$peak_condition[which(histone$Geneid %in% peaks$V4)] <- "In peaks"
  bar_outpeak <- quantile(histone$row_avg[which(histone$peak_condition=="Out of peaks" & histone$row_avg !=0)], probs = c(0.75))
  bar_inpeak <- quantile(histone$row_avg[which(histone$peak_condition=="In peaks" & histone$row_avg !=0)], probs = c(0.25))
  bar <- c(histone$row_avg[which(histone$peak_condition=="Out of peaks" & histone$row_avg > bar_outpeak)],histone$row_avg[which(histone$peak_condition=="In peaks" & histone$row_avg < bar_inpeak & histone$row_avg >0)])
  bar <- median(bar)
  p_list[[1]] <- ggplot(histone[which(histone$row_avg !=0),], aes(x = peak_condition, y = log2(row_avg), fill = peak_condition)) + 
    geom_boxplot() + 
    labs(x = NULL,
         y = "log2(RPKM)") + 
    theme_minimal()+
    ggtitle(tissue_label_change(tissue),paste(antibody,"log2(RPKM)",capitalize(age_label)))+
    geom_hline(yintercept = log2(bar_outpeak), linetype = "dashed", color = "blue")+
    geom_hline(yintercept = log2(bar_inpeak), linetype = "dashed", color = "blue")+
    geom_hline(yintercept = log2(bar), linetype = "dashed", color = "red")
  # bar <- max(histone$row_avg[which(histone$peak_condition=="Out of peaks")])
 
  p_list[[2]] <- ggplot(histone[which(histone$row_avg!=0),], aes(x = log2(row_avg), fill = peak_condition)) +  
    geom_histogram(aes(y = ..density.. * 100),bins = 30, position = "identity", alpha = 0.6) +  
    labs(x = "log2(RPKM)", y = "Percentage") +  
    ggtitle(paste(tissue_label_change(tissue),antibody,"log2(RPKM)",capitalize(age_label)))+
    theme_minimal() +  
    scale_fill_brewer(palette = "Set1") +
    ylim(0,100) +
    geom_vline(xintercept = log2(bar), linetype = "dashed", color = "red")+
    geom_vline(xintercept = log2(bar_outpeak), linetype = "dashed", color = "blue")+
    geom_vline(xintercept = log2(bar_inpeak), linetype = "dashed", color = "blue")
  
  histone <- histone[,c("Chr","Start","End","row_avg")]
  histone$Start <- histone$Start+1
  histone <- as.data.table(histone)
  setDT(histone)
  setkey(histone,Chr,Start,End)
  
  resolution <- "200000"
  HiC_search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  HiC_search_table <- HiC_search_table[which(HiC_search_table$tissue==tissue),]
  bed <- read.table(paste0("data/samples/HiC/",tissue,"/raw_matrix/",HiC_search_table$sample_name[1],"_",resolution,"_abs.bed"))
  bed <- bed[which(bed$V1 %in% paste0("chr",c(1:19,"X"))),]
  bed$V5 <- "NA"
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
  bed$V2 <- bed$V2+1
  bed <- as.data.table(bed)
  setDT(bed)
  setkey(bed,V1,V2,V3)
  overlaps <- foverlaps(histone, bed, type = "any", nomatch = 0L)  
  overlaps_summary <- overlaps %>%  
    group_by(V4) %>%  
    summarise(Median_RPKM = median(row_avg), .groups = "drop")  
  overlaps_summary <- merge(bed,overlaps_summary,by="V4")
  
  re <-  read.table(paste0("data/samples/HiC/",tissue,"/differential_analysis/",tissue,"_",resolution,".FDR"))
  re$Significant <- "Stable"
  re$Significant[which((re$V5 < 0.01 & re$V6 < 0.01) & re$V4 < 0)] <- "Down"
  re$Significant[which((re$V5 < 0.01 & re$V6 < 0.01) & re$V4 > 0)] <- "Up"
  re_sig <- re[which(re$Significant!="Stable"),]
  
  re_sig <- re_sig[which(abs(re_sig$V2 - re_sig$V3)>4),]
  re_sig$V1 <- paste0("chr",re_sig$V1)
  re_sig$V1[which(re_sig$V1=="chr20")] <- "chrX"

  re_sig_df <- re_sig
  re_sig_df$label1 <- paste(re_sig_df$V1,re_sig_df$V2,sep = "-")
  re_sig_df$label2 <- paste(re_sig_df$V1,re_sig_df$V3,sep = "-")
  overlaps_summary$label <- paste(overlaps_summary$V1,overlaps_summary$V5,sep = "-")
  overlaps_summary$region <- paste(overlaps_summary$V1,overlaps_summary$V2,overlaps_summary$V3,sep = "-")
  re_sig_df <- merge(re_sig_df,overlaps_summary,by.x="label1",by.y="label")
  colnames(re_sig_df)[which(colnames(re_sig_df)=="region")] <- "region1"
  re_sig_df <- merge(re_sig_df,overlaps_summary,by.x="label2",by.y="label")
  colnames(re_sig_df)[which(colnames(re_sig_df)=="region")] <- "region2"
  colnames(re_sig_df)[which(colnames(re_sig_df)=="Median_RPKM.x")] <- "region1_Median_RPKM"
  colnames(re_sig_df)[which(colnames(re_sig_df)=="Median_RPKM.y")] <- "region2_Median_RPKM"
  
  df <- re_sig_df[,c("region1","region2","Significant")]
  df <- df %>%  
    separate(region1, into = c("chr1", "x1", "x2"), sep = "-", convert = TRUE)  
  df <- df %>%  
    separate(region2, into = c("chr2", "y1", "y2"), sep = "-", convert = TRUE)  
  write.table(df[which(df$Significant=="Up"),c(1:6)],paste0("data/samples/HiC/",tissue,"/differential_analysis/Wang_output_",resolution,"_increase.bedpe"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)
  write.table(df[which(df$Significant=="Down"),c(1:6)],paste0("data/samples/HiC/",tissue,"/differential_analysis/Wang_output_",resolution,"_decrease.bedpe"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)
  
  
  re_sig_reverse <- data.frame(V1=re_sig$V1,V2=re_sig$V3, V3=re_sig$V2, V4=re_sig$V4,V5=re_sig$V5,V6=re_sig$V6,Significant=re_sig$Significant)
  re_sig <- rbind(re_sig,re_sig_reverse)
  re_sig$label1 <- paste(re_sig$V1,re_sig$V2,sep = "-")
  re_sig$label2 <- paste(re_sig$V1,re_sig$V3,sep = "-")
  re_sig <- merge(re_sig,overlaps_summary,by.x="label1",by.y="label")
  colnames(re_sig)[which(colnames(re_sig)=="region")] <- "region1"
  re_sig <- merge(re_sig,overlaps_summary,by.x="label2",by.y="label")
  colnames(re_sig)[which(colnames(re_sig)=="region")] <- "region2"
  colnames(re_sig)[which(colnames(re_sig)=="Median_RPKM.x")] <- "region1_Median_RPKM"
  colnames(re_sig)[which(colnames(re_sig)=="Median_RPKM.y")] <- "region2_Median_RPKM"
  color <- setNames(c("red","blue"),c("Up","Down"))
  re_sig_plot <- sample_n(re_sig,min(50000,nrow(re_sig)))  
  p_list[[3]] <- ggplot(data=re_sig_plot, aes(log2(region1_Median_RPKM+1), log2(region2_Median_RPKM+1),color=Significant))+
    geom_point(size=1,alpha=0.2) +
    scale_color_manual(values = color)+
    geom_vline(xintercept = log2(bar_outpeak+1), linetype = "dashed", color = "black") +  
    geom_hline(yintercept = log2(bar_outpeak+1), linetype = "dashed", color = "black") + 
    theme_bw()+xlab("Enrichment in anchor1")+ylab("Enrichment in anchor2")+
    ggtitle(paste0(antibody," in ",capitalize(age_label))) + 
    theme(  
      plot.title = element_text(size=15, hjust=0.5), 
      axis.title = element_text(size=15),            
      axis.text = element_text(size=10),             
      legend.title = element_text(size=12),         
      legend.text = element_text(size=10)           
    ) +
    ylim(0,max(log2(re_sig_plot$region1_Median_RPKM +1))) + xlim(0,max(log2(re_sig_plot$region1_Median_RPKM +1)))
  
  type1 <- re_sig[which(re_sig$Significant != "Stable" & re_sig$region1_Median_RPKM > bar_outpeak & re_sig$region2_Median_RPKM > bar_outpeak),]
  type2 <- re_sig[which(re_sig$Significant != "Stable" & re_sig$region1_Median_RPKM > bar_outpeak & re_sig$region2_Median_RPKM < bar_outpeak),]
  type3 <- re_sig[which(re_sig$Significant != "Stable" & re_sig$region1_Median_RPKM < bar_outpeak & re_sig$region2_Median_RPKM < bar_outpeak),]
  
  type1 <- as.data.frame(table(type1$Significant))
  type1$condition <- "type1"
  type1$percent <- type1$Freq/sum(type1$Freq)*100
  type2 <- as.data.frame(table(type2$Significant))
  type2$condition <- "type2"
  type2$percent <- type2$Freq/sum(type2$Freq)*100
  type3 <- as.data.frame(table(type3$Significant))
  type3$condition <- "type3"
  type3$percent <- type3$Freq/sum(type3$Freq)*100
  to_plot <- rbind(type1,type2,type3)
  to_plot$position <- 100
  to_plot$position[which(to_plot$Var1 == "Up")] <- to_plot$percent[which(to_plot$Var1 == "Up")]
  color <- setNames(c("red","blue"),c("Up","Down"))
  
  p_list[[4]] <- ggplot(to_plot, aes(x = condition, y = percent, fill = Var1)) +  
    geom_bar(stat = 'identity',color="white") +   
    theme_minimal() +   
    scale_fill_manual(values = color) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 20),legend.title = element_blank()) +
    ylab("Proportion")+
    ggtitle(paste0(tissue_label_change(tissue)," ",antibody),paste0("Three ",capitalize(age_label)," types"))+
    geom_text(data = to_plot,   
              aes(label = Freq, y = position),   
              color = "black", size = 5, vjust = 0.5)
  return(p_list)
}

peak_rpkm_plot <- list()
peak_rpkm_distribution_plot <- list()
three_type_plot <- list()
dot_plot <- list()
tissues <- c("brain","CB","stomach","colon","lung","liver","thymus","heart","bonemarrow","kidney")
age <- "3m"
if(age=="3m"){
  age_label <- "Young"
}else{
  age_label <- "Old"
}
for(tissue in tissues){
  for(antibody in c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me1","H3K4me3")){
    p_list <- interaction_histone_RPKM_three_type(tissue,antibody,age)
    peak_rpkm_plot[[antibody]] <- p_list[[1]]
    peak_rpkm_distribution_plot[[antibody]] <- p_list[[2]]
    dot_plot[[antibody]] <- p_list[[3]]
    three_type_plot[[antibody]] <- p_list[[4]]
  }
  dir.create(paste0("result/HiC/",tissue,"/differential_analysis/with_histone/"),recursive = T,showWarnings = F)
  peak_rpkm_combined_plot <- plot_a_list(peak_rpkm_plot,2,3)
  ggsave(paste0("result/HiC/",tissue,"/differential_analysis/with_histone/",tissue,"_peak_rpkm_critical_point_",age_label,".png"),peak_rpkm_combined_plot,width = 10,height = 8,type="cairo")
  dot_combined_plot <-  plot_a_list(dot_plot,2,3)
  ggsave(paste0("result/HiC/",tissue,"/differential_analysis/with_histone/",tissue,"_histone_signal_rpkm_interaction_change_",age_label,".png"),dot_combined_plot,width = 12,height = 8,type="cairo")
  three_type_combined_plot <- plot_a_list(three_type_plot,2,3)
  ggsave(paste0("result/HiC/",tissue,"/differential_analysis/with_histone/",tissue,"_three_types_proportion_",age_label,".png"),three_type_combined_plot,width = 12,height = 12,type="cairo")
  peak_rpkm_distribution_combined_plot <- plot_a_list(peak_rpkm_distribution_plot,2,3)
  ggsave(paste0("result/HiC/",tissue,"/differential_analysis/with_histone/",tissue,"_peak_rpkm_distribution_",age_label,".png"),peak_rpkm_distribution_combined_plot,width = 12,height = 6,type="cairo")
}



interaction_histone_RPKM_logFC_three_type <- function(tissue,antibody){
  p_list <- list()
  # if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
  #   bin_size <- "10kb"
  # }else{
  #   bin_size <- "1kb"
  # }
  bin_size <- "200kb"
  histone <- read.table(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins.counts"),header = T)
  counts <- histone[,c(7:ncol(histone))]
  rownames(counts)= histone$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff.csv")
  search_table <- search_table[which(search_table$tissue==tissue &  search_table$antibody==antibody),]
  young <- counts[,c(search_table$sample_name[which(search_table$age=="3m")])]
  old <- counts[,c(search_table$sample_name[which(search_table$age=="24m")])]
  lengths_kb <- histone$Length/1000
  young_rpkm_matrix <- matrix(0, nrow = nrow(young), ncol = ncol(young))  
  old_rpkm_matrix <- matrix(0,nrow = nrow(old), ncol = ncol(old))
  rownames(young_rpkm_matrix) <- rownames(young)  
  colnames(young_rpkm_matrix) <- colnames(young) 
  rownames(old_rpkm_matrix) <- rownames(old)  
  colnames(old_rpkm_matrix) <- colnames(old) 
  
  total_young_mapped_reads <- colSums(young)
  for (i in 1:ncol(young)) {  
    young_rpkm_matrix[, i] <- (young[, i] / lengths_kb) / (total_young_mapped_reads[i] / 1e6)  
  } 
  total_old_mapped_reads <- colSums(old)
  for (i in 1:ncol(old)) {  
    old_rpkm_matrix[, i] <- (old[, i] / lengths_kb) / (total_old_mapped_reads[i] / 1e6)  
  } 
  young_rpkm_matrix <- as.data.frame(young_rpkm_matrix)
  old_rpkm_matrix <- as.data.frame(old_rpkm_matrix)
  
  logFC <- as.data.frame(log2(rowMeans(old_rpkm_matrix)/rowMeans(young_rpkm_matrix)))
  histone <- cbind(histone[,c(1:6)],logFC)
  colnames(histone)[7] <- "log2FC"
  histone$Start <- histone$Start+1
  histone <- as.data.table(histone)
  setDT(histone)
  setkey(histone,Chr,Start,End)
  
  resolution <- "200000"
  HiC_search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  HiC_search_table <- HiC_search_table[which(HiC_search_table$tissue==tissue),]
  bed <- read.table(paste0("data/samples/HiC/",tissue,"/raw_matrix/",HiC_search_table$sample_name[1],"_",resolution,"_abs.bed"))
  bed <- bed[which(bed$V1 %in% paste0("chr",c(1:19,"X"))),]
  bed$V5 <- "NA"
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
  bed$V2 <- bed$V2+1
  bed <- as.data.table(bed)
  setDT(bed)
  setkey(bed,V1,V2,V3)
  
  overlaps <- foverlaps(histone, bed, type = "any", nomatch = 0L)  
  
  overlaps_summary <- overlaps %>%  
    group_by(V4) %>%  
    summarise(Median_logFC = median(log2FC), .groups = "drop")  
  overlaps_summary <- merge(bed,overlaps_summary,by="V4")
  overlaps_summary <- overlaps_summary[!is.na(overlaps_summary$Median_logFC),]
  
  re <-  read.table(paste0("data/samples/HiC/",tissue,"/differential_analysis/",tissue,"_",resolution,".FDR"))
  re$Significant <- "Stable"
  re$Significant[which((re$V5 < 0.01 & re$V6 < 0.01) & re$V4 < 0)] <- "Down"
  re$Significant[which((re$V5 < 0.01 & re$V6 < 0.01) & re$V4 > 0)] <- "Up"
  re_sig <- re[which(re$Significant!="Stable"),]
  
  re_sig <- re_sig[which(abs(re_sig$V2 - re_sig$V3)>4),]
  re_sig$V1 <- paste0("chr",re_sig$V1)
  re_sig$V1[which(re_sig$V1=="chr20")] <- "chrX"
  re_sig_reverse <- data.frame(V1=re_sig$V1,V2=re_sig$V3, V3=re_sig$V2, V4=re_sig$V4,V5=re_sig$V5,V6=re_sig$V6,Significant=re_sig$Significant)
  re_sig <- rbind(re_sig,re_sig_reverse)
  
  re_sig$label1 <- paste(re_sig$V1,re_sig$V2,sep = "-")
  re_sig$label2 <- paste(re_sig$V1,re_sig$V3,sep = "-")
  
  overlaps_summary$label <- paste(overlaps_summary$V1,overlaps_summary$V5,sep = "-")
  re_sig <- merge(re_sig,overlaps_summary,by.x="label1",by.y="label")
  colnames(re_sig)[which(colnames(re_sig)=="region")] <- "region1"
  re_sig <- merge(re_sig,overlaps_summary,by.x="label2",by.y="label")
  colnames(re_sig)[which(colnames(re_sig)=="region")] <- "region2"
  colnames(re_sig)[which(colnames(re_sig)=="Median_logFC.x")] <- "region1_median_logFC"
  colnames(re_sig)[which(colnames(re_sig)=="Median_logFC.y")] <- "region2_median_logFC"
  color <- setNames(c("red","blue"),c("Up","Down"))
  re_sig_plot <- re_sig
  p_list[[1]] <- ggplot(data=re_sig_plot, aes(region1_median_logFC, region2_median_logFC,color=Significant))+
    geom_point(size=1,alpha=0.2) +
    scale_color_manual(values = color)+
    # geom_vline(xintercept = 1.5, linetype = "dashed", color = "black") +  
    # geom_hline(yintercept = 1.5, linetype = "dashed", color = "black") + 
    theme_bw()+xlab("log2(Fold change) in anchor1")+ylab("log2(Fold change) in anchor2")+
    ggtitle(paste0(tissue_label_change(tissue)," ",antibody," log2(Fold change)")) + 
    theme(  
      plot.title = element_text(size=15, hjust=0.5), 
      axis.title = element_text(size=15),            
      axis.text = element_text(size=10),             
      legend.title = element_text(size=12),         
      legend.text = element_text(size=10)           
    ) +
    coord_cartesian(xlim = c(-max(abs(re_sig_plot$region1_median_logFC)), max(abs(re_sig_plot$region1_median_logFC))),  
                    ylim = c(-max(abs(re_sig_plot$region2_median_logFC)), max(abs(re_sig_plot$region2_median_logFC))))+
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") 
  
  type1 <- re_sig[which(re_sig$Significant != "Stable" & re_sig$region1_median_logFC>0 & re_sig$region2_median_logFC>0),]
  type2 <- re_sig[which(re_sig$Significant != "Stable" & re_sig$region1_median_logFC>0 & re_sig$region2_median_logFC<0),]
  type3 <- re_sig[which(re_sig$Significant != "Stable" & re_sig$region1_median_logFC<0 & re_sig$region2_median_logFC<0),]
  type1 <- as.data.frame(table(type1$Significant))
  type1$condition <- "type1"
  type1$percent <- type1$Freq/sum(type1$Freq)*100
  type2 <- as.data.frame(table(type2$Significant))
  type2$condition <- "type2"
  type2$percent <- type2$Freq/sum(type2$Freq)*100
  type3 <- as.data.frame(table(type3$Significant))
  type3$condition <- "type3"
  type3$percent <- type3$Freq/sum(type3$Freq)*100
  to_plot <- rbind(type1,type2,type3)
  color <- setNames(c("red","blue"),c("Up","Down"))
  to_plot$position <- 100
  to_plot$position[which(to_plot$Var1 == "Up")] <- to_plot$percent[which(to_plot$Var1 == "Up")]
  p_list[[2]] <- ggplot(to_plot, aes(x = condition, y = percent, fill = Var1)) +  
    geom_bar(stat = 'identity',color="white") +   
    theme_minimal() +   
    scale_fill_manual(values = color) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 20),legend.title = element_blank()) +
    ylab("Proportion")+
    ggtitle(paste0(tissue_label_change(tissue)," ",antibody),"Three types")+
    geom_text(data = to_plot,   
              aes(label = Freq, y = position),   
              color = "black", size = 5, vjust = 0.5)
  return(p_list)
  }

plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}
tissues <- c("brain","stomach","colon","kidney","liver","bonemarrow","thymus")
for(tissue in tissues){
  three_type_plot <- list()
  dot_plot <- list()
  for(antibody in c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me3","H3K4me1")){
    p_list <- interaction_histone_RPKM_logFC_three_type(tissue,antibody)
    dot_plot[[antibody]] <- p_list[[1]]
    three_type_plot[[antibody]] <- p_list[[2]]
  }
  dot_combined_plot <-  plot_a_list(dot_plot,2,3)
  ggsave(paste0("result/HiC/",tissue,"/differential_analysis/with_histone/",tissue,"_histone_signal_200kb_rpkm_logFC_interaction_change.png"),dot_combined_plot,width = 12,height = 8,type="cairo")
  three_type_combined_plot <- plot_a_list(three_type_plot,2,3)
  ggsave(paste0("result/HiC/",tissue,"/differential_analysis/with_histone/",tissue,"_three_types_proportion_200kb_rpkm_logFC_interaction_change.png"),three_type_combined_plot,width = 12,height = 12,type="cairo")
}

