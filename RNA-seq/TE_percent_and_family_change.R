rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(tidyr) 
library(dplyr)
library(stringr)
library(reshape2)
tissue <- "mammarygland"
search_table <- read.csv("data/samples/all/RNA_search_table.csv")
colnames(search_table)[3] <- "sample"
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
tissue_label_change <- function(tissue){
  if(tissue=="FC"){
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
    }else if (tissue_label == "Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
} 
TE_percent <- function(tissue){
  p_list <- list()
  tab <- read.delim(paste0("data/samples/RNA/",tissue,"/TEcount/combined.cntTable"),row.names = 1)
  counts <- tab
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SRR[0-9]+).*"
  colnames(counts) <- gsub(pattern, "\\1", colnames(counts))
  y= DGEList(counts=counts)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  CPMs <-  cpm(y, log = F)
  ens_sums <- list()  
  non_ens_sums <- list()  
  for(col in colnames(CPMs)){
    ens_sum <- sum(CPMs[grep("^ENSMUSG", rownames(CPMs)), col])  
    non_ens_sum <- sum(CPMs[!grepl("^ENSMUSG", rownames(CPMs)), col])  
    ens_sums[[col]] <- ens_sum  
    non_ens_sums[[col]] <- non_ens_sum
  }
  ens_sums <- as.data.frame(ens_sums)
  non_ens_sums <- as.data.frame(non_ens_sums)
  rownames(ens_sums) <- "gene_sum"
  rownames(non_ens_sums) <- "TE_sum"
  to_plot <- rbind(ens_sums,non_ens_sums)
  to_plot <- as.data.frame(t(to_plot))
  to_plot$total_sum <- to_plot$gene_sum + to_plot$TE_sum  
  to_plot$gene_sum_ratio <- to_plot$gene_sum / to_plot$total_sum  
  to_plot$TE_sum_ratio <- to_plot$TE_sum / to_plot$total_sum  
  
  to_plot$sample <- rownames(to_plot)
  to_plot <- merge(to_plot, search_table, by="sample")  
  to_plot$label <- paste0(to_plot$sample,"-",to_plot$age,"-",to_plot$mouse_ID)
  to_plot$age <- factor(to_plot$age, levels=c("3m","24m"))
  to_plot <- to_plot[order(to_plot$age),]
  to_plot$label <- factor(to_plot$label, levels = unique(to_plot$label))
  to_plot_long <- to_plot %>%  
    select(sample,label, gene_sum_ratio, TE_sum_ratio) %>%  
    pivot_longer(cols = c("gene_sum_ratio", "TE_sum_ratio"),   
                 names_to = "Category",   
                 values_to = "Ratio")  
  to_plot_long$Ratio <- to_plot_long$Ratio * 100  
  
  to_plot_long$position <- 100
  to_plot_long$position[which(to_plot_long$Category=="TE_sum_ratio")] <- to_plot_long$Ratio[which(to_plot_long$Category=="TE_sum_ratio")]
  to_plot_long$Category_label <- "Gene"
  to_plot_long$Category_label[which(to_plot_long$Category=="TE_sum_ratio")] <- "TE"  
  color <- read.table("data/samples/7_distinct_color.txt")
  color <- setNames(color$V1[c(4,7)],unique(to_plot_long$Category))
  p_list[[1]] <-ggplot(to_plot_long, aes(x = label, y = Ratio, fill = Category)) +  
                      geom_bar(stat = "identity") +  
                      labs(x = "Sample", y = "Percentage", fill = "Category") +  
                      theme_minimal() +  
                      scale_fill_manual(values=color,labels = to_plot_long$Category_label)+
                      geom_text(data = subset(to_plot_long, Category == "TE_sum_ratio"),   
                                aes(label = paste0(round(Ratio, 2),"%"), y = position),   
                                color = "white", size = 5, vjust = 0.5) + 
                      theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
                      ggtitle(tissue_label_change(tissue)) +
                      labs(x = "", y = "Percentage", fill = "Category")
  
  TE_CPMs <- as.data.frame(CPMs[!grepl("^ENSMUSG", rownames(CPMs)), ])
  TE_CPMs$family <- rownames(TE_CPMs)
  to_plot <- melt(TE_CPMs)
  split_names <- strsplit(TE_CPMs$family, ":")  
  split_df <- do.call(rbind, split_names) 
  to_plot <- cbind(to_plot, split_df)
  to_plot <- to_plot[,-1]
  colnames(to_plot)[1] <- "sample"
  colnames(to_plot)[c(3,5)] <- c("TE","family")
  to_plot <- merge(to_plot,search_table, by = "sample")
  to_plot$age <- factor(to_plot$age, levels = c("3m","24m"))
  to_plot <- to_plot[order(to_plot$age),]
  to_plot$sample <- factor(to_plot$sample,levels = unique(to_plot$sample))
  to_plot$label <- paste0(to_plot$sample,"-",to_plot$age,"-",to_plot$mouse_ID)
  to_plot$label <-  factor(to_plot$label,levels = unique(to_plot$label))
  color <- read.table("data/samples/7_distinct_color.txt")
  color <- setNames(color$V1,unique(to_plot$label))
  p_list[[2]] <- ggplot(to_plot, aes(x = family, y = value, fill=label)) +  
    geom_boxplot(outlier.shape = NA) +  
    scale_fill_manual(values=color)+
    labs(title = tissue_label_change(tissue),  
         x = "Category",  
         y = "CPM") +  
    theme_minimal() +
    theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))
  
  
  logCPMs <-  cpm(y, log = T)
  TE_CPMs <- as.data.frame(logCPMs[!grepl("^ENSMUSG", rownames(logCPMs)), ])
  TE_CPMs$family <- rownames(TE_CPMs)
  to_plot <- melt(TE_CPMs)
  split_names <- strsplit(TE_CPMs$family, ":")  
  split_df <- do.call(rbind, split_names) 
  to_plot <- cbind(to_plot, split_df)
  to_plot <- to_plot[,-1]
  colnames(to_plot)[1] <- "sample"
  colnames(to_plot)[c(3,5)] <- c("TE","family")
  to_plot <- merge(to_plot,search_table, by = "sample")
  to_plot$age <- factor(to_plot$age, levels = c("3m","24m"))
  to_plot <- to_plot[order(to_plot$age),]
  to_plot$sample <- factor(to_plot$sample,levels = unique(to_plot$sample))
  to_plot$label <- paste0(to_plot$sample,"-",to_plot$age,"-",to_plot$mouse_ID)
  to_plot$label <-  factor(to_plot$label,levels = unique(to_plot$label))
  color <- read.table("data/samples/7_distinct_color.txt")
  color <- setNames(color$V1,unique(to_plot$label))
  p_list[[3]] <- ggplot(to_plot, aes(x = family, y = value, fill=label)) +  
                        geom_boxplot(outlier.shape = NA) +  
                        scale_fill_manual(values=color)+
                        labs(title = tissue_label_change(tissue),  
                             x = "Category",  
                             y = "log2(CPM)") +  
                        theme_minimal() +
                        theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p_list)
  }

tissues <- sort(c("skin","CB","spleen","heart","bladder","tongue","uterus","aorta","thymus","stomach","Hip","FC","BAT","iWAT","muscle","bonemarrow","lung","kidney","liver","testis","colon","cecum","ileum","jejunum","pancreas"))
percent_plot <- list()
CPM_plot <- list()
logCPM_plot <- list()
for(i in c(1:length(tissues))){
  p_list <- TE_percent(tissues[i])
  percent_plot[[i]] <- p_list[[1]]
  CPM_plot[[i]] <- p_list[[2]]
  logCPM_plot[[i]] <- p_list[[3]]
}
percent_plot_combined <- plot_a_list(percent_plot, no_of_rows = 5,no_of_cols = 5)
CPM_plot_combined <- plot_a_list(CPM_plot, no_of_rows = 5,no_of_cols = 5)
logCPM_plot_combined <- plot_a_list(logCPM_plot, no_of_rows = 5,no_of_cols = 5)
ggsave("result/RNA/TE/TE_percent_all_tissues.png",percent_plot_combined, width = 30,height = 25, type="cairo")
ggsave("result/RNA/TE/TE_family_change_CPM_all_tissues.png",CPM_plot_combined, width = 50,height = 40, type="cairo",limitsize = FALSE)
ggsave("result/RNA/TE/TE_family_change_logCPM_all_tissues.png",logCPM_plot_combined,width = 50,height = 40, type="cairo",limitsize = FALSE)

# search_table <- read.csv("data/public_data/GSE132040/GSE132040_MACA_Bulk_metadata.csv")
# colnames(search_table)[11] <- "sample"
# tab <- read.delim(paste0("data/public_data/GSE132040/Pancreas/TEcount/combined.cntTable"),row.names = 1)
# 
# to_plot$label <- paste0(to_plot$sample,"-",to_plot$characteristics..age)
# to_plot_long$characteristics..age <- as.numeric(to_plot_long$characteristics..age)  
# to_plot_long <- to_plot_long[order(to_plot_long$characteristics..age), ]
# levels=to_plot_long$label
# to_plot_long$label <- factor(to_plot_long$label,levels=unique(levels))
