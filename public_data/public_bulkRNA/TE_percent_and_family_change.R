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
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}

search_table <- read.csv("data/public_data/GSE132040/GSE132040_MACA_Bulk_metadata.csv")
colnames(search_table)[11] <- "sample"

tab <- read.delim(paste0("data/public_data/GSE132040/te_counts.txt"),row.names = 1)
counts <- tab
pattern <- ".*te\\.(SRR[0-9]+).*"
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
to_plot$source.name <- sub("_(\\w+)_\\d*|_(\\w+)$", "\\1",  to_plot$source.name)  

to_plot$label <- paste0(to_plot$sample,"-",to_plot$characteristics..age)
to_plot$characteristics..age <- as.numeric(to_plot$characteristics..age)  
to_plot <- to_plot[which(to_plot$characteristics..sex=="m"),]
tissues <- unique(to_plot$source.name)
p_list <- list()
for(i in c(1:length(tissues))){
  t_to_plot <- to_plot[which(to_plot$source.name==tissues[i]),]
  t_to_plot <- t_to_plot[which(t_to_plot$characteristics..age %in% c("3","24","27")),]
  t_to_plot <- t_to_plot[order(t_to_plot$characteristics..age), ]
  t_to_plot$label <- factor(t_to_plot$label, levels = unique(t_to_plot$label))
  t_to_plot_long <- t_to_plot %>%  
    select(sample,label, gene_sum_ratio, TE_sum_ratio) %>%  
    pivot_longer(cols = c("gene_sum_ratio", "TE_sum_ratio"),   
                 names_to = "Category",   
                 values_to = "Ratio")  
    t_to_plot_long$Ratio <- t_to_plot_long$Ratio * 100  
    t_to_plot_long$position <- 100
    t_to_plot_long$position[which(t_to_plot_long$Category=="TE_sum_ratio")] <- t_to_plot_long$Ratio[which(t_to_plot_long$Category=="TE_sum_ratio")]
    t_to_plot_long$Category_label <- "Gene"
    t_to_plot_long$Category_label[which(t_to_plot_long$Category=="TE_sum_ratio")] <- "TE" 
    color <- read.table("data/samples/7_distinct_color.txt")
    color <- setNames(color$V1[c(4,7)],unique(t_to_plot_long$Category))
    p_list[[i]] <- ggplot(t_to_plot_long, aes(x = label, y = Ratio, fill = Category)) +  
      geom_bar(stat = "identity") +  
      labs(x = "Sample", y = "Percentage", fill = "Category") +  
      theme_minimal() +  
      scale_fill_manual(values=color,labels = t_to_plot_long$Category_label)+
      geom_text(data = subset(t_to_plot_long, Category == "TE_sum_ratio"),   
                aes(label = paste0(round(Ratio, 2),"%"), y = position),   
                color = "white", size = 5, vjust = 0.5) + 
      theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
      ggtitle(tissues[i]) +
      labs(x = "", y = "Percentage", fill = "Category")
}
combined_plot <- plot_a_list(p_list,no_of_rows = 6,no_of_cols = 3)
ggsave("result/RNA/TE/public_data_TE_percent.png",combined_plot,width = 45,height = 36,limitsize = FALSE,type="cairo")






TE_CPMs <- as.data.frame(CPMs[!grepl("^ENSMUSG", rownames(CPMs)), ])
TE_CPMs$family <- rownames(TE_CPMs)
p_list <- list()
to_plot <- melt(TE_CPMs)
split_names <- strsplit(TE_CPMs$family, ":")  
split_df <- do.call(rbind, split_names) 
to_plot <- cbind(to_plot, split_df)
to_plot <- to_plot[,-1]
colnames(to_plot)[1] <- "sample"
colnames(to_plot)[c(3,5)] <- c("TE","family")
to_plot <- merge(to_plot,search_table, by = "sample")
to_plot$label <- paste0(to_plot$sample,"-",to_plot$characteristics..age)
to_plot$characteristics..age <- as.numeric(to_plot$characteristics..age)  
to_plot$source.name <- sub("_(\\w+)_\\d*|_(\\w+)$", "\\1",  to_plot$source.name)  
p_list <- list()
for(i in c(1:length(tissues))){
  t_to_plot <- to_plot[which(to_plot$source.name==tissues[i]),]
  t_to_plot <- t_to_plot[which(t_to_plot$characteristics..age %in% c("3","24","27")),]
  t_to_plot <- t_to_plot[order(t_to_plot$characteristics..age), ]
  t_to_plot$label <- factor(t_to_plot$label, levels = unique(t_to_plot$label))
  color <- read.table("data/samples/7_distinct_color.txt")
  t_to_plot$age <- t_to_plot$characteristics..age
  color <- data.frame(age=unique(t_to_plot$age),color=color$V1[1:length(unique(t_to_plot$age))])
  t_to_plot <- merge(t_to_plot,color,by="age")
  color <- unique(t_to_plot[,c("label","color")])
  color <- setNames(color$color,color$label)
  p_list[[i]] <- ggplot(t_to_plot, aes(x = family, y = value, fill=label)) +  
    geom_boxplot(outlier.shape = NA) +  
    scale_fill_manual(values=color)+
    labs(title = tissues[i],  
         x = "Category",  
         y = "CPM") +  
    theme_minimal() +
    ylim(0,5000)+
    theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))
  }
combined_plot <- plot_a_list(p_list, no_of_rows = 6,no_of_cols = 3)
ggsave("result/RNA/TE/public_data_TE_family_change_CPM_all_tissues.png.png",combined_plot,width = 45,height = 36,limitsize = FALSE,type="cairo")


