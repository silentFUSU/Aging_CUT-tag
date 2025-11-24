rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(bitmapType="cairo")  
library(grid)
library(stringr)
library(dplyr)
library(ggplot2)
library(data.table)
library(edgeR)
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
tissues <-  c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
              "thymus","skin","bladder","bonemarrow","Hip","heart",
              "muscle","jejunum","uterus","ovary","liver","tongue",
              "cecum","colon","testis","stomach","pancreas","iWAT","ileum")
antibodys <-c("H3K9me3","H3K27me3","H3K36me3","H3K27ac","H3K4me1","H3K4me3","ATAC")
for(kmeans in c("kmeans1","kmeans3","kmeans2","kmeans4")){
  print(kmeans)
  tissue_summary <- data.frame()
  print(nrow(tissue_summary))
  for(antibody in antibodys){
    print(antibody)
    for(tissue in tissues){ 
      print(tissue)
      if(antibody=="ATAC"){
        tab <- read.table(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_TE_all.counts"),header = T)
        tab_summary <- read.table(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_TE_all.counts.summary"),header = T)
        search_table <- read.csv("data/samples/all/ATAC_search_table_batch.csv")
      }else{
        tab <- read.table(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_TE_all.counts"),header = T)
        tab_summary <- read.table(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_TE_all.counts.summary"),header = T)
        search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
      }
      tab$Start[which(tab$Geneid=="chr9:119078473-119078571")] <- 119078473
      tab$End[which(tab$Geneid=="chr9:119078473-119078571")] <- 119078571
      regions <- read.table(paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/bed/",kmeans,"_uinon_recursion_peaks.bed"))
      regions <- as.data.table(regions)
      setDT(regions)
      setkey(regions,V1,V2,V3)
      TE <- tab[,c("Geneid","Chr","Start","End")]
      TE$Start <- as.numeric(TE$Start)
      TE$End <- as.numeric(TE$End)
      setDT(TE)
      setkey(TE,Chr,Start,End)
      overlaps <- foverlaps(regions, TE, type = "any", nomatch = 0L)  
      tab <- tab[which(tab$Geneid %in% overlaps$Geneid),]
      if(tissue %in% c("mammarygland","ovary","uterus")){
        tab <- tab[which(tab$Chr %in% paste0("chr",c(1:19,"X"))),]
      }
      t_search_table <- search_table[which(search_table$tissue==tissue & search_table$antibody==antibody),]
      rownames(tab) <- tab$Geneid
      counts = tab[,c(7:ncol(tab))]
      pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
      colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
      colnames(tab_summary) <- gsub(pattern,"\\1",colnames(tab_summary))
      counts <- counts[,t_search_table$sample_name]
      tab_summary <- tab_summary[,t_search_table$sample_name]
      
      keep = which(rowSums(cpm(counts)>1)>=2)
      counts <- counts[keep,]
      tab <- tab[rownames(counts),]
      t_search_table <- t_search_table[which(t_search_table$sample_name %in% colnames(counts)),]
      counts <- counts[,t_search_table$sample_name]
      tab_summary <- tab_summary[-2,t_search_table$sample_name]
      
      length_kb <- tab$Length / 1000  
      total_reads <- colSums(tab_summary)
      total_reads_million <- total_reads / 1e6  
      for (i in c(1:ncol(counts))) {  
        counts[[i]] <- (counts[[i]] / (length_kb * total_reads_million[i]))  
      }  
      
      RPKM <- counts
      
      young_cols <- RPKM[, t_search_table$age=="3m"]
      young_cols$rowmeans <- rowMeans(young_cols)
      old_cols <- RPKM[, t_search_table$age=="24m"]
      old_cols$rowmeans <- rowMeans(old_cols)
      
      young_cols$label <- rownames(young_cols)
      colnames(young_cols)[which(colnames(young_cols)=="rowmeans")] <- "young"
      young_cols <- young_cols[,c("label","young")]
      
      old_cols$label <- rownames(old_cols)
      colnames(old_cols)[which(colnames(old_cols)=="rowmeans")] <- "old"
      old_cols <- old_cols[,c("label","old")]
      t_tissue_summary <- merge(young_cols,old_cols,by="label")
      t_tissue_summary$logFC <- log2(t_tissue_summary$old/t_tissue_summary$young)
      t_tissue_summary$logFC[is.na(t_tissue_summary$logFC)] <- 0
      max_value <- max(t_tissue_summary$logFC[is.finite(t_tissue_summary$logFC)], na.rm = TRUE)
      min_value <- min(t_tissue_summary$logFC[is.finite(t_tissue_summary$logFC)], na.rm = TRUE)
      
      t_tissue_summary$logFC[t_tissue_summary$logFC == Inf] <- max_value
      t_tissue_summary$logFC[t_tissue_summary$logFC == -Inf] <- min_value
      # t_tissue_summary <- data.frame(tissue=tissue_label_change(tissue),logFC=median(t_tissue_summary$logFC,na.rm = T),antibody=antibody)
      t_tissue_summary$tissue <- tissue_label_change(tissue)
      t_tissue_summary$antibody <- antibody
      tissue_summary <- rbind(tissue_summary,t_tissue_summary)
    }  
  }
  
  to_plot <- tissue_summary
  tissue_order <- c("Lung","Cerebellum","BAT","Muscle","Heart","Aorta","Skin","Kidney","Hippocampus","Cortex","Liver","Tongue","Uterus","Testis","Bladder","Ovary",
                    "Colon","Stomach","Thymus","Cecum","Jejunum","Pancreas","Bone Marrow","Ileum","Spleen","iWAT","Mammary Gland")
  
  p_value_summary <- data.frame(matrix(ncol = length(antibodys), nrow = 27))
  rownames(p_value_summary) <- tissue_order
  colnames(p_value_summary) <- antibodys
  for(tissue in tissue_order){
    for(antibody in antibodys){
      t_df <- to_plot[which(to_plot$tissue==tissue & to_plot$antibody==antibody),]
      test <- wilcox.test(t_df$old, t_df$young,paired = T)
      p_value_summary[tissue,antibody] <- test$p.value
    }
  }
  mark_significance <- function(p_value) {
    if (is.na(p_value)) {
      return(NA)
    } else if (p_value < 0.001) {
      return("***")
    } else if (p_value < 0.01) {
      return("**")
    } else if (p_value < 0.05) {
      return("*")
    } else {
      return(NA)
    }
  }
  p_value_summary <- as.data.frame(lapply(p_value_summary, function(column) {
    sapply(column, mark_significance)
  }))
  rownames(p_value_summary) <- tissue_order
  p_value_summary$tissue <- rownames(p_value_summary)
  p_value_long <- reshape2::melt(p_value_summary,id.vars = "tissue")
  colnames(p_value_long)[c(2,3)] <- c("antibody","Label")
  to_plot_avg <- to_plot %>%
    group_by(antibody,tissue) %>%
    summarise(logFC = median(logFC, na.rm = TRUE))
  to_plot_avg$tissue <- factor(to_plot_avg$tissue,levels=tissue_order)
  
  merged_data <- merge(to_plot_avg, p_value_long, by = c("tissue", "antibody"), all.x = TRUE)
  merged_data$logFC[which(merged_data$logFC > 1)] <- 1
  merged_data$logFC[which(merged_data$logFC < -1)] <- -1
  merged_data$antibody <- factor(merged_data$antibody,levels = c("H3K9me3","H3K27me3","H3K36me3","H3K27ac","H3K4me1","H3K4me3","ATAC"))
  p <- ggplot(merged_data, aes(x = antibody, y = tissue, fill = logFC)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",limits = c(-1, 1), midpoint = 0) +
    theme_minimal() +
    geom_text(aes(label = Label), color = "black", size = 4, na.rm = TRUE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p
  ggsave(paste0("result/Sup_figures/TE_",kmeans,"_histone_change.pdf"),p,width = 6,height = 8)
}


# 
# to_plot <- reshape2::dcast(to_plot,tissue~antibody,value.var = "logFC")
# rownames(to_plot) <- to_plot$tissue  
# to_plot <- to_plot[,-1]
# color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
# breaks <- c(seq(-1, -0.29, length.out = 40), seq(-0.3, 0.3, length.out = 20), seq(0.31, 1, length.out = 40))
# tissue_order <- c("Lung","Cerebellum","BAT","Muscle","Heart","Aorta","Skin","Kidney","Hippocampus","Cortex","Liver","Tongue","Uterus","Testis","Bladder","Ovary",
#                   "Colon","Stomach","Thymus","Cecum","Jejunum","Pancreas","Bone Marrow","Ileum","Spleen","iWAT","Mammary Gland")
# to_plot <- to_plot[tissue_order,c("H3K9me3","H3K27me3","H3K36me3","H3K27ac","H3K4me1","H3K4me3","ATAC")]
# pheatmap::pheatmap(to_plot,breaks = breaks,color = color_palette,cluster_rows = F,cluster_cols = F,filename = "result/Sup_figures/TE_kmeans1_histone_change.pdf",width = 6,height = 8)

