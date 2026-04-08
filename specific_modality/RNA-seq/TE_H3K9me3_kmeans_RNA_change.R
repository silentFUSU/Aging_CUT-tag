rm(list=ls()) 
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)
library(rtracklayer)
library(gridExtra)
library(grid)  
library(GenomicRanges)
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggsignif)
library(edgeR)
options(bitmapType="cairo")  

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

gtf <- import("~/ref_data/TE_reference/mm10_rmsk_TE.gtf", format = "gtf")
family_data <- as.data.frame(gtf[,c("gene_id","transcript_id","family_id","class_id")])
family_data <- as.data.table(family_data)
setDT(family_data)
setkey(family_data,seqnames,start,end)
kmean <- "kmeans1"
regions <- read.table(paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/bed/",kmean,"_uinon_recursion_peaks.bed"))
regions <- as.data.table(regions)
setDT(regions)
setkey(regions,V1,V2,V3)
overlaps <- foverlaps(family_data, regions, type = "any", nomatch = 0L)  
overlaps$label <- paste(overlaps$transcript_id,overlaps$gene_id,overlaps$family_id,overlaps$class_id,sep=":")
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")

TE_info <- as.data.frame(gtf[,c("gene_id","transcript_id","family_id","class_id")])
TE_info$label <- paste(TE_info$transcript_id,TE_info$gene_id,TE_info$family_id,TE_info$class_id,sep=":")
TE_info <- TE_info[,c("label","width","gene_id")]
tissue_summary <- data.frame()
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/RNA_search_table.csv")
  search_table <- search_table[which(search_table$tissue_label==tissue),]
  tab <- read.table(paste0("data/samples/RNA/",tissue,"/TElocal/combined.cntTable"),header = T,row.names = 1)  
  tab_summary <- read.table(paste0("data/samples/RNA/",tissue,"/counts/",tissue,"_H3K9me3_peaks.counts.summary"),header = T)
  
  tab <- tab[which(rownames(tab)%in% overlaps$label),]
  counts <- tab
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SRR[0-9]+|HM[0-9]+).*"
  colnames(counts) <- gsub(pattern, "\\1", colnames(counts))  
  colnames(tab_summary) <- gsub(pattern,"\\1",colnames(tab_summary))
  counts <- counts[,search_table$sample_name]
  tab_summary <- tab_summary[-2,search_table$sample_name]
  
  counts <- merge(counts,TE_info,by.x="row.names",by.y="label")
  rownames(counts) <- counts$Row.names
  counts <- counts[,-1]
  
  length <- counts %>%
    group_by(gene_id) %>%
    summarise(length = sum(width))
  length <- as.data.frame(length)
  
  counts_summary <- counts %>%
    group_by(gene_id) %>%
    summarise(
      across(.cols = 1:ncol(tab), sum, .names = "{.col}")
    )
  
  counts_summary <- as.data.frame(counts_summary)
  rownames(counts_summary) <- counts_summary$gene_id
  counts_summary <- counts_summary[,-1]
  keep = which(rowSums(cpm(counts_summary)>0)>=2)
  counts_summary = counts_summary[keep,]
  
  rownames(length) <- length$gene_id
  length <- length[rownames(counts_summary),]
  length_kb <- length$length / 1000  
  total_reads <- colSums(tab_summary)
  total_reads_million <- total_reads / 1e6  
  for (i in c(1:ncol(counts_summary))) {  
    counts_summary[[i]] <- (counts_summary[[i]] / (length_kb * total_reads_million[i]))  
  }  
  RPKM <- counts_summary
  
  young_cols <- RPKM[, search_table$sample_name[which(search_table$age=="3m")]]
  young_cols$rowmeans <- rowMeans(young_cols)
  old_cols <- RPKM[, search_table$sample_name[which(search_table$age=="24m")]]
  old_cols$rowmeans <- rowMeans(old_cols)
  
  young_cols$label <- rownames(young_cols)
  young_cols$age <- "young"
  young_cols <- young_cols[,c("label","rowmeans","age")]
  
  old_cols$label <- rownames(old_cols)
  old_cols$age <- "old"
  old_cols <- old_cols[,c("label","rowmeans","age")]
  t_tissue_summary <- rbind(young_cols,old_cols)
  t_tissue_summary$tissue <- tissue_label_change(tissue)
  tissue_summary <- rbind(tissue_summary,t_tissue_summary)
  }

t.test(tissue_summary$rowmeans[which(tissue_summary$age=="young")],tissue_summary$rowmeans[which(tissue_summary$age=="old")])
to_plot <- tissue_summary
TE_info <- as.data.frame(family_data[,c("gene_id","family_id","class_id")])
TE_info <-  unique(TE_info)
to_plot <- merge(to_plot,TE_info[,c("gene_id","family_id","class_id")],by.x="label",by.y="gene_id")

to_plot_avg <- to_plot %>%
  group_by(class_id,tissue, age) %>%
  summarise(rowmeans = mean(rowmeans, na.rm = TRUE))
to_plot_avg <- to_plot_avg[which(to_plot_avg$class_id %in% c("DNA","LINE","LTR","RC","RNA","Satellite","SINE")),]
to_plot_avg$age <- factor(to_plot_avg$age,levels=c("young","old"))
ggplot(to_plot_avg, aes(x = class_id, y = rowmeans,fill=age)) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal()+  
  scale_fill_brewer(palette = "Pastel1") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12,face = "bold", color = "black"),  
    axis.text.y = element_text(size = 12,face = "bold", color = "black"),  
    axis.title.x = element_text(size = 14,face = "bold", color = "black"), 
    axis.title.y = element_text(size = 14,face = "bold", color = "black"), 
    legend.text = element_text(size = 12)
  )+ylim(0,1)

to_plot_avg <- to_plot %>%
  group_by(class_id,family_id,tissue, age) %>%
  summarise(rowmeans = mean(rowmeans, na.rm = TRUE))
to_plot_avg <- to_plot_avg[which(to_plot_avg$class_id %in% c("DNA","LINE","LTR","RC","RNA","Satellite","SINE")),]
to_plot_avg$age <- factor(to_plot_avg$age,levels=c("young","old"))
ggplot(to_plot_avg, aes(x = family_id, y = rowmeans,fill=age)) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal()+  
  scale_fill_brewer(palette = "Pastel1") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12,face = "bold", color = "black"),  
    axis.text.y = element_text(size = 12,face = "bold", color = "black"),  
    axis.title.x = element_text(size = 14,face = "bold", color = "black"), 
    axis.title.y = element_text(size = 14,face = "bold", color = "black"), 
    legend.text = element_text(size = 12)
  )+ylim(0,1)

#### heatmap
gtf <- import("~/ref_data/TE_reference/mm10_rmsk_TE.gtf", format = "gtf")
family_data <- as.data.frame(gtf[,c("gene_id","transcript_id","family_id","class_id")])
family_data <- as.data.table(family_data)
setDT(family_data)
setkey(family_data,seqnames,start,end)
kmean <- "kmeans4"
regions <- read.table(paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/bed/",kmean,"_uinon_recursion_peaks.bed"))
regions <- as.data.table(regions)
setDT(regions)
setkey(regions,V1,V2,V3)
overlaps <- foverlaps(family_data, regions, type = "any", nomatch = 0L)  
overlaps$label <- paste(overlaps$transcript_id,overlaps$gene_id,overlaps$family_id,overlaps$class_id,sep=":")
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
TE_info <- as.data.frame(gtf[,c("gene_id","transcript_id","family_id","class_id")])
TE_info$label <- paste(TE_info$transcript_id,TE_info$gene_id,TE_info$family_id,TE_info$class_id,sep=":")
TE_info <- TE_info[,c("label","width","gene_id")]
tissue_summary <- data.frame()
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/RNA_search_table.csv")
  search_table <- search_table[which(search_table$tissue_label==tissue),]
  tab <- read.table(paste0("data/samples/RNA/",tissue,"/TElocal/combined.cntTable"),header = T,row.names = 1)  
  tab_summary <- read.table(paste0("data/samples/RNA/",tissue,"/counts/",tissue,"_H3K9me3_peaks.counts.summary"),header = T)
  
  tab <- tab[which(rownames(tab)%in% overlaps$label),]
  counts <- tab
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SRR[0-9]+|HM[0-9]+).*"
  colnames(counts) <- gsub(pattern, "\\1", colnames(counts))  
  colnames(tab_summary) <- gsub(pattern,"\\1",colnames(tab_summary))
  counts <- counts[,search_table$sample_name]
  tab_summary <- tab_summary[-2,search_table$sample_name]
  
  counts <- merge(counts,TE_info,by.x="row.names",by.y="label")
  rownames(counts) <- counts$Row.names
  counts <- counts[,-1]
  
  length <- counts %>%
    group_by(gene_id) %>%
    summarise(length = sum(width))
  length <- as.data.frame(length)
  
  counts_summary <- counts %>%
    group_by(gene_id) %>%
    summarise(
      across(.cols = 1:ncol(tab), sum, .names = "{.col}")
    )
  
  counts_summary <- as.data.frame(counts_summary)
  rownames(counts_summary) <- counts_summary$gene_id
  counts_summary <- counts_summary[,-1]
  keep = which(rowSums(cpm(counts_summary)>0)>=2)
  counts_summary = counts_summary[keep,]
  
  rownames(length) <- length$gene_id
  length <- length[rownames(counts_summary),]
  length_kb <- length$length / 1000  
  total_reads <- colSums(tab_summary)
  total_reads_million <- total_reads / 1e6  
  for (i in c(1:ncol(counts_summary))) {  
    counts_summary[[i]] <- (counts_summary[[i]] / (length_kb * total_reads_million[i]))  
  }  
  RPKM <- counts_summary
  
  young_cols <- RPKM[, search_table$sample_name[which(search_table$age=="3m")]]
  young_cols$young <- rowMeans(young_cols)
  
  old_cols <- RPKM[, search_table$sample_name[which(search_table$age=="24m")]]
  old_cols$old <- rowMeans(old_cols)
  
  young_cols$label <- rownames(young_cols)
  young_cols <- young_cols[,c("label","young")]
  
  old_cols$label <- rownames(old_cols)
  old_cols <- old_cols[,c("label","old")]
  
  
  t_tissue_summary <- merge(young_cols,old_cols,by="label")
  t_tissue_summary$logFC <- log2(t_tissue_summary$old / t_tissue_summary$young)
  finite_values <- t_tissue_summary$logFC[is.finite(t_tissue_summary$logFC)]
  min_value <- min(finite_values)
  max_value <- max(finite_values) 
  t_tissue_summary$logFC[t_tissue_summary$logFC == -Inf] <- min_value
  t_tissue_summary$logFC[t_tissue_summary$logFC == Inf] <- max_value
  t_tissue_summary$tissue <- tissue_label_change(tissue)
  tissue_summary <- rbind(tissue_summary,t_tissue_summary)
}

to_plot <- tissue_summary
TE_info <- as.data.frame(family_data[,c("gene_id","family_id","class_id")])
TE_info <-  unique(TE_info)
to_plot <- merge(to_plot,TE_info[,c("gene_id","family_id","class_id")],by.x="label",by.y="gene_id")
to_plot <- to_plot[which(to_plot$class_id %in% c("DNA","LINE","LTR","Satellite","SINE")),]

## class_id
tissue_order <- c("Lung","Cerebellum","BAT","Muscle","Heart","Aorta","Skin","Kidney","Hippocampus","Cortex","Liver","Tongue","Uterus","Testis","Bladder","Ovary",
                  "Colon","Stomach","Thymus","Cecum","Jejunum","Pancreas","Bone Marrow","Ileum","Spleen","iWAT","Mammary Gland")

p_value_summary <- data.frame(matrix(ncol = length(unique(to_plot$class_id)), nrow = 27))
rownames(p_value_summary) <- tissue_order
colnames(p_value_summary) <- sort(unique(to_plot$class_id))
for(tissue in tissue_order){
  for(class_id in c("DNA","LINE","LTR","Satellite","SINE")){
    t_df <- to_plot[which(to_plot$tissue==tissue & to_plot$class_id==class_id),]
    test <- wilcox.test(t_df$old, t_df$young,paired = T)
    p_value_summary[tissue,class_id] <- test$p.value
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
colnames(p_value_long)[c(2,3)] <- c("class_id","Label")
to_plot_avg <- to_plot %>%
  group_by(class_id,tissue) %>%
  summarise(logFC = median(logFC, na.rm = TRUE))
to_plot_avg$tissue <- factor(to_plot_avg$tissue,levels=tissue_order)

merged_data <- merge(to_plot_avg, p_value_long, by = c("tissue", "class_id"), all.x = TRUE)
merged_data$logFC[which(merged_data$logFC > 1)] <- 1
merged_data$logFC[which(merged_data$logFC < -1)] <- -1
p <- ggplot(merged_data, aes(x = class_id, y = tissue, fill = logFC)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",limits = c(-1, 1), midpoint = 0) +
  theme_minimal() +
  ggtitle(kmean)+
  geom_text(aes(label = Label), color = "black", size = 4, na.rm = TRUE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
ggsave(paste0("result/Sup_figures/TE_",kmean,"_RNA_change.pdf"),p,width = 6,height = 8)
## family_id
tissue_order <- c("Lung","Cerebellum","BAT","Muscle","Heart","Aorta","Skin","Kidney","Hippocampus","Cortex","Liver","Tongue","Uterus","Testis","Bladder","Ovary",
                  "Colon","Stomach","Thymus","Cecum","Jejunum","Pancreas","Bone Marrow","Ileum","Spleen","iWAT","Mammary Gland")

p_value_summary <- data.frame(matrix(ncol = length(unique(to_plot$family_id)), nrow = 27))
rownames(p_value_summary) <- tissue_order
colnames(p_value_summary) <- sort(unique(to_plot$family_id))
for(tissue in tissue_order){
  for(family_id in sort(unique(to_plot$family_id))){
    t_df <- to_plot[which(to_plot$tissue==tissue & to_plot$family_id==family_id),]
    if(nrow(t_df)>0){  
      test <- wilcox.test(t_df$old, t_df$young,paired = T)
      p_value_summary[tissue,family_id] <- test$p.value
      }
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
colnames(p_value_long)[c(2,3)] <- c("family_id","Label")
to_plot_avg <- to_plot %>%
  group_by(family_id,tissue) %>%
  summarise(logFC = median(logFC, na.rm = TRUE))
to_plot_avg$tissue <- factor(to_plot_avg$tissue,levels=tissue_order)

merged_data <- merge(to_plot_avg, p_value_long, by = c("tissue", "family_id"), all.x = TRUE)
merged_data$logFC[which(merged_data$logFC > 1)] <- 1
merged_data$logFC[which(merged_data$logFC < -1)] <- -1
p <- ggplot(merged_data, aes(x = family_id, y = tissue, fill = logFC)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",limits = c(-1, 1), midpoint = 0) +
  theme_minimal() +
  ggtitle(kmean)+
  geom_text(aes(label = Label), color = "black", size = 4, na.rm = TRUE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
