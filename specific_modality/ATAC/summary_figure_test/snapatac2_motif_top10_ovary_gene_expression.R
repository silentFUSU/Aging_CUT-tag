rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(bitmapType="cairo")  
library(grid)
library(ggnewscale)
options(scipen=0)
extract_before_bracket <- function(s) {  
  parts <- strsplit(s, "\\(")[[1]]  
  return(parts[1])  
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
    }else if(tissue_label=="Mammarygland"){
      tissue_label <- "Mammary Gland"
    }else if(tissue_label=="Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
}
top <- 10
conditions <- c("up","down")
tissues <-  sort(c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
                   "thymus","skin","bladder","bonemarrow","Hip","heart",
                   "muscle","jejunum","uterus","ovary","liver","tongue",
                   "cecum","colon","testis","stomach","pancreas","iWAT","ileum"))
data_path <- "data/samples/ATAC/all/ATAC/snapatac2_macs/strict_stable_peaks_summits_spm3/"
p_list <-list()
tissue_p_value <- data.frame()
for(tissue in tissues){
  tissue_summary <- list()
  max_value <- 0
  for(condition in conditions){
    if(file.exists(paste0(data_path,condition,"/enrichment_results_",tissue,"_summits_spm3.bed.csv"))){
      motif <- read.csv(paste0(data_path,condition,"/enrichment_results_",tissue,"_summits_spm3.bed.csv"))
      motif <- motif[which(motif$log2.fold.change. > 0 & motif$adjusted.p.value < 0.05),]
      
      if(nrow(motif) > 0){
        min_non_zero <- min(motif$adjusted.p.value[motif$adjusted.p.value != 0])
        max_larger_zero_without_Inf <- max(motif$log2.fold.change.[is.finite(motif$log2.fold.change.)])
        motif$id <- ifelse(
          grepl("\\(.*\\)", motif$id),
          sub(".*\\((.*?)\\).*", "\\1", motif$id), 
          sub("^M\\d+_2\\.00\\s*", "", motif$id) 
        )
        
        motif <- motif[order(motif$adjusted.p.value, -motif$log2.fold.change.),]
        # motif <- motif[c(1:min(nrow(motif),top)),]
        motif$adjusted.p.value[motif$adjusted.p.value == 0] <- min_non_zero
        motif$log2.fold.change.[is.infinite(motif$log2.fold.change.)] <- max_larger_zero_without_Inf
        motif$log10fdr <- -log10(motif$adjusted.p.value) 
        if(max(motif$log10fdr) > max_value){
          max_value <- max(motif$log10fdr)
        }
        if(condition == "down"){
          motif$log2.fold.change. <- -motif$log2.fold.change.
        }
      }
      tissue_summary[[condition]] <- motif
    }
  }
  rank <- merge(tissue_summary[["up"]][,c("id","log2.fold.change.","log10fdr")],tissue_summary[["down"]][,c("id","log2.fold.change.","log10fdr")],by="id",all=T)
  rank[is.na(rank)] <- 0 
  rank$diff <- abs(abs(rank$log2.fold.change..x)-abs(rank$log2.fold.change..y))
  rank <- rank[order(rank$diff,decreasing = T),]
  opening <- rank[which(abs(rank$log2.fold.change..x) > abs(rank$log2.fold.change..y)),]
  opening <- opening[order(-opening$log10fdr.x,-opening$diff),]
  
  closing <- rank[which(abs(rank$log2.fold.change..x) < abs(rank$log2.fold.change..y)),]
  closing <- closing[order(-closing$log10fdr.y,-closing$diff),]
  
  to_plot <- rbind(head(opening,10),head(closing,10))
  to_plot <- rank[grep("Hoxd", rank$id),]
  to_plot$log2FC <- ifelse(abs(to_plot$log2.fold.change..x) > abs(to_plot$log2.fold.change..y), to_plot$log2.fold.change..x, to_plot$log2.fold.change..y)
  to_plot <- to_plot[order(to_plot$log2FC,decreasing = T),]
  to_plot$condition <- "up"
  to_plot$condition[which(to_plot$log2FC <0)] <- "down"
  color <- setNames(c("#FBB4AEFF","#B3CDE3FF"),c("up","down"))
  p <- ggplot(to_plot, aes(x = reorder(id, log2FC), y = log2FC,fill=condition)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = color)+
    coord_flip() +
    labs(x = "Motif", y = "Log2 Fold Change",title = "Uterus motif enrichment") +
    theme_minimal() +  
    theme(  
      axis.title.x = element_text(size = 14),  
      axis.title.y = element_text(size = 14),  
      axis.text.x = element_text(size = 14,angle = 45,hjust = 1),  
      axis.text.y = element_text(size = 14),  
      plot.title = element_text(size = 16, face = "bold")
    )
  ggsave("result/figures/Uterus_Hoxd_motif_enrichment.pdf",p,width = 5,height = 6)
}

tissue <- "ovary"
gene <- "Sox"
tab <- read.table(paste0("data/samples/RNA/",tissue,"/combined-chrM.counts"),header = T)
rownames(tab) <- tab$Geneid
counts = tab[,c(7:ncol(tab))]
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+|HM[0-9]+).*"
colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
search_table <- read.csv("data/samples/all/RNA_search_table.csv")
t_search_table <- search_table[which(search_table$tissue==tissue_label_change(tissue)),]
counts <- counts[,t_search_table$sample_name]
CPM <- as.data.frame(edgeR::cpm(counts))
CPM <- CPM[grep(gene, rownames(CPM)),]
young_cols <- CPM[,  t_search_table$sample_name[which(t_search_table$age=="3m")]]
young_cols$rowmeans <- rowMeans(young_cols)
young_cols$age <- "young"
young_cols$Geneid <- rownames(young_cols)
old_cols <- CPM[, t_search_table$sample_name[which(t_search_table$age=="24m")]]
old_cols$rowmeans <- rowMeans(old_cols)
old_cols$age <- "old"
old_cols$Geneid <- rownames(old_cols)
RNA_CPM <- rbind(young_cols[,c("Geneid","rowmeans","age")],old_cols[,c("Geneid","rowmeans","age")])

RNA <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
colnames(RNA)[1] <- "Geneid"
RNA <- RNA[grep(gene, RNA$Geneid),c("Geneid","logFC","logCPM","Significant") ]
RNA <- RNA[order(RNA$logFC),]
RNA$Geneid <- factor(RNA$Geneid,levels=rev(RNA$Geneid))
ggplot(RNA, aes(x = logFC, y = Geneid, fill = Significant)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Up" = "red", "Stable" = "gray", "Down" = "blue")) + 
  labs(title = "RNA logFC Bar Plot", x = "logFC", y = "RNA") +
  theme_bw() +
  ylab("") +
  xlab("log2(Fold change)")+
  theme(text = element_text(size = 13),
        axis.text.y = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 18, face = "bold"))+
  ggtitle(tissue_label_change(tissue))+
  scale_y_discrete(labels = function(labels) {
    labels[grepl("^Empty", labels)] <- ""
    labels
  })

RNA_CPM$Geneid <- factor(RNA_CPM$Geneid,levels=rev(RNA$Geneid))
RNA_CPM2 <- RNA_CPM[!is.na(RNA_CPM$Geneid),]
p <- ggplot(RNA_CPM2, aes(x = rowmeans, y = Geneid, fill = age)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("young" = "#FBB4AEFF", "old"="#B3CDE3FF")) + 
  ylab("") +
  xlab("CPM")+
  theme_minimal() +  
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 14,angle = 45,hjust = 1),  
    axis.text.y = element_text(size = 14),  
    plot.title = element_text(size = 16, face = "bold")
  )+
  ggtitle(paste0(tissue_label_change(tissue)," Gene expression"))+
  scale_y_discrete(labels = function(labels) {
    labels[grepl("^Empty", labels)] <- ""
    labels
  })
print(p)
ggsave("result/figures/ovary_Sox_gene_expression.pdf",p,width = 5,height = 6)
