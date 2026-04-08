rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
library(ggrepel)
library(tidyr)
library(dplyr)
library(tidyverse)
library(plotly)
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
load("data/samples/GRN/grn_union_tissue.rdata")
load("data/samples/GRN/grn_union_skin.rdata")
grn_tissue[["skin"]] <- grn_union
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")

### raw value calculate
summary <- data.frame()
for(tissue in tissues){
  
  search_table <- read.csv("data/samples/all/RNA_search_table.csv")
  search_table <- search_table[which(search_table$tissue_label==tissue),]
  
  pagerank <- read.table("data/samples/GRN/TF_pagerank_sample_new.txt")
  pagerank <- pagerank[,search_table$sample_name]
  young_pagerank <- pagerank[,search_table$sample_name[which(search_table$age=="3m")]]
  young_pagerank$young <- rowMeans(young_pagerank)
  
  old_pagerank <- pagerank[,search_table$sample_name[which(search_table$age=="24m")]]
  old_pagerank$old <- rowMeans(old_pagerank)
  pagerank <- merge(young_pagerank[,"young",drop=F],old_pagerank[,"old",drop=F],by="row.names")
  rownames(pagerank) <- pagerank$Row.names
  pagerank <- pagerank[,-1]
  pagerank <- pagerank[rowSums(pagerank) != 0, ]
  pagerank$logFC <- log2(pagerank$old/pagerank$young)
  logfc_values <- pagerank$logFC
  logfc_filtered <- logfc_values[is.finite(logfc_values)]
  max_logfc <- max(logfc_filtered)
  min_logfc <- min(logfc_filtered)
  pagerank$logFC[which(pagerank$logFC== Inf)] <- max_logfc
  pagerank$logFC[which(pagerank$logFC== -Inf)] <- min_logfc
  
  
  df <- grn_tissue[[tissue]]
  df <- merge(df,pagerank[,c("logFC"),drop=F],by.x="TF",by.y="row.names")
  
  # get TF pagerank pvalue and p.adj
  pagerank_limma <- read.table(paste0("data/samples/GRN/TF_pagerank_limma_remove_zero_row_log/",tissue,"_TF_diff.txt"))
  df <- merge(df,pagerank_limma[,c("P.Value","adj.P.Val"),drop=F],by.x="TF",by.y="row.names")
  
  #keep TF with expression change significantly
  gene <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
  gene <- gene[which(gene$Significant !="Stable"),]
  df <- df[which(df$TF %in% gene$X),]
  
  colnames(df)[which(colnames(df)%in% c("logFC","P.Value","adj.P.Val"))] <- c("pagerank_logFC","P.Value","adj.P.Val")
  
  if(nrow(df) > 0){
    df$gene_logFC <- log2(df$gene_old/df$gene_young)
    logfc_values <- df$gene_logFC
    logfc_filtered <- logfc_values[is.finite(logfc_values)]
    max_logfc <- max(logfc_filtered)
    min_logfc <- min(logfc_filtered)
    df$gene_logFC[which(df$gene_logFC== Inf)] <- max_logfc
    df$gene_logFC[which(df$gene_logFC== -Inf)] <- min_logfc
    
    # get number of TF-gene pairs  
    TF_counts <- as.data.frame(table(df$TF))
    df <- merge(df,TF_counts,by.x="TF",by.y="Var1")
    df$TF <- paste0(tissue_label_change(tissue),"-",df$TF)
    
    
    t_summary <- df[,c("TF","gene_logFC","Freq","pagerank_logFC","P.Value","adj.P.Val")]
    summary <- rbind(summary,t_summary)
  }
}

mean_logFC_summary <- summary %>%
  filter(!is.infinite(gene_logFC)) %>%  # 过滤掉无穷大值（如 Inf 和 -Inf）
  group_by(TF,pagerank_logFC,Freq,P.Value,adj.P.Val) %>%
  filter(n() >= 50) %>%  # 仅保留至少有50行数据的TF
  summarise(mean_logFC = mean(gene_logFC, na.rm = TRUE))  # 计算平均数

mean_logFC_summary <- mean_logFC_summary[order(mean_logFC_summary$mean_logFC,decreasing = T),]

to_plot <- mean_logFC_summary

# keep padj < 0.05
to_plot_sig <- to_plot[which(to_plot$adj.P.Val < 0.05),]
to_plot_sig <- to_plot_sig[order(to_plot_sig$pagerank_logFC,decreasing = T),]

# highest_points <- 
# lowest_points <- to_plot_sig[(nrow(to_plot_sig)-9):nrow(to_plot_sig),]
annotate_points <- to_plot_sig[which(to_plot_sig$TF %in% c("Ovary-Esr2","Ovary-Nr5a2","Aorta-Dbp","BAT-Nr5a2","Colon-Klf8","Mammary Gland-Lef1","Ileum-Hoxb9","Uterus-Irf4","Uterus-Spi1","Ovary-Irf4","Pancreas-Egr2")),]

p <- ggplot(to_plot[which(to_plot$adj.P.Val <0.05),], aes(x = pagerank_logFC, y = mean_logFC)) +
  geom_point(aes(size = Freq),color = "black",alpha=0.3) +
  geom_point(data=annotate_points,aes(x = pagerank_logFC, y = mean_logFC,size = Freq),color = "red",alpha=0.3) +
  labs(x = "raw value logFC", y = "TF related gene logFC") +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank()  # 去掉次要网格线
  )+
  ylim(-0.5,0.5)+
  xlim(-4,4)+
  geom_text_repel(data = annotate_points, aes(label = TF), color = "red", size = 4,max.overlaps = 50)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "blue")
p
ggsave("result/figures/correlation_pagerank_with_gene_padj.pdf",p,width = 8,height = 7)
