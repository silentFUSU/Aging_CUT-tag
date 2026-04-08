rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(stringr)
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
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
      tissue_label <- "Mammary gland"
    }else if(tissue_label=="Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
}
gene <- "RLTR6B_Mm"


gene_summary_plot <- function(gene){
  summary <- data.frame()
  for(tissue in tissues){
    df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_TE_change_filter_bar.csv"))
    matching_rows <- grep(paste0(gene,":"), df$X)
    df <- df[matching_rows,]
    if(nrow(df)>0){
      df <- df[,c("logFC","fdr","Significant")]
      df$tissue <- tissue_label_change(tissue)
      
      if(nrow(summary)==0){
        summary <- df
      }else{
        summary <- rbind(summary,df)
      }
    }
  }
  color <- setNames(c("red","grey","blue"),c("TE-Up","Stable","TE-Down"))
  p <- ggplot(summary,mapping = aes(x=logFC,y=tissue,fill = Significant))+
    geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+ylab("")+
    scale_fill_manual(values=color)+
    theme(text = element_text(size = 13))+ 
    ggtitle(paste0(gene," RNA")) +
    geom_text(data = summary[which(summary$logFC<0),],aes(label = round(logFC,3)), position = position_dodge2(width = 0.9), hjust = 0.5, size = 5) +
    geom_text(data = summary[which(summary$logFC>0),],aes(label = round(logFC,3)), position = position_dodge2(width = 0.9), hjust = 0.5, size = 5)
  print(p)
}

# gene_summary_plot <- function(gene){
#   summary <- data.frame()
#   for(tissue in tissues){
#     df <- read.delim(paste0("data/samples/RNA/",tissue,"/TEcount/combined.cntTable"),row.names = 1)
#     pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SRR[0-9]+|HM[0-9]+).*"
#     colnames(df) <- gsub(pattern, "\\1", colnames(df))
#     CPM <- as.data.frame(cpm(df))
#     matching_rows <- grep(paste0(gene,":"), rownames(CPM))
#     CPM <- CPM[matching_rows, ]
#     search_table <- read.csv("data/samples/all/RNA_search_table.csv")
#     search_table <- search_table[which(search_table$sample_name %in% colnames(df)),]
#     search_table$sample_name <- factor(search_table$sample_name, levels = colnames(df))
#     search_table <- search_table[order(search_table$sample_name),]
#     young_samples <- CPM[,search_table$sample_name[which(search_table$age=="3m")]]
#     young_samples <- colSums(young_samples)
#     young_samples <- mean(young_samples)
#     old_samples <- CPM[,search_table$sample_name[which(search_table$age=="24m")]]
#     old_samples <- colSums(old_samples)
#     old_samples <- mean(old_samples)
#     t_summary <- data.frame(tissue=tissue_label_change(tissue),logFC=log2(old_samples/young_samples))
#     summary <- rbind(summary,t_summary)
#   }
#   summary <- summary[order(summary$logFC),]
#   summary$tissue <- factor(summary$tissue,levels = summary$tissue)
#   p <- ggplot(summary,mapping = aes(x=logFC,y=tissue))+
#     geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+ylab("")+
#     theme(text = element_text(size = 13))+ 
#     ggtitle(paste0(gene," RNA")) +
#     geom_text(data = summary[which(summary$logFC<0),],aes(label = round(logFC,3)), position = position_dodge2(width = 0.9), hjust = 0.5, size = 5) +
#     geom_text(data = summary[which(summary$logFC>0),],aes(label = round(logFC,3)), position = position_dodge2(width = 0.9), hjust = 0.5, size = 5)
#   print(p)
# }