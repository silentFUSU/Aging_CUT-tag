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

gene <- "Gata"
gene_summary_plot <- function(gene,tissue){
  df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
  df <- df[grep("Gata", df$X),]
  df <- df[,c("X","logFC","fdr","Significant")]
  df <- df[order(df$logFC),]
  df$X <- factor(df$X,levels=rev(df$X))
  color <- setNames(c("red","grey","blue"),c("Up","Stable","Down"))
  p <- ggplot(df,mapping = aes(x=logFC,y=X,fill = Significant))+
    geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+ylab("")+
    scale_fill_manual(values=color)+
    theme(text = element_text(size = 13))+ 
    ggtitle(paste0(tissue_label_change(tissue)," RNA")) +
    geom_text(data = df[which(df$logFC<0),],aes(label = round(logFC,3)), position = position_dodge2(width = 0.9), hjust = 0.5, size = 5) +
    geom_text(data = df[which(df$logFC>0),],aes(label = round(logFC,3)), position = position_dodge2(width = 0.9), hjust = 0.5, size = 5)
  print(p)
  genes_order <- df$X
  search_table <- read.csv("data/samples/all/RNA_search_table.csv")
  search_table$age[which(search_table$age=="3m")] <- "young"
  search_table$age[which(search_table$age=="24m")] <- "old"
  
  search_table <- search_table[which(search_table$tissue == tissue_label_change(tissue)),]
  df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
  df <- df[grep("Gata", df$X),]
  rownames(df) <- df$X
  if(nrow(df)>0){
    if(tissue == "skin"){
      df <- df[,which(colnames(df) %in% paste0(search_table$sample_name,".",search_table$mouse_ID,".",search_table$age))]
    }else{
      df <- df[,which(colnames(df) %in% paste0(search_table$sample_name,".",search_table$mouse_ID,".",search_table$age))]
    }
    df <- as.data.frame(t(df))
    df$sample <- rownames(df)
    df <- reshape2::melt(df)
    colnames(df)[3] <- "CPM"
    df$sample <- sapply(strsplit(df$sample, "\\."), `[`, 1)
    search_table <- search_table[which(search_table$sample_name %in% df$sample),]
    df <- merge(df,search_table[,c("sample_name","age")],by.x="sample",by.y="sample_name")

  }

  df$age <- factor(df$age, levels = c("young","old"))
  df$variable <- factor(df$variable,levels=rev(genes_order))
  p <- ggplot(df,aes(x=CPM,y=variable,color = age))+    
    geom_point(size = 2, alpha = 0.7)+
    geom_text(aes(label = sample), vjust = -1, size = 3) +
    scale_fill_brewer(palette="Set3")+
    ggtitle(paste0(tissue_label_change(tissue)," ",gene," CPM"))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+xlab("CPM")+labs(fill = "", color = "") +ylab(NULL)
  print(p)
  return(p)
  } 

gene_summary_plot(gene,tissues)
