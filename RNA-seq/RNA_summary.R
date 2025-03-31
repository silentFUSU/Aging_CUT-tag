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
summary <- data.frame()
gene <- "C4b"
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene.csv"))
  df <- df[which(df$X==gene),]
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
color <- setNames(c("red","grey","blue"),c("Up","Stable","Down"))
ggplot(summary,mapping = aes(x=logFC,y=tissue,fill = Significant))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+ylab("")+
  scale_fill_manual(values=color)+
  theme(text = element_text(size = 13))+ 
  ggtitle(paste0(gene," RNA")) +
  geom_text(data = summary[which(summary$logFC<0),],aes(label = round(logFC,3)), position = position_dodge2(width = 0.9), hjust = 0.5, size = 5) +
  geom_text(data = summary[which(summary$logFC>0),],aes(label = round(logFC,3)), position = position_dodge2(width = 0.9), hjust = 0.5, size = 5)


summary <- data.frame()
gene <- "Setdb1"
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/RNA_search_table.csv")
  search_table$age[which(search_table$age=="3m")] <- "young"
  search_table$age[which(search_table$age=="24m")] <- "old"
  
  search_table <- search_table[which(search_table$tissue == tissue_label_change(tissue)),]
  df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene.csv"))
  df <- df[which(df$X==gene),]
  if(nrow(df)>0){
    df <- df[,which(colnames(df) %in% paste0(search_table$sample_name,".",search_table$mouse_ID,".",search_table$age))]
    df <- as.data.frame(t(df))
    colnames(df)[1] <- "CPM"
    df$tissue <- tissue_label_change(tissue)
    rownames(df) <- sapply(strsplit(rownames(df), "\\."), `[`, 1)
    search_table <- search_table[which(search_table$sample_name%in% rownames(df)),]
    df$sample_name <- search_table$sample_name
    df <- merge(df,search_table[,c("sample_name","age")])
    if(nrow(summary)==0){
      summary <- df
    }else{
      summary <- rbind(summary,df)
    }
  }
}
summary$age <- factor(summary$age, levels = c("young","old"))
ggplot(summary,aes(x=tissue,y=CPM,color = age))+    
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.7)+
  geom_text(aes(label = sample_name), position = position_jitter(width = 0.2), vjust = -1, size = 3) +
  scale_fill_brewer(palette="Set3")+
  ggtitle(paste0(gene," CPM"))+
  theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")+labs(fill = "", color = "") +ylab("CPM")
