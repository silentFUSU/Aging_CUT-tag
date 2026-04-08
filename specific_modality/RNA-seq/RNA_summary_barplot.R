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

gene <- "Htra1"
tissue <- "ovary"
search_table <- read.csv("data/samples/all/RNA_search_table.csv")
search_table$age[which(search_table$age=="3m")] <- "young"
search_table$age[which(search_table$age=="24m")] <- "old"

search_table <- search_table[which(search_table$tissue == tissue_label_change(tissue)),]
df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
df <- df[which(df$X==gene),]
if(nrow(df)>0){
  if(tissue == "skin"){
    df <- df[,which(colnames(df) %in% paste0(search_table$sample_name,".",search_table$mouse_ID,".",search_table$age))]
  }else{
    df <- df[,which(colnames(df) %in% paste0(search_table$sample_name,".",search_table$mouse_ID,".",search_table$age))]
  }
  df <- as.data.frame(t(df))
  colnames(df)[1] <- "CPM"
  df$tissue <- tissue_label_change(tissue)
  rownames(df) <- sapply(strsplit(rownames(df), "\\."), `[`, 1)
  search_table <- search_table[which(search_table$sample_name%in% rownames(df)),]
  df$sample_name <- search_table$sample_name
  df <- merge(df,search_table[,c("sample_name","age")])
}
to_plot<- df %>%
  group_by(age) %>%
  summarise(mean_CPM = mean(CPM, na.rm = TRUE))
to_plot$age <- factor(to_plot$age,levels=c("young","old"))

p <- ggplot(to_plot, aes(x = age, y = mean_CPM,fill=age)) +
  geom_bar(stat = "identity",color="black") +
  theme_bw() +
  labs(title = paste0(tissue,"-",gene),
       x = "Age",
       y = "Mean CPM")
ggsave(paste0("result/Sup_figures/ovary_",gene,"_gene_expression.pdf"),p,width = 4,height = 5)
