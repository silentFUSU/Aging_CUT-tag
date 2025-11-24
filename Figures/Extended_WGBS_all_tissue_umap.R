rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(Polychrome)
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
tissues <- c("mammarygland","lung","liver","kidney","ileum","Hip","skin","bonemarrow","jejunum","colon","ovary","CB","BAT","thymus","testis","heart","stomach","muscle","bladder","aorta","tongue","spleen","pancreas","brain","cecum","uterus","iWAT")
bin_size <- "10kb"

for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  t_df <- read.csv(paste0("data/samples/WGBS/",tissue,"/compress2bin/",bin_size,"_bins_all_depth.csv"))
  t_df <- reshape2::dcast(t_df,label ~ sample, value.var = "percent")
  if(i == 1){
    df <- t_df
  }else{
    df <- merge(df,t_df,by="label")
  }
}
search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
rownames(df) <- df$label
df <- df[,-1]
df <- na.omit(df) 
df_pca <- df
variances <- apply(df_pca, 1, var)
high_var_features <- names(sort(variances, decreasing = TRUE))[1:20000]
selected_df_pca <- df_pca[high_var_features, ]
pca_result <- prcomp(t(selected_df_pca), center = TRUE, scale. = TRUE)
selected_pcs <- pca_result$x[, 1:30]
umap_result <- umap::umap(selected_pcs,random_state=42)
umap_df <- data.frame(UMAP1 = umap_result$layout[, 1], UMAP2 = umap_result$layout[, 2], names = colnames(df_pca))
umap_df <- merge(umap_df,search_table[,c("sample_name","tissue","age")],by.x="names",by.y="sample_name")
umap_df$tissue <- sapply(umap_df$tissue,tissue_label_change)

umap_df$age[which(umap_df$age=="3M")] <- "young"
umap_df$age[which(umap_df$age=="24M")] <- "old"
umap_df$age <- factor(umap_df$age,levels=c("young","old"))
colours <- read.table("data/samples/30_distinct_color.txt")
colours <- setNames(colours$V1,sort(unique(umap_df$tissue)))
p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2,color=tissue,shape=age)) + 
  scale_color_manual(values = colours) +
  geom_point(size=4) +
  theme_bw()+
  theme(text = element_text(size = 20))+
  ggtitle("WGBS")

ggsave("result/Sup_figures/WGBS_all_tissues_UMAP.pdf",p,width = 10,height = 6)
