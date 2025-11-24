rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
library(ggplot2)
set.seed(1)
tissue <- "lung"
resolution <- "50000"
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
    }
  }
  return(tissue_label)
}

tissues <- c("brain","CB","kidney","liver","lung","bonemarrow","colon","heart","Hip","mammarygland","stomach","thymus","skin","muscle","cecum","ileum","pancreas","spleen")
df_pca <- data.frame()
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")  
  search_table <- search_table[which(search_table$tissue==tissue),]
  for(i in c(1:length(search_table$sample_name))){
    sample <- search_table$sample_name[i]
    df <- read.table(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_",resolution,".PC1.txt"))
    df <- df[which(df$V2 %in% c(paste0("chr",c(1:19,"X")))),]
    df <- df[,c(1,6)]
    colnames(df)[2] <- sample
    if(nrow(df_pca)==0){
      df_pca <- df
    }else{
      df_pca <- merge(df_pca,df,by="V1")
    }
  }
}
search_table <- read.csv("data/samples/all/HiC_search_table.csv")  
rownames(df_pca) <- df_pca[,1]
df_pca <- df_pca[,-1]
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

to_plot <- umap_df
tissues_label <- c("mammarygland","lung","liver","kidney","ileum","Hip","skin","bonemarrow","jejunum","colon","ovary","CB","BAT","thymus","testis","heart","stomach","muscle","bladder","aorta","tongue","spleen","pancreas","brain","cecum","uterus","iWAT")
colours <- read.table("data/samples/30_distinct_color.txt")
colours <- setNames(colours$V1,sort(sapply(tissues_label,tissue_label_change)))

p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2,color=tissue,shape=age)) + 
  scale_color_manual(values = colours) +
  geom_point(size=4) +
  theme_bw()+
  theme(text = element_text(size = 20))+
  ggtitle("HiC")
ggsave("result/Sup_figures/HiC_all_tissues_UMAP.pdf",p,width = 8,height = 6)

