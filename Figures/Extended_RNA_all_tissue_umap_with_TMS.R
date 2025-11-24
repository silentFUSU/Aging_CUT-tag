rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
library(limma)
library(ggrepel)
tab = read.csv("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/GSE132040/GSE132040_190214_A00111_0269_AHH3J3DSXX_190214_A00111_0270_BHHMFWDSXX.csv")
tab <- tab[-c(54353:54357),]
tab <- as.data.frame(t(tab))
colnames(tab) <- tab[1,]
tab <- tab[-1,]
new_row_name <- sub("\\..*", "", rownames(tab)) 
tab$Sample.name <- new_row_name
table = read.delim("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/GSE132040/GSE132040_MACA_Bulk_metadata.csv",sep = ',')
tab <- merge(tab,table[,c(1,3,5,7)],by="Sample.name")
tab <- tab[!grepl("^NA", tab$source.name), ] 
rownames(tab) <- tab$Sample.name
tab <- tab[,-1]
tab <- tab[which(tab$characteristics..sex=="m"),]
tab <- tab[which(tab$characteristics..age %in% c("3","24")),]
count <- as.data.frame(t(tab[,-c(54353:54355)])) 

tab_ourdata <- read.table("data/samples/RNA/all_tissues_combined-chrM.counts",header = T)
rownames(tab_ourdata) <- tab_ourdata$Geneid
tab_ourdata <- tab_ourdata[,-1]
colnames <- colnames(tab_ourdata)[6:ncol(tab_ourdata)]
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|HM[0-9]+).*"
new_colnames <- gsub(pattern, "\\1", colnames)
sorted_index <- order(new_colnames)
order_colnames <- new_colnames[sorted_index] 
colnames(tab_ourdata)[6:ncol(tab_ourdata)] <- new_colnames
counts <- tab_ourdata[,order_colnames] 

count_merge <- merge(count, counts,by = 'row.names')
rownames(count_merge) <- count_merge$Row.names
count_merge <- count_merge[,-1]
count_merge[] <- lapply(count_merge, as.numeric) 
df_pca <- as.data.frame(edgeR::cpm(count_merge))
variances <- apply(df_pca, 1, var)
high_var_features <- names(sort(variances, decreasing = TRUE))[1:2000]
selected_df_pca <- df_pca[high_var_features, ]
pca_result <- prcomp(t(selected_df_pca), center = TRUE, scale. = TRUE)
selected_pcs <- pca_result$x[, 1:30]
umap_result <- umap::umap(selected_pcs,random_state=42)
umap_df <- data.frame(UMAP1 = umap_result$layout[, 1], UMAP2 = umap_result$layout[, 2], names = colnames(df_pca))


group_info <- tab[,c(54353:54355)]
group_info$source.name <- sub("_(\\w+)_\\d*|_(\\w+)$", "\\1",  group_info$source.name)  
group_info$source.name <- paste0("TMS-",group_info$source.name)
group <- read.csv("data/samples/all/RNA_search_table.csv",sep = ',')
group <- group[which(group$sample_name %in% colnames(counts)),]
group <- group[order(group$sample_name),]
tissue <- group$tissue
age <- group$age


new_row <- data.frame(source.name=tissue,characteristics..age=age,characteristics..sex="m")
rownames(new_row) <- group$sample_name
group_info <- rbind(group_info,new_row)

umap_df <- merge(umap_df,group_info[,c("source.name","characteristics..age")],by.x="names",by.y="row.names")
umap_df$tissue <-gsub("^TMS-", "", umap_df$source.name)
umap_df$condition <- ifelse(startsWith(umap_df$source.name, "TMS"), "TMS", "Own")

search_table <- read.csv("data/samples/all/RNA_search_table.csv")
color <- read.table("data/samples/30_distinct_color.txt")
tissues <- sort(unique(search_table$tissue))
tissues <- tissues[-which(tissues=="MEF")]
tissues <- c(tissues,sort(unique(umap_df$tissue[which(! umap_df$tissue %in% tissues)])))
color <- setNames(c(color$V1,"#FF5733","#33FF57","#5733FF","#FF33A1","#33A1FF","#A1FF33"),tissues)
condition_shape <- c("TMS" = 0, "Own" = 16)


colnames(umap_df)[5] <- "age"
umap_df$age[which(umap_df$age %in% c("3", "3m"))] <- "young"
umap_df$age[which(umap_df$age %in% c("24","24m"))] <- "old"
condition_shape <- c("TMS" = 0, "Own" = 16)
umap_df$shape_combo <- interaction(umap_df$age, umap_df$condition, sep = "_")
condition_shape <- c("young_Own" = 19, "old_Own" = 17, "young_TMS" = 21, "old_TMS" = 24)
unique_tissue_df <- umap_df[!duplicated(umap_df$source.name), ]

p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2,color=tissue,shape=shape_combo, alpha = condition)) + 
  scale_color_manual(values = color) +
  scale_shape_manual(values = condition_shape) +
  scale_alpha_manual(values = c("TMS" = 0.5, "Own" = 1)) + 
  geom_point(size=3) +
  theme_bw()+
  geom_text_repel(data = unique_tissue_df, aes(label = source.name),
                  size = 5,color="black", hjust = 0.5, vjust = -1) +
  theme(text = element_text(size = 20))+
  ggtitle("RNA with TMS data")
ggsave("result/Sup_figures/RNA_with_TMS_all_tissues_UMAP.pdf",p,width = 15,height = 12)
