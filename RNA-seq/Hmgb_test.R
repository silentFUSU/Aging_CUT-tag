rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
tissues <- c("skin","CB","spleen","heart","bladder","tongue","uterus","aorta","thymus","stomach","Hip","FC","BAT","iWAT","muscle","bonemarrow","lung","kidney","liver","testis","colon","cecum","ileum","jejunum")
df <- data.frame(Gene = as.character(),
                 logFC = as.numeric(),
                 Significant = as.character(),
                 tissue = as.character())
for(tissue in tissues){
  t_df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_nodup.csv"))
  t_df <- t_df[which(t_df$X %in% c("Hmgb1","Hmgb2")),c("X","logFC","Significant")]
  colnames(t_df)[1] <- "Gene"
  t_df$tissue <- tissue
  df <- rbind(df,t_df)
}
