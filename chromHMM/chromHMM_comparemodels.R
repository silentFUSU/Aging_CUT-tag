rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
df <- read.delim("result/all/ChromHMM/all_tissues/comparedir/compare_to_25_state_models.txt")
rownames(df) <- paste0("State",df$State)
df <- df[,-1]
state_mean_value <- colMeans(df)
state_mean_value <- data.frame(model = paste0("model_",c(2:25)),state_mean_value=state_mean_value, state=2:25)
state_mean_value$model <- factor(state_mean_value$model,levels = paste0("model_",c(2:25)))

ggplot(state_mean_value, aes(x = model, y = state_mean_value)) +  
  geom_point() + 
  labs(x = NULL, y = NULL)  +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("ChromHMM state similarity")
