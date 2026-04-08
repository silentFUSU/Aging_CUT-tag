rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(RColorBrewer)  
library(stringr)
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
    }else if(tissue_label == "Mammarygland"){
      tissue_label <- "Mammary Gland"
    }else if(tissue_label == "Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
} 
resolution=10000
resolution_k <- paste0(as.character(resolution/1000),"k")
tissues <- c("brain","CB","kidney", "liver", "lung", "bonemarrow", "colon", "heart", "Hip", "mammarygland", "stomach", "thymus")
short_summary <- data.frame()
long_summary <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/HiC/",tissue,"/distance_contact/all_samples_distance_contact_log2_",resolution_k,".csv"))
  short <- df[which(df$dist_log2>=13 & df$dist_log2 <=20),]
  long <- df[which(df$dist_log2>=21),]
  short <- short %>%
    group_by(sample) %>%
    summarise(total_prob = sum(prob))
  long <- long %>%
    group_by(sample) %>%
    summarise(total_prob = sum(prob))
  short <- as.data.frame(short)
  long <- as.data.frame(long)
  short_summary <- rbind(short_summary,short)
  long_summary <- rbind(long_summary,long)
}

quality_control <- read.csv("data/samples/all/HiC_Quality_control.csv")
quality_control$VP_cis.VP_unique <- as.numeric(gsub("%", "", quality_control$VP_cis.VP_unique))
quality_control$VP_cis.VP_unique <- as.numeric(quality_control$VP_cis.VP_unique)
quality_control <- quality_control %>%
  mutate(TissueName = sapply(TissueName, tissue_label_change))

short_summary <- merge(short_summary,quality_control,by.x="sample",by.y="SampleID")
long_summary <- merge(long_summary,quality_control,by.x="sample",by.y="SampleID")
color <- read.table("data/samples/20_distinct_color.txt")
color <- setNames(color$V1,sort(unique(quality_control$TissueName)))
ggplot(short_summary, aes(x =VP_cis.VP_unique , y =  total_prob , color = TissueName,shape=Age)) +
  geom_point(size = 3) +  
  scale_color_manual(values=color)+
  labs(
    title = "Short distance",
    x = "VP_cis/VP_unique",
    y = "Percentage of short-range interactions"
  ) + 
  theme_minimal()
cortest <- cor.test(short_summary$total_prob,short_summary$VP_cis.VP_unique)

ggplot(long_summary, aes(x =VP_cis.VP_unique , y =  total_prob , color = TissueName,shape=Age)) +
  geom_point(size = 3) +  
  scale_color_manual(values=color)+
  labs(
    title = "Long distance",
    x = "VP_cis/VP_unique",
    y = "Percentage of long-range interactions"
  ) + 
  theme_minimal()



