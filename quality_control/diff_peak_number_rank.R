rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(dplyr)
diff_peak_number <- read.csv("data/samples/all/diff_bins_number.csv")
diff_peak_number$tissue[which(diff_peak_number$tissue=="brain")] <- "Cortex"
diff_peak_number$tissue[which(diff_peak_number$tissue=="Hip")] <- "Hippocampus"
diff_peak_number$tissue[which(diff_peak_number$tissue=="CB")] <- "Cerebellum"
diff_peak_number <- diff_peak_number %>% mutate(tissue = str_to_title(tissue))  
diff_peak_number$tissue[which(diff_peak_number$tissue=="Bonemarrow")] <- "Bone Marrow"

antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K4me3","H3K4me1","H3K27ac")
increase <- data.frame(Var1 = as.character(),
                       Rank = as.numeric(),
                       tissue = as.character(),
                       antibodys = as.character())

for(i in c(1:length(antibodys))){
  antibody <- antibodys[i]
  df <- diff_peak_number[which(diff_peak_number$antibody==antibody & diff_peak_number$Var1=="Up"),]
  max <- max(df$Freq)
  df$Rank <- df$Freq/max
  increase <- rbind(increase,df[,c("Var1","Rank","tissue","antibody")])
  }

increase_rank <- increase %>%
  group_by(tissue) %>%
  summarize(mean_rank = mean(Rank))

decrease <- data.frame(Var1 = as.character(),
                       Rank = as.numeric(),
                       tissue = as.character(),
                       antibodys = as.character())
for(i in c(1:length(antibodys))){
  antibody <- antibodys[i]
  df <- diff_peak_number[which(diff_peak_number$antibody==antibody & diff_peak_number$Var1=="Down"),]
  max <- max(df$Freq)
  df$Rank <- df$Freq/max
  decrease <- rbind(decrease,df[,c("Var1","Rank","tissue","antibody")])
}
decrease_rank <- decrease %>%
  group_by(tissue) %>%
  summarize(mean_rank = mean(Rank))

palette <- read.table("data/samples/30_distinct_color.txt")
palette <- palette$V1
tissues <-unique(diff_peak_number$tissue)
palette <- setNames(palette,tissues)

increase_rank <- arrange(increase_rank,mean_rank)
increase_rank$tissue <- factor(increase_rank$tissue,levels = increase_rank$tissue)
increase_p <- ggplot(increase_rank,mapping = aes(x=mean_rank,y=tissue,fill = tissue))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+ylab("")+xlab("Mean Rank")+ggtitle(paste0("Increase"))+
  theme(text = element_text(size = 18))+ scale_fill_manual(values = palette)+ xlim(0,1) +guides(fill= guide_legend(title = ""))

decrease_rank <- arrange(decrease_rank,mean_rank)
decrease_rank$tissue <- factor(decrease_rank$tissue,levels = decrease_rank$tissue)
decrease_p <- ggplot(decrease_rank,mapping = aes(x=mean_rank,y=tissue,fill = tissue))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+ylab("")+xlab("Mean Rank")+ggtitle(paste0("Decrease"))+
  theme(text = element_text(size = 18))+ scale_fill_manual(values = palette)+ xlim(0,1) +guides(fill= guide_legend(title = ""))

increase_p + decrease_p

combined_change <- rbind(increase, decrease)
combined_change <- combined_change  %>%
  group_by(tissue) %>%
  summarize(mean_rank = mean(Rank))

combined_change <- arrange(combined_change,mean_rank)
combined_change$tissue <- factor(combined_change$tissue,levels = combined_change$tissue)
ggplot(combined_change,mapping = aes(x=mean_rank,y=tissue,fill = tissue))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+ylab("")+xlab("Mean Rank")+
  theme(text = element_text(size = 18))+ scale_fill_manual(values = palette)+ xlim(0,1) +guides(fill= guide_legend(title = ""))
