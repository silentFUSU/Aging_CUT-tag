rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(edgeR)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
library(limma)
library(data.table)
df <- read.csv("data/samples/MEF_EZH2_inhibit/MEF_Inhibit_poliferation.csv")
colnames(df)[1] <- "condition"

to_plot <- df[,1,drop=F]
to_plot$day0 <- 0
to_plot$day1 <- to_plot$day0 + log2(df$day1/df$day0)
to_plot$day3 <- to_plot$day1 + log2(df$day3/df$day1)
to_plot$day5 <- to_plot$day3 + log2(df$day5/df$day3)
to_plot$day7 <- to_plot$day5 + log2(df$day7/df$day5)

df_long <- to_plot %>%
  pivot_longer(cols = starts_with("day"),
               names_to = "Day", values_to = "value") %>%
  mutate(Day = as.numeric(sub("day", "", Day))) 
sum_df <- df_long %>%
  group_by(condition, Day) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value,   na.rm = TRUE),
    n    = sum(!is.na(value)),
    se   = sd / sqrt(n),
    .groups = "drop"
  )
p <- ggplot(sum_df, aes(x = Day, y = mean, color = condition, group = condition)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.6) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.4, linewidth = 0.6) +
  scale_x_continuous(breaks = sort(unique(sum_df$Day))) +
  labs(x = "Day", y = "Value", color = "Condition") +
  theme_classic()
ggsave("result/figures/MEF_inhibit_poliferation.pdf",p,width = 8,height = 6)
