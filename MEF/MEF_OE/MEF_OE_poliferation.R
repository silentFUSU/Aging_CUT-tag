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

df <- read.csv("data/samples/MEF_OE/MEF_OE_poliferation.csv")
df <- df[,-8]
to_plot <- df[,1:2]
to_plot$Day0 <- 0
to_plot$Day3 <- to_plot$Day0 + df$PD.D3
to_plot$Day6 <- to_plot$Day3 + df$PD.D6
to_plot$Day9 <- to_plot$Day6 + df$PD.D9
to_plot$Day12 <- to_plot$Day9 + df$PD.D12
to_plot$Day15 <- to_plot$Day12 + df$PD.D15

df_long <- to_plot %>%
  pivot_longer(cols = starts_with("Day"),
               names_to = "Day", values_to = "value") %>%
  mutate(Day = as.numeric(sub("Day", "", Day)))  # Day1 -> 1
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
  theme_bw()
ggsave("result/figures/MEF_OE_poliferation.pdf",p,width = 8,height = 4)
