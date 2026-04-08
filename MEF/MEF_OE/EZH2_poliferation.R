rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)

df <- read.csv("data/samples/MEF_OE/MEF_OE_poliferation.csv")
df <- df[1:3,1:4]

df2 <- read.csv("data/samples/MEF_OE/EZH2_poliferation.csv")
df2$PD.D3 <- log2(df2$PD.D3/0.2)
df2$PD.D6 <- log2(df2$PD.D6/0.3)

df <- rbind(df,df2)
df <- df[,-5]

to_plot <- df[,1:2]
to_plot$Day0 <- 0
to_plot$Day3 <- to_plot$Day0 + df$PD.D3
to_plot$Day6 <- to_plot$Day3 + df$PD.D6


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
p
ggsave("result/figures/MEF_OE_EZH2_poliferation.pdf",p,width = 6,height = 4)
