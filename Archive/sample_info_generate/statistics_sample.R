setwd("~/projects/Aging_CUT_Tag/")
library(dplyr)
library(ggplot2)
summary <- data.frame()
for(i in 1:5){
  df <- readxl::read_excel("data/samples/all/mouse_information.xlsx",sheet = i)
  df$sex <- "male"  
  df$sex[which(df$organ %in% c("Mammary Gland","Ovary","Uterus"))] <- "female"
  df <- df[,c("organ","experiment","sex")]
  summary <- rbind(summary,df)
  }

sex <- summary %>%
  count(sex, name = "n") %>%
  arrange(sex) %>%
  mutate(p = n / sum(n),
         ymax = cumsum(p),
         ymin = lag(ymax, default = 0))
exp <- summary %>%
  count(experiment, name = "n") %>%
  arrange(experiment) %>%
  mutate(p = n / sum(n),
         ymax = cumsum(p),
         ymin = lag(ymax, default = 0))

gap <- 0.3   # 两圈间隔
w_in  <- 1.2 # 内环厚度（大一点）
w_out <- 0.6 # 外环厚度（细一点）

r_in  <- c(0, w_in)
r_out <- c(w_in + gap, w_in + gap + w_out)

p <- ggplot() +
  geom_rect(data = exp,
            aes(xmin = r_in[1], xmax = r_in[2], ymin = ymin, ymax = ymax, fill = experiment),
            color = "white", linewidth = 0.3) +
  geom_rect(data = sex,
            aes(xmin = r_out[1], xmax = r_out[2], ymin = ymin, ymax = ymax, fill = sex),
            color = "white", linewidth = 0.3) +
  coord_polar(theta = "y") +
  xlim(0, r_out[2]) +
  theme_void() +
  theme(legend.position = "right")
ggsave("~/projects/Aging_CUT_Tag/result/figures/sample_statistics.pdf",p,width = 8,height =6)
sex <- sex[,c(1:3)]
sex$sex[which(sex$sex=="male")] <- "Male"
sex$sex[which(sex$sex=="female")] <- "Female"
sex[,3] <- sex[,3] *100
sex <- sex[c(2,1),]
colnames(sex) <- c("condition","# of samples","% of samples")
saveRDS(sex,"~/projects/database_web/data/sex_summary.rds")
exp <- exp[,c(1:3)]
exp[,3] <- exp[,3] *100
colnames(exp) <- c("condition","# of samples","% of samples")
saveRDS(exp,"~/projects/database_web/data/exp_summary.rds")
organ <- summary %>%
  count(organ, name = "n") %>%
  arrange(organ) %>%
  mutate(p = n / sum(n),
         ymax = cumsum(p),
         ymin = lag(ymax, default = 0))
