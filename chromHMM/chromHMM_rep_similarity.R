rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(corrplot)

model <- read.delim(paste0("result/all/ChromHMM/all_tissues_normal_chr/15_all_tissues/emissions_15.txt"))
dictionary <- list("1"=1, "2"=2, "3"=3,
                   "4"=4, "5"=5, "6"=7,
                   "7"=8, "8"=6, "9"=9,
                   "10"=10,"11"=11,"12"=15,
                   "13"=12,"14"=13,"15"=14)
keys <- names(dictionary)
values <- unlist(dictionary)
model$State..Emission.order. <- values[match(model$State..Emission.order., keys)]
model <- model[order(model$State..Emission.order.),]
rownames(model) <- model$State..Emission.order.
model <- model[,-1]
model <- model[,c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me1","H3K4me3")]
model <- as.data.frame(t(model))

model_rep1 <- read.delim(paste0("result/all/ChromHMM/all_tissues_normal_chr_split_rep/15_all_tissues_rep1/emissions_15.txt"))
model_rep1$State..Emission.order. <- paste0("rep1_",model_rep1$State..Emission.order.)
rownames(model_rep1) <- model_rep1$State..Emission.order.
model_rep1 <- model_rep1[,-1]
model_rep1 <- as.data.frame(t(model_rep1))
model_rep1 <- model_rep1[rownames(model),]

correlation_matrix <- cor(model_rep1, model)
selected_columns <- apply(correlation_matrix, 1, function(x) {
  max_index <- which.max(x)
  return(colnames(correlation_matrix)[max_index])
})
selected_columns <- as.data.frame(selected_columns)
selected_columns[3,1] <- 2
model_rep1 <- read.delim(paste0("result/all/ChromHMM/all_tissues_normal_chr_split_rep/15_all_tissues_rep1/emissions_15.txt"))
model_rep1$state <- as.numeric(selected_columns$selected_columns)
model_rep1<- model_rep1[order(model_rep1$state),]
rownames(model_rep1) <- model_rep1$state
model_rep1 <- model_rep1[,-c(1,8)]
model_rep1 <- model_rep1[,c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me1","H3K4me3")]
breaks <- c(seq(0, 0.8, length.out = 50))
color_palette <- colorRampPalette(c("white", "#defcf9","#4589C8FF"))(50) 
# pheatmap::pheatmap(model_rep1,cluster_cols = F,cluster_rows = F,color = color_palette,breaks = breaks,filename = "result/Sup_figures/chromHMM_15_rep1.pdf",width = 4,height = 6)

model_rep2 <- read.delim(paste0("result/all/ChromHMM/all_tissues_normal_chr_split_rep/15_all_tissues_rep2/emissions_15.txt"))
model_rep2$State..Emission.order. <- paste0("rep2_",model_rep2$State..Emission.order.)
rownames(model_rep2) <- model_rep2$State..Emission.order.
model_rep2 <- model_rep2[,-1]
model_rep2 <- as.data.frame(t(model_rep2))
model_rep2 <- model_rep2[rownames(model),]

correlation_matrix <- cor(model_rep2, model)
selected_columns <- apply(correlation_matrix, 1, function(x) {
  max_index <- which.max(x)
  return(colnames(correlation_matrix)[max_index])
})
selected_columns <- as.data.frame(selected_columns)
selected_columns[1,1] <- 2
model_rep2 <- read.delim(paste0("result/all/ChromHMM/all_tissues_normal_chr_split_rep/15_all_tissues_rep2/emissions_15.txt"))
model_rep2$state <- as.numeric(selected_columns$selected_columns)
model_rep2<- model_rep2[order(model_rep2$state),]
rownames(model_rep2) <- model_rep2$state
model_rep2 <- model_rep2[,-c(1,8)]
model_rep2 <- model_rep2[,c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me1","H3K4me3")]
breaks <- c(seq(0, 0.8, length.out = 50))
color_palette <- colorRampPalette(c("white", "#defcf9","#4589C8FF"))(50) 
# pheatmap::pheatmap(model_rep2,cluster_cols = F,cluster_rows = F,color = color_palette,breaks = breaks,filename = "result/Sup_figures/chromHMM_15_rep2.pdf",width = 4,height = 6)

model_rep1$label <- paste0("state",rownames(model_rep1))
model_rep2$label <- paste0("state",rownames(model_rep2))
model_rep1 <- reshape2::melt(model_rep1)
model_rep2 <- reshape2::melt(model_rep2)
model_rep1$label <- paste0(model_rep1$label,"-",model_rep1$variable)
model_rep2$label <- paste0(model_rep2$label,"-",model_rep2$variable)

to_plot <- merge(model_rep1[,c("label","value")],model_rep2[,c("label","value")],by="label")
to_plot <- to_plot %>%
  separate(label, into = c("state", "histone"), sep = "-")
to_plot$label <- paste0(to_plot$state,"-",to_plot$histone)
color <- read.table("data/samples/20_distinct_color.txt")
color <- setNames(color$V1,paste0("state",c(1:15)))
to_plot$state <- factor(to_plot$state,levels = paste0("state",c(1:15)))
p <- ggplot(data = to_plot, aes(x = value.x, y = value.y,color=state)) +
  scale_color_manual(values = color)+
  geom_point(size=2) +
  scale_x_continuous(
    trans = 'log10',
    breaks = c(0.0001,  0.001, 0.01, 0.1, 1),
    labels = c("0.0001",  "0.001", "0.01", "0.1", "1"),
    limits = c(0.0001, 1)
  ) +
  scale_y_continuous(
    trans = 'log10',
    breaks = c(0.0001,  0.001, 0.01, 0.1, 1),
    labels = c("0.0001",  "0.001", "0.01", "0.1", "1"),
    limits = c(0.0001, 1)
  ) +
  labs(x = "Replicate (1) [Emission Probability]", y = "Replicate (2) [Emission Probability]") +
  theme_bw()
ggsave("result/Sup_figures/chromHMM_15_rep_compare_state.pdf",p,width = 5,height = 4)
p <- ggplot(data = to_plot, aes(x = value.x, y = value.y,shape=histone)) +
  geom_point(size=2) +
  scale_x_continuous(
    trans = 'log10',
    breaks = c(0.0001,  0.001, 0.01, 0.1, 1),
    labels = c("0.0001",  "0.001", "0.01", "0.1", "1"),
    limits = c(0.0001, 1)
  ) +
  scale_y_continuous(
    trans = 'log10',
    breaks = c(0.0001,  0.001, 0.01, 0.1, 1),
    labels = c("0.0001",  "0.001", "0.01", "0.1", "1"),
    limits = c(0.0001, 1)
  ) +
  labs(x = "Replicate (1) [Emission Probability]", y = "Replicate (2) [Emission Probability]") +
  theme_bw()
ggsave("result/Sup_figures/chromHMM_15_rep_compare_histone.pdf",p,width = 5,height = 4)

cor_test <- cor.test(to_plot$value.x,to_plot$value.y,method = "spearman")

