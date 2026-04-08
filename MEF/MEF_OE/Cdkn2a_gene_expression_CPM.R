rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
search_table <- read.csv("data/samples/all/RNA_search_table.csv")

tab <- read.table("data/samples/RNA/MEF_OE_RNA/all_merge.counts",header = T)
new_names <- sub("^.*\\.bam\\.([A-Za-z]+\\d+\\.\\d+)(?:[_.].*)?$", "\\1", colnames(tab)[7:ncol(tab)])
new_names <- gsub("\\.", "-", new_names)
colnames(tab)[7:ncol(tab)] <- new_names
rownames(tab) <- tab$Geneid
counts <- tab[,c(7:ncol(tab))]
cpm <- edgeR::cpm(counts)
df1 <- as.data.frame(cpm["Cdkn2a",,drop=F])
df1 <- df1[,c("HM038-1","HM038-2","HM038-3","HM038-4","HM038-5","HM038-6","HM038-7","HM038-8","HM040-7","HM040-8")]

to_plot <- reshape2::melt(df1)
to_plot <- merge(to_plot,search_table[,c("sample_name","tissue")],by.x="variable",by.y="sample_name")
df_sum <- to_plot %>%
  group_by(tissue) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    se   = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
    .groups = "drop"
  )
to_plot$variable <- factor(to_plot$variable,levels=to_plot$variable)
p <- ggplot(df_sum, aes(x = mean, y = tissue)) +
  geom_col(width = 0.7, fill = "steelblue") +
  geom_errorbar(aes(xmin = mean - se, xmax = mean + se),
                width = 0.2, linewidth = 0.6) +
  theme_bw() +
  labs(x = NULL, y = NULL)+xlim(0,100)
p
ggsave("result/figures/MEF_OE_Cdkn2a_RNA.pdf",p,height = 6,width = 8)

df1 <- as.data.frame(cpm["Cdkn2a",,drop=F])
df1 <- df1[,c("HM040-1","HM040-2","HM040-3","HM040-4","HM040-5","HM040-6")]

to_plot <- reshape2::melt(df1)
to_plot <- merge(to_plot,search_table[,c("sample_name","tissue")],by.x="variable",by.y="sample_name")
df_sum <- to_plot %>%
  group_by(tissue) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    se   = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
    .groups = "drop"
  )
to_plot$variable <- factor(to_plot$variable,levels=to_plot$variable)
p <- ggplot(df_sum, aes(x = tissue, y = mean)) +
  geom_col(width = 0.7, fill = "steelblue") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.2, linewidth = 0.6) +
  theme_bw() +
  labs(x = "tissue", y = "value")+ylim(0,100)
p
ggsave("result/figures/MEF_OE_EZH2_Cdkn2a_RNA.pdf",p,height = 8,width = 6)
