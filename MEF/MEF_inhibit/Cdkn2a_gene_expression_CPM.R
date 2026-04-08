rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
df1 <- read.table("data/samples/RNA/MEF/combined-chrM.counts",header = T)
colnames <- colnames(df1)[6:length(df1)]
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|XM[0-9]+|DYQ[0-9]+|TK[0-9]+).*"
new_colnames <- gsub(pattern, "\\1", colnames)
colnames(df1)[6:length(df1)] <- new_colnames
rownames(df1) <- df1$Geneid
df1 <- df1[,c(7:ncol(df1))]
df1 <- edgeR::cpm(df1)
df1 <- as.data.frame(df1["Cdkn2a",,drop=F])

# df1 <- data.frame(sample="replicative senescence",
#                     young=mean(df1$TK2141,df1$TK2142),
#                     middle=mean(df1$TK2143,df1$TK2144),
#                     old=mean(df1$TK2145,df1$TK2146))

df2 <- read.table("data/samples/RNA/MEF_EZH2_inhibit/combined-chrM.counts",header = T)
colnames <- colnames(df2)[6:length(df2)]
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|XM[0-9]+|DYQ[0-9]+|TK[0-9]+).*"
new_colnames <- gsub(pattern, "\\1", colnames)
colnames(df2)[6:length(df2)] <- new_colnames
rownames(df2) <- df2$Geneid
df2 <- df2[,c(7:ncol(df2))]
df2 <- edgeR::cpm(df2)
df2 <- as.data.frame(df2["Cdkn2a",,drop=F])
# df2 <- data.frame(sample="EZH2 inhibition",
#                   young=mean(df2$XM0021,df2$XM0022),
#                   middle=NA,
#                   old=mean(df2$XM0024,df2$XM0023))
df <- cbind(df1,df2)
search_table <- read.csv("data/samples/all/RNA_search_table.csv")
to_plot <- reshape2::melt(df)
to_plot <- merge(to_plot,search_table[,c("sample_name","tissue")],by.x="variable",by.y="sample_name")
to_plot <- to_plot[c(1,3,5,7,9,11,13:16),]
to_plot$tissue <- c("P2","P2","P6","P6","P10","P10","DMSO","DMSO","GSK126","GSK126")
df_sum <- to_plot %>%
  group_by(tissue) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    se   = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
    .groups = "drop"
  )
df_sum$tissue <- factor(df_sum$tissue,levels=rev(c("P2","P6","P10","DMSO","GSK126")))
p <- ggplot(df_sum, aes(x = mean, y = tissue)) +
  geom_col(width = 0.7, fill = "steelblue") +
  geom_errorbar(aes(xmin = mean - se, xmax = mean + se),
                width = 0.2, linewidth = 0.6) +
  theme_bw() +
  labs(x = NULL, y = NULL)+xlim(0,300)
p

ggsave("result/figures/MEF_Cdkn2a_expression.pdf",p,width = 6,height = 8)





