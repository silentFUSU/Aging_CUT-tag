rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
diff1 <- read.csv("data/samples/RNA/MEF/diff_expression_gene_change_filter_bar.csv")
diff1 <- diff1[which(diff1$X == "Cdkn2a"),]
diff1 <- diff1[,1:5]
diff1 <- data.frame(sample="replicative senescence",
                    young=mean(diff1$TK2141.p2.young,diff1$TK2142.p2.young),
                    old=mean(diff1$TK2145.p10.old,diff1$TK2146.p10.old))

diff2 <- read.csv("data/samples/RNA/MEF_EZH2_inhibit/diff_expression_gene_change_filter_bar.csv")
diff2 <- diff2[which(diff2$X == "Cdkn2a"),]
diff2 <- diff2[,c(1:5)]
diff2 <- data.frame(sample="EZH2 inhibition",
                    young=mean(diff2$XM0021.DMSO.young,diff2$XM0022.DMSO.young),
                    old=mean(diff2$XM0023.GSK126.old,diff2$XM0024.GSK126.old))
to_plot <- rbind(diff1,diff2)

to_plot <- reshape2::melt(to_plot)
to_plot$variable <- as.character(to_plot$variable)
to_plot$variable <- c("P2","Normal","P10","Inhibition")
to_plot$variable <- factor(to_plot$variable,levels=c("P2","Normal","P10","Inhibition"))
p <- ggplot(to_plot, aes(x = sample, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Sample", y = "Value", fill = "Variable") +
  theme_bw() +
  ggtitle("Stacked Barplot of Value by Sample and Variable") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("result/figures/MEF_Cdkn2a_expression.pdf",p,width = 6,height = 8)
