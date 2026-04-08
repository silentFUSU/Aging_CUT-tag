rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(edgeR)
library(ggplot2)
library(clusterProfiler)
library(stringr)
diff1 <- read.csv("data/samples/RNA/MEF_mid_age/diff_expression_gene_change_filter_bar.csv")
# diff1 <- diff1[order(diff1$fdr),]
# diff1 <- diff1[1:1000,]
# 
# down_sig <- subset(diff1, Significant == "Down")
# up_sig <- subset(diff1, Significant == "Up")
# down_sig_sorted <- down_sig[order(down_sig$fdr),]
# up_sig_sorted <- up_sig[order(up_sig$fdr),]
# down_top_2000 <- head(down_sig_sorted, 1000)
# up_top_2000 <- head(up_sig_sorted, 1000)
# diff1$Significant <- "Stable"
# diff1$Significant[which(diff1$X %in% down_top_2000$X)] <- "Down"
# diff1$Significant[which(diff1$X %in% up_top_2000$X)] <- "Up"

diff2 <- read.csv("data/samples/RNA/MEF_EZH2_inhibit/diff_expression_gene_change_filter_bar.csv")
# down_sig <- subset(diff2, Significant == "Down")
# up_sig <- subset(diff2, Significant == "Up")
# down_sig_sorted <- down_sig[order(down_sig$fdr),]
# up_sig_sorted <- up_sig[order(up_sig$fdr),]
# down_top_2000 <- head(down_sig_sorted, 1000)
# up_top_2000 <- head(up_sig_sorted, 1000)
# diff2$Significant <- "Stable"
# diff2$Significant[which(diff2$X %in% down_top_2000$X)] <- "Down"
# diff2$Significant[which(diff2$X %in% up_top_2000$X)] <- "Up"

colnames(diff1)[which(colnames(diff1)=="logFC")] <- "MEF"
colnames(diff2)[which(colnames(diff2)=="logFC")] <- "MEF_EZH2_inhibit"

to_plot <- merge(diff1[,c("X","MEF","Significant")],diff2[,c("X","MEF_EZH2_inhibit","Significant")],by="X",all=T)
to_plot$Significant.x[is.na(to_plot$Significant.x)] <- "Stable"
to_plot$Significant.y[is.na(to_plot$Significant.y)] <- "Stable"
to_plot$condition <- "Stable"
to_plot$condition[which(to_plot$Significant.x=="Up" & to_plot$Significant.y=="Up")] <- "Up"
to_plot$condition[which(to_plot$Significant.x=="Down" & to_plot$Significant.y=="Down")] <- "Down"
to_plot$condition[which(to_plot$Significant.x=="Stable" & to_plot$Significant.y=="Stable")] <- "Stable"
to_plot$condition[which(to_plot$Significant.x=="Up" & to_plot$Significant.y=="Down")] <- "Inconsistent"
to_plot$condition[which(to_plot$Significant.x=="Down" & to_plot$Significant.y=="Up")] <- "Inconsistent"
color <- setNames(c("#e64b35","#3c5488","gray","#00a087"),c("Up","Down","Stable","Inconsistent"))

x_range <- range(to_plot$MEF, na.rm = TRUE)  
y_range <- range(to_plot$MEF_EZH2_inhibit, na.rm = TRUE)  
x_pos_right <- x_range[2] * 0.9    
x_pos_left <- x_range[1] * 0.9   
y_pos_top <- y_range[2] * 0.9    
y_pos_bottom <- y_range[1] * 0.9 
p <- ggplot(to_plot[which(to_plot$condition=="Stable"),], aes(x = MEF, y = MEF_EZH2_inhibit,color=condition)) +
  geom_point(alpha=0.1) +
  geom_point(data = to_plot[which(to_plot$condition!="Stable"),], aes(x = MEF, y = MEF_EZH2_inhibit,color=condition)) +
  scale_color_manual(values=color)+
  labs(title = "Scatter Plot of MEF vs MEF_EZH2_inhibit",
       x = "MEF",
       y = "MEF_EZH2_inhibit") +
  theme_bw()+
  annotate("text", label = paste0(nrow(to_plot[which(to_plot$Significant.x=="Up" & to_plot$Significant.y=="Up"),])),  
           x = x_pos_right, y = y_pos_top, colour = "#00b8a9", size = 5) +  
  annotate("text", label = paste0(nrow(to_plot[which(to_plot$Significant.x=="Down" & to_plot$Significant.y=="Down"),])),  
           x = x_pos_left, y = y_pos_bottom, colour = "#ff9a00", size = 5) +  
  annotate("text", label = paste0(nrow(to_plot[which(to_plot$Significant.x=="Down" & to_plot$Significant.y=="Up"),])),  
           x = x_pos_left, y = y_pos_top, colour = "#f6416c", size = 5) +  
  annotate("text", label = paste0(nrow(to_plot[which(to_plot$Significant.x=="Up" & to_plot$Significant.y=="Down"),])),  
           x = x_pos_right, y = y_pos_bottom, colour = "#48466d", size = 5) 
ggsave("result/figures/MEF_p10_p6_vs_MEF_EZH2_inhibit_RNA.pdf",p,width = 7,height = 6)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
common_increased <- to_plot$X[which(to_plot$Significant.x=="Up" & to_plot$Significant.y=="Up")]
genelist_up <- bitr(common_increased,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_up_GO <- enrichGO( genelist_up$ENTREZID,#GO富集分析
                            OrgDb = GO_database,
                            keyType = "ENTREZID",#设定读取的gene ID类型
                            ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                            pvalueCutoff = 0.05,#设定p值阈值
                            qvalueCutoff = 0.05,#设定q值阈值
                            readable = T)
barplot(genelist_up_GO,label_format = 50,showCategory = 20)

common_decreased <- to_plot$X[which(to_plot$Significant.x=="Down" & to_plot$Significant.y=="Down")]
genelist_down <- bitr(common_decreased,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_down_GO <- enrichGO( genelist_down$ENTREZID,#GO富集分析
                            OrgDb = GO_database,
                            keyType = "ENTREZID",#设定读取的gene ID类型
                            ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                            pvalueCutoff = 0.05,#设定p值阈值
                            qvalueCutoff = 0.05,#设定q值阈值
                            readable = T)
barplot(genelist_down_GO,label_format = 50,showCategory = 20)


second <- to_plot$X[which(to_plot$Significant.x=="Down" & to_plot$Significant.y=="Up")]
genelist <- bitr(second,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_GO <- enrichGO( genelist$ENTREZID,#GO富集分析
                            OrgDb = GO_database,
                            keyType = "ENTREZID",#设定读取的gene ID类型
                            ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                            pvalueCutoff = 0.05,#设定p值阈值
                            qvalueCutoff = 0.05,#设定q值阈值
                            readable = T)
barplot(genelist_GO,label_format = 50,showCategory = 20)



fourth <- to_plot$X[which(to_plot$Significant.x=="Up" & to_plot$Significant.y=="Down")]
genelist <- bitr(fourth,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_GO <- enrichGO( genelist$ENTREZID,#GO富集分析
                         OrgDb = GO_database,
                         keyType = "ENTREZID",#设定读取的gene ID类型
                         ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                         pvalueCutoff = 0.05,#设定p值阈值
                         qvalueCutoff = 0.05,#设定q值阈值
                         readable = T)
barplot(genelist_GO,label_format = 50,showCategory = 20)





