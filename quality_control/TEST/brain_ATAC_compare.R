rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyr)
library(edgeR)
library(corrplot)
tissue <- "brain"
antibody <- "ATAC"
bin_size <- "1kb"
tab = read.delim(paste0("data/samples/ATAC/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins.counts"),skip=1)
counts = tab[,c(7:ncol(tab))]
rownames(counts)= tab$Geneid
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))

counts_list <- list(previous = counts[,c(1:4)], last = counts[,c(5:8)], all = counts)
conditions <- c("previous","last","all")
out_list <- list()
for(condition in conditions){
  search_table <- read.csv("data/samples/all/ATAC_search_table.csv")
  search_table <- search_table[which(search_table$sample_name %in% colnames(counts_list[[condition]])),]
  search_table$sample_name <- factor(search_table$sample_name, levels = colnames(counts_list[[condition]]))
  search_table <- search_table[order(search_table$sample_name),]
  age <- search_table$age
  mouse_ID <- search_table$mouse_ID
  age[which(age=="3m")] <- "young"
  age[which(age=="24m")] <- "old"
  colnames(counts_list[[condition]]) <- paste0(colnames(counts_list[[condition]]),"-",age,"-",mouse_ID)
  y= DGEList(counts=counts_list[[condition]],group=age)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  y$samples$year <- age
  y$samples$year <- factor(y$samples$year,c("young","old"))
  y <- calcNormFactors(y)
  design <- model.matrix(~year, y$samples)
  y<-estimateCommonDisp(y)
  y<-estimateGLMTagwiseDisp(y,design)
  fit_tag = glmFit(y,design)
  lrt = glmLRT(fit_tag, coef = 2)
  t_tab<-tab[keep,]
  out_list[[condition]] <- cbind(t_tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
                                 "PValue.old-young"=lrt$table$PValue,"FDR.old-young"= p.adjust(lrt$table$PValue,method="BH"),
                                 "LogFC.old-young"=lrt$table$logFC)
  out_list[[condition]]$Significant <- ifelse(out_list[[condition]]$`FDR.old-young` < 0.05 & abs(out_list[[condition]]$`LogFC.old-young`) >= 0, 
                            ifelse(out_list[[condition]]$`LogFC.old-young` > 0, "Up", "Down"), "Stable")
}
df <- out_list[["previous"]][-which(out_list[["previous"]]$Geneid %in% out_list[["last"]]$Geneid[which(out_list[["last"]]$`FDR.old-young`<0.05)]),]

selected_out_list <- lapply(out_list, function(df) {  
  df[, c("Geneid","LogFC.old-young","FDR.old-young"), drop = FALSE]
})  
compare_df <- Reduce(function(x, y) merge(x, y, by="Geneid"), selected_out_list) 
# compare_df <- compare_df[which(compare_df$`FDR.old-young.x` < 0.01 & compare_df$`LogFC.old-young.x` <0),]
colnames(compare_df)[c(2:7)] <- c("previous","fdr.previous","last","fdr.last","all","fdr.all")
to_plot <- compare_df
ggplot()+
  geom_point(data=to_plot, mapping=aes(previous,last),color = "grey",alpha=0.5) +  
  geom_point(data=to_plot[which(to_plot$previous>0 & to_plot$last>0),], mapping=aes(previous,last),color = "#00b8a9") +
  geom_point(data=to_plot[which(to_plot$previous<0 & to_plot$last<0),], mapping=aes(previous,last),color = "#ff9a00") +
  labs(x="previous",
       y="last") +
  theme_bw()+theme(text = element_text(size = 18))+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ggtitle(paste0("Cortex ATAC previous and last"))+
  annotate("text",label = paste0(nrow(to_plot[which(to_plot$previous>0 & to_plot$last>0),])),x = Inf, y = Inf,hjust = 1.1, vjust = 1.2,colour="#00b8a9",size=5)+
  annotate("text",label = paste0(nrow(to_plot[which(to_plot$previous<0 & to_plot$last>0),])),x = -Inf, y = Inf,hjust = -0.1, vjust = 1.2,colour="#f6416c",size=5)+
  annotate("text",label = paste0(nrow(to_plot[which(to_plot$previous<0 & to_plot$last<0),])),x = -Inf, y = -Inf,hjust = -0.1, vjust = -1.2,colour="#ff9a00",size=5)+
  annotate("text",label = paste0(nrow(to_plot[which(to_plot$previous>0 & to_plot$last<0),])),x = Inf, y = -Inf,hjust = 1.1, vjust = -1.2,colour="#48466d",size=5)

other_antibodys <- c("H3K4me1","H3K4me3","H3K27ac")
for(antibody in other_antibodys){
  df <- read.csv(paste0("data/samples/brain/",antibody,"/",antibody,"_1kb_bins_diff.csv"))
  df <- df[,c("Geneid","LogFC.old.young","FDR.old.young")]
  colnames(df)[2] <- antibody
  colnames(df)[3] <- paste0("fdr.",antibody)
  compare_df <- merge(compare_df,df,by="Geneid")
  }

to_plot <- compare_df[,c(1,2,4,6,8,10,12)]
rownames(to_plot) <- compare_df$Geneid
to_plot <- to_plot[,-1]
to_plot_cor <- cor(to_plot)
to_plot_cor[upper.tri(to_plot_cor)] <- NA  
pheatmap::pheatmap(to_plot_cor,cluster_rows = F,cluster_cols = F,breaks = seq(-1, 1, length.out = 100),display_numbers = TRUE, )

conditions <- c("previous","last")
for(condition in conditions){
  for(antibody in other_antibodys){
    to_plot <- compare_df[,c(condition,paste0("fdr.",condition),antibody,paste0("fdr.",antibody))]
    to_plot <- to_plot[which(to_plot[,2] < 0.05),]
    colnames(to_plot)[c(1,3)] <- c("ATAC","antibody")
    p <- ggplot()+
      geom_point(data=to_plot, mapping=aes(ATAC,antibody),color = "grey",alpha=0.5) +  
      geom_point(data=to_plot[which(to_plot$ATAC>0 & to_plot$antibody>0),], mapping=aes(ATAC,antibody),color = "#00b8a9") +
      geom_point(data=to_plot[which(to_plot$ATAC<0 & to_plot$antibody<0),], mapping=aes(ATAC,antibody),color = "#ff9a00") +
      labs(x=paste0(condition," ATAC"),
           y=antibody) +
      theme_bw()+theme(text = element_text(size = 18))+
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") + 
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ggtitle(paste0("Cortex ",condition," ATAC and ",antibody))+
      annotate("text",label = paste0(nrow(to_plot[which(to_plot$ATAC>0 & to_plot$antibody>0),])),x = Inf, y = Inf,hjust = 1.1, vjust = 1.2,colour="#00b8a9",size=5)+
      annotate("text",label = paste0(nrow(to_plot[which(to_plot$ATAC<0 & to_plot$antibody>0),])),x = -Inf, y = Inf,hjust = -0.1, vjust = 1.2,colour="#f6416c",size=5)+
      annotate("text",label = paste0(nrow(to_plot[which(to_plot$ATAC<0 & to_plot$antibody<0),])),x = -Inf, y = -Inf,hjust = -0.1, vjust = -1.2,colour="#ff9a00",size=5)+
      annotate("text",label = paste0(nrow(to_plot[which(to_plot$ATAC>0 & to_plot$antibody<0),])),x = Inf, y = -Inf,hjust = 1.1, vjust = -1.2,colour="#48466d",size=5)
    print(p)
    }
}
