rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(edgeR)
tissue <- "CB"
antibody <- "H3K4me1"
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,axis_titles = "collect")
}

unusually_high_value_check <- function(tissue, antibody){
  if(antibody %in% c("H3K27ac","H3K4me3","H3K4me1","ATAC")){
    df <- read.table(paste0("result/all/QC/FRiP/union/",tissue,"_",antibody,"_macs_young_old_narrowpeak_rm_blacklist.counts"),header = T)
  }else{
    df <- read.table(paste0("result/all/QC/FRiP/union/",tissue,"_",antibody,"_young_old_merge-W1000-G3000-E100_rm_blacklist.counts"),header = T)
  }
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(df)[7:ncol(df)] <- gsub(pattern, "\\1",  colnames(df)[7:ncol(df)])
  rownames(df) <- df$Geneid
  counts <- df[,7:ncol(df)]
  y <- DGEList(counts = df[,7:ncol(df)], genes = data.frame(Length=df$Length))
  y <- calcNormFactors(y)  
  rpkm <- as.data.frame(rpkm(y,gene.length = df$Length))
  percent <- as.data.frame(apply(counts,2,function(x) x/sum(x)))
  samples <- colnames(counts)
  rpkm$label <- rownames(rpkm)
  percent$label <- rownames(percent)
  p_list <- list()
  for(i in c(1:length(samples))){
    to_plot <- merge(rpkm[,c(samples[i],"label")],percent[,c(samples[i],"label")],by="label")
    colnames(to_plot)[2:3] <- c("rpkm","percent")
    # par(cex.lab = 1.5, cex.axis = 1.2)
    # if(antibody %in% c("H3K27ac","H3K4me3","H3K4me1","ATAC")){
    #   smoothScatter(to_plot[,2] ~ to_plot[,3]*100,xlab = "percent",ylab = "rpkm", main = samples[i],pch = 10, col = "red")
    # }else{
    #   smoothScatter(to_plot[,2] ~ to_plot[,3]*100,bandwidth = 0.01,xlab = "percent",ylab = "rpkm", main = samples[i],pch = 10, col = "red")
    # }
    # 

    # write.csv(to_plot$label[which(to_plot$rpkm >50)],"data/raw_data/tmp.csv",row.names = F)
    to_plot_sort <- to_plot %>%  
      arrange(desc(rpkm)) %>%             # 按value列进行排序  
      mutate(cumulative_percent = cumsum(percent))
    to_plot_sort$rank <- c(0:(nrow(to_plot_sort)-1))
    p_list[[i]] <- ggplot(to_plot_sort, aes(x = rank, y = cumulative_percent)) +
      geom_point() +
      labs(title = samples[i],
           x = "rpkm Rank",
           y = "cumulative_percent") +
      geom_vline(xintercept = 10, linetype = "dashed", color = "red") +
      theme_bw()
    
  }
  
  combined_plot <- plot_a_list(p_list,2,2)
}
