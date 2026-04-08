rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
library(ggrepel)
library(grid)
library(ggsignif)
library(data.table)
library(dplyr)
library(GenomeInfoDb)
library("GenomicRanges")
library(genomation)
tissue_label_change <- function(tissue){
  if(tissue=="brain"){
    tissue_label <- "Cortex"
  }else if(tissue == "Hip"){
    tissue_label <- "Hippocampus"
  }else if(tissue == "CB"){
    tissue_label <- "Cerebellum"
  }else{
    tissue_label <- str_to_title(tissue)
    if(tissue_label == "Bonemarrow"){
      tissue_label <- "Bone Marrow"
    }else if(tissue_label == "Bat"){
      tissue_label <- "BAT"
    }else if(tissue_label=="Mammarygland"){
      tissue_label <- "Mammary Gland"
    }else if(tissue_label=="Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
}
tissue <- "CB"
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
                  "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT","ileum")) 
tissue_summary <- data.frame()
p_list <- list()
for(tissue in sort(tissues)){
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  tab <- read.table(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_young_old_merge-W5000-G10000-E100.counts"),header = T)
  peak_regions <- as.data.table(tab[,c(1:4)])
  setDT(peak_regions)
  setkey(peak_regions,Chr,Start,End)
  blacklist <- read.table("~/ref_data/mm10-blacklist.v2.bed",sep = "\t")
  blacklist <- as.data.table(blacklist)
  setDT(blacklist)
  setkey(blacklist,V1,V2,V3)
  overlaps <- foverlaps(peak_regions, blacklist, type = "any", nomatch = 0L)
  tab <- tab[which(!tab$Geneid %in% overlaps$Geneid),]
  
  tab_summary <- read.table(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_young_old_merge-W5000-G10000-E100.counts.summary"),header = T)
  
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+|DYQ[0-9]+|HM[0-9]+).*"
  colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
  colnames(tab_summary) <- gsub(pattern,"\\1",colnames(tab_summary))
  
  search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
  counts <- counts[,search_table$sample_name]
  tab_summary <- tab_summary[-2,search_table$sample_name]
  total_reads <- colSums(tab_summary)
  length <- as.numeric(tab$Length)
  rpkm <- sweep(counts,2,total_reads,"/")
  rpkm <- sweep(rpkm,1,length,"/") * 1000000000
  rpkm$mean_all <- rowMeans(rpkm)
  rpkm_young <- rpkm[,search_table$sample_name[which(search_table$age=="3m")]]
  rpkm_old <- rpkm[,search_table$sample_name[which(search_table$age=="24m")]]
  rpkm_young$mean_young <- rowMeans(rpkm_young)
  rpkm_old$mean_old <- rowMeans(rpkm_old)
  rpkm_mean_summary <- merge(rpkm_young[,"mean_young",drop=F],rpkm_old[,"mean_old",drop=F],by="row.names")
  rpkm_mean_summary <- merge(rpkm_mean_summary,rpkm[,"mean_all",drop=F],by.x="Row.names",by.y="row.names")
  
  rpkm_mean_summary$log2FC <- log2(rpkm_mean_summary$mean_old/rpkm_mean_summary$mean_young)
  colnames(rpkm_mean_summary)[1] <- "Geneid"
  rpkm_mean_summary$tissue <- tissue_label_change(tissue)
  
  to_plot <- rpkm_mean_summary
  to_plot$condition <- NA
  to_plot$condition[which(to_plot$log2FC < 0)] <- "Down"
  to_plot$condition[which(to_plot$log2FC > 0)] <- "Up"
  to_plot$condition <- factor(to_plot$condition,levels=c("Up","Down"))
  to_plot$quadrant <- "first"
  to_plot$quadrant[which(log2(to_plot$mean_young)<0 & to_plot$log2FC>0)] <- "second"
  to_plot$quadrant[which(log2(to_plot$mean_young)<0 & to_plot$log2FC<0)] <- "third"
  to_plot$quadrant[which(log2(to_plot$mean_young)>0 & to_plot$log2FC<0)] <- "fourth"
  to_plot$quadrant <- factor(to_plot$quadrant,levels=c("first","second","third","fourth"))
  
  x_range <- range(log2(to_plot$mean_young), na.rm = TRUE)  
  y_range <- range(to_plot$log2FC, na.rm = TRUE)  
  x_pos_right <- x_range[2] * 0.9    
  x_pos_left <- x_range[1] * 0.9   
  y_pos_top <- y_range[2] * 0.9    
  y_pos_bottom <- y_range[1] * 0.9 
  ### log2 RPKM log2(O/Y)
  ggplot(to_plot,aes(x=log2(mean_young),y=log2FC,color=condition))+    
    geom_jitter(size = 3, alpha = 0.7)+
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # 添加水平线
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # 添加竖直线
    geom_hline(yintercept = 1, color = "red") +  # 添加水平线
    geom_hline(yintercept = -1, color = "red") +
    theme_bw()+theme(text = element_text(size = 18))+
    xlab("log2(RPKM)")+
    ylab("log2(O/Y)")+
    ggtitle(tissue_label_change(tissue))+
    labs(fill = "", color = "")+
    annotate("text", label = paste0(nrow(to_plot[which(log2(to_plot$mean_young)>0 & to_plot$log2FC>0),])),  
             x = x_pos_right, y = y_pos_top, colour = "#00b8a9", size = 5) +  
    annotate("text", label = paste0(nrow(to_plot[which(log2(to_plot$mean_young)<0 & to_plot$log2FC<0),])),  
             x = x_pos_left, y = y_pos_bottom, colour = "#ff9a00", size = 5) +  
    annotate("text", label = paste0(nrow(to_plot[which(log2(to_plot$mean_young)<0 & to_plot$log2FC>0),])),  
             x = x_pos_left, y = y_pos_top, colour = "#f6416c", size = 5) +  
    annotate("text", label = paste0(nrow(to_plot[which(log2(to_plot$mean_young)>0 & to_plot$log2FC<0),])),  
             x = x_pos_right, y = y_pos_bottom, colour = "#48466d", size = 5) 
  
  ## EZH2 SUZ12
  EZH2 <- fread("data/public_data/EZH2_SUZ12_E14/peaks/macs_narrowpeak_input/EZH2_comp.bedGraph")
  colnames(EZH2)[4] <- "EZH2"
  
  SUZ12 <- fread("data/public_data/EZH2_SUZ12_E14/peaks/macs_narrowpeak_input/SUZ12_comp.bedGraph")
  colnames(SUZ12)[4] <- "SUZ12"
  
  regions <- as.data.table(tab[,c(1:4)])
  setDT(regions)
  setkey(regions,Chr,Start,End)
  
  setDT(EZH2)
  setkey(EZH2,V1,V2,V3)
  
  
  setDT(SUZ12)
  setkey(SUZ12,V1,V2,V3)
  
  overlaps <- foverlaps(EZH2,regions, type = "any", nomatch = 0L) 
  overlaps<- overlaps[, .(EZH2 = mean(EZH2)), by = .(V1, Start, End, Geneid)]
  overlaps <- as.data.table(overlaps[,c("V1","Start","End","Geneid","EZH2")])
  
  setDT(overlaps)
  setkey(overlaps,V1,Start,End)
  
  overlaps_SUZ12 <- foverlaps(SUZ12,overlaps, type = "any", nomatch = 0L)  
  overlaps_SUZ12<- overlaps_SUZ12[, .(SUZ12 = mean(SUZ12)), by = .(V1, Start, End, Geneid, EZH2)]
  
  overlaps_SUZ12$score <- (overlaps_SUZ12$SUZ12 + overlaps_SUZ12$EZH2)/2
  
  overlaps_SUZ12 <- as.data.frame(overlaps_SUZ12)
  overlaps_SUZ12 <- overlaps_SUZ12[,c("V1","Start","End","Geneid","score")]
  
  to_plot <- merge(to_plot,overlaps_SUZ12[,c("Geneid","score")],by="Geneid")
  to_plot$score_label <- to_plot$score
  to_plot$score_label[which(to_plot$score_label > 10)] <- 10
  ggplot(to_plot,aes(x=log2(mean_young),y=log2FC,color=score_label))+    
    geom_point(size = 1, alpha = 0.7)+
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # 添加水平线
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # 添加竖直线
    geom_hline(yintercept = 1, color = "red") +  # 添加水平线
    geom_hline(yintercept = -1, color = "red") +
    theme_bw()+theme(text = element_text(size = 18))+
    xlab("log2(RPKM)")+
    ylab("log2(O/Y)")+
    ggtitle(tissue_label_change(tissue))+
    labs(fill = "", color = "")+
    scale_color_gradientn(colors = c("blue", "yellow", "red"),limits = c(0, 10)) +
    annotate("text", label = paste0(nrow(to_plot[which(log2(to_plot$mean_young)>0 & to_plot$log2FC>0),])),  
             x = x_pos_right, y = y_pos_top, colour = "#00b8a9", size = 5) +  
    annotate("text", label = paste0(nrow(to_plot[which(log2(to_plot$mean_young)<0 & to_plot$log2FC<0),])),  
             x = x_pos_left, y = y_pos_bottom, colour = "#ff9a00", size = 5) +  
    annotate("text", label = paste0(nrow(to_plot[which(log2(to_plot$mean_young)<0 & to_plot$log2FC>0),])),  
             x = x_pos_left, y = y_pos_top, colour = "#f6416c", size = 5) +  
    annotate("text", label = paste0(nrow(to_plot[which(log2(to_plot$mean_young)>0 & to_plot$log2FC<0),])),  
             x = x_pos_right, y = y_pos_bottom, colour = "#48466d", size = 5) 
  
  p_list[[tissue]] <- ggplot(to_plot, aes(x = quadrant, y = score,fill=condition)) +
    geom_boxplot(outliers = F) +
    ggtitle(tissue_label_change(tissue))+
    # scale_fill_manual(values = color) +
    labs(x = NULL, y = "EZH2 SUZ12 binding levels") +
    theme_bw()+ 
    theme(
      plot.title = element_text(size = 16, face = "bold"),           # 图标题字体
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),  # x轴标签字体和旋转设置
      axis.text.y = element_text(size = 14),                        # y轴标签字体
      axis.title.x = element_text(size = 16),                       # x轴标题字体
      axis.title.y = element_text(size = 16),                       # y轴标题字体
      legend.text = element_text(size = 14),                        # 图例文本字体
      legend.title = element_text(size = 16)                        # 图例标题字体
    )
}
tissue_order <- c("BAT","aorta","iWAT","muscle","brain","Hip","bladder","thymus","tongue","stomach","pancreas",
                  "ovary","lung","testis","mammarygland","bonemarrow","uterus","heart","kidney","ileum","jejunum",
                  "skin","liver","colon","cecum","spleen","CB")
p_list_order <- p_list[tissue_order]
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
combined_plot <- plot_a_list(p_list_order, 4, 7)
ggsave("tmp.png",combined_plot,width = 35,height = 20,type="cairo")


