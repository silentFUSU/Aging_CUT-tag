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

tissue <- "lung"
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
                  "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT","ileum")) 
tissue_summary <- data.frame()
state_num <- 15
for(tissue in tissues){
  file_dir <- paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/")  
  files_to_read <- list.files(path = file_dir, pattern = paste0(tissue,"_young[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
  file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  t_search_table <- search_table[which(search_table$tissue==tissue),]
  chromHMM_young <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
  chromHMM_young <- chromHMM_young[which(chromHMM_young$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  chromHMM_young$V2 <- chromHMM_young$V2+1
  bin <- chromHMM_young
  if(tissue %in% c("ovary","uterus","mammarygland")){
    bin <- bin[-which(bin$V1=="chrY"),]
    bin$V1 <- factor(bin$V1,levels=paste0("chr",c(1:19,"X")))
  }else{
    bin$V1 <- factor(bin$V1,levels=paste0("chr",c(1:19,"X","Y")))
  }
  bin <- as.data.table(bin)
  setDT(bin)
  setkey(bin,V1,V2,V3)
  summary <- data.frame()
  for(sample in t_search_table$sample_name){
    df <- fread(paste0("data/samples/WGBS/",tissue,"/bdg/",sample,"_CpG.bdg"),sep = "\t")
    setDT(df)
    setkey(df,V1,V2,V3)  
    overlaps <- foverlaps(df, bin, type = "any", nomatch = 0L)  
    
    result <- overlaps[, .(V4_sum = sum(i.V4), V5_sum = sum(V5)), by = V4]
    result <- as.data.frame(result)
    result$methylation <- result$V4_sum/result$V5_sum
    result <- result[,c("V4","methylation")]
    colnames(result) <- c("label",sample)
    if(nrow(summary)==0){
      summary <- result
    }else{
      summary <- merge(summary,result,by="label")
    }
  }
  young_samples <- t_search_table$sample_name[which(t_search_table$age=="3M")]
  old_samples <- t_search_table$sample_name[which(t_search_table$age=="24M")]
  combinations <- as.data.frame(expand.grid(young = young_samples, old = old_samples))
  compare_summary <- data.frame()
  for(i in c(1:nrow(combinations))){
    young_samples <- summary[,c("label",as.character(combinations$young[i]))]
    colnames(young_samples)[2] <- "young"
    old_samples <- summary[,c("label",as.character(combinations$old[i]))]
    colnames(old_samples)[2] <- "old"
    compare <- merge(young_samples,old_samples,by="label")
    compare$delta  <- (compare$old - compare$young)
    compare <- compare[,c("label","delta")]
    colnames(compare)[2] <- paste0(as.character(combinations$old[i]),"-",as.character(combinations$young[i]))
    if(nrow(compare_summary)==0){
      compare_summary <- compare
    }else{
      compare_summary <- merge(compare_summary,compare,by="label")
    }
  }
  t_tissue_summary <- compare_summary
  colnames(t_tissue_summary)[2:ncol(t_tissue_summary)] <- paste0(tissue_label_change(tissue),"-",colnames(t_tissue_summary)[2:ncol(t_tissue_summary)])
  if(nrow(tissue_summary)==0){
    tissue_summary <- t_tissue_summary
  }else{
    tissue_summary <- merge(tissue_summary,t_tissue_summary,by="label",all=T)
  }
}
# write.csv(tissue_summary,"data/samples/WGBS/all_tissues_delta_in_cross_comparison_chromHMM_state.csv")
tissue_summary <- read.csv("data/samples/WGBS/all_tissues_delta_in_cross_comparison_chromHMM_state.csv",row.names = 1)

H3K9me3_tissue_order_label <- c() 
annotation_col <- data.frame()
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  t_search_table <- search_table[which(search_table$tissue==tissue),]
  young_samples <- t_search_table$sample_name[which(t_search_table$age=="3M")]
  old_samples <- t_search_table$sample_name[which(t_search_table$age=="24M")]
  combinations <- as.data.frame(expand.grid(young = young_samples, old = old_samples))
  if(tissue=="bonemarrow"){
    H3K9me3_tissue_order_label <- c(H3K9me3_tissue_order_label, paste0("Bone.Marrow",".",paste0(combinations$old,".",combinations$young)))  
    t_annotation_col <- data.frame(tissue=tissue_label_change(tissue),sample=paste0("Bone.Marrow",".",paste0(combinations$old,".",combinations$young)))  
  }else if(tissue=="mammarygland"){
    H3K9me3_tissue_order_label <- c(H3K9me3_tissue_order_label, paste0("Mammary.Gland",".",paste0(combinations$old,".",combinations$young)))  
    t_annotation_col <- data.frame(tissue=tissue_label_change(tissue),sample=paste0("Mammary.Gland",".",paste0(combinations$old,".",combinations$young)))  
  }
  else{
    H3K9me3_tissue_order_label <- c(H3K9me3_tissue_order_label, paste0(tissue_label_change(tissue),".",paste0(combinations$old,".",combinations$young)))  
    t_annotation_col <- data.frame(tissue=tissue_label_change(tissue),sample=paste0(tissue_label_change(tissue),".",paste0(combinations$old,".",combinations$young)))  
  }
  annotation_col <- rbind(annotation_col,t_annotation_col)
}
rownames(annotation_col) <- annotation_col$sample
annotation_col <- annotation_col[,c("tissue"),drop = F]

annotation_row <- tissue_summary[,"label",drop=F]
colnames(annotation_row)[1] <- "state"
rownames(annotation_row) <- annotation_row$state

to_plot <- tissue_summary
rownames(to_plot) <- to_plot$label
to_plot$label <- factor(to_plot$label,levels=paste0("E",1:state_num))
to_plot <- to_plot[order(to_plot$label),]
to_plot <- to_plot[,-1]
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1[1:27],sapply(sort(tissues), tissue_label_change, USE.NAMES = FALSE))
color_row <- read.table("data/samples/20_distinct_color.txt")
color_row <- setNames(color_row$V1[1:state_num],paste0("E",1:state_num))

annotation_color <- list(state=color_row,tissue=color)

color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- c(seq(-0.5, -0.11, length.out = 40), seq(-0.1, 0.1, length.out = 20), seq(0.11, 0.5, length.out = 40))
pheatmap::pheatmap(to_plot,cluster_rows = F,cluster_cols = F,show_rownames = T,breaks = breaks, color = color_palette,annotation_row = annotation_row,annotation_col = annotation_col,annotation_colors = annotation_color,main = "Whole genome 200Kb bins CpG methylation Delta(Old - Young)")

rownames(tissue_summary) <- tissue_summary$label
tissue_mean_summary <- data.frame() 
annotation <- annotation_col
annotation$sample <- rownames(annotation)
for(tissue in tissues){
  t_annotation <- annotation[which(annotation$tissue==tissue_label_change(tissue)),]
  t_tissue_mean_summary <- tissue_summary[,t_annotation$sample]
  t_tissue_mean_summary$mean_delta <- rowMeans(t_tissue_mean_summary) 
  t_tissue_mean_summary$label <- rownames(t_tissue_mean_summary)
  t_tissue_mean_summary <- t_tissue_mean_summary[,c("label","mean_delta")]
  colnames(t_tissue_mean_summary)[2] <- tissue_label_change(tissue)
  if(nrow(tissue_mean_summary)==0){
    tissue_mean_summary <- t_tissue_mean_summary  
  }else{
    tissue_mean_summary <- merge(tissue_mean_summary,t_tissue_mean_summary,by="label")
  }
}
to_plot_mean <- tissue_mean_summary
rownames(to_plot_mean) <- to_plot_mean$label
to_plot_mean <- to_plot_mean[,-1]
breaks <- c(seq(-0.1, -0.03, length.out = 40), seq(-0.02, 0.02, length.out = 20), seq(0.03, 0.1, length.out = 40))
pheatmap::pheatmap(to_plot_mean,cluster_rows = T,cluster_cols = T,show_rownames = T,breaks = breaks, color = color_palette,annotation_row = annotation_row,annotation_colors = annotation_color,main = "Whole genome chromHMM state CpG methylation Delta(Old - Young)")








