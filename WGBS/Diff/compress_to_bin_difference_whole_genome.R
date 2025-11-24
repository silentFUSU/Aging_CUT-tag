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
bin_size <- "200kb"
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  t_search_table <- search_table[which(search_table$tissue==tissue),]
  bin <- read.table(paste0("~/ref_data/mm10_",bin_size,"_bins.bed"))
  if(tissue %in% c("ovary","uterus","mammarygland")){
    bin <- bin[-which(bin$V1=="chrY"),]
    bin$V1 <- factor(bin$V1,levels=paste0("chr",c(1:19,"X")))
    bin <- bin[order(bin$V1),]
    bin$V4 <- paste0("bin",1:nrow(bin))
  }else{
    bin$V1 <- factor(bin$V1,levels=paste0("chr",c(1:19,"X","Y")))
    bin <- bin[order(bin$V1),]
    bin$V4 <- paste0("bin",1:nrow(bin))
  }
  bin$V2 <- bin$V2+1
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
  young_summary <- summary[,c("label",t_search_table$sample_name[which(t_search_table$age=="3M")])]
  old_summary <- summary[,c("label",t_search_table$sample_name[which(t_search_table$age=="24M")])]
  
  young_summary$young_methylation <- rowMeans(young_summary[,-1])
  old_summary$old_methylation <- rowMeans(old_summary[,-1])
  t_tissue_summary <- merge(young_summary,old_summary,by="label")
  t_tissue_summary$delta <- (t_tissue_summary$old_methylation - t_tissue_summary$young_methylation)
  t_tissue_summary <- t_tissue_summary[,c("label","delta")]
  colnames(t_tissue_summary)[2] <- tissue_label_change(tissue)
  if(nrow(tissue_summary)==0){
    tissue_summary <- t_tissue_summary
  }else{
    tissue_summary <- merge(tissue_summary,t_tissue_summary,by="label",all=T)
  }
}
# write.csv(tissue_summary,"data/samples/WGBS/all_tissues_delta_in_200kb_bins.csv")
bin_size <- "200kb"
bin <- read.table(paste0("~/ref_data/mm10_",bin_size,"_bins.bed"))
bin$V1 <- factor(bin$V1,levels=paste0("chr",c(1:19,"X","Y")))
bin <- bin[order(bin$V1),]
bin$V4 <- paste0("bin",1:nrow(bin))

tissue_summary <- read.csv("data/samples/WGBS/all_tissues_delta_in_200kb_bins.csv")
tissue_summary <- tissue_summary[,-1]
tissue_summary$label <- factor(tissue_summary$label,levels=bin$V4)
tissue_summary <- tissue_summary[order(tissue_summary$label),]
rownames(tissue_summary) <- as.character(tissue_summary$label)
tissue_summary <- tissue_summary[,-1]
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- c(seq(-0.5, -0.11, length.out = 40), seq(-0.1, 0.1, length.out = 20), seq(0.11, 0.5, length.out = 40))
annotation <- bin[,c("V1","V4")]
rownames(annotation) <- annotation$V4
annotation <- annotation[,c("V1"),drop=F]
colnames(annotation) <- "Chromosome"
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1[1:21],paste0("chr",c(1:19,"X","Y")))
annotation_color <- list(Chromosome=color)
pheatmap::pheatmap(tissue_summary,cluster_rows = F,show_rownames = F,breaks = breaks, color = color_palette,annotation_row = annotation,annotation_colors = annotation_color,main = "Whole genome 200Kb bins CpG methylation Delta")

H3K9me3_tissue_order <- c("Lung","Cerebellum","BAT","Muscle","Heart","Aorta","Skin","Kidney","Hippocampus","Cortex",
                          "Liver","Tongue","Uterus","Testis","Bladder","Ovary","Colon","Stomach","Thymus","Cecum","Jejunum",
                          "Pancreas","Bone Marrow","Ileum","Spleen","iWAT","Mammary Gland")
to_plot <- tissue_summary
colnames(to_plot)[which(colnames(to_plot)=="IWAT")] <- "iWAT"
colnames(to_plot)[which(colnames(to_plot)=="Bone.Marrow")] <- "Bone Marrow"
colnames(to_plot)[which(colnames(to_plot)=="Mammary.Gland")] <- "Mammary Gland"
to_plot_H3K9me3_order <- to_plot[,H3K9me3_tissue_order]
pheatmap::pheatmap(to_plot_H3K9me3_order,cluster_rows = F,cluster_cols = F,annotation_row = annotation,breaks = breaks,color = color_palette,show_rownames = F,main = "CpG methylation")
