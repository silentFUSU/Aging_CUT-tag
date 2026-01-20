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
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","iWAT","jejunum","kidney","liver",
                  "lung","mammarygland","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus"))


PMD_HMD_region <- read.table("data/public_data/PMD_coordinates_mm10.bed")
PMD_HMD_region$V2 <- PMD_HMD_region$V2 + 1
PMD_HMD_region$label <- paste0("bin",c(1:nrow(PMD_HMD_region)))
PMD_HMD_region <- PMD_HMD_region[,c("V1","V2","V3","V5","label")]
PMD_HMD_region$V5[is.na(PMD_HMD_region$V5)] <- "other"
PMD_HMD_region <- as.data.table(PMD_HMD_region)
setDT(PMD_HMD_region)
setkey(PMD_HMD_region,V1,V2,V3)
tissue_regions_summary <- list()
tissue_overall_summary <- data.frame()
for(tissue in tissues){
  print(tissue)
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  t_tissue_regions_summary <- data.frame()
  t_tissue_overall_summary <- data.frame()
  for(sample in search_table$sample_name){
    df <- fread(paste0("data/samples/WGBS/",tissue,"/bdg/",sample,"_CpG.bdg"),sep = "\t")
    setDT(df)
    setkey(df,V1,V2,V3)  
    overlaps <- foverlaps(df, PMD_HMD_region, type = "any", nomatch = 0L)  
    
    sample_regions_summary <- overlaps[, .(V4_sum = sum(V4), V5_sum = sum(i.V5)),by=label]
    sample_regions_summary <- as.data.frame(sample_regions_summary)
    sample_regions_summary$methylation <- sample_regions_summary$V4_sum/sample_regions_summary$V5_sum*100
    sample_regions_summary <- sample_regions_summary[,c("label","methylation")]
    colnames(sample_regions_summary)[2] <- sample 
    if(nrow(t_tissue_regions_summary) ==0){
      t_tissue_regions_summary <- sample_regions_summary
    }else{
      t_tissue_regions_summary <- merge(t_tissue_regions_summary,sample_regions_summary,by="label")
    }
    
    sample_overall_summary <- overlaps[, .(V4_sum = sum(V4), V5_sum = sum(i.V5)),by=V5]
    sample_overall_summary <- as.data.frame(sample_overall_summary)
    sample_overall_summary$methylation <- sample_overall_summary$V4_sum/sample_overall_summary$V5_sum*100
    sample_overall_summary$sample <- sample
    sample_overall_summary$tissue <- tissue_label_change(tissue)
    t_tissue_overall_summary <- rbind(t_tissue_overall_summary,sample_overall_summary)
  }
  tissue_regions_summary[[tissue_label_change(tissue)]] <- t_tissue_regions_summary
  tissue_overall_summary <- rbind(tissue_overall_summary,t_tissue_overall_summary)
}

# saveRDS(tissue_regions_summary,"data/samples/WGBS/all/DNA_methylation_change_in_PMD_HMD_regions_summary.rds")
# write.csv(tissue_overall_summary,"data/samples/WGBS/all/DNA_methylation_change_in_PMD_HMD_overall_summary.csv")

tissue_regions_summary <- readRDS("data/samples/WGBS/all/DNA_methylation_change_in_PMD_HMD_regions_summary.rds")
to_plot <- data.frame()
for(tissue in tissues){
  df <- tissue_regions_summary[[tissue_label_change(tissue)]]
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  young_df <- df[,c("label",search_table$sample_name[which(search_table$age=="3M")])]  
  young_df$young <- rowMeans(young_df[,-1])
  old_df <- df[,c("label",search_table$sample_name[which(search_table$age=="24M")])] 
  old_df$old <- rowMeans(old_df[,-1])
  
  df <- merge(young_df[,c("label","young")],old_df[,c("label","old")],by="label")
  df$delta <- df$old - df$young
  df <- df[,c("label","delta")]
  colnames(df) <- c("label",tissue_label_change(tissue))
  
  if(nrow(to_plot)==0){
    to_plot <- df
  }else{
    to_plot <- merge(to_plot,df,by="label")
  }
}
rownames(to_plot) <- to_plot$label
to_plot <- to_plot[,-1]

annotation <- read.table("data/public_data/PMD_coordinates_mm10.bed")
annotation$V2 <- annotation$V2 + 1
annotation$label <- paste0("bin",c(1:nrow(annotation)))
annotation <- annotation[,c("label","V5")]
colnames(annotation)[2] <- "condition"
annotation <- annotation[which(annotation$label %in% rownames(to_plot)),]
bin_order <- annotation[order(annotation$condition),]
to_plot <- to_plot[bin_order$label,]
rownames(annotation) <- annotation$label
annotation <- annotation[,-1,drop=F]
annotation_color <- list(condition=setNames(c("#e64b35","#3c5488"),c("PMD","HMD")))

breaks <- c(seq(-5, -2.1, length.out = 40), seq(-2, 2, length.out = 20), seq(2.1, 5, length.out = 40))
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
tissue_order <- c("Mammary Gland","Cecum","Thymus","Uterus","iWAT","Stomach","Skin","Spleen",
                  "Muscle","Bone Marrow","Liver","Ileum","Testis","Cortex","Jejunum","Tongue",
                  "Hippocampus","Colon","Bladder","Aorta","Cerebellum","Lung","Heart","Kidney","BAT","Ovary","Pancreas")
to_plot <- to_plot[,tissue_order]
pheatmap::pheatmap(to_plot,breaks = breaks,color = color_palette,annotation_row = annotation,
                   annotation_colors = annotation_color,show_rownames = F,cluster_rows = F,cluster_cols = F,
                   filename = "result/figures/DNAm_change_in_PMD_HMD.pdf",height = 8,width = 6)



