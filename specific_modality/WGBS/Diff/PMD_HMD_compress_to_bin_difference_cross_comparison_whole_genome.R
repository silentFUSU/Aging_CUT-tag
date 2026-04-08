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

tissue <- "mammarygland"
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","iWAT","jejunum","kidney","liver",
                  "lung","mammarygland","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus"))
# tissues <- c("muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus")
tissue_summary <- data.frame()
PMD_HMD_region <- read.table("data/public_data/PMD_coordinates_mm10.bed")
PMD_HMD_region$V2 <- PMD_HMD_region$V2 + 1
PMD_HMD_region$label <- paste0("bin",c(1:nrow(PMD_HMD_region)))
PMD_HMD_region <- PMD_HMD_region[,c("V1","V2","V3","V5","label")]
PMD_HMD_region$V5[is.na(PMD_HMD_region$V5)] <- "other"
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  t_search_table <- search_table[which(search_table$tissue==tissue),]
  bin <- PMD_HMD_region
  if(tissue %in% c("ovary","uterus","mammarygland")){
    bin <- bin[which(bin$V1 %in% paste0("chr",c(1:19,"X"))),]
    bin$V1 <- factor(bin$V1,levels=paste0("chr",c(1:19,"X")))
    bin <- bin[order(bin$V1),]
  }else{
    bin$V1 <- factor(bin$V1,levels=paste0("chr",c(1:19,"X","Y")))
    bin <- bin[order(bin$V1),]
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
    
    result <- overlaps[, .(V4_sum = sum(V4), V5_sum = sum(i.V5)), by = label]
    result <- as.data.frame(result)
    result$methylation <- result$V4_sum/result$V5_sum
    result <- result[,c("label","methylation")]
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
# write.csv(tissue_summary,"data/samples/WGBS/all_tissues_delta_in_100kb_bins_PMD_HMD_cross_comparison.csv")
tissue_summary <- read.csv("data/samples/WGBS/all_tissues_delta_in_100kb_bins_PMD_HMD_cross_comparison.csv",row.names = 1)
annotation_col <- data.frame()
tissue_label <- c()
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  t_search_table <- search_table[which(search_table$tissue==tissue),]
  young_samples <- t_search_table$sample_name[which(t_search_table$age=="3M")]
  old_samples <- t_search_table$sample_name[which(t_search_table$age=="24M")]
  combinations <- as.data.frame(expand.grid(young = young_samples, old = old_samples))
  if(tissue=="bonemarrow"){
    tissue_label <- c(tissue_label, paste0("Bone.Marrow",".",paste0(combinations$old,".",combinations$young)))  
    t_annotation_col <- data.frame(tissue=tissue_label_change(tissue),sample=paste0("Bone.Marrow",".",paste0(combinations$old,".",combinations$young)))  
  }else if(tissue=="mammarygland"){
    tissue_label <- c(tissue_label, paste0("Mammary.Gland",".",paste0(combinations$old,".",combinations$young)))  
    t_annotation_col <- data.frame(tissue=tissue_label_change(tissue),sample=paste0("Mammary.Gland",".",paste0(combinations$old,".",combinations$young)))  
  }
  else{
    tissue_label <- c(tissue_label, paste0(tissue_label_change(tissue),".",paste0(combinations$old,".",combinations$young)))  
    t_annotation_col <- data.frame(tissue=tissue_label_change(tissue),sample=paste0(tissue_label_change(tissue),".",paste0(combinations$old,".",combinations$young)))  
  }
  annotation_col <- rbind(annotation_col,t_annotation_col)
}

rownames(tissue_summary) <- tissue_summary$label
tissue_mean_summary <- data.frame() 
annotation <- annotation_col
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
write.csv(tissue_mean_summary,"data/samples/WGBS/all_tissues_delta_in_100kb_bins_PMD_HMD_cross_comparison_means.csv")

to_plot <- tissue_mean_summary[which(tissue_mean_summary$label %in% PMD_HMD_region$label[which(PMD_HMD_region$V1=="chr2")]),]
rownames(to_plot) <- to_plot$label
to_plot$label <- factor(to_plot$label,levels=PMD_HMD_region$label)
to_plot <- to_plot[order(to_plot$label),]
to_plot <- to_plot[,-1]
annotation_row <- PMD_HMD_region[which(PMD_HMD_region$V1=="chr2"),]
annotation_row <- annotation_row[,c("label","V5")]
rownames(annotation_row) <- annotation_row$label
annotation_row <- annotation_row[,c("V5"),drop=F]
colnames(annotation_row) <-"condition"
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- c(seq(-0.2, -0.06, length.out = 40), seq(-0.05, 0.05, length.out = 20), seq(0.06, 0.2, length.out = 40))
pheatmap::pheatmap(to_plot,cluster_rows =T,cluster_cols = T,annotation_row = annotation_row,breaks = breaks,color = color_palette,show_rownames = F,main = "Chromosome 2 CpG methylation Delta(Old - Young)")

H3K9me3_tissue_order <- c("Lung","Cerebellum","BAT","Muscle","Heart","Aorta","Skin","Kidney","Hippocampus","Cortex",
                          "Liver","Tongue","Testis","Bladder","Pancreas","Cecum","Spleen","Stomach","Colon","Bone Marrow","Jejunum",
                          "iWAT","Thymus","Ileum")
to_plot <- to_plot[,H3K9me3_tissue_order]
pheatmap::pheatmap(to_plot,cluster_rows =T,cluster_cols = F,annotation_row = annotation_row,breaks = breaks,color = color_palette,show_rownames = F,main = "Chromosome 1 CpG methylation Delta(Old - Young)")

PMD_to_plot <- tissue_mean_summary[which(tissue_mean_summary$label %in% PMD_HMD_region$label[which(PMD_HMD_region$V5=="PMD")]),]
PMD_to_plot$condition <- "PMD"
HMD_to_plot <- tissue_mean_summary[which(tissue_mean_summary$label %in% PMD_HMD_region$label[which(PMD_HMD_region$V5=="HMD")]),]
HMD_to_plot$condition <- "HMD"

to_plot <- rbind(PMD_to_plot,HMD_to_plot)
to_plot <- reshape2::melt(to_plot)
to_plot$tissue_condition <- "Large changes"
to_plot$tissue_condition[which(to_plot$variable %in% c("Ileum","Bladder","Testis","Tongue","Cecum","Pancreas","Stomach","Bone Marrow","Colon","Spleen","Thymus","Jejunum","iWAT"))] <- "Small changes"

small_change_test <- wilcox.test(to_plot$value[which(to_plot$condition=="HMD" & to_plot$tissue_condition=="Small changes")],
                            to_plot$value[which(to_plot$condition=="PMD" & to_plot$tissue_condition=="Small changes")])
large_change_test <- wilcox.test(to_plot$value[which(to_plot$condition=="HMD" & to_plot$tissue_condition=="Large changes")],
                            to_plot$value[which(to_plot$condition=="PMD" & to_plot$tissue_condition=="Large changes")])



ggplot(to_plot, aes(x = tissue_condition, y = value,fill=condition)) +
  geom_boxplot(outliers =F) +
  labs(title = "Whole genome CpG methylation Delta(Old - Young)",
       x = NULL,
       y = "Delta(Old - Young)") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 14),  # 调整图例标题大小
    legend.text = element_text(size = 12),   # 调整图例项的字符大小
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12)  # 调整x轴文本角度和大小
  )+
  geom_signif(comparisons = list(c("Large changes", "Small changes")),
              textsize = 4,test = "t.test",
              map_signif_level = TRUE,
              y_position = c(0.07),
              tip_length = c(1/100))

#### per tissue PMD HMD compare
p_list <- list()
for(tissue in H3K9me3_tissue_order){
  # t_tissue_mean_summary <- tissue_mean_summary[,c("label",tissue_label_change(tissue))]  
  t_tissue_mean_summary <- tissue_mean_summary[,c("label",tissue)]
  t_tissue_mean_summary <- merge(t_tissue_mean_summary,PMD_HMD_region[,c("label","V5")],by="label")
  to_plot <- t_tissue_mean_summary[which(t_tissue_mean_summary$V5 %in% c("PMD","HMD")),]
  colnames(to_plot)[2:3] <- c("value","condition")
  p_list[[tissue]] <- ggplot(to_plot, aes(x = condition, y = value,fill=condition)) +
    geom_boxplot(outliers =F) +
    # labs(title = tissue_label_change(tissue),
    labs(title = tissue,
         x = NULL,
         y = "Delta(Old - Young)") +
    theme_minimal() +
    theme(
      legend.title = element_text(size = 14),  # 调整图例标题大小
      legend.text = element_text(size = 12),   # 调整图例项的字符大小
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12)  # 调整x轴文本角度和大小
    )+
    geom_signif(comparisons = list(c("PMD", "HMD")),
                textsize = 4,test = "t.test",
                map_signif_level = TRUE,
                y_position = 0.10,
                tip_length = c(1/100))
  }

plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}
combined_plot <- plot_a_list(p_list,no_of_rows = 4,no_of_cols = 6)
ggsave("result/WGBS/all_tissue_PMD_HMD_delta_change.png",combined_plot,width = 12,height = 18,type="cairo")


