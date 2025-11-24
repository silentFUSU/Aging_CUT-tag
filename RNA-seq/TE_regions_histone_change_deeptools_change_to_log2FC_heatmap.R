rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(deepToolsDownstream)
library(ggplot2)
library(patchwork)
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
    }else if(tissue_label == "Mammarygland"){
      tissue_label <- "Mammary Gland"
    }else if(tissue_label == "Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
} 

antibodys <- c("H3K9me3","H3K27me3","H3K36me3","H3K27ac","H3K4me1","H3K4me3")
tissues <- sort(c("skin","CB","spleen","heart","bladder","tongue","uterus","aorta","thymus","stomach","Hip","brain","BAT","iWAT","muscle","bonemarrow","lung","kidney","liver","testis","colon","cecum","ileum","jejunum","ovary","mammarygland","pancreas"))
class <- "ERV1_kmeans1"
##### dot plot
p_list <- list()
for(tissue in tissues){
  to_plot <- data.frame()
  for(antibody in antibodys){
    se <- importCount(paste0("result/RNA/TE/TE_with_histone/",antibody,"/matrix/",tissue,"_",antibody,"_change_in_",class,"_regions.mat.gz"))
    
    data_list <- se@assays@data@listData
    mean_list <- lapply(data_list, function(df) {
      colMeans(df, na.rm = TRUE)
    })
    mean_df <- as.data.frame(do.call(cbind, mean_list))
    pattern <- "(LLX[0-9]+|CKJ[0-9]+|DYQ[0-9]+|HM[0-9]+|HJC[0-9]+|HJC_[0-9]+|SZJ[0-9]+|NTY[0-9]+).*"
    colnames(mean_df) <- gsub(pattern, "\\1",colnames(mean_df))
    
    search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
    search_table <- search_table[which(search_table$tissue==tissue & search_table$antibody==antibody),]
    search_table <- search_table[which(search_table$sample_name %in% colnames(mean_df)),]
    mean_df_young <- mean_df[,search_table$sample_name[which(search_table$age=="3m")]]
    mean_df_young$young_mean <- rowMeans(mean_df_young)
    
    mean_df_old <- mean_df[,search_table$sample_name[which(search_table$age=="24m")]]
    mean_df_old$old_mean <- rowMeans(mean_df_old)
    
    t_to_plot <- merge(mean_df_young[,"young_mean",drop=F],mean_df_old[,"old_mean",drop=F],by="row.names")
    t_to_plot$Row.names <- factor(t_to_plot$Row.names,levels = rownames(mean_df))
    t_to_plot <- t_to_plot[order(t_to_plot$Row.names),]
    t_to_plot$log2FC <- log2(t_to_plot$old_mean/t_to_plot$young_mean)
    t_to_plot$histone <- antibody
    to_plot <- rbind(to_plot,t_to_plot)
  }
  color <- read.table("data/samples/7_distinct_color.txt")
  color$V1[1] <- "blue"
  color <- setNames(c(color$V1),antibodys)
  colnames(to_plot)[1] <- "distance"
  p_list[[tissue]] <- ggplot(to_plot, aes(x = distance, y = log2FC,color=histone)) +
    geom_point() +
    scale_color_manual(values = color)+
    labs(x = "Distance",
         y = "log2(Fold change)") +
    ggtitle(paste0(tissue_label_change(tissue)," histone modifiaction change in ",class))+
    theme_minimal()+ 
    scale_x_discrete(
      breaks = c("B5", "B995"), 
      labels = c("Start", "End")  
    ) +
    ylim(-3,3) +
    theme(
      text = element_text(size = 14), 
      axis.text.x = element_text(angle = 45, hjust = 1) )
}

H3K9me3_rank <- data.frame(
  tissue = c("lung", "CB", "BAT", "muscle", "heart", "aorta", "skin", 
             "kidney", "Hip", "brain", "liver", "tongue", "testis", 
             "bladder", "pancreas", "cecum", "spleen", "stomach", "colon", 
             "bonemarrow", "jejunum", "iWAT", "thymus", "ileum"),
  Median_Log2_FC = c(-0.36, -0.33, -0.32, -0.19, -0.19, -0.18, -0.18, -0.17, -0.15,
                     -0.10, -0.10, -0.09, -0.08, -0.06, -0.03, -0.03, -0.02, 0.01, 
                     0.02, 0.03, 0.04, 0.05, 0.05, 0.09),
  H3K9me3_rank = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 
                   21, 22, 23, 24)
)
p_list_final <- p_list[H3K9me3_rank$tissue]    
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}
combined_plot <- plot_a_list(p_list_final,no_of_rows = 4,no_of_cols = 6)
ggsave(paste0("result/RNA/TE/TE_with_histone/all_tissue_histone_change_in_",class,".png"),combined_plot,width = 40,height = 30,limitsize = FALSE)

### per tissue dot plot
color <- read.table("data/samples/7_distinct_color.txt")
color$V1[1] <- "blue"
color <- setNames(c(color$V1),antibodys)
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols, main_title) {
  patchwork::wrap_plots(
    master_list_with_plots, 
    nrow = no_of_rows, 
    ncol = no_of_cols
  ) + 
    plot_annotation(title = main_title, theme = theme(plot.title = element_text(hjust = 0.5, size = 16)))
}
for(tissue in tissues){
  p_list <- list()
  for(antibody in antibodys){
    se <- importCount(paste0("result/RNA/TE/TE_with_histone/",antibody,"/matrix/",tissue,"_",antibody,"_change_in_",class,"_regions.mat.gz"))
    data_list <- se@assays@data@listData
    mean_list <- lapply(data_list, function(df) {
      colMeans(df, na.rm = TRUE)
    })
    mean_df <- as.data.frame(do.call(cbind, mean_list))
    pattern <- "(LLX[0-9]+|CKJ[0-9]+|DYQ[0-9]+|HM[0-9]+|HJC[0-9]+|HJC_[0-9]+|SZJ[0-9]+|NTY[0-9]+).*"
    colnames(mean_df) <- gsub(pattern, "\\1",colnames(mean_df))
    
    search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
    search_table <- search_table[which(search_table$tissue==tissue & search_table$antibody==antibody),]
    search_table <- search_table[which(search_table$sample_name %in% colnames(mean_df)),]
    mean_df_young <- mean_df[,search_table$sample_name[which(search_table$age=="3m")]]
    mean_df_young$young_mean <- rowMeans(mean_df_young)
    
    mean_df_old <- mean_df[,search_table$sample_name[which(search_table$age=="24m")]]
    mean_df_old$old_mean <- rowMeans(mean_df_old)
    
    t_to_plot <- merge(mean_df_young[,"young_mean",drop=F],mean_df_old[,"old_mean",drop=F],by="row.names")
    t_to_plot$Row.names <- factor(t_to_plot$Row.names,levels = rownames(mean_df))
    t_to_plot <- t_to_plot[order(t_to_plot$Row.names),]
    t_to_plot$log2FC <- log2(t_to_plot$old_mean/t_to_plot$young_mean)
    t_to_plot$histone <- antibody
    colnames(t_to_plot)[1] <- "distance"
    p_list[[antibody]] <- ggplot(t_to_plot, aes(x = distance, y = log2FC,color=histone, group = histone)) +
      geom_point() +
      geom_line() + 
      scale_color_manual(values = color)+
      labs(x = "Distance",
           y = "log2(Fold change)") +
      ggtitle(paste0(antibody))+
      theme_minimal()+ 
      scale_x_discrete(
        breaks = c("B5", "B995"), 
        labels = c("Start", "End")  
      ) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red")
      theme(
        text = element_text(size = 14), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
  }
  combined_plot <- plot_a_list(p_list,no_of_rows = 2,no_of_cols = 3,main_title = paste0(tissue_label_change(tissue)," histone modification changed in ",class))
  ggsave(paste0("result/RNA/TE/TE_with_histone/per_tissue/log2FC/",tissue,"_histone_change_in_",class,"_log2FC.png"),combined_plot,width = 18,height = 8,limitsize = FALSE)

}
# heatmap
H3K9me3_rank <- data.frame(
  tissue = c("Lung", "Cerebellum", "BAT", "Muscle", "Heart", "Aorta", "Skin", 
             "Kidney", "Hippocampus", "Cortex", "Liver", "Tongue", "Testis", 
             "Bladder", "Pancreas", "Cecum", "Spleen", "Stomach", "Colon", 
             "Bone Marrow", "Jejunum", "iWAT", "Thymus", "Ileum"),
  Median_Log2_FC = c(-0.36, -0.33, -0.32, -0.19, -0.19, -0.18, -0.18, -0.17, -0.15,
                     -0.10, -0.10, -0.09, -0.08, -0.06, -0.03, -0.03, -0.02, 0.01, 
                     0.02, 0.03, 0.04, 0.05, 0.05, 0.09),
  H3K9me3_rank = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 
                   21, 22, 23, 24)
)
to_plot_summary <- data.frame()
for(tissue in tissues){
  to_plot <- data.frame()
  for(antibody in antibodys){
    se <- importCount(paste0("result/RNA/TE/TE_with_histone/",antibody,"/matrix/",tissue,"_",antibody,"_change_in_",class,"_regions.mat.gz"))
    data_list <- se@assays@data@listData
    mean_list <- lapply(data_list, function(df) {
      colMeans(df, na.rm = TRUE)
    })
    mean_df <- as.data.frame(do.call(cbind, mean_list))
    pattern <- "(LLX[0-9]+|CKJ[0-9]+|DYQ[0-9]+|HM[0-9]+|HJC[0-9]+|HJC_[0-9]+|SZJ[0-9]+|NTY[0-9]+).*"
    colnames(mean_df) <- gsub(pattern, "\\1",colnames(mean_df))
    
    search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
    search_table <- search_table[which(search_table$tissue==tissue & search_table$antibody==antibody),]
    search_table <- search_table[which(search_table$sample_name %in% colnames(mean_df)),]
    mean_df_young <- mean_df[,search_table$sample_name[which(search_table$age=="3m")]]
    mean_df_young$young_mean <- rowMeans(mean_df_young)
    
    mean_df_old <- mean_df[,search_table$sample_name[which(search_table$age=="24m")]]
    mean_df_old$old_mean <- rowMeans(mean_df_old)
    
    t_to_plot <- merge(mean_df_young[,"young_mean",drop=F],mean_df_old[,"old_mean",drop=F],by="row.names")
    t_to_plot$Row.names <- factor(t_to_plot$Row.names,levels = rownames(mean_df))
    t_to_plot <- t_to_plot[order(t_to_plot$Row.names),]
    t_to_plot$log2FC <- log2(t_to_plot$old_mean/t_to_plot$young_mean)
    t_to_plot <- t_to_plot[str_starts(t_to_plot$Row.names, "B"), ]
    t_to_plot <- data.frame(histone=antibody,log2FC=median(t_to_plot$log2FC))
    to_plot <- rbind(to_plot,t_to_plot)
  }
  colnames(to_plot)[2] <- tissue_label_change(tissue)
  if(nrow(to_plot_summary)==0){
    to_plot_summary <- to_plot
  }else{
    to_plot_summary <- merge(to_plot_summary,to_plot,by="histone")
  }
}

to_plot_summary_heatmap <- to_plot_summary
rownames(to_plot_summary_heatmap) <- to_plot_summary$histone 
to_plot_summary_heatmap <- to_plot_summary_heatmap[,-1]
to_plot_summary_heatmap <- as.data.frame(t(to_plot_summary_heatmap))
breaks <- c(seq(-1, -0.29, length.out = 40), seq(-0.3, 0.3, length.out = 20), seq(0.31, 1, length.out = 40))
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
pheatmap::pheatmap(to_plot_summary_heatmap,breaks = breaks,color = color_palette)

to_plot_summary_heatmap <- to_plot_summary_heatmap[H3K9me3_rank$tissue,]
pheatmap::pheatmap(to_plot_summary_heatmap,breaks = breaks,color = color_palette,cluster_rows = F)
