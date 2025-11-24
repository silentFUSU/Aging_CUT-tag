rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
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


tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
tissue_summary <- data.frame()
for(tissue in tissues){
  annotation <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation_add_stable.csv",row.names = 1)
  split_names <- strsplit(annotation$label, "[:-]")
  annotation_df <- data.frame(
    chr = sapply(split_names, "[", 1),
    start = sapply(split_names, "[", 2),
    end = sapply(split_names, "[", 3),
    cluster = annotation$cluster
  )
  random_regions1 <- read.table("data/samples/all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY.bed")
  random_regions1$cluster <- "random whole genome"
  colnames(random_regions1)[1:3] <- c("chr","start","end")
  annotation_df <- rbind(annotation_df,random_regions1)
  
  random_regions2<- read.table("data/samples/all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY_out_recursion_peaks.bed")
  random_regions2$cluster <- "random out of peak"
  colnames(random_regions2)[1:3] <- c("chr","start","end")
  annotation_df <- rbind(annotation_df,random_regions2)
  if(tissue %in% c("mammarygland","ovary","uterus")){
    annotation_df <- annotation_df[which(annotation_df$chr %in% paste0("chr",c(1:19,"X"))),]
  }
  annotation_df$start <- as.numeric(annotation_df$start)
  annotation_df$end <- as.numeric(annotation_df$end)
  annotation_df$label <- paste0(annotation_df$chr,":",annotation_df$start,"-",annotation_df$end)
  annotation_df <- as.data.table(annotation_df)
  setDT(annotation_df)
  setkey(annotation_df,chr,start,end)
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  t_search_table <- search_table[which(search_table$tissue==tissue & search_table$age=="3M"),]
  summary <- data.frame()
  for(sample in t_search_table$sample_name){
    df <- fread(paste0("data/samples/WGBS/",tissue,"/bdg/",sample,"_CpG.bdg"),sep = "\t")
    setDT(df)
    setkey(df,V1,V2,V3)  
    overlaps <- foverlaps(df, annotation_df, type = "any", nomatch = 0L)  
    
    result <- overlaps[, .(V4_sum = sum(V4), V5_sum = sum(V5)), by = label]
    result <- as.data.frame(result)
    result$methylation <- result$V4_sum/result$V5_sum * 100
    result <- result[,c("label","methylation")]
    colnames(result) <- c("label",sample)
    if(nrow(summary)==0){
      summary <- result
    }else{
      summary <- merge(summary,result,by="label")
    }
  }
  summary$mean <- rowMeans(summary[,-1])
  summary <- summary[,c("label","mean")]
  colnames(summary)[2] <- tissue_label_change(tissue)
  if(nrow(tissue_summary)==0){
    tissue_summary <- summary
  }else{
    tissue_summary <- merge(tissue_summary,summary,by="label",all=T)
  }
}

to_plot <- tissue_summary
annotation <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation_add_stable.csv",row.names = 1)
# annotation <- annotation %>%
#   mutate(cluster = ifelse(grepl("chrY", label) & cluster == "kmeans1", "kmeans1-chrY", cluster))

random_regions1 <- read.table("data/samples/all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY.bed")
random_regions1$cluster <- "random whole genome"
random_regions1$label <- paste0(random_regions1$V1,":",random_regions1$V2,"-",random_regions1$V3)
rownames(random_regions1) <- random_regions1$label
random_regions1 <- random_regions1[,c("label","cluster"),drop=F]
annotation <- rbind(annotation,random_regions1)
random_regions2<- read.table("data/samples/all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY_out_recursion_peaks.bed")
random_regions2$cluster <- "random out of peak"
random_regions2$label <- paste0(random_regions2$V1,":",random_regions2$V2,"-",random_regions2$V3)
rownames(random_regions2) <- random_regions2$label
random_regions2 <- random_regions2[,c("label","cluster"),drop=F]
annotation <- rbind(annotation,random_regions2)

to_plot <- merge(to_plot,annotation,by="label")
to_plot$cluster <- factor(to_plot$cluster,levels=c("kmeans1-chrY",paste0("kmeans",1:4),"Stable","random whole genome","random out of peak"))
to_plot <- to_plot[order(to_plot$cluster),]
write.csv(to_plot,"data/samples/all/H3K9me3/recursion_peaks_diff_table/feature_table_combined_chrY/kmeans_median_DNA_methylation_detail.csv")


rownames(to_plot) <- to_plot$label
to_plot <- to_plot[,-c(1,ncol(to_plot))]

rownames(annotation) <- annotation$label
annotation <- annotation[,-1,drop=F]
breaks <- c(seq(60, 100, length.out = 100))
color_palette <- colorRampPalette(c("white", "red"))(100)  
tissues_order <- c("Kidney","Muscle","Skin","Bladder","Stomach","Heart","Hippocampus","Uterus","Liver","Aorta","Testis","Cortex","Tongue","Cerebellum","BAT","Lung","Mammary Gland",
                   "Pancreas","Bone Marrow","iWAT","Cecum","Colon","Jejunum","Spleen","Thymus","Ileum","Ovary")
to_plot <- to_plot[,tissues_order]
annotation_color <- list(cluster=setNames(c("yellow","#f6416c", "#f8f3d4", "#ffde7d", "#00b8a9","grey","blue","purple"),c("kmeans1-chrY",paste0("kmeans",1:4),"Stable","random whole genome","random out of peak")))
pheatmap::pheatmap(to_plot,color=color_palette,breaks=breaks,cluster_cols = F,cluster_rows = F,show_rownames = F,annotation_row = annotation,annotation_colors = annotation_color)

## box plot
annotation <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation_add_stable.csv",row.names = 1)
# annotation <- annotation %>%
#   mutate(cluster = ifelse(grepl("chrY", label) & cluster == "kmeans1", "kmeans1-chrY", cluster))
random_regions1 <- read.table("data/samples/all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY.bed")
random_regions1$cluster <- "random whole genome"
random_regions1$label <- paste0(random_regions1$V1,":",random_regions1$V2,"-",random_regions1$V3)
rownames(random_regions1) <- random_regions1$label
random_regions1 <- random_regions1[,c("label","cluster"),drop=F]
annotation <- rbind(annotation,random_regions1)
random_regions2<- read.table("data/samples/all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY_out_recursion_peaks.bed")
random_regions2$cluster <- "random out of peak"
random_regions2$label <- paste0(random_regions2$V1,":",random_regions2$V2,"-",random_regions2$V3)
rownames(random_regions2) <- random_regions2$label
random_regions2 <- random_regions2[,c("label","cluster"),drop=F]
annotation <- rbind(annotation,random_regions2)

to_plot_box <- tissue_summary
to_plot_box <- merge(to_plot_box,annotation,by="label")
rownames(to_plot_box) <- to_plot_box$label
to_plot_box <- to_plot_box[,-1]
to_plot_box <- reshape2::melt(to_plot_box)
to_plot_box$cluster <- factor(to_plot_box$cluster,levels=c("kmeans1-chrY",paste0("kmeans",1:4),"Stable","random whole genome","random out of peak"))
ggplot(to_plot_box, aes(x = cluster, y = value, fill = cluster)) +  
  geom_boxplot(outliers = F) + 
  theme_minimal() +   
  ylab("DNA methylation")+
  ggtitle("Young sample DNA methylation") +
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 14,angle = 45,hjust = 1),  
    axis.text.y = element_text(size = 14),  
    plot.title = element_text(size = 16, face = "bold")
  ) 
result <- to_plot_box %>%
  group_by(cluster) %>%
  summarise(median_value = median(value, na.rm = TRUE))
write.csv(result,paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/feature_table_combined_chrY/kmeans_median_DNA_methylation.csv"))

tissue_order <- c("Kidney","Muscle","Skin","Bladder","Stomach","Heart","Hippocampus","Uterus","Liver","Aorta","Testis","Cortex","Tongue","Cerebellum","BAT","Lung",
                  "Mammary Gland","Pancreas","Bone Marrow","iWAT","Cecum","Colon","Jejunum","Spleen","Thymus","Ileum","Ovary")
to_plot_box$variable <- factor(to_plot_box$variable,levels=tissue_order)
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(tissue_order))
ggplot(to_plot_box, aes(x = cluster, y = value, fill = variable)) +  
  geom_boxplot(outliers = F) + 
  scale_fill_manual(values = color) +
  theme_minimal() +   
  ylab("DNA methylation")+
  ggtitle("Young") +
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 14,angle = 45,hjust = 1),  
    axis.text.y = element_text(size = 14),  
    plot.title = element_text(size = 16, face = "bold")
  ) 

to_plot_box <- tissue_summary
to_plot_box <- to_plot_box[!grepl("chrY", to_plot_box$label), ]
to_plot_box <- merge(to_plot_box,annotation,by="label")
rownames(to_plot_box) <- to_plot_box$label
to_plot_box <- to_plot_box[,-1]
to_plot_box <- reshape2::melt(to_plot_box)
ggplot(to_plot_box, aes(x = cluster, y = value, fill = cluster)) +  
  geom_boxplot(outliers = F) + 
  theme_minimal() +   
  ylab("DNA methylation")+
  ggtitle("Young")

tissue_order <- c("Kidney","Muscle","Skin","Bladder","Stomach","Heart","Hippocampus","Uterus","Liver","Aorta","Testis","Cortex","Tongue","Cerebellum","BAT","Lung",
                  "Mammary Gland","Pancreas","Bone Marrow","iWAT","Cecum","Colon","Jejunum","Spleen","Thymus","Ileum","Ovary")
to_plot_box$variable <- factor(to_plot_box$variable,levels=tissue_order)
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(tissue_order))
ggplot(to_plot_box, aes(x = cluster, y = value, fill = variable)) +  
  geom_boxplot(outliers = F) + 
  scale_fill_manual(values = color) +
  theme_minimal() +   
  ylab("DNA methylation")+
  ggtitle("Young")

### per age
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
ages <- c("3M","24M")

tissue_summary <- data.frame()
for(age in ages){
  for(tissue in tissues){
    annotation <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation_add_stable.csv",row.names = 1)
    split_names <- strsplit(annotation$label, "[:-]")
    annotation_df <- data.frame(
      chr = sapply(split_names, "[", 1),
      start = sapply(split_names, "[", 2),
      end = sapply(split_names, "[", 3),
      cluster = annotation$cluster
    )
    random_regions1 <- read.table("data/samples/all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY.bed")
    random_regions1$cluster <- "random whole genome"
    colnames(random_regions1)[1:3] <- c("chr","start","end")
    annotation_df <- rbind(annotation_df,random_regions1)
    
    random_regions2<- read.table("data/samples/all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY_out_recursion_peaks.bed")
    random_regions2$cluster <- "random out of peak"
    colnames(random_regions2)[1:3] <- c("chr","start","end")
    annotation_df <- rbind(annotation_df,random_regions2)
    if(tissue %in% c("mammarygland","ovary","uterus")){
      annotation_df <- annotation_df[which(annotation_df$chr %in% paste0("chr",c(1:19,"X"))),]
    }
    annotation_df$start <- as.numeric(annotation_df$start)
    annotation_df$end <- as.numeric(annotation_df$end)
    annotation_df$label <- paste0(annotation_df$chr,":",annotation_df$start,"-",annotation_df$end)
    annotation_df <- as.data.table(annotation_df)
    setDT(annotation_df)
    setkey(annotation_df,chr,start,end)
    search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
    t_search_table <- search_table[which(search_table$tissue==tissue & search_table$age==age),]
    summary <- data.frame()
    for(sample in t_search_table$sample_name){
      df <- fread(paste0("data/samples/WGBS/",tissue,"/bdg/",sample,"_CpG.bdg"),sep = "\t")
      setDT(df)
      setkey(df,V1,V2,V3)  
      overlaps <- foverlaps(df, annotation_df, type = "any", nomatch = 0L)  
      
      result <- overlaps[, .(V4_sum = sum(V4), V5_sum = sum(V5)), by = label]
      result <- as.data.frame(result)
      result$methylation <- result$V4_sum/result$V5_sum * 100
      result <- result[,c("label","methylation")]
      colnames(result) <- c("label",sample)
      if(nrow(summary)==0){
        summary <- result
      }else{
        summary <- merge(summary,result,by="label")
      }
    }
    summary$mean <- rowMeans(summary[,-1])
    summary <- summary[,c("label","mean")]
    summary$tissue <- tissue_label_change(tissue)
    summary$age <- age
    tissue_summary <- rbind(tissue_summary,summary)
  }
}
result <- tissue_summary %>%
  group_by(label,age) %>%
  summarise(median_methylation = median(mean, na.rm = TRUE))
annotation <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation_add_stable.csv",row.names = 1)
# annotation <- annotation %>%
#   mutate(cluster = ifelse(grepl("chrY", label) & cluster == "kmeans1", "kmeans1-chrY", cluster))
random_regions1 <- read.table("data/samples/all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY.bed")
random_regions1$cluster <- "random whole genome"
random_regions1$label <- paste0(random_regions1$V1,":",random_regions1$V2,"-",random_regions1$V3)
rownames(random_regions1) <- random_regions1$label
random_regions1 <- random_regions1[,c("label","cluster"),drop=F]
annotation <- rbind(annotation,random_regions1)
random_regions2<- read.table("data/samples/all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY_out_recursion_peaks.bed")
random_regions2$cluster <- "random out of peak"
random_regions2$label <- paste0(random_regions2$V1,":",random_regions2$V2,"-",random_regions2$V3)
rownames(random_regions2) <- random_regions2$label
random_regions2 <- random_regions2[,c("label","cluster"),drop=F]
annotation <- rbind(annotation,random_regions2)

result <- merge(result,annotation,by="label")
result$cluster <- factor(result$cluster,levels=c("kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"))
result$age <- factor(result$age, levels = c("3M","24M"))
color <- setNames(c("#f39b7f","#4dbbd5"),c("3M","24M"))
p <- ggplot(result, aes(x = cluster, y =median_methylation, fill = age)) +  
  geom_boxplot(alpha = 0.7,outliers = F) +  
  scale_fill_manual(values = color) +
  labs(  
    title = paste0("Distribution of DNA methylation"),  
    x = "Cluster",  
    y="RPKM"
  ) + 
  theme_bw() +  
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 14),  
    axis.text.y = element_text(size = 14),  
    plot.title = element_text(size = 16, face = "bold")
  ) 
ggsave(paste0("result/Sup_figures/H3K9me3_kmeans_DNA_methylation_per_age.pdf"),p,height = 4,width = 8)
