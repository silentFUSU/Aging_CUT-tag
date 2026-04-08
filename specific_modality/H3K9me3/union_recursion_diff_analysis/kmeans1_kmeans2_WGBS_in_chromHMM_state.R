rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(corrplot)
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
state_num <- 15

kmeans_df <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation_add_stable.csv",row.names = 1)
split_chr <- strsplit(as.character(kmeans_df$label), ":")  
chr_column <- sapply(split_chr, `[[`, 1)  
split_start_end <- strsplit(sapply(split_chr, `[[`, 2), "-")  
start_column <- sapply(split_start_end, `[[`, 1)  
end_column <- sapply(split_start_end, `[[`, 2)  
kmeans_df <- data.frame(chr = chr_column,start = start_column, end = end_column, cluster=kmeans_df$cluster)
kmeans_df$cluster[which(kmeans_df$chr == "chrY" & kmeans_df$cluster=="kmeans1")] <- "kmeans1-chrY"
kmeans_df$start <- as.numeric(kmeans_df$start)
kmeans_df$end <- as.numeric(kmeans_df$end)
kmeans_df <- as.data.table(kmeans_df)
setDT(kmeans_df)
setkey(kmeans_df,chr,start,end)
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver","ileum",
                  "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT"))

kmeans_summary <- data.frame()
for(kmean in c("kmeans1","kmeans2","kmeans1-chrY")){
  print(kmean)
  summary_df <- data.frame()
    for(tissue in tissues){
      print(tissue)
      if(kmean == "kmeans1-chrY" & tissue %in% c("mammarygland","ovary","uterus")){
        print("female tissues without chrY")
      }else{
        file_dir <- paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/")  
        files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_young[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
        file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
        chromHMM_young <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
        if(tissue %in% c("mammarygland","uterus","ovary")){
          chromHMM_young <- chromHMM_young[which(chromHMM_young$V1 %in% paste0("chr",c(1:19,"X"))),]
        }else{
          chromHMM_young <- chromHMM_young[which(chromHMM_young$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
        }
        
        chromHMM_young$V2 <- chromHMM_young$V2+1
        chromHMM_young <- data.table(chromHMM_young)
        setDT(chromHMM_young) 
        setkey(chromHMM_young, V1, V2, V3) 
        if(tissue %in% c("mammarygland","uterus","ovary")){
          overlaps_young <- foverlaps(kmeans_df[cluster == kmean & chr %in% paste0("chr",c(1:19,"X"))] , chromHMM_young, type = "any", nomatch = 0L)  
        }else{
          overlaps_young <- foverlaps(kmeans_df[cluster == kmean & chr %in% paste0("chr",c(1:19,"X","Y"))] , chromHMM_young, type = "any", nomatch = 0L)  
        }
        
        overlaps_young <- as.data.table(overlaps_young)
        setDT(overlaps_young)
        setkey(overlaps_young,chr,V2,V3)
        search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
        search_table <- search_table[which(search_table$tissue==tissue & search_table$age=="3M"),]
        summary <- data.frame()
        for(sample in search_table$sample_name){
          df <- fread(paste0("data/samples/WGBS/",tissue,"/bdg/",sample,"_CpG.bdg"),sep = "\t")
          setDT(df)
          setkey(df,V1,V2,V3)  
          overlaps <- foverlaps(df,overlaps_young, type = "any", nomatch = 0L)  
          overlaps[, label := paste(V1, V2, V3, V4, sep = "-")]
          result <- overlaps[, .(V4_sum = sum(i.V4), V5_sum = sum(V5)), by = label]
          result <- as.data.frame(result)
          result$methylation <- result$V4_sum/result$V5_sum*100
          result <- result[,c("label","methylation")]
          colnames(result)[2] <- sample
          if(nrow(summary)==0){
            summary <- result
          }else{
            summary <- merge(summary,result,by="label")
          }
        }
        summary$mean <- rowMeans(summary[,-1])
        summary <- summary[,c("label","mean")]
        summary$tissue <- tissue_label_change(tissue)
        summary$kmeans <- kmean
        summary$state <- sub(".*-(E\\d+)$", "\\1", summary$label)
        kmeans_summary <- rbind(kmeans_summary,summary)
    }
  }
}
kmeans_summary <- read.csv("data/samples/WGBS/all_tissues_DNA_methylation_in_H3K9me3_recursion_peaks_chromHMM_state.csv",row.names = 1)

kmeans_summary$state <- factor(kmeans_summary$state,levels = paste0("E",1:state_num))
ggplot(kmeans_summary, aes(x = state, y =mean, fill =kmeans)) +  
  geom_boxplot(alpha = 0.7,outliers = F)+
  labs(  
    title = paste0("Distribution of Young samples DNA methylation"),  
    x = "Cluster",  
    # y = "log2(mean_RPKM)"  
    y="mean_RPKM"
  ) +  
  theme_minimal() +  
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 14),  
    axis.text.y = element_text(size = 14),  
    plot.title = element_text(size = 16, face = "bold")
  ) 

# percentage <- as.data.frame(table(kmeans_summary$state))
# percentage$percentage <- percentage$Freq/sum(percentage$Freq) *100


#### delta heatmap
result <- kmeans_summary %>%
  group_by(tissue, kmeans, state) %>% 
  summarise(average_mean = mean(mean, na.rm = TRUE)) 
result <- as.data.frame(result)
to_plot <- result %>%
  group_by(tissue, state) %>% 
  filter(kmeans %in% c("kmeans1", "kmeans2")) %>%
  summarise(difference = diff(average_mean[kmeans %in% c("kmeans1", "kmeans2")])) 
to_plot <- as.data.frame(to_plot)

to_plot_heatmap <- reshape2::dcast(to_plot, tissue ~ state, value.var = "difference")
rownames(to_plot_heatmap) <- to_plot_heatmap$tissue
to_plot_heatmap <- to_plot_heatmap[,-1]
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- c(seq(-20, -6, length.out = 40), seq(-5, 5, length.out = 20), seq(6, 20, length.out = 40))
pheatmap::pheatmap(to_plot_heatmap,breaks = breaks,color = color_palette)
