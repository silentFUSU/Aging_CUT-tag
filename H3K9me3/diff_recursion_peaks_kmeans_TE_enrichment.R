rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)
library(rtracklayer)
library(gridExtra)
library(grid)  
options(bitmapType="cairo")  
kmeans <- paste0("kmeans",c(1:4))
TE <- fread("~/ref_data/TE_reference/mm10_TE_metadata_20250116.bed")
setDT(TE)
setkey(TE,V1,V2,V3)
gtf <- import("~/ref_data/TE_reference/mm10_rmsk_TE.gtf", format = "gtf")
family_data <- as.data.frame(gtf[,c("gene_id","transcript_id","family_id","class_id")])

##### all class compare
summary_percentage <- data.frame()
summary_enrichment <- data.frame()
Classifications <-c("gene_id","family_id","class_id")
Classification <- "family_id"
for(kmean in kmeans){
  df <- read.table(paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/bed/",kmean,"_uinon_recursion_peaks.bed"))
  # df <- df[which(df$V1 %in% paste0("chr",c(1:19,"X"))),]
  df$label <- paste0(df$V1,":",df$V2,"-",df$V3)
  length <- data.frame(label=df$label,length=df$V3-df$V2+1)
  df <- as.data.table(df)
  setDT(df)
  setkey(df,V1,V2,V3)
  overlaps <- foverlaps(df, TE, type = "any", nomatch = 0L)  
  overlaps <- merge(overlaps,family_data[,c("transcript_id","gene_id","family_id","class_id")],by.x="V4",by.y="transcript_id")
  overlaps <- as.data.frame(overlaps)
  overlaps$TE_length <- overlaps$V3 - overlaps$V2 +1
  colnames(overlaps)[which(colnames(overlaps)==Classification)] <- "classification"
  overlaps$classification <- paste0(overlaps$classification,"-",overlaps$class_id)
  t_summary_percentage <- overlaps %>%
    group_by(label,classification) %>%
    summarize(total_TE_length = sum(TE_length, na.rm = TRUE))
  t_summary_percentage <- reshape2::dcast(t_summary_percentage, classification ~ label, value.var = "total_TE_length")
  t_summary_percentage[is.na(t_summary_percentage)] <- 0
  t_summary_percentage <- reshape2::melt(t_summary_percentage)
  colnames(t_summary_percentage)[2] <- "label"
  colnames(t_summary_percentage)[3] <- "total_TE_length"
  t_summary_percentage$condition <- kmean
  t_summary_percentage <- merge(t_summary_percentage,length,by="label")
  t_summary_percentage$TE_length_percentage <- t_summary_percentage$total_TE_length/t_summary_percentage$length*100
  
  
  t_summary_enrichment <- overlaps %>%
    group_by(label, classification) %>%
    summarise(count = n(), .groups = "drop")
  t_summary_enrichment <- reshape2::dcast(t_summary_enrichment, classification ~ label, value.var = "count")
  t_summary_enrichment[is.na(t_summary_enrichment)] <- 0
  t_summary_enrichment <- reshape2::melt(t_summary_enrichment)
  colnames(t_summary_enrichment)[2] <- "label"
  colnames(t_summary_enrichment)[3] <- "count"
  t_summary_enrichment$condition <- kmean
  t_summary_enrichment <- merge(t_summary_enrichment,length,by="label")
  t_summary_enrichment$enrichment <- t_summary_enrichment$count/t_summary_enrichment$length*1000
  summary_percentage <- rbind(summary_percentage,t_summary_percentage[,c("condition","classification","TE_length_percentage")])
  summary_enrichment <- rbind(summary_enrichment,t_summary_enrichment[,c("condition","classification","enrichment")])
}

# annotation <- family_data[,c(Classification,"class_id")]
# annotation <- unique(annotation)
# rownames(annotation) <- annotation[,1]

# summary_percentage_median <- summary_percentage %>%
#   group_by(condition,classification) %>%
#   summarise(median_TE_length_percentage = median(TE_length_percentage))

annotation <- family_data[,c(Classification,"class_id")]
annotation <- unique(annotation)
rownames(annotation) <- paste0(annotation[,1],"-",annotation$class_id)
classes <- as.data.frame(table(family_data$class_id))
color <- read.table("data/samples/20_distinct_color.txt")
color <- setNames(color$V1[1:nrow(classes)],classes$Var1)
annotation <- annotation[,-1,drop=F]
summary_percentage_mean <- summary_percentage %>%
  group_by(condition,classification) %>%
  summarise(mean_TE_length_percentage = mean(TE_length_percentage))
to_plot <- reshape2::dcast(summary_percentage_mean, classification ~ condition, value.var = "mean_TE_length_percentage")
rownames(to_plot) <- to_plot$classification
to_plot <- to_plot[,-1]
to_plot[is.na(to_plot)] <- 0
p <- pheatmap::pheatmap(to_plot,main = paste0(Classification," TE length percentage"),annotation_colors = list(class_id=color),scale = "row",cluster_cols = F,annotation_row = annotation)0
order <- p$gtable$grobs[[5]]$label

p_list <- list()
for(kmean in paste0("kmeans",c(1:4))){
  bar_to_plot <- to_plot[order,kmean,drop=F]
  bar_to_plot$TE <- rownames(bar_to_plot)
  bar_to_plot$TE <- factor(bar_to_plot$TE,levels=rev(order))
  bar_to_plot <- merge(bar_to_plot,annotation,by.x="TE",by.y="row.names")
  colnames(bar_to_plot)[2] <- "kmean"
  if(kmean == "kmeans1"){
    p_list[[kmean]] <- ggplot(bar_to_plot,mapping = aes(x=kmean,y=TE,fill= class_id))+
      scale_fill_manual(values = color)+
      geom_bar(stat = "identity", position = position_dodge2())+
      theme_bw()+xlim(0,35)+
      theme(  
        text = element_text(size = 10),  
        axis.text.y = element_text(size = 9), 
        axis.title.y = element_blank(),  
        axis.ticks.y = element_blank(),  
        legend.position = "none") +
    xlab(kmean)
  }else{
    p_list[[kmean]] <- ggplot(bar_to_plot,mapping = aes(x=kmean,y=TE,fill= class_id))+
      scale_fill_manual(values = color)+
      geom_bar(stat = "identity", position = position_dodge2())+
      theme_bw()+xlim(0,35)+
      theme(  
        text = element_text(size = 10),  
        axis.text.y = element_text(size = 9), 
        axis.title.y = element_blank(),  
        axis.ticks.y = element_blank(),  
        legend.position = "none") +
      theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),legend.position = "none")+
      xlab(kmean)
  }
}
combined_plot <- arrangeGrob(  
  grobs = p_list,  
  ncol = length(p_list),  
  widths = c(1.4,rep(1,(length(p_list)-1)))
)  
grid.draw(combined_plot) 



# summary_enrichment_median <- summary_enrichment %>%
#   group_by(condition,classification) %>%
#   summarise(median_enrichment = median(enrichment))
summary_enrichment_mean <- summary_enrichment %>%
  group_by(condition,classification) %>%
  summarise(mean_enrichment = mean(enrichment))
to_plot <- reshape2::dcast(summary_enrichment_mean, classification ~ condition, value.var = "mean_enrichment")
rownames(to_plot) <- to_plot$classification
to_plot <- to_plot[,-1]
to_plot[is.na(to_plot)] <- 0
annotation <- annotation[,-1,drop=F]
pheatmap::pheatmap(to_plot,annotation_colors = list(class_id=color),main = paste0(Classification," TE enrichment"),scale = "row",cluster_cols = F,annotation_row = annotation)

#### ERVK ERV1
family_data <- family_data[which(family_data$family_id %in% c("ERVK")),]
Classification <- "gene_id"
summary_percentage <- data.frame()
summary_enrichment <- data.frame()
for(kmean in kmeans){
  df <- read.table(paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/bed/",kmean,"_uinon_recursion_peaks.bed"))
  # df <- df[which(df$V1 %in% paste0("chr",c(1:19,"X"))),]
  df$label <- paste0(df$V1,":",df$V2,"-",df$V3)
  length <- data.frame(label=df$label,length=df$V3-df$V2+1)
  df <- as.data.table(df)
  setDT(df)
  setkey(df,V1,V2,V3)
  overlaps <- foverlaps(df, TE, type = "any", nomatch = 0L)  
  overlaps <- merge(overlaps,family_data[,c("transcript_id","gene_id","family_id","class_id")],by.x="V4",by.y="transcript_id")
  overlaps <- as.data.frame(overlaps)
  overlaps$TE_length <- overlaps$V3 - overlaps$V2 +1
  colnames(overlaps)[which(colnames(overlaps)==Classification)] <- "classification"
  overlaps$classification <- paste0(overlaps$classification,"-",overlaps$family_id)
  t_summary_percentage <- overlaps %>%
    group_by(label,classification) %>%
    summarize(total_TE_length = sum(TE_length, na.rm = TRUE))
  t_summary_percentage <- reshape2::dcast(t_summary_percentage, classification ~ label, value.var = "total_TE_length")
  t_summary_percentage[is.na(t_summary_percentage)] <- 0
  t_summary_percentage <- reshape2::melt(t_summary_percentage)
  colnames(t_summary_percentage)[2] <- "label"
  colnames(t_summary_percentage)[3] <- "total_TE_length"
  t_summary_percentage$condition <- kmean
  t_summary_percentage <- merge(t_summary_percentage,length,by="label")
  t_summary_percentage$TE_length_percentage <- t_summary_percentage$total_TE_length/t_summary_percentage$length*100
  
  
  t_summary_enrichment <- overlaps %>%
    group_by(label, classification) %>%
    summarise(count = n(), .groups = "drop")
  t_summary_enrichment <- reshape2::dcast(t_summary_enrichment, classification ~ label, value.var = "count")
  t_summary_enrichment[is.na(t_summary_enrichment)] <- 0
  t_summary_enrichment <- reshape2::melt(t_summary_enrichment)
  colnames(t_summary_enrichment)[2] <- "label"
  colnames(t_summary_enrichment)[3] <- "count"
  t_summary_enrichment$condition <- kmean
  t_summary_enrichment <- merge(t_summary_enrichment,length,by="label")
  t_summary_enrichment$enrichment <- t_summary_enrichment$count/t_summary_enrichment$length*1000
  summary_percentage <- rbind(summary_percentage,t_summary_percentage[,c("condition","classification","TE_length_percentage")])
  summary_enrichment <- rbind(summary_enrichment,t_summary_enrichment[,c("condition","classification","enrichment")])
}
# annotation <- family_data[,c(Classification,"family_id")]
# annotation <- unique(annotation)
# rownames(annotation) <- paste0(annotation[,1],"-",annotation$family_id)

# gene_id_table <- as.data.frame(table(summary_enrichment$classification))
# gene_id_table <- gene_id_table[which(gene_id_table$Freq >=100),]
# summary_enrichment <- summary_enrichment[which(summary_enrichment$classification %in% gene_id_table$Var1),]
# classes <- as.data.frame(table(family_data$class_id))
# color <- read.table("data/samples/20_distinct_color.txt")
# color <- setNames(color$V1[1:nrow(classes)],classes$Var1)
# summary_enrichment_median <- summary_enrichment %>%
#   group_by(condition,classification) %>%
#   summarise(median_enrichment = median(enrichment,na.rm = T))
summary_enrichment_mean <- summary_enrichment %>%
  group_by(condition,classification) %>%
  summarise(mean_enrichment = mean(enrichment,na.rm = T))
to_plot <- reshape2::dcast(summary_enrichment_mean, classification ~ condition, value.var = "mean_enrichment")

summary_percentage_mean <- summary_percentage %>%
  group_by(condition,classification) %>%
  summarise(mean_percentage = mean(TE_length_percentage,na.rm = T))
to_plot <- reshape2::dcast(summary_percentage_mean, classification ~ condition, value.var = "mean_percentage")

to_plot <- to_plot[order(to_plot$kmeans1,decreasing = T),]
to_plot <- to_plot[c(1:50),]
rownames(to_plot) <- to_plot$classification
to_plot <- to_plot[,-1]
to_plot[is.na(to_plot)] <- 0
# annotation <- annotation[,-1,drop=F]
# to_plot <- to_plot[apply(to_plot, 1, function(x) sum(x == 0) < 2), ]


pheatmap::pheatmap(to_plot,main = paste0(Classification," TE percentage"),
                   scale = "row",cluster_cols = F,cluster_rows = F,
                   # filename = "result/all/diff/H3K9me3/uinon_recursion_peaks_ERV1_TE_percentage.png",
                   width = 10,height = 20)

p_list <- list()
for(kmean in paste0("kmeans",c(1:4))){
  bar_to_plot <- to_plot[,kmean,drop=F]
  bar_to_plot$TE <- rownames(bar_to_plot)
  colnames(bar_to_plot)[1] <- "kmean"
  bar_to_plot$TE <- factor(bar_to_plot$TE,levels=rev(bar_to_plot$TE))
  if(kmean == "kmeans1"){
    p_list[[kmean]] <- ggplot(bar_to_plot,mapping = aes(x=kmean,y=TE))+
      geom_bar(stat = "identity", position = position_dodge2())+
      theme_bw()+xlim(0,5)+
      theme(  
        text = element_text(size = 10),  
        axis.text.y = element_text(size = 9), 
        axis.title.y = element_blank(),  
        axis.ticks.y = element_blank(),  
        legend.position = "none") +
      xlab(kmean)
  }else{
    p_list[[kmean]] <- ggplot(bar_to_plot,mapping = aes(x=kmean,y=TE))+
      geom_bar(stat = "identity", position = position_dodge2())+
      theme_bw()+xlim(0,5)+
      theme(  
        text = element_text(size = 10),  
        axis.text.y = element_text(size = 9), 
        axis.title.y = element_blank(),  
        axis.ticks.y = element_blank(),  
        legend.position = "none") +
      theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),legend.position = "none")+
      xlab(kmean)
  }
}
combined_plot <- arrangeGrob(  
  grobs = p_list,  
  ncol = length(p_list),  
  widths = c(1.4,rep(1,(length(p_list)-1)))
)  
grid.draw(combined_plot) 
