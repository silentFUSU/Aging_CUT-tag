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
library(GenomicRanges)
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggsignif)
options(bitmapType="cairo")  

gtf <- import("~/ref_data/TE_reference/mm10_rmsk_TE.gtf", format = "gtf")
family_data <- as.data.frame(gtf[,c("gene_id","transcript_id","family_id","class_id")])
TE <- family_data[,c(1:3,7)]
colnames(TE)[1:3] <- c("V1","V2","V3")
TE <- as.data.table(TE)
setDT(TE)
setkey(TE,V1,V2,V3)

regions <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation_add_stable.csv",row.names = 1)
split_chr <- strsplit(as.character(regions$label), ":")  
chr_column <- sapply(split_chr, `[[`, 1)  
split_start_end <- strsplit(sapply(split_chr, `[[`, 2), "-")  
start_column <- sapply(split_start_end, `[[`, 1)  
end_column <- sapply(split_start_end, `[[`, 2)  
regions <- data.frame(chr = chr_column,start = start_column, end = end_column, cluster=regions$cluster)
regions$start <- as.numeric(regions$start)
regions$end <- as.numeric(regions$end)

random_regions1 <- read.table("data/samples/all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY.bed")
random_regions1$cluster <- "random whole genome"
colnames(random_regions1)[1:3] <- c("chr","start","end")
regions <- rbind(regions,random_regions1)

random_regions2 <- read.table("data/samples/all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY_out_recursion_peaks.bed")
random_regions2$cluster <- "random out of peak"
colnames(random_regions2)[1:3] <- c("chr","start","end")
regions <- rbind(regions,random_regions2)
# regions$cluster[which(regions$chr=="chrY" & regions$cluster=="kmeans1")] <- "kmeans1-chrY"

### all class compare
summary_percentage <- data.frame()
Classifications <-c("gene_id","family_id","class_id")
Classification <- "family_id"

for(cluster in sort(unique(regions$cluster))){
  df <- regions[which(regions$cluster==cluster),]
  # df <- df[which(df$chr %in% paste0("chr",c(1:19,"X"))),]
  df$label <- paste0(df$chr,":",df$start,"-",df$end)
  length <- data.frame(label=df$label,length=df$end-df$start+1)
  df <- as.data.table(df)
  setDT(df)
  setkey(df,chr,start,end)
  overlaps <- foverlaps(df, TE, type = "any", nomatch = 0L) 
  overlaps <- as.data.frame(overlaps)
  overlaps <- merge(overlaps,family_data[,c("transcript_id","gene_id","family_id","class_id")],by="transcript_id")
  overlaps$TE_length <- pmin(overlaps$V3,overlaps$end) - pmax(overlaps$V2,overlaps$start) +1
  colnames(overlaps)[which(colnames(overlaps)==Classification)] <- "classification"
  if(Classification == "family_id"){
    overlaps$classification <- paste0(overlaps$classification,"-",overlaps$class_id)
  }else if(Classification == "gene_id"){
    overlaps$classification <- paste0(overlaps$classification,"-",overlaps$family_id,"-",overlaps$class_id)
  }
  t_summary_percentage <- overlaps %>%
    group_by(label,classification) %>%
    summarize(total_TE_length = sum(TE_length, na.rm = TRUE))
  t_summary_percentage <- reshape2::dcast(t_summary_percentage, classification ~ label, value.var = "total_TE_length")
  t_summary_percentage[is.na(t_summary_percentage)] <- 0
  t_summary_percentage <- reshape2::melt(t_summary_percentage)
  colnames(t_summary_percentage)[2] <- "label"
  colnames(t_summary_percentage)[3] <- "total_TE_length"
  t_summary_percentage <- merge(t_summary_percentage,length,by="label")
  class <- unique(t_summary_percentage$classification[!is.na(t_summary_percentage$classification)])
  missed_peak <- df$label[which(! df$label %in% t_summary_percentage$label)]
  if(length(missed_peak) != 0){
    missed_peak_length <- length
    rownames(missed_peak_length) <- missed_peak_length$label
    missed_peak_length <- missed_peak_length[missed_peak,]
    missed_peak_df <- data.frame(label=rep(missed_peak,length(class)),classification=rep(class,each = length(missed_peak)),total_TE_length=0,length=rep(missed_peak_length$length,length(class)))
    t_summary_percentage <- rbind(t_summary_percentage,missed_peak_df)
  }
  t_summary_percentage$condition <- cluster
  t_summary_percentage$total_TE_length[is.na(t_summary_percentage$total_TE_length)] <- 0
  t_summary_percentage$TE_length_percentage <- t_summary_percentage$total_TE_length/t_summary_percentage$length*100
  summary_percentage <- rbind(summary_percentage,t_summary_percentage[,c("condition","classification","TE_length_percentage")])
  }

if(Classification=="class_id"){
  to_plot <- summary_percentage[which(summary_percentage$classification=="LTR"),]
  write.csv(to_plot,paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/feature_table_combined_chrY/kmeans_median_LTR_median_TE_length_percentage_detail.csv"))
  result <- to_plot %>%
    group_by(condition) %>%
    summarise(TE_length_percentage = median(TE_length_percentage, na.rm = TRUE))
  write.csv(result,paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/feature_table_combined_chrY/kmeans_median_LTR_median_TE_length_percentage.csv"))
  
  to_plot$condition <- factor(to_plot$condition, levels=c("kmeans1-chrY","kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"))
  p <- ggplot(to_plot, aes(x = condition, y = TE_length_percentage,fill=condition)) +
    geom_boxplot() +
    labs(x = NULL, y = "TE_length_percentage", title = "ERV1-LTR") +
    theme_minimal()+
    theme(
      axis.text.x = element_text(size = 12,face = "bold", color = "black", angle = 45, hjust = 1),  
      axis.text.y = element_text(size = 12,face = "bold", color = "black"),  
      axis.title.x = element_text(size = 14,face = "bold", color = "black"), 
      axis.title.y = element_text(size = 14,face = "bold", color = "black"), 
      legend.text = element_text(size = 12),  
      panel.background = element_blank(),  
      panel.grid.major = element_line(size = 0.1, linetype = 'solid', color = "grey"),
      panel.grid.minor = element_line(size = 0.1, linetype = 'solid', color = "lightgrey"),
      plot.title = element_text(size = 16, face = "bold")
    )+
    geom_signif(
      comparisons = list(c("kmeans1", "kmeans2"), c("kmeans1", "kmeans3"),c("kmeans1", "kmeans4"),c("kmeans1", "Stable"),c("kmeans1", "random out of peak"),c("kmeans1", "random whole genome")),
      textsize = 4,test = "t.test",
      map_signif_level = TRUE,
      y_position = c(40, 42, 44, 46, 48, 50,52),
      tip_length = c(1/100,1/100,1/100,1/100,1/100,1/100)
    )
  # geom_signif(
  #   comparisons = list(c("kmeans1", "kmeans2"), c("kmeans1", "kmeans3"),c("kmeans1", "kmeans4"),c("kmeans1", "random out of peak"),c("kmeans1", "random whole genome")),
  #   textsize = 4,test = "t.test",
  #   map_signif_level = TRUE,
  #   y_position = c(15, 17, 19, 21, 23, 25),
  #   tip_length = c(1/100,1/100,1/100,1/100,1/100)
  # )
  # ggsave("result/figures/LTR_H3K9me3_kmeans_box_plot.pdf",p,width = 10,height = 4)
}

if(Classification=="class_id"){
  annotation <- family_data[,c("class_id"),drop=F]
}else if(Classification=="family_id"){
  annotation <- family_data[,c(Classification,"class_id")]
}else{
  annotation <- family_data[,c(Classification,"family_id","class_id")]
}
annotation <- unique(annotation)
if(Classification == "family_id"){
  rownames(annotation) <- paste0(annotation[,1],"-",annotation$class_id)
  annotation <- annotation[,-1,drop=F]
}else if(Classification == "gene_id"){
  rownames(annotation) <- paste0(annotation[,1],"-",annotation$family_id,"-",annotation$class_id)
  annotation <- annotation[,-c(1:2),drop=F]
}else{
  rownames(annotation) <- annotation[,1]
}
classes <- as.data.frame(table(family_data$class_id))
color <- read.table("data/samples/20_distinct_color.txt")
color <- setNames(color$V1[1:nrow(classes)],classes$Var1)

summary_percentage_mean <- summary_percentage %>%
  group_by(condition,classification) %>%
  summarise(mean_TE_length_percentage = mean(TE_length_percentage))
to_plot <- reshape2::dcast(summary_percentage_mean, classification ~ condition, value.var = "mean_TE_length_percentage")
rownames(to_plot) <- to_plot$classification
to_plot <- to_plot[,-1]
to_plot[is.na(to_plot)] <- 0
color_palette <- colorRampPalette(c("blue", "white", "red"))(100) 
to_plot <- to_plot[!grepl("\\?|Unknown", rownames(to_plot)), ]
annotation <- annotation[!grepl("\\?|Unknown", rownames(annotation)),,drop=F]
to_plot <- to_plot[,c("kmeans1-chrY","kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome")]
to_plot <- to_plot[order(to_plot$kmeans1,decreasing = T),]
# write.csv(to_plot,"data/samples/all/H3K9me3/recursion_peaks_diff_table/gene_id_level_TE_length_percentage_chrY.csv")
# p <- pheatmap::pheatmap(to_plot,main = paste0(Classification," TE length percentage"),color = color_palette,annotation_colors = list(class_id=color),scale = "row",cluster_cols = F,cluster_rows = T,annotation_row = annotation,border_color = "black",width = 6,height = 8,filename = "result/Sup_figures/family_level_TE_length_percentage.pdf")
p <- pheatmap::pheatmap(to_plot,main = paste0(Classification," TE length percentage"),color = color_palette,annotation_colors = list(class_id=color),scale = "row",cluster_cols = F,cluster_rows = T,annotation_row = annotation,border_color = "black")

order <- p$gtable$grobs[[5]]$label

p_list <- list()
for(kmean in c("kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome")){
  bar_to_plot <- to_plot[order,kmean,drop=F]
  bar_to_plot$TE <- rownames(bar_to_plot)
  bar_to_plot$TE <- factor(bar_to_plot$TE,levels=rev(order))
  bar_to_plot <- merge(bar_to_plot,annotation,by.x="TE",by.y="row.names")
  colnames(bar_to_plot)[2] <- "kmean"
  if(kmean == "kmeans1"){
    p_list[[kmean]] <- ggplot(bar_to_plot,mapping = aes(x=kmean,y=TE,fill= class_id))+
      scale_fill_manual(values = color)+
      geom_bar(stat = "identity", position = position_dodge2())+
      theme_bw()+xlim(0,40)+
      theme(  
        text = element_text(size = 10),  
        axis.text.y = element_text(size = 9), 
        legend.position = "none") +
      xlab(kmean)
  }else{
    p_list[[kmean]] <- ggplot(bar_to_plot,mapping = aes(x=kmean,y=TE,fill= class_id))+
      scale_fill_manual(values = color)+
      geom_bar(stat = "identity", position = position_dodge2())+
      theme_bw()+xlim(0,40)+
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
ggsave("result/Sup_figures/family_level_TE_length_percentage_bar_plot.pdf",combined_plot,width = 24,height = 8)
