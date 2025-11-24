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
# write.csv(tissue_summary,"data/samples/WGBS/all_tissues_delta_in_200kb_bins_cross_comparison.csv")
tissue_summary <- read.csv("data/samples/WGBS/all_tissues_delta_in_200kb_bins_cross_comparison.csv",row.names = 1)
H3K9me3_tissue_order <- c("lung","CB","BAT","muscle","heart","aorta","skin","kidney","Hip","brain",
                          "liver","tongue","uterus","testis","bladder","ovary","colon","stomach","thymus","cecum","jejunum",
                          "pancreas","bonemarrow","ileum","spleen","iWAT","mammarygland")
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

bin_size <- "200kb"
bin <- read.table(paste0("~/ref_data/mm10_",bin_size,"_bins.bed"))
bin$V1 <- factor(bin$V1,levels=paste0("chr",c(1:19,"X","Y")))
bin <- bin[order(bin$V1),]
bin$V4 <- paste0("bin",1:nrow(bin))
annotation_row <- bin[,c("V1","V4")]
rownames(annotation_row) <- annotation_row$V4
annotation_row <- annotation_row[,c("V1"),drop=F]
colnames(annotation_row) <- "Chromosome"

color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sapply(sort(tissues), tissue_label_change, USE.NAMES = FALSE))

color_row <- read.table("data/samples/30_distinct_color.txt")
color_row <- setNames(color_row$V1[1:21],paste0("chr",c(1:19,"X","Y")))
annotation_color <- list(Chromosome=color_row,tissue=color)

color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- c(seq(-0.5, -0.11, length.out = 40), seq(-0.1, 0.1, length.out = 20), seq(0.11, 0.5, length.out = 40))

tissue_summary$label <- factor(tissue_summary$label,levels=bin$V4)
tissue_summary <- tissue_summary[order(tissue_summary$label),]
rownames(tissue_summary) <- as.character(tissue_summary$label)
tissue_summary <- tissue_summary[,-1]
# to_plot_H3K9me3_order <- tissue_summary[,H3K9me3_tissue_order_label]
# tissue_summary <- tissue_summary[,]
H3K9me3_tissue_order <- c("Lung","Cerebellum","BAT","Muscle","Heart","Aorta","Skin","Kidney","Hippocampus","Cortex",
                          "Liver","Tongue","Uterus","Testis","Bladder","Ovary","Colon","Stomach","Thymus","Cecum","Jejunum",
                          "Pancreas","Bone Marrow","Ileum","Spleen","iWAT","Mammary Gland")
annotation_col$tissue <- factor(annotation_col$tissue,levels = H3K9me3_tissue_order)
annotation_col <- annotation_col[order(annotation_col$tissue),,drop=F]
tissue_summary <- tissue_summary[,rownames(annotation_col)]
pheatmap::pheatmap(tissue_summary,cluster_rows = F,cluster_cols = F,show_rownames = F,breaks = breaks, color = color_palette,annotation_row = annotation_row,annotation_col = annotation_col,annotation_colors = annotation_color,main = "Whole genome 200Kb bins CpG methylation Delta(Old - Young)")

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

tissue_mean_summary$label <- factor(tissue_mean_summary$label,levels=bin$V4)
tissue_mean_summary <- tissue_mean_summary[order(tissue_mean_summary$label),]
rownames(tissue_mean_summary) <- as.character(tissue_mean_summary$label)
tissue_mean_summary <- tissue_mean_summary[,-1]
pheatmap::pheatmap(tissue_mean_summary,cluster_rows = F,cluster_cols = F,show_rownames = F,breaks = breaks, color = color_palette,annotation_row = annotation_row,annotation_col = annotation_col,annotation_colors = annotation_color,main = "Whole genome 200Kb bins CpG methylation Delta(Old - Young)")

breaks <- c(seq(-0.2, -0.06, length.out = 40), seq(-0.05, 0.05, length.out = 20), seq(0.06, 0.2, length.out = 40))
to_plot <- tissue_mean_summary[,sapply(H3K9me3_tissue_order, tissue_label_change)]
pheatmap::pheatmap(to_plot,cluster_rows = F,cluster_cols = F,show_rownames = F,breaks = breaks, color = color_palette,annotation_row = annotation_row,annotation_colors = annotation_color,main = "Whole genome 200Kb bins CpG methylation Delta(Old - Young)")

to_plot_H3K9me3_order_long <- to_plot
to_plot_H3K9me3_order_long$label <- rownames(to_plot_H3K9me3_order_long)
to_plot_H3K9me3_order_long <- reshape2::melt(to_plot_H3K9me3_order_long)
to_plot_H3K9me3_order_long$variable <- factor(to_plot_H3K9me3_order_long$variable,levels = sapply(H3K9me3_tissue_order, tissue_label_change))
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(colnames(to_plot)))
ggplot(to_plot_H3K9me3_order_long, aes(x = variable, y = value,fill=variable)) +
  geom_boxplot(outliers = F) +
  scale_fill_manual(values = color, name = "Tissue") +
  labs(title = "CpG methylation Delta (Old - Young)",
       x = NULL,
       y = "Delta") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 14),  # 调整图例标题大小
    legend.text = element_text(size = 12),   # 调整图例项的字符大小
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12)  # 调整x轴文本角度和大小
  )+
  ylim(-0.2,0.2)

to_plot_H3K9me3_order_long$condition <- "Large changes"
to_plot_H3K9me3_order_long$condition[which(to_plot_H3K9me3_order_long$variable %in% c("Ileum","Bladder","Testis","Tongue","Cecum","Pancreas","Stomach","Bone Marrow","Colon","Spleen","Thymus","Jejunum","iWAT"))] <- "Small changes"

ggplot(to_plot_H3K9me3_order_long, aes(x = condition, y = value,fill=condition)) +
  geom_boxplot(outliers =F) +
  # scale_fill_manual(values = color, name = "Tissue") +
  labs(title = "CpG methylation log2(Fold change) in H3K9me3 peaks",
       x = NULL,
       y = "log2(old/young)") +
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

###### ssGSEA
counts <- data.frame() 
search_table <- read.csv("data/samples/all/RNA_search_table.csv")
for(tissue in tissues){
  df <- read.table(paste0("data/samples/RNA/",tissue,"/combined-chrM.counts"),header = T)
  df <- df[,c(1,7:ncol(df))]
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|HM[0-9]+).*"
  colnames(df)[-1] <- gsub(pattern, "\\1",colnames(df)[-1])
  if(tissue=="mammarygland"){
    df <- df[,c("Geneid",search_table$sample_name[which(search_table$tissue=="Mammary gland")])]
  }else{
    df <- df[,c("Geneid",search_table$sample_name[which(search_table$tissue==tissue_label_change(tissue))])]
  }
  if(nrow(counts)==0){
    counts <- df
  }else{
    counts <- merge(counts,df,by="Geneid")
  }
}
rownames(counts) <- counts$Geneid
counts <- counts[,-1]

mitotic_nuclear_division <- read.csv("data/public_data/GO_term_summary_0140014.csv")
mitotic_nuclear_division <- unique(mitotic_nuclear_division$Symbol)
description <- "mitotic nuclear division"
target_genes <- mitotic_nuclear_division
search_table <- read.csv("data/samples/all/RNA_search_table.csv")
genelist <- list(score=target_genes)
counts_matrix <- as.matrix(counts)
re <- gsva(counts_matrix,genelist , method="ssgsea",ssgsea.norm=TRUE) 
re <- as.data.frame(t(re))
to_plot <- merge(re,search_table,by.x="row.names",by.y="sample_name")
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(unique(as.character(to_plot$tissue))))

average_scores <-aggregate(score ~ tissue, data = to_plot, FUN = mean)
average_scores <- average_scores[order(average_scores$score),]
to_plot$tissue <- factor(to_plot$tissue,levels = average_scores$tissue)
to_plot$age[which(to_plot$age=="3m")] <- "Young"
to_plot$age[which(to_plot$age=="24m")] <- "Old"
to_plot$age <- factor(to_plot$age,levels=c("Young","Old"))
ggplot(to_plot,aes(x=tissue,y=score,color = tissue,shape=age))+    
  geom_jitter(position = position_jitter(width = 0.2), size = 3, alpha = 0.7)+
  scale_color_manual(values = color)+
  ggtitle(paste0(description))+
  theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("Tissues")+labs(fill = "", color = "") 

medians <- to_plot_H3K9me3_order_long %>%
  group_by(variable) %>%
  summarise(median_value = median(value, na.rm = TRUE))
medians <- medians[order(medians$median_value),]
medians$rank <- c(1:nrow(medians))
medians <- as.data.frame(medians)
to_plot <- merge(to_plot,medians,by.x="tissue",by.y="variable")

colnames(to_plot)[ncol(to_plot)] <- "histone"
cor_test <- cor.test(to_plot$score,to_plot$histone,method="spearman")
average_scores <-aggregate(score ~ tissue, data = to_plot, FUN = mean)
average_scores <- average_scores[order(average_scores$score),]
average_scores <- merge(average_scores,medians,by="tissue",by.x="tissue",by.y="variable")
colnames(average_scores)[ncol(average_scores)] <- "histone"
cor_test <- cor.test(average_scores$score,average_scores$histone,method="spearman")
to_plot$histone <- factor(to_plot$histone,levels = c(24:1))

ggplot(to_plot,aes(x=histone,y=score,color = tissue,shape=age))+    
  geom_jitter(position = position_jitter(width = 0.2), size = 3, alpha = 0.7)+
  scale_color_manual(values = color)+
  ggtitle(paste0(description," with CpG DNA methylation"))+
  theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("Rank")+labs(fill = "", color = "")  




