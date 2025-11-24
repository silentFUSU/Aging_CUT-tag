rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(factoextra)
library(cluster)
library(umap)
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
increase <- c()
decrease <- c()
stable <- c()
non_significant <- c()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_macs_young_old_narrowpeak_summits_spm3_all_tissues_merge_diff_after_remove_batch_effect.csv"))
  rownames(df) <- df$Geneid
  t_increase <- rownames(df)[which(df$FDR.old.young < 0.05 & df$LogFC.old.young > 0)]
  t_decrease <- rownames(df)[which(df$FDR.old.young < 0.05 & df$LogFC.old.young < 0)]
  t_stable <- rownames(df)[which(df$FDR.old.young > 0.9 & abs(df$LogFC.old.young) < 0.05)]
  t_non_significant <- rownames(df)[which(df$FDR.old.young > 0.05)]
  t_non_significant <- t_non_significant[which(! t_non_significant %in% t_stable)]
  increase <- c(increase,t_increase)
  decrease <- c(decrease,t_decrease)
  stable <- c(stable,t_stable)
  non_significant <- c(non_significant,t_non_significant)
}

increase <- unique(increase)
decrease <- unique(decrease)
non_significant <- unique(non_significant)
stable <- unique(stable)
length(increase)+length(decrease)+length(non_significant)+length(stable)

non_significant <- non_significant[which(! non_significant %in% increase)]
non_significant <- non_significant[which(! non_significant %in% decrease)]
non_significant <- non_significant[which(! non_significant %in% stable)]


highly_stable <- stable[which(! stable %in% increase)]
highly_stable <- highly_stable[which(! highly_stable %in% decrease)]

Robustly_opening <- increase[which(! increase %in% decrease)]
Robustly_opening <- Robustly_opening[which(! Robustly_opening %in% stable)]

Robustly_closing<- decrease[which(! decrease %in% increase)]
Robustly_closing <- Robustly_closing[which(! Robustly_closing %in% stable)]

Secondary_opening <- increase[which(! increase %in% decrease)]
Secondary_opening <- Secondary_opening[which(! Secondary_opening %in% Robustly_opening)]

Secondary_closing<- decrease[which(! decrease %in% increase)]
Secondary_closing <- Secondary_closing[which(! Secondary_closing %in% Robustly_closing)]

Discordant <- intersect(increase,decrease)

length(highly_stable)+length(Robustly_opening)+length(Robustly_closing)+length(Secondary_opening)+length(Secondary_closing)+length(Discordant)+length(non_significant)
peaks <- read.table("data/samples/ATAC/all/ATAC/macs2_summit_all_tissues_merge/spm3/ATAC_macs_young_old_narrowpeak_summits_spm3.bed")
peaks <- peaks[which(!peaks$V4 %in% c(highly_stable,Robustly_opening,Robustly_closing,Secondary_opening,Secondary_closing,Secondary_closing,Discordant,non_significant)),]
non_significant <- c(non_significant,peaks$V4)
annotation <- data.frame(label=c(Robustly_opening,Robustly_closing,Secondary_opening,Secondary_closing,Discordant,highly_stable,non_significant),
                        condition=c(rep("Robustly_opening",length(Robustly_opening)),
                        rep("Robustly_closing",length(Robustly_closing)),
                        rep("Secondary_opening",length(Secondary_opening)),
                        rep("Secondary_closing",length(Secondary_closing)),
                        rep("Discordant",length(Discordant)),
                        rep("highly_stable",length(highly_stable)),
                        rep("non_significant",length(non_significant)))
                        )
result <- as.data.frame(table(annotation$condition))
result$proportion <- result$Freq/sum(result$Freq)*100

color <- read.table("data/samples/7_distinct_color.txt")
color <- setNames(color$V1,result$Var1)
result$condition <- "diff_peak"

result$Var1 <- factor(result$Var1,levels = c("Robustly_opening","Robustly_closing","Secondary_opening","Secondary_closing","Discordant","highly_stable","non_significant"))
p <- ggplot(result, aes(x = condition, y = proportion, fill = Var1)) +  
  geom_bar(stat = 'identity',color="white") +   
  theme_bw() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Proportion")
ggsave("result/Sup_figures/ATAC_merged_peak_diff_condition_proportion.pdf",p,height = 6,width = 4)

diff_summary <- data.frame()
diff_num <- data.frame()
for(tissue in tissues){
  df <- read.table(paste0("data/samples/ATAC/ATAC_peak_from_LMJ/common_peak_set_DAR/",tissue,"_DARs.txt"))
  t_diff_num <- data.frame(tissue=tissue_label_change(tissue),nrow(df[which(df$FDR < 0.05),]))  
  diff_num <- rbind(diff_num,t_diff_num)
  df$label <- rownames(df)
  df <- df[,c("label","logFC")]
  colnames(df)[2] <- tissue_label_change(tissue)
  
  if(nrow(diff_summary) == 0){
    diff_summary <- df
  }else{
    diff_summary <- merge(diff_summary,df,by="label",all=T)    
  }
}
rownames(diff_summary) <- diff_summary$label
diff_summary <- merge(diff_summary,annotation,by="label")
diff_summary$condition <- factor(diff_summary$condition,levels=c("Robustly_opening","Robustly_closing","Secondary_opening","Secondary_closing","Discordant","highly_stable","non_significant"))
diff_summary$row_sum <- apply(diff_summary[, -c(1, ncol(diff_summary))], 1, function(row) {
  sum(abs(row), na.rm = TRUE)
})
diff_summary <- diff_summary[order(diff_summary$condition,diff_summary$row_sum),]
rownames(diff_summary) <- diff_summary$label
colnames(diff_num)[2] <- "counts"
diff_num <- diff_num[order(diff_num$counts,decreasing = T),]
to_plot <- diff_summary[,diff_num$tissue]
to_plot[is.na(to_plot)]<- 0
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- c(seq(-1, -0.6, length.out = 40), seq(-0.59, 0.59, length.out = 20), seq(0.6, 1, length.out = 40))
annotation_row <- annotation
rownames(annotation_row) <- annotation_row$label
annotation_row <- annotation_row[,-1,drop=F]
annotation_row$condition <- factor(annotation_row$condition,levels=c("Robustly_opening","Robustly_closing","Secondary_opening","Secondary_closing","Discordant","highly_stable","non_significant"))
annotation_row_col <- list(condition=color)
pheatmap::pheatmap(to_plot,cluster_rows = F,cluster_cols = F,show_rownames = F,breaks = breaks, color = color_palette,annotation_row = annotation_row,annotation_colors = annotation_row_col)

conditions <- c("Robustly_opening","Robustly_closing","Secondary_opening","Secondary_closing","Discordant","highly_stable","non_significant")
for(tissue in tissues){
  df <- read.table(paste0("data/samples/ATAC/ATAC_peak_from_LMJ/tissue_peak_set_DAR/",tissue,"_DARs.txt"))
  peaks <- read.table(paste0("data/samples/ATAC/ATAC_peak_from_LMJ/tissue_peak_set/",tissue,".bed"))
  rownames(peaks) <- paste0(peaks$V1,":",peaks$V2,"-",peaks$V3)
  df <- merge(df,peaks,by="row.names")
  for(condition in conditions){
    dir.create(paste0("data/samples/ATAC/ATAC_peak_from_LMJ/common_peak_set_DAR/",condition),showWarnings = F,recursive = T)
    df_condition <- df[which(df$Row.names %in% annotation$label[which(annotation$condition==condition)]),]    
    write.table(df_condition,paste0("data/samples/ATAC/ATAC_peak_from_LMJ/common_peak_set_DAR/",condition,"/",tissue,".bed"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)
  }
}

