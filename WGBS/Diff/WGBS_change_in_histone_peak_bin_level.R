rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(stringr)
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
antibody <- "H3K9me3"
delta_summary <- data.frame()
test_summary <- data.frame()
for(tissue in tissues){
  if(antibody %in% c("H3K9me3","H3K27me3","H3K36me3")){
    bin_size <- "10kb"
    bin_region <- read.table("~/ref_data/mm10_10kb_bins.bed")
    peak <- read.table(paste0("data/samples/",tissue,"/H3K9me3/bed/H3K9me3_young_old_merge-W5000-G10000-E100.bed"))
  }else{
    bin_size <- "1kb"
    bin_region <- read.table("~/ref_data/mm10_1kb_bins.bed")
  }
  bin_region <- as.data.table(bin_region)
  setDT(bin_region)
  setkey(bin_region,V1,V2,V3)
  peak <- as.data.table(peak)
  setDT(peak)
  setkey(peak,V1,V2,V3)
  overlaps <- as.data.frame(foverlaps(peak,bin_region, type = "any", nomatch = 0L))
  
  df <- read.csv(paste0("data/samples/WGBS/",tissue,"/compress2bin/",bin_size,"_bins_all_depth.csv"))
  # df <- df[which(df$total_V5>15),]
  colnames(df)[6] <- "sample_name"
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  df <- merge(df,search_table[,c(3,4,5)],by="sample_name")  
  df$age[which(df$age=="3M")] <- "young"
  df$age[which(df$age=="24M")] <- "old"
  df$age <- factor(df$age,levels=c("young","old"))
  df$percent <- 100*(df$total_V4/df$total_V5)
  
  result <- df %>% 
    group_by(age, label) %>%
    summarise(avg_percent = mean(percent, na.rm = TRUE))
  young_result <- result[which(result$age=="young"),]
  colnames(young_result)[3] <- "young_methylation"
  old_result <- result[which(result$age=="old"),]
  colnames(old_result)[3] <- "old_methylation"
  
  result <- merge(young_result[,2:3],old_result[,2:3],by="label")
  result$delta <- result$old_methylation - result$young_methylation
  result_to_plot <- result
  result_to_plot$condition <- "peak region"
  result_to_plot$condition[which(! result_to_plot$label %in% overlaps$V4)] <- "outside region"
  result_to_plot$condition <- factor(result_to_plot$condition, levels = c("peak region","outside region"))
  result_to_plot$tissue <- tissue_label_change(tissue)
  delta_summary <- rbind(delta_summary,result_to_plot) 
  t <- t.test(result_to_plot$delta[which(result_to_plot$condition=="peak region")],result_to_plot$delta[which(result_to_plot$condition=="outside region")])
  t_test_summary <- data.frame(peak_region=t$estimate[[1]],out_region=t$estimate[[2]],p_value=t$p.value,tissue=tissue_label_change(tissue))
  test_summary <- rbind(test_summary,t_test_summary)
  # ggplot(result_to_plot, aes(x = condition, y = delta,fill=condition)) +  
  #   scale_fill_brewer(palette = "Pastel1") +
  #   geom_boxplot(outliers = F) +  
  #   theme_minimal() +
  #   theme(text = element_text(size = 20)) +
  #   guides(fill = FALSE)+
  #   labs(title = paste0(tissue_label_change(tissue)), x = NULL, y = "CG% delta")+
  #   annotate("text", x = Inf, y = -Inf, label = paste("p-value =",  format(t$p.value, scientific = TRUE, digits = 3)  ),   
  #            hjust = 1.1, vjust = -1.1, size = 5, colour = "red")
}
tissues_order <- c("Mammary Gland","Cecum","Thymus","Uterus","iWAT","Stomach","Skin","Spleen","Muscle","Bone Marrow","Liver","Ileum","Testis",
                   "Cortex","Jejunum","Tongue","Hippocampus","Colon","Bladder","Aorta","Cerebellum","Lung","Heart","Kidney","BAT","Ovary","Pancreas")

to_plot <- delta_summary
to_plot$tissue <- factor(to_plot$tissue,levels = tissues_order)
ggplot(to_plot, aes(x = tissue, y = delta,fill=condition)) +
  scale_fill_brewer(palette = "Pastel1") +
  geom_boxplot(outliers = F) +
  theme_minimal() +
  theme(text = element_text(size = 20)) +
  labs(x = NULL, y = "CG% delta")+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) 

test_summary_filter <- test_summary[which(abs(test_summary$peak_region) > abs(test_summary$out_region)),]
test_summary_to_plot <- test_summary
test_summary_to_plot$peak_region <- abs(test_summary_to_plot$peak_region)
test_summary_to_plot$out_region <- abs(test_summary_to_plot$out_region)
test_summary_to_plot <- test_summary_to_plot[,-3]
test_summary_to_plot <- reshape2::melt(test_summary_to_plot)
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(tissues_order))
p <- ggplot(test_summary_to_plot, aes(x = variable, y = value,fill=variable)) +
  geom_boxplot(outliers = F) +
  scale_fill_brewer(palette = "Pastel1") +
  geom_point(aes(color = tissue), size = 2) + 
  scale_color_manual(values = color) +
  labs(x = NULL, y = "abs(Delta)")+
  guides(fill = F,color=F)+ 
  theme_minimal() +  
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 14,angle = 45,hjust = 1),  
    axis.text.y = element_text(size = 14),  
    plot.title = element_text(size = 16, face = "bold")
  )

t.test(test_summary_to_plot$value[which(test_summary_to_plot$variable=="peak_region")],test_summary_to_plot$value[which(test_summary_to_plot$variable=="out_region")])
ggsave("result/figures/DNA_methylation_change_in_H3K9me3_peaks_boxplot2.pdf",p,width = 2,height = 4)


summary <- read.csv("data/samples/WGBS/CG_manual.csv")
summary$tissue <- sapply(summary$tissue,tissue_label_change)
summary <- merge(summary,search_table[,c("sample_name","age")],by.x="sample",by.y="sample_name")
result <- summary %>%
  group_by(tissue, age) %>%
  summarise(mean_CG = mean(CG, na.rm = TRUE), .groups = 'drop')
result$mean_CG <- result$mean_CG *100
result_young <- result[which(result$age=="3M"),]
result_old <- result[which(result$age=="24M"),]
result <- merge(result_young,result_old,by="tissue")
result$delta <- result$mean_CG.y - result$mean_CG.x
tissues_filter <- result$tissue[which(abs(result$delta)>1)]
test_summary_to_plot2 <- test_summary_to_plot[which(test_summary_to_plot$tissue %in% tissues_filter),]
p <- ggplot(test_summary_to_plot2, aes(x = variable, y = value,fill=variable)) +
  geom_boxplot(outliers = F) +
  scale_fill_brewer(palette = "Pastel1") +
  geom_point(aes(color = tissue), size = 2) + 
  scale_color_manual(values = color) +
  labs(x = NULL, y = "abs(Delta)")+
  guides(fill = F,color=F)+ 
  theme_minimal() +  
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 14,angle = 45,hjust = 1),  
    axis.text.y = element_text(size = 14),  
    plot.title = element_text(size = 16, face = "bold")
  )
t.test(test_summary_to_plot2$value[which(test_summary_to_plot2$variable=="peak_region")],test_summary_to_plot2$value[which(test_summary_to_plot2$variable=="out_region")])
df1 <- test_summary_to_plot2[which(test_summary_to_plot2$variable=="peak_region"),]
df2 <- test_summary_to_plot2[which(test_summary_to_plot2$variable=="out_region"),]
df <- merge(df1,df2,by="tissue")
wilcox.test(df$value.x,df$value.y,paired = T)
ggsave("result/figures/DNA_methylation_change_in_H3K9me3_peaks_boxplot_delta_larger_1.pdf",p,width = 2,height = 4)




