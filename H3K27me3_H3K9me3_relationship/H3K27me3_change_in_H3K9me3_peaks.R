rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(data.table)
H3K27me3_tissues_label <- c("Bladder","Kidney","Liver","Pancreas","Stomach",
                            "Bone Marrow","Cecum","Colon","Ileum","Jejunum","iWAT","Mammary Gland","Spleen","Thymus",
                            "Aorta","BAT","Cerebellum","Cortex","Hippocampus","Heart","Lung","Muscle","Ovary","Skin","Tetis","Tongue","Uterus")
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}
tissue <- "CB"
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
#########opposite_change_region change in H3K27me3 or H3K9me3 peaks region
tissues <- sort(c("mammarygland","BAT","CB","lung","kidney","aorta","brain",
                  "spleen","thymus","skin","bladder","bonemarrow","Hip",
                  "muscle","iWAT","jejunum","uterus","ovary","liver"))
to_plot <- data.frame()
for(tissue in tissues){
  opposite_change_region <- read.table(paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/10kb_all_significant_second_quadrant_after_remove_batch_effect.bed"))
  H3K27me3_peak_region <- read.table(paste0("data/samples/",tissue,"/H3K27me3/bed/H3K27me3_10kb_in_young_merge-W1000-G3000-E100.bed"))
  opposite_change_region$condition <- "Out of H3K27me3 peaks"
  opposite_change_region$condition[which(opposite_change_region$V4 %in% H3K27me3_peak_region$V4)] <- "In H3K27me3 peaks"
  t_to_plot <- as.data.frame(table(opposite_change_region$condition))
  t_to_plot$percent <- t_to_plot$Freq/sum(t_to_plot$Freq) *100
  t_to_plot$tissue <- tissue_label_change(tissue)
  if(nrow(to_plot)==0){
    to_plot <- t_to_plot
  }else{
    to_plot <- rbind(to_plot,t_to_plot)
  }
}
to_plot$position <- 100
to_plot$position[which(to_plot$Var1=="Out of H3K27me3 peaks")] <- to_plot[which(to_plot$Var1=="Out of H3K27me3 peaks"),"percent"]
to_plot_peaks <- to_plot[which(to_plot$Var1=="In H3K27me3 peaks"),]
to_plot_peaks <- to_plot_peaks[order(to_plot_peaks$percent),]
to_plot$tissue <- factor(to_plot$tissue,levels=to_plot_peaks$tissue)
total_label <- to_plot %>%  
  group_by(tissue) %>%  
  summarise(Freq_sum = sum(Freq))  
total_label$position <- 0
color <- setNames(c("#009980","#838B8B"),c("In H3K27me3 peaks","Out of H3K27me3 peaks"))
ggplot(to_plot, aes(x = tissue, y = percent, fill = Var1)) +  
  geom_bar(stat = 'identity',colour = "white") +   
  theme_minimal() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  geom_text(data = to_plot,   
            aes(label = Freq, y = position),   
            color = "black", size = 5, vjust = 0.5) +
  ylab("Proportion")+
  ggtitle("H3K27me3 increase and H3K9me3 decrease bins")


############### H3K27me3 change whether in H3K9me3 peaks
tissues <- sort(c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
             "thymus","skin","bladder","bonemarrow","Hip","heart",
             "muscle","jejunum","uterus","ovary","liver","tongue",
             "cecum","colon","testis","stomach","pancreas","iWAT","ileum"))
condition <- "Up"
if(condition=="Up"){
  condition_label <- "increase"
}else{
  condition_label <- "decrease" 
}
to_plot <- data.frame()
for(tissue in tissues){
  H3K27me3_change_region <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  H3K27me3_change_region <- H3K27me3_change_region[which(H3K27me3_change_region$Significant==condition),]
  H3K9me3_peak_region <- read.table(paste0("data/samples/",tissue,"/H3K9me3/bed/H3K9me3_10kb_in_young_merge-W1000-G3000-E100.bed"))
  H3K27me3_change_region$condition <- "Out of H3K9me3 peaks"
  H3K27me3_change_region$condition[which(H3K27me3_change_region$Geneid %in% H3K9me3_peak_region$V4)] <- "In H3K9me3 peaks"
  
  t_to_plot <- as.data.frame(table(H3K27me3_change_region$condition))
  t_to_plot$percent <- t_to_plot$Freq/sum(t_to_plot$Freq) *100
  t_to_plot$tissue <- tissue_label_change(tissue)
  if(nrow(to_plot)==0){
    to_plot <- t_to_plot
  }else{
    to_plot <- rbind(to_plot,t_to_plot)
  }
}
to_plot_peaks <- to_plot[which(to_plot$Var1=="In H3K9me3 peaks"),]
to_plot_peaks <- to_plot_peaks[order(to_plot_peaks$percent),]
to_plot$position <- 100
to_plot$position[which(to_plot$Var1=="Out of H3K9me3 peaks")] <- to_plot[which(to_plot$Var1=="Out of H3K9me3 peaks"),"percent"]
# to_plot$tissue <- factor(to_plot$tissue, levels=to_plot_peaks$tissue)
color <- setNames(c("#009980","#838B8B"),c("In H3K9me3 peaks","Out of H3K9me3 peaks"))
ggplot(to_plot, aes(x = tissue, y = percent, fill = Var1)) +  
  geom_bar(stat = 'identity',colour = "white") +   
  theme_minimal() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Proportion")+
  geom_text(data = to_plot,   
            aes(label = Freq, y = position),   
            color = "black", size = 5, vjust = 0.5) +
  ggtitle(paste0("H3K27me3 ",condition_label," regions"))

p_list <- list()
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff.csv")
  search_table <- search_table[which(search_table$tissue==tissue & search_table$antibody=="H3K27me3"),]
  search_table$age[which(search_table$age=="3m")] <- "young"
  search_table$age[which(search_table$age=="24m")] <- "old"
  H3K27me3_change_region <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_10kb_bins_diff.csv"))
  H3K27me3_change_region <- H3K27me3_change_region[which(H3K27me3_change_region$Significant=="Up"),]
  H3K9me3_peak_region <- read.table(paste0("data/samples/",tissue,"/H3K9me3/bed/H3K9me3_10kb_in_young_merge-W1000-G3000-E100.bed"))
  H3K27me3_change_region$condition <- "Out of H3K9me3 peaks"
  H3K27me3_change_region$condition[which(H3K27me3_change_region$Geneid %in% H3K9me3_peak_region$V4)] <- "In H3K9me3 peaks"
  H3K27me3_change_region <- H3K27me3_change_region[,c(paste(search_table$sample_name,search_table$age,search_table$mouse_ID,sep = "."),"condition","Geneid")]
  colnames(H3K27me3_change_region)[which(colnames(H3K27me3_change_region)!="condition")] <- sapply(strsplit(colnames(H3K27me3_change_region)[which(colnames(H3K27me3_change_region)!="condition")], "\\."), `[`, 1)
  H3K27me3_change_region <- reshape2::melt(H3K27me3_change_region)
  H3K27me3_change_region <- merge(H3K27me3_change_region,search_table,by.x="variable",by.y="sample_name")
  H3K27me3_change_region$age <- factor(H3K27me3_change_region$age, levels=c("young","old"))
  p_list[[tissue]] <- ggplot(H3K27me3_change_region,aes(x=condition,y=log2(value),fill=age))+
    geom_boxplot()+
    xlab(NULL)+
    ylab("H3K27me3 log2(CPM)")+
    theme_bw()+
    theme(axis.title.x = element_blank(),
          text = element_text(size = 13),legend.title = element_blank())+
    ggtitle(tissue_label_change(tissue)," H3K27me3 increase regions")
}
combined_plot <- plot_a_list(p_list, 4, 7)
ggsave("result/all/H3K27me3_H3K9me3/H3K27me3_increase_H3K9me3_decrease_region_H3K27me3_in_H3K9me3_peaks.png",combined_plot,width = 35,height = 20,type="cairo")

ggplot(H3K27me3_change_region,aes(x=condition,y=log2(value),fill=age))+
  geom_boxplot()+
  xlab(NULL)+
  ylab("H3K27me3 log2(CPM)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        text = element_text(size = 13),legend.title = element_blank()) +
  ggtitle(tissue_label_change(tissue))


############### H3K27me3 increase in H3K9me3 peaks, How H3K9me3 change
tissues <- sort(c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
                  "thymus","skin","bladder","bonemarrow","Hip","heart",
                  "muscle","jejunum","uterus","ovary","liver","tongue",
                  "cecum","colon","testis","stomach","pancreas","iWAT","ileum"))
to_plot <- data.frame()
for(tissue in tissues){
  H3K27me3_increase_region <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  H3K27me3_increase_region <- H3K27me3_increase_region[which(H3K27me3_increase_region$Significant=="Up"),]
  H3K9me3_peak_region <- read.table(paste0("data/samples/",tissue,"/H3K9me3/bed/H3K9me3_10kb_in_young_merge-W1000-G3000-E100.bed"))
  H3K27me3_increase_region$condition <- "Out of H3K9me3 peaks"
  H3K27me3_increase_region$condition[which(H3K27me3_increase_region$Geneid %in% H3K9me3_peak_region$V4)] <- "In H3K9me3 peaks"
  
  H3K9me3 <- read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  H3K9me3 <- H3K9me3[which(H3K9me3$Geneid %in% H3K27me3_increase_region$Geneid[which(H3K27me3_increase_region$condition=="In H3K9me3 peaks")]),]
  if(nrow(H3K9me3[which(H3K9me3$Significant!="Stable" ),])>100){
    # H3K9me3 <- H3K9me3[which(H3K9me3$Significant !="Stable"),]
    t_to_plot <- as.data.frame(table(H3K9me3$Significant))
    t_to_plot$percent <- t_to_plot$Freq/sum(t_to_plot$Freq) *100
    t_to_plot$tissue <- tissue_label_change(tissue)
    if(nrow(to_plot)==0){
      to_plot <- t_to_plot
    }else{
      to_plot <- rbind(to_plot,t_to_plot)
    }
  }
}
to_plot$Var1 <- as.character(to_plot$Var1)
to_plot$Var1[which(to_plot$Var1=="Down")] <- "H3K9me3 decrease"
to_plot$Var1[which(to_plot$Var1=="Up")] <- "H3K9me3 increase"
to_plot$Var1[which(to_plot$Var1=="Stable")] <- "H3K9me3 stable"
color <- setNames(c("#009980","#E69900","#838B8B"),c("H3K9me3 decrease","H3K9me3 increase","H3K9me3 stable"))
to_plot$Var1 <- factor(to_plot$Var1, levels=c("H3K9me3 decrease","H3K9me3 increase","H3K9me3 stable"))
freq_sum <- to_plot %>%  
  group_by(tissue) %>%  
  summarise(Freq_sum = sum(Freq))  
ggplot(to_plot, aes(x = tissue, y = percent, fill = Var1)) +  
  geom_bar(stat = 'identity',colour = "white") +   
  theme_minimal() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Proportion")+
  geom_text(data = freq_sum, aes(x = tissue, y = -5, label = Freq_sum),   
            inherit.aes = FALSE,  # 避免继承主 ggplot 中的 aes 映射  
            vjust = 1, size = 5, color = "black")  +
  ggtitle("H3K27me3 increase in H3K9me3 peak regions")

############### H3K27me3 increase, H3K9me3 whether decrease
tissues <- c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
             "thymus","skin","bladder","bonemarrow","Hip","heart",
             "muscle","jejunum","uterus","ovary","liver","tongue",
             "cecum","colon","testis","stomach","pancreas","iWAT","ileum")
to_plot <- data.frame()
for(tissue in tissues){
  H3K27me3_increase_region <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  H3K27me3_increase_region <- H3K27me3_increase_region[which(H3K27me3_increase_region$Significant=="Up"),]
  H3K9me3_decrease_region <- read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  H3K9me3_decrease_region <- H3K9me3_decrease_region[which(H3K9me3_decrease_region$Significant=="Down"),]
  H3K27me3_increase_region$H3K9me3_condition <- "Other"
  H3K27me3_increase_region$H3K9me3_condition[which(H3K27me3_increase_region$Geneid %in% H3K9me3_decrease_region$Geneid)] <- "H3K9me3 decrease"
  t_to_plot <- as.data.frame(table(H3K27me3_increase_region$H3K9me3_condition))
  t_to_plot$percent <- t_to_plot$Freq/sum(t_to_plot$Freq) *100
  t_to_plot$tissue <- tissue_label_change(tissue)
  if(nrow(to_plot)==0){
    to_plot <- t_to_plot
  }else{
    to_plot <- rbind(to_plot,t_to_plot)
  }
}
color <- setNames(c("#009980","#838B8B"),c("H3K9me3 decrease","Other"))
ggplot(to_plot, aes(x = tissue, y = percent, fill = Var1)) +  
  geom_bar(stat = 'identity',colour = "white") +   
  theme_minimal() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Proportion")+
  ggtitle("H3K27me3 increase regions")

############### H3K9me3 decrease whether in H3K27me3 peaks
tissues <- sort(c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
                  "thymus","skin","bladder","bonemarrow","Hip","heart",
                  "muscle","jejunum","uterus","ovary","liver","tongue",
                  "cecum","colon","testis","stomach","pancreas","iWAT","ileum"))
to_plot <- data.frame()
for(tissue in tissues){
  H3K9me3_decrease_region <- read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_10kb_bins_diff.csv"))
  H3K9me3_decrease_region <- H3K9me3_decrease_region[which(H3K9me3_decrease_region$Significant=="Down"),]
  if(nrow(H3K9me3_decrease_region) > 100) {
    H3K27me3_peak_region <- read.table(paste0("data/samples/",tissue,"/H3K27me3/bed/H3K27me3_10kb_in_young_merge-W1000-G3000-E100.bed"))
    H3K9me3_decrease_region$condition <- "Out of H3K27me3 peaks"
    H3K9me3_decrease_region$condition[which(H3K9me3_decrease_region$Geneid %in% H3K27me3_peak_region$V4)] <- "In H3K27me3 peaks"
    t_to_plot <- as.data.frame(table(H3K9me3_decrease_region$condition))
    t_to_plot$percent <- t_to_plot$Freq/sum(t_to_plot$Freq) *100
    t_to_plot$tissue <- tissue_label_change(tissue)
    if(nrow(to_plot)==0){
      to_plot <- t_to_plot
    }else{
      to_plot <- rbind(to_plot,t_to_plot)
    }
  }
}
color <- setNames(c("#009980","#838B8B"),c("In H3K27me3 peaks","Out of H3K27me3 peaks"))
ggplot(to_plot, aes(x = tissue, y = percent, fill = Var1)) +  
  geom_bar(stat = 'identity',colour = "white") +   
  theme_minimal() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Proportion")+
  ggtitle("H3K9me3 decrease regions")

############### H3K9me3 decrease, H3K27me3 whether increase
tissues <- c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
             "thymus","skin","bladder","bonemarrow","Hip","heart",
             "muscle","jejunum","uterus","ovary","liver","tongue",
             "cecum","colon","testis","stomach","pancreas","iWAT","ileum")
to_plot <- data.frame()
for(tissue in tissues){
  H3K27me3_increase_region <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  H3K27me3_increase_region <- H3K27me3_increase_region[which(H3K27me3_increase_region$Significant=="Up"),]
  H3K9me3_decrease_region <- read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  H3K9me3_decrease_region <- H3K9me3_decrease_region[which(H3K9me3_decrease_region$Significant=="Down"),]
  if(nrow(H3K9me3_decrease_region) > 100) {
    H3K9me3_decrease_region$H3K27me3_condition <- "Other"
    H3K9me3_decrease_region$H3K27me3_condition[which(H3K9me3_decrease_region$Geneid %in% H3K27me3_increase_region$Geneid)] <- "H3K27me3 increase"
    t_to_plot <- as.data.frame(table(H3K9me3_decrease_region$H3K27me3_condition))
    t_to_plot$percent <- t_to_plot$Freq/sum(t_to_plot$Freq) *100
    t_to_plot$tissue <- tissue_label_change(tissue)
    if(nrow(to_plot)==0){
      to_plot <- t_to_plot
    }else{
      to_plot <- rbind(to_plot,t_to_plot)
    }
  }
}
color <- setNames(c("#009980","#838B8B"),c("H3K27me3 increase","Other"))
ggplot(to_plot, aes(x = tissue, y = percent, fill = Var1)) +  
  geom_bar(stat = 'identity',colour = "white") +   
  theme_minimal() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Proportion")+
  ggtitle("H3K9me3 decrease regions")

############### H3K9me3 peak region H3K27me3 condition
tissues <- c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
             "thymus","skin","bladder","bonemarrow","Hip","heart",
             "muscle","jejunum","uterus","ovary","liver","tongue",
             "cecum","colon","testis","stomach","pancreas","iWAT","ileum")
to_plot <- data.frame()
for(tissue in tissues){
  H3K27me3 <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  H3K9me3_peak_region <- read.table(paste0("data/samples/",tissue,"/H3K9me3/bed/H3K9me3_10kb_in_young_merge-W1000-G3000-E100_compress.bed"))
  H3K9me3_peak_region$condition <- "H3K27me3 stable"
  H3K9me3_peak_region$condition[which(H3K9me3_peak_region$V4 %in% H3K27me3$Geneid[which(H3K27me3$Significant=="Up")])] <- "H3K27me3 increase"
  H3K9me3_peak_region$condition[which(H3K9me3_peak_region$V4 %in% H3K27me3$Geneid[which(H3K27me3$Significant=="Down")])] <- "H3K27me3 decrease"

  H3K9me3_peak_region <- H3K9me3_peak_region[-which(H3K9me3_peak_region$condition == "H3K27me3 stable"),]  
  t_to_plot <- as.data.frame(table(H3K9me3_peak_region$condition))
  t_to_plot$percent <- t_to_plot$Freq/sum(t_to_plot$Freq) * 100
  t_to_plot$tissue <- tissue_label_change(tissue)
  if(nrow(to_plot)==0){
    to_plot <- t_to_plot
  }else{
    to_plot <- rbind(to_plot,t_to_plot)
  }
}
color <- setNames(c("#009980","#E69900","#838B8B"),c("H3K27me3 decrease","H3K27me3 increase","H3K27me3 stable"))
freq_sum <- to_plot %>%  
  group_by(tissue) %>%  
  summarise(Freq_sum = sum(Freq))  
to_plot_order <- to_plot[which(to_plot$Var1=="H3K27me3 increase"),]
to_plot_order$tissue <-as.character(to_plot_order$tissue)
to_plot_order <- to_plot_order[order(to_plot_order$percent),]
to_plot$tissue <- factor(to_plot$tissue, levels= to_plot_order$tissue)
ggplot(to_plot, aes(x = tissue, y = percent, fill = Var1)) +  
  geom_bar(stat = 'identity',colour = "white") +   
  theme_minimal() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Proportion")+
  geom_text(data = freq_sum, aes(x = tissue, y = -5, label = Freq_sum),   
            inherit.aes = FALSE,  # 避免继承主 ggplot 中的 aes 映射  
            vjust = 1, size = 4, color = "black")  +
  ggtitle("H3K9me3 peak regions")


############### H3K9me3 peak region H3K27me3 condition peak level
tissues <- c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
             "thymus","skin","bladder","bonemarrow","Hip","heart",
             "muscle","jejunum","uterus","ovary","liver","tongue",
             "cecum","colon","testis","stomach","pancreas","iWAT","ileum")
to_plot <- data.frame()
for(tissue in tissues){
  H3K27me3 <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  H3K9me3_peak_region <- read.table(paste0("data/samples/",tissue,"/H3K9me3/bed/H3K9me3_young_merge-W1000-G3000-E100_compress.bed"))
  H3K9me3_peak_region$V2 <- H3K9me3_peak_region$V2+1
  H3K27me3$Start <- H3K27me3$Start + 1
  H3K27me3 <- H3K27me3[,c("Chr","Start","End","Significant")]
  
  H3K27me3 <- as.data.table(H3K27me3)
  H3K9me3_peak_region <- as.data.table(H3K9me3_peak_region)
  setDT(H3K27me3)
  setDT(H3K9me3_peak_region)
  setkey(H3K27me3,Chr,Start,End)
  setkey(H3K9me3_peak_region,V1,V2,V3)
  overlap <- foverlaps(H3K27me3, H3K9me3_peak_region, type = "any", nomatch = 0L)  
  result <- overlap %>%  
    group_by(V4) %>%  
    summarise(  
      total = n(),  
      up_count = sum(Significant == "Up"),  
      down_count = sum(Significant == "Down"),  
      up_ratio = up_count / total,  
      down_ratio = down_count / total  
    ) %>%  
    mutate(  
      Classification = case_when(  
        up_ratio > 0.3 ~ "H3K27me3 increase",  
        down_ratio > 0.3 ~ "H3K27me3 decrease",  
        TRUE ~ "H3K27me3 stable"  
      )  
    ) %>%  
    select(V4, Classification)  
  # H3K9me3_peak_region <- H3K9me3_peak_region[-which(H3K9me3_peak_region$condition == "H3K27me3 stable"),]  
  t_to_plot <- as.data.frame(table(result$Classification))
  t_to_plot$percent <- t_to_plot$Freq/sum(t_to_plot$Freq) * 100
  t_to_plot$tissue <- tissue_label_change(tissue)
  if(nrow(to_plot)==0){
    to_plot <- t_to_plot
  }else{
    to_plot <- rbind(to_plot,t_to_plot)
  }
}
color <- setNames(c("#009980","#E69900","#838B8B"),c("H3K27me3 decrease","H3K27me3 increase","H3K27me3 stable"))
freq_sum <- to_plot %>%  
  group_by(tissue) %>%  
  summarise(Freq_sum = sum(Freq))  
to_plot_order <- to_plot[which(to_plot$Var1=="H3K27me3 increase"),]
to_plot_order$tissue <-as.character(to_plot_order$tissue)
to_plot_order <- to_plot_order[order(to_plot_order$percent),]
to_plot$tissue <- factor(to_plot$tissue, levels= to_plot_order$tissue)
ggplot(to_plot, aes(x = tissue, y = percent, fill = Var1)) +  
  geom_bar(stat = 'identity',colour = "white") +   
  theme_minimal() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Proportion")+
  geom_text(data = freq_sum, aes(x = tissue, y = -5, label = Freq_sum),   
            inherit.aes = FALSE,  # 避免继承主 ggplot 中的 aes 映射  
            vjust = 1, size = 4, color = "black")  +
  ggtitle("H3K9me3 peak regions")
