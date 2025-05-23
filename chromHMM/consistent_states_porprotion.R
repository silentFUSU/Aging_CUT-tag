rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(ggalluvial)
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver","ileum",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT"))
state_num <- "11"
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}
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
p_list <- list()
for(tissue in tissues){
  young1 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_young1_",state_num,"_segments_1k.bed"),header = F)
  young2 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_young2_",state_num,"_segments_1k.bed"),header = F)
  young1$label <- paste(young1$V1,young1$V2,young1$V3,young1$V4,sep = "-")
  young2$label <- paste(young2$V1,young2$V2,young2$V3,young2$V4,sep = "-")
  young <- young1[which(young1$label %in% intersect(young1$label,young2$label)),]
  
  old1 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_old1_",state_num,"_segments_1k.bed"),header = F)
  old2 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_old2_",state_num,"_segments_1k.bed"),header = F)
  old1$label <- paste(old1$V1,old1$V2,old1$V3,old1$V4,sep = "-")
  old2$label <- paste(old2$V1,old2$V2,old2$V3,old2$V4,sep = "-")
  old <- old1[which(old1$label %in% intersect(old1$label,old2$label)),]
  
  to_plot <- data.frame(age=c("young","young","old","old"),
                        condition=c("Consistent","Inconsistent","Consistent","Inconsistent"),
                        proportion=c(nrow(young)/nrow(young1),(nrow(young1)-nrow(young))/nrow(young1),nrow(old)/nrow(old1),(nrow(old1)-nrow(old))/nrow(old1)))
  to_plot$proportion <- to_plot$proportion*100
  color <- setNames(c("#009980","#838B8B"),c("Consistent","Inconsistent"))
  to_plot$condition <- factor(to_plot$condition,levels=c("Inconsistent","Consistent"))
  to_plot$position <- 100
  to_plot$position[which(to_plot$condition=="Consistent")] <- to_plot$proportion[which(to_plot$condition=="Consistent")]
  to_plot$label <- paste0(round(to_plot$proportion,2),"%")
  to_plot$age <- factor(to_plot$age,levels=c("young","old"))
  p_list[[tissue]] <- ggplot(to_plot, aes(x = age, y = proportion, fill = condition)) +  
    geom_bar(stat = 'identity',colour = "white") +   
    theme_minimal() +   
    scale_fill_manual(values = color) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 15),legend.title = element_blank()) +
    ylab("Proportion")+
    geom_text(data = to_plot,   
              aes(label = label, y = position),   
              color = "black", size = 5, vjust = 0.5)+
    ggtitle(tissue_label_change(tissue),"ChromHMM state consistent proportion")
}

combined_plot <- plot_a_list(p_list, 4, 7)
ggsave("result/all/ChromHMM/all_tissues/11_all_tissues/consistent_state_proportion.png",combined_plot,width = 30,height = 20,type="cairo")

p_list <- list()
for(tissue in tissues){
  young1 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_young1_",state_num,"_segments_1k.bed"),header = F)
  young2 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_young2_",state_num,"_segments_1k.bed"),header = F)
  young1$label <- paste(young1$V1,young1$V2,young1$V3,young1$V4,sep = "-")
  young2$label <- paste(young2$V1,young2$V2,young2$V3,young2$V4,sep = "-")
  young <- young1[which(young1$label %in% intersect(young1$label,young2$label)),]
  
  old1 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_old1_",state_num,"_segments_1k.bed"),header = F)
  old2 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_old2_",state_num,"_segments_1k.bed"),header = F)
  old1$label <- paste(old1$V1,old1$V2,old1$V3,old1$V4,sep = "-")
  old2$label <- paste(old2$V1,old2$V2,old2$V3,old2$V4,sep = "-")
  old <- old1[which(old1$label %in% intersect(old1$label,old2$label)),]
  
  young_summary <- as.data.frame(table(young$V4))
  young_summary$age <- "young"
  young_summary$percent <- young_summary$Freq/sum(young_summary$Freq)*100
  
  old_summary <- as.data.frame(table(old$V4))
  old_summary$age <- "old"
  old_summary$percent <- old_summary$Freq/sum(old_summary$Freq)*100
  
  to_plot <- rbind(young_summary,old_summary)
  to_plot$age <- factor(to_plot$age, levels = c("young","old"))
  color <- read.csv("data/samples/20_distinct_color.txt",header = F)
  color <- setNames(color$V1,paste0("E",1:state_num))
  to_plot$Var1 <- factor(to_plot$Var1,levels=paste0("E",1:state_num))
  p_list[[tissue]] <- ggplot(to_plot, aes(x = age, y = percent, fill = Var1)) +  
    geom_bar(stat = 'identity',colour = "white") +   
    theme_minimal() +   
    scale_fill_manual(values = color) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 15),legend.title = element_blank()) +
    ylab("Proportion")+
    ggtitle(tissue_label_change(tissue),"ChromHMM consistent state proportion")
}
combined_plot <- plot_a_list(p_list, 4, 7)
ggsave("result/all/ChromHMM/all_tissues/11_all_tissues/consistent_state_proportion_of_each_state.png",combined_plot,width = 35,height = 20,type="cairo")

####inconsistent state composition
p_list <- list()
for(tissue in tissues){
  young1 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_young1_",state_num,"_segments_1k.bed"),header = F)
  young2 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_young2_",state_num,"_segments_1k.bed"),header = F)
  young1$label <- paste(young1$V1,young1$V2,young1$V3,young1$V4,sep = "-")
  young2$label <- paste(young2$V1,young2$V2,young2$V3,young2$V4,sep = "-")
  young <- young1[which(young1$label %in% intersect(young1$label,young2$label)),]
  
  old1 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_old1_",state_num,"_segments_1k.bed"),header = F)
  old2 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_old2_",state_num,"_segments_1k.bed"),header = F)
  old1$label <- paste(old1$V1,old1$V2,old1$V3,old1$V4,sep = "-")
  old2$label <- paste(old2$V1,old2$V2,old2$V3,old2$V4,sep = "-")
  old <- old1[which(old1$label %in% intersect(old1$label,old2$label)),]
  
  young1_inconsistent <- young1[-which(young1$label %in% young$label),]
  young2_inconsistent <- young2[-which(young2$label %in% young$label),]
  young1_inconsistent$label <- paste(young1_inconsistent$V1,young1_inconsistent$V2,young1_inconsistent$V3,sep = "-")
  young2_inconsistent$label <- paste(young2_inconsistent$V1,young2_inconsistent$V2,young2_inconsistent$V3,sep = "-")  
  colnames(young1_inconsistent)[4] <- "young1_state"
  colnames(young2_inconsistent)[4] <- "young2_state"
  young_inconsistent <- merge(young1_inconsistent[,c("label","young1_state")],young2_inconsistent[,c("label","young2_state")],by="label")
  young_inconsistent$young1_state <- factor(young_inconsistent$young1_state, levels=c(paste0("E",1:state_num)))
  young_inconsistent$young2_state <- factor(young_inconsistent$young2_state, levels=c(paste0("E",1:state_num)))
  
  old1_inconsistent <- old1[-which(old1$label %in% old$label),]
  old2_inconsistent <- old2[-which(old2$label %in% old$label),]
  old1_inconsistent$label <- paste(old1_inconsistent$V1,old1_inconsistent$V2,old1_inconsistent$V3,sep = "-")
  old2_inconsistent$label <- paste(old2_inconsistent$V1,old2_inconsistent$V2,old2_inconsistent$V3,sep = "-")
  old_inconsistent$old1_state <- factor(old_inconsistent$old1_state, levels=c(paste0("E",1:state_num)))
  old_inconsistent$old2_state <- factor(old_inconsistent$old2_state, levels=c(paste0("E",1:state_num)))
  colnames(old1_inconsistent)[4] <- "old1_state"
  colnames(old2_inconsistent)[4] <- "old2_state"
  old_inconsistent <- merge(old1_inconsistent[,c("label","old1_state")],old2_inconsistent[,c("label","old2_state")],by="label")
  
  to_plot <- young_inconsistent %>%
    group_by(young1_state, young2_state) %>%
    summarize(freq = n())
  to_plot <- to_plot %>%
    group_by(young1_state) %>%
    mutate(percent_freq = freq/sum(freq) *100)
  to_plot <- reshape2::dcast(to_plot, formula = young1_state~young2_state, value.var = "percent_freq") 
  labels <- paste0("E",1:state_num)
  rownames(to_plot) <- to_plot$young1_state
  to_plot <- to_plot[,-1]
  to_plot[is.na(to_plot)] <- 0
  for(i in 1:nrow(to_plot)){
    to_plot[i,i] <- NA
  }
  pheatmap::pheatmap(to_plot,
                     cluster_rows = F,cluster_cols = F, 
                     breaks = seq(0, 100, length.out = 101),
                     display_numbers = T,labels_row = labels,
                     labels_col = labels,fontsize = 10,na_col = "white",
                     main = paste0(tissue_label_change(tissue)," Young inconsistent states"))
  
  to_plot <- old_inconsistent %>%
    group_by(old1_state, old2_state) %>%
    summarize(freq = n())
  to_plot <- to_plot %>%
    group_by(old1_state) %>%
    mutate(percent_freq = freq/sum(freq) *100)
  to_plot <- reshape2::dcast(to_plot, formula = old1_state~old2_state, value.var = "percent_freq") 
  labels <- paste0("E",1:state_num)
  rownames(to_plot) <- to_plot$old1_state
  to_plot <- to_plot[,-1]
  to_plot[is.na(to_plot)] <- 0
  for(i in 1:nrow(to_plot)){
    to_plot[i,i] <- NA
  }
  pheatmap::pheatmap(to_plot,
                     cluster_rows = F,cluster_cols = F, 
                     breaks = seq(0, 100, length.out = 101),
                     display_numbers = T,labels_row = labels,
                     labels_col = labels,fontsize = 10,na_col = "white",
                     main = paste0(tissue_label_change(tissue)," Old inconsistent states"))
  
  # to_plot <- young_inconsistent %>%  
  #   group_by(young1_state, young2_state) %>%  
  #   summarize(freq = n())  
  # to_plot$young1_state <- factor(to_plot$young1_state,levels=paste0("E",1:state_num))
  # to_plot$young2_state <- factor(to_plot$young2_state,levels=paste0("E",1:state_num))
  # p_list[[paste0(tissue,"_young")]] <- ggplot(to_plot, aes(axis1 = young1_state, axis2 = young2_state, y = freq)) +  
  #   geom_alluvium(aes(fill = young2_state)) +  
  #   geom_stratum() +  
  #   geom_text(stat = "stratum", aes(label = after_stat(stratum))) +  
  #   theme_minimal() +  
  #   labs(y = "Count", x = "State Transition", 
  #        fill = "Young State")+
  #   ggtitle(paste0(tissue_label_change(tissue)," young inconsistent regions"))
  
  # to_plot <- old_inconsistent %>%  
  #   group_by(old1_state, old2_state) %>%  
  #   summarize(freq = n())  
  # to_plot$old1_state <- factor(to_plot$old1_state,levels=paste0("E",1:state_num))
  # to_plot$old2_state <- factor(to_plot$old2_state,levels=paste0("E",1:state_num))
  # p_list[[paste0(tissue,"_old")]] <- ggplot(to_plot, aes(axis1 = old1_state, axis2 = old2_state, y = freq)) +  
  #   geom_alluvium(aes(fill = old2_state)) +  
  #   geom_stratum() +  
  #   geom_text(stat = "stratum", aes(label = after_stat(stratum))) +  
  #   theme_minimal() +  
  #   labs(y = "Count", x = "State Transition", 
  #        fill = "old State")+
  #   ggtitle(paste0(tissue_label_change(tissue)," old inconsistent regions"))
}

combined_plot <- plot_a_list(p_list,no_of_cols = 4,no_of_rows = length(p_list)/4)
ggsave("result/all/ChromHMM/all_tissues/11_all_tissues/inconsistent_state_proportion_of_each_state.png",combined_plot,width = 20,height = length(p_list)/4*7,limitsize = FALSE,type="cairo")
