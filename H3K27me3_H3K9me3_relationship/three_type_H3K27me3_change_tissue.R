rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(VennDiagram) 
H3K27me3_tissues_label <- list(type1=c("bladder","kidney","liver","pancreas","stomach"),
                            type2=c("bonemarrow","cecum","colon","ileum","jejunum","iWAT","mammarygland","spleen","thymus"),
                            type3=c("aorta","BAT","CB","brain","Hip","heart","lung","muscle","ovary","skin","testis","tongue","uterus"))
common_change_region <- list()
types <- c("type1","type2","type3")
for(type in types){
  tissues <- H3K27me3_tissues_label[[type]]
  bin <- data.frame()
  for(tissue in tissues){
    df <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_10kb_bins_diff_after_remove_batch_effect.csv"))
    df <- df[which(df$Significant!="Stable"),c("Geneid","Chr","Start","End","Significant")]
    if(nrow(df) >0){
      df$tissue <- tissue
      df <- unique(df)
      bin <- rbind(bin,df)
    }
  }
  bin_file <- read.table(paste0("/storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_10kb_bins.bed"))
  colnames(bin_file)[4]<-"Geneid"
  
  increase <- bin[which(bin$Significant=="Up"),]
  increase_count <- increase %>%   
    count(Geneid)
  increase_tissue <- increase %>%   
    group_by(Geneid) %>%   
    summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
  increase_count <- merge(increase_count,increase_tissue,by="Geneid")
  increase_count <- merge(increase_count,bin_file,by="Geneid")
  colnames(increase_count)[4:6] <- c("chr","start","end")
  
  decrease <- bin[which(bin$Significant=="Down"),]
  decrease_count <- decrease %>%   
    count(Geneid)
  decrease_tissue <- decrease %>%   
    group_by(Geneid) %>%   
    summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
  decrease_count <- merge(decrease_count,decrease_tissue,by="Geneid")
  decrease_count <- merge(decrease_count,bin_file,by="Geneid")
  
  increase_count <- increase_count[which(increase_count$n >= (0.8 * length(tissues))),]
  decrease_count <- decrease_count[which(decrease_count$n >= (0.8 * length(tissues))),]
  common_change_region[[type]] <- list(increase=increase_count,decrease=decrease_count)
}

condition_counts <- data.frame()
for(type in types){
  t_condition_counts <- data.frame(Increase=nrow(common_change_region[[type]][["increase"]]), Decrease=nrow(common_change_region[[type]][["decrease"]]), condition=type)    
  condition_counts <- rbind(condition_counts,t_condition_counts)
}
to_plot <- reshape2::melt(condition_counts)
ggplot(to_plot, aes(x = condition, y = value, fill = variable)) +  
  geom_bar(stat = 'identity') +   
  theme_minimal() +   
  scale_fill_brewer(palette = "Pastel1") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Counts")+
  ggtitle("H3K27me3 common change bins")

to_plot <- condition_counts
to_plot_rowsum <- rowSums(to_plot[,c(1:2)])
to_plot$Increase <- to_plot$Increase/to_plot_rowsum*100
to_plot$Decrease <- to_plot$Decrease/to_plot_rowsum*100
to_plot <- reshape2::melt(to_plot)
ggplot(to_plot, aes(x = condition, y = value, fill = variable)) +  
  geom_bar(stat = 'identity') +   
  theme_minimal() +   
  scale_fill_brewer(palette = "Pastel1") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Counts")+
  ggtitle("H3K27me3 common change bins")

venn.plot <- draw.triple.venn(  
  area1 = nrow(common_change_region[["type1"]][["increase"]]),  
  area2 = nrow(common_change_region[["type2"]][["increase"]]),  
  area3 = nrow(common_change_region[["type3"]][["increase"]]),  
  n12 = length(intersect(common_change_region[["type1"]][["increase"]]$Geneid, common_change_region[["type2"]][["increase"]]$Geneid)),  
  n23 = length(intersect(common_change_region[["type2"]][["increase"]]$Geneid, common_change_region[["type3"]][["increase"]]$Geneid)),  
  n13 = length(intersect(common_change_region[["type1"]][["increase"]]$Geneid, common_change_region[["type3"]][["increase"]]$Geneid)),  
  n123 = length(Reduce(intersect, list(common_change_region[["type1"]][["increase"]]$Geneid, common_change_region[["type2"]][["increase"]]$Geneid, common_change_region[["type3"]][["increase"]]$Geneid))),  
  category = c("Type1", "Type2", "Type3"),  
  fill = c("red", "green", "blue"),  
  lty = "dashed",  
  cex = 2,  
  cat.cex = 2,  
  cat.col = c("red", "green", "blue")
) 
venn.plot <- draw.triple.venn(  
  area1 = nrow(common_change_region[["type1"]][["decrease"]]),  
  area2 = nrow(common_change_region[["type2"]][["decrease"]]),  
  area3 = nrow(common_change_region[["type3"]][["decrease"]]),  
  n12 = length(intersect(common_change_region[["type1"]][["decrease"]]$Geneid, common_change_region[["type2"]][["decrease"]]$Geneid)),  
  n23 = length(intersect(common_change_region[["type2"]][["decrease"]]$Geneid, common_change_region[["type3"]][["decrease"]]$Geneid)),  
  n13 = length(intersect(common_change_region[["type1"]][["decrease"]]$Geneid, common_change_region[["type3"]][["decrease"]]$Geneid)),  
  n123 = length(Reduce(intersect, list(common_change_region[["type1"]][["decrease"]]$Geneid, common_change_region[["type2"]][["decrease"]]$Geneid, common_change_region[["type3"]][["decrease"]]$Geneid))),  
  category = c("Type1", "Type2", "Type3"),  
  fill = c("red", "green", "blue"),  
  lty = "dashed",  
  cex = 2,  
  cat.cex = 2,  
  cat.col = c("red", "green", "blue")
) 

########################unique
common_change_region_unique <- common_change_region
for(i in c(1:length(types))){
  type1=types[i]
  for(j in c(1:length(types))){
    if(j != i){
      type2=types[j]
      print(paste0(type1," ",type2))
      for(condition in c("increase","decrease")){
        print(condition)
        df1 <- common_change_region_unique[[type1]][[condition]]
        print(paste0("df1 ",nrow(df1)))
        df2 <- common_change_region[[type2]][[condition]]
        df1 <- df1[-which(df1$Geneid %in% df2$Geneid),]
        print(paste0("df1 ",nrow(df1)))
        common_change_region_unique[[type1]][[condition]] <- df1
      }
    }
  }
}

condition_counts <- data.frame()
for(type in types){
  t_condition_counts <- data.frame(Increase=nrow(common_change_region_unique[[type]][["increase"]]), Decrease=nrow(common_change_region_unique[[type]][["decrease"]]), condition=type)    
  condition_counts <- rbind(condition_counts,t_condition_counts)
}
to_plot <- reshape2::melt(condition_counts)
ggplot(to_plot, aes(x = condition, y = value, fill = variable)) +  
  geom_bar(stat = 'identity') +   
  theme_minimal() +   
  scale_fill_brewer(palette = "Pastel1") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Counts")+
  ggtitle("H3K27me3 common change bins")

to_plot <- condition_counts
to_plot_rowsum <- rowSums(to_plot[,c(1:2)])
to_plot$Increase <- to_plot$Increase/to_plot_rowsum*100
to_plot$Decrease <- to_plot$Decrease/to_plot_rowsum*100
to_plot <- reshape2::melt(to_plot)
ggplot(to_plot, aes(x = condition, y = value, fill = variable)) +  
  geom_bar(stat = 'identity') +   
  theme_minimal() +   
  scale_fill_brewer(palette = "Pastel1") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Counts")+
  ggtitle("H3K27me3 common change bins")


