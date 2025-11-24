rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(bitmapType="cairo")  
library(grid)
library(stringr)
library(dplyr)
library(ggplot2)
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
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
                  "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT","ileum")) 
tissue_summary <- data.frame()
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  t_search_table <- search_table[which(search_table$tissue==tissue),]
  H3K27me3_peaks <- read.table(paste0("data/samples/all/H3K27me3/bed/H3K27me3_young_merge-W5000-G10000-E100_peak_merged_exist_in_larger_27_tissues.bed"))
  if(tissue %in% c("ovary","mammarygland","uterus")){
    H3K27me3_peaks <- H3K27me3_peaks[which(H3K27me3_peaks$V1 %in% paste0("chr",c(1:19,"X"))),]
  }
  setDT(H3K27me3_peaks)
  setkey(H3K27me3_peaks,V1,V2,V3)
  summary <- data.frame()
  for(sample in t_search_table$sample_name){
    df <- fread(paste0("data/samples/WGBS/",tissue,"/bdg/",sample,"_CpG.bdg"),sep = "\t")
    setDT(df)
    setkey(df,V1,V2,V3)  
    overlaps <- foverlaps(df, H3K27me3_peaks, type = "any", nomatch = 0L)  
    result <- data.frame(tissue=tissue_label_change(tissue),sample=sample,methylation=sum(overlaps$V4)/sum(overlaps$V5)*100)
    summary <- rbind(summary,result)

  }
  young_summary <- data.frame(tissue=tissue_label_change(tissue),young=mean(summary$methylation[which(summary$sample %in% t_search_table$sample_name[which(t_search_table$age=="3M")])]))
  old_summary <-  data.frame(tissue=tissue_label_change(tissue),old=mean(summary$methylation[which(summary$sample %in% t_search_table$sample_name[which(t_search_table$age=="24M")])]))
  t_tissue_summary <- merge(young_summary,old_summary,by="tissue")
  t_tissue_summary$delta <- t_tissue_summary$old - t_tissue_summary$young
  t_tissue_summary <- t_tissue_summary[,c("tissue","delta")]

  tissue_summary <- rbind(tissue_summary,t_tissue_summary)

}

tissue_summary <- readRDS("tmp_WGBS_in_H3K27me3.rds")
tissue_order <- c("Mammary Gland","Cecum","Thymus","Uterus","iWAT","Stomach","Skin","Spleen",
                  "Muscle","Bone Marrow","Liver","Ileum","Testis","Cortex","Jejunum","Tongue",
                  "Hippocampus","Colon","Bladder","Aorta","Cerebellum","Lung","Heart","Kidney","BAT","Ovary","Pancreas")

to_plot <- tissue_summary
to_plot$condition <- "Up"
to_plot$condition[which(to_plot$delta < 0 )] <- "Down"
to_plot$tissue <- factor(to_plot$tissue,levels=tissue_order)
color <- setNames(c("#f39b7f","#4dbbd5"),c("Up","Down"))

ggplot(to_plot, aes(x = tissue, y = delta, fill = condition)) +  
  geom_bar(stat = 'identity') +   
  theme_bw() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Delta")+
  ylim(-6,6) +
  ggtitle(NULL) +
  guides(fill = FALSE) +
  geom_hline(yintercept = c(-1, 1), color = "black", linetype = "dashed")
2