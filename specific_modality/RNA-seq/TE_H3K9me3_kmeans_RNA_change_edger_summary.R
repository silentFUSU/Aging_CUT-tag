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
library(edgeR)
library(ggbreak) 
options(bitmapType="cairo")  
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
tissue_summary <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_TE_local_H3K9me3_kmeans_change_filter_bar.csv"))
  df <- df %>%
    mutate(condition = sub(".*_", "", X))
  df <- df[,c("condition","Significant")]
  count_df <- df %>%
    group_by(condition, Significant) %>%
    summarise(count = n()) %>%
    arrange(desc(count))
  count_df <- count_df[which(count_df$Significant!="Stable"),]
  count_df$tissue <- tissue_label_change(tissue)
  tissue_summary <- rbind(tissue_summary,count_df)  
}
tissue_order <- c("Lung","Cerebellum","BAT","Muscle","Heart","Aorta","Skin","Kidney","Hippocampus","Cortex","Liver","Tongue","Uterus","Testis","Bladder","Ovary",
                  "Colon","Stomach","Thymus","Cecum","Jejunum","Pancreas","Bone Marrow","Ileum","Spleen","iWAT","Mammary Gland")
conditions <- c("kmeans1", "kmeans2", "kmeans3", "kmeans4", "Stable","other")
significance <- c("Up", "Down")
all_combinations <- expand_grid(
  tissue =tissue_order,
  condition = conditions,
  Significant = significance
)
tissue_summary<- all_combinations %>%
  left_join(tissue_summary, by = c("tissue", "condition", "Significant"))

color <- setNames(c("#e64b35","#3c5488"),c("Up","Down"))
for(kmeans in conditions){
  to_plot <- tissue_summary[which(tissue_summary$condition==kmeans),]
  to_plot$tissue <- factor(to_plot$tissue,levels=tissue_order)
  to_plot <- data.frame(to_plot)
  to_plot$count[which(to_plot$Significant=="Down")] <- -to_plot$count[which(to_plot$Significant=="Down")]
  p <- ggplot(to_plot, aes(x = count, y = tissue, fill = Significant)) +
    geom_bar(stat = "identity") +  
    labs(x = "Count" , y = "Tissue") +  
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12,face = "bold", color = "black"),  
      axis.text.y = element_text(size = 12,face = "bold", color = "black"),  
      axis.title.x = element_text(size = 14,face = "bold", color = "black"), 
      axis.title.y = element_text(size = 14,face = "bold", color = "black"), 
      legend.text = element_text(size = 12),  
      panel.background = element_blank(),  
      panel.grid.major = element_line(size = 0.1, linetype = 'solid', color = "grey"),
      panel.grid.minor = element_line(size = 0.1, linetype = 'solid', color = "lightgrey"),
      panel.border = element_rect(color = "black", fill = NA, size = 1) 
    ) + 
    ggtitle(kmeans)+
    scale_fill_manual(values = color) +
    geom_vline(xintercept = 0, color = "white") +
    scale_x_continuous(limits = c(-300, 300),
                       breaks = seq(-300, 300, by = 50), 
                       labels = function(x) format(abs(x), scientific = FALSE))   
  ggsave(paste0("result/Sup_figures/",kmeans,"_TE_change_count.pdf"),p,width = 6,height = 8)
  }
gtf <- import("~/ref_data/TE_reference/mm10_rmsk_TE.gtf", format = "gtf")
TE_info <- as.data.frame(gtf[,c("gene_id","class_id")])
TE_info <- as.data.frame(TE_info[,c("gene_id","class_id")])
TE_info <- TE_info[!duplicated(TE_info),]

for(kmeans in conditions){
  to_plot_summary <- data.frame()
  for(condition in c("Up","Down")){
    tissue_summary <- data.frame()
    for(tissue in tissues){
      df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_TE_local_H3K9me3_kmeans_change_filter_bar.csv"))
      df <- df %>%
        mutate(
          gene_id = str_extract(X, ".*(?=_.+$)"), 
          condition = str_extract(X, "(?<=_)[^_]+$")
        )
      df <- df[,c("gene_id","Significant","condition")]
      df <- df[which(df$Significant !="Stable" & df$condition==kmeans & df$Significant==condition),]
      if(nrow(df) > 0){
        df <- merge(df,TE_info,by="gene_id")
        t_tissue_summary <- as.data.frame(table(df$class_id))
        t_tissue_summary$tissue <- tissue_label_change(tissue)
        tissue_summary <- rbind(tissue_summary,t_tissue_summary)
      }
    }
    tissue_summary <- reshape2::dcast(tissue_summary,tissue~Var1,value.var = "Freq")
    tissue_summary <- tissue_summary[,which(colnames(tissue_summary) %in%c("tissue","DNA","LINE","LTR","RC","RNA","Satellite","SINE"))]
    tissue_summary[is.na(tissue_summary)] <- 0
    to_plot <- reshape2::melt(tissue_summary)
    to_plot <- to_plot %>%
      group_by(variable) %>%             
      summarize(count = sum(value, na.rm = TRUE)) 
    to_plot$conidtion <- condition
    to_plot_summary <- rbind(to_plot_summary,to_plot)
  }
  color <- setNames(c("#f39b7f","#3c5488","#e64b35","#00a087","#4dbbd5","#fdc010","#bf904a"),c("DNA","LINE","LTR","RC","RNA","Satellite","SINE"))
  to_plot_summary$conidtion <- factor(to_plot_summary$conidtion,levels = c("Up","Down"))
  to_plot_summary$variable <- factor(to_plot_summary$variable,levels = c("DNA","LINE","LTR","RC","RNA","Satellite","SINE"))
  p <- ggplot(to_plot_summary, aes(x = conidtion, y = count, fill = variable)) +
    geom_bar(stat = "identity") +  
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12,face = "bold", color = "black"),  
      axis.text.y = element_text(size = 12,face = "bold", color = "black"),  
      axis.title.x = element_text(size = 14,face = "bold", color = "black"), 
      axis.title.y = element_text(size = 14,face = "bold", color = "black"), 
      legend.text = element_text(size = 12),  
      panel.background = element_blank(),  
      panel.grid.major = element_line(size = 0.1, linetype = 'solid', color = "grey"),
      panel.grid.minor = element_line(size = 0.1, linetype = 'solid', color = "lightgrey"),
      panel.border = element_rect(color = "black", fill = NA, size = 1) 
    ) + 
    scale_fill_manual(values=color)+
    ggtitle(kmeans)+
    ylim(0,800)
  ggsave(paste0("result/Sup_figures/",kmeans,"_TE_change_count_class_level.pdf"),p,width = 4,height = 6)   
  
  to_plot_summary$variable2 <- as.character(to_plot_summary$variable)
  to_plot_summary$variable2[which(!to_plot_summary$variable %in% c("LINE","LTR"))] <- "Other"
  to_plot_summary$variable2 <- factor(to_plot_summary$variable2,levels = c("LINE","LTR","Other"))
  p <- ggplot(to_plot_summary, aes(x = conidtion, y = count, fill = variable2)) +
    geom_bar(stat = "identity") +  
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12,face = "bold", color = "black"),  
      axis.text.y = element_text(size = 12,face = "bold", color = "black"),  
      axis.title.x = element_text(size = 14,face = "bold", color = "black"), 
      axis.title.y = element_text(size = 14,face = "bold", color = "black"), 
      legend.text = element_text(size = 12),  
      panel.background = element_blank(),  
      panel.grid.major = element_line(size = 0.1, linetype = 'solid', color = "grey"),
      panel.grid.minor = element_line(size = 0.1, linetype = 'solid', color = "lightgrey"),
      panel.border = element_rect(color = "black", fill = NA, size = 1) 
    ) + 
    scale_fill_manual(values=color)+
    ggtitle(kmeans)+
    ylim(0,800)
  ggsave(paste0("result/Sup_figures/",kmeans,"_TE_change_count_class_level_only_show_LTR_LINE.pdf"),p,width = 4,height = 6)  
  }



