rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
library(limma)
library(readr) 
library(ggrepel)
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
search_table_CUTTAG <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
search_table_ATAC <- read.csv("data/samples/all/ATAC_search_table_batch.csv")
search_table_RNA <- read.csv("data/samples/all/RNA_search_table.csv")
search_table <- rbind(search_table_CUTTAG[,c("tissue","antibody","sample_name","age")],search_table_ATAC[,c("tissue","antibody","sample_name","age")])
search_table <- rbind(search_table,search_table_RNA[,c("tissue","antibody","sample_name","age")])
search_table$tissue <- sapply(search_table$tissue,tissue_label_change)

qc_df <- read.csv("data/samples/all/CUTTAG_ATAC_RNA_depth_map.csv")

qc_df <- merge(search_table,qc_df,by.x="sample_name",by.y="sample")

# qc_df$depth <- qc_df$depth*2*150/1000000000
qc_df$depth <- qc_df$depth/1000000
to_plot <- qc_df
to_plot$antibody <- factor(to_plot$antibody,levels=c("H3K9me3","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K27ac","ATAC","RNA"))
color <- read.table("data/samples/20_distinct_color.txt")
color <- setNames(color$V1,c("H3K9me3","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K27ac","ATAC","RNA"))
p <- ggplot(to_plot, aes(x = antibody, y = depth)) +
  geom_violin(fill="gray") +
  geom_boxplot(fill="white",width = 0.2,) +
  # scale_fill_manual(values = color) +
  labs(y = "#total reads (million)") +
  theme_minimal() +
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
  ylim(0,125)
ggsave("result/Sup_figures/histone_ATAC_RNA_depth.pdf",p,height=6,width = 10)

to_plot$map <- to_plot$map *100
p <- ggplot(to_plot, aes(x = antibody, y = map)) +
  geom_violin(fill="gray") +
  geom_boxplot(fill="white",width = 0.2,) +
  # scale_fill_manual(values = color) +
  labs(y = "Mapping%") +
  theme_minimal()+
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
  ylim(0,100)
ggsave("result/Sup_figures/histone_ATAC_RNA_map.pdf",p,height=6,width = 10)

#### WGBS
search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
qc_df <- read.csv("data/samples/all/WGBS_depth_map.csv")
qc_df <- merge(search_table,qc_df,by.x="sample_name",by.y="sample")

# qc_df$depth <- qc_df$depth*2*150/1000000000
qc_df$depth <- qc_df$depth/1000000
to_plot <- qc_df
# to_plot$antibody <- factor(to_plot$antibody,levels=c("H3K9me3","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K27ac","ATAC","RNA"))
# color <- read.table("data/samples/20_distinct_color.txt")
# color <- setNames(color$V1,c("H3K9me3","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K27ac","ATAC","RNA","WGBS"))
p <- ggplot(to_plot, aes(x = antibody, y = depth)) +
  geom_violin(fill="gray") +
  geom_boxplot(fill="white",width = 0.2,) +
  # scale_fill_manual(values = color) +
  labs(y = "# total reads (million)") +
  theme_minimal()+
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
  ylim(0,500)
ggsave("result/Sup_figures/WGBS_depth.pdf",p,height=6,width = 3)

qc_df <- qc_df %>%
  mutate(
    map = parse_number(map),
    Lambda.DNA.coversion = parse_number(Lambda.DNA.coversion) 
  )
to_plot <- qc_df
p <- ggplot(to_plot, aes(x = antibody, y = map)) +
  geom_violin(fill= "gray") +
  geom_boxplot(fill="white",width = 0.2,) +
  # scale_fill_manual(values = color) +
  labs(y = "Mapping%") +
  theme_minimal()+
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
  ylim(0,100)
ggsave("result/Sup_figures/WGBS_mapping.pdf",p,height=6,width = 3)

p <- ggplot(to_plot, aes(x = antibody, y = Lambda.DNA.coversion)) +
  geom_violin(fill="gray") +
  geom_boxplot(fill="white",width = 0.2,) +
  # scale_fill_manual(values = color) +
  labs(y = "Lambda DNA conversion %") +
  theme_minimal()+
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
  ylim(95,100)
ggsave("result/Sup_figures/WGBS_lambada_DNA.pdf",p,height=6,width = 3)

#### HiC
search_table <- read.csv("data/samples/all/HiC_search_table.csv")
qc_df <- read.csv("data/samples/all/HiC_depth_map.csv")
qc_df <- merge(search_table,qc_df,by.x="sample_name",by.y="sample")
qc_df <- qc_df %>%
  mutate(
    map = parse_number(map),
    VP_cis_20K = parse_number(VP_cis_20K) 
  )

# qc_df$depth <- qc_df$depth*2*150/1000000000
qc_df$depth <- qc_df$depth/1000000
to_plot <- qc_df
# color <- read.table("data/samples/20_distinct_color.txt")
# color <- setNames(color$V1,c("H3K9me3","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K27ac","ATAC","RNA","WGBS","HiC"))
p <- ggplot(to_plot, aes(x = antibody, y = depth)) +
  geom_violin(fill="gray") +
  geom_boxplot(fill="white",width = 0.2,) +
  # scale_fill_manual(values = color) +
  labs(y = "# total reads (million)") +
  theme_minimal()+
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
  ylim(0,1500)
ggsave("result/Sup_figures/HiC_depth.pdf",p,height=6,width = 3)


p <- ggplot(to_plot, aes(x = antibody, y = map)) +
  geom_violin(fill="gray") +
  geom_boxplot(fill="white",width = 0.2,) +
  # scale_fill_manual(values = color) +
  labs(y = "Mapping%") +
  theme_minimal()+
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
  ylim(0,100)
ggsave("result/Sup_figures/HiC_mapping.pdf",p,height=6,width = 3)

p <- ggplot(to_plot, aes(x = antibody, y = VP_cis_20K)) +
  geom_violin(fill="gray") +
  geom_boxplot(fill="white",width = 0.2,) +
  # scale_fill_manual(values = color) +
  labs(y = "Long distance cis%") +
  theme_minimal()+
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
  ylim(0,100)
ggsave("result/Sup_figures/HiC_long_distance_cis.pdf",p,height=6,width = 3)
