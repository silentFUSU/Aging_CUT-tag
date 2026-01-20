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
library(DSS)
library(patchwork)

tissues <- c("liver","lung","kidney","ileum","Hip","mammarygland","skin","bonemarrow",
             "jejunum","colon","ovary","CB","BAT","thymus","testis","stomach","heart",
             "muscle","bladder","aorta","tongue","spleen","pancreas","brain",
             "cecum","uterus","iWAT")

state_num <- 15
### WGBS
WGBS_tissue_summary <- data.frame()
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  search_table <- search_table[which(search_table$tissue == tissue),]
  search_table$age <- factor(search_table$age, c("3M","24M"))
  search_table <- search_table[order(search_table$age),]
  file_dir <- paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/split_1k/")  
  files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_young[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
  file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
  chromHMM_young <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
  chromHMM_young <- chromHMM_young[which(chromHMM_young$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  chromHMM_young$V2 <- chromHMM_young$V2+1
  chromHMM_young <- data.table(chromHMM_young)
  
  file_dir <- paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/split_1k/")  
  files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_old[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
  file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
  chromHMM_old <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
  chromHMM_old <- chromHMM_old[which(chromHMM_old$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  chromHMM_old$V2 <- chromHMM_old$V2+1
  chromHMM_old <- data.table(chromHMM_old)
  df_list <- list()
  depth_threshold <- 5
  sample_summary <- data.frame()
  for(i in c(1:nrow(search_table))){
    df_list[[i]] <- fread(paste0("data/samples/WGBS/",tissue,"/bdg/",search_table$sample_name[i],"_CpG.bdg"),sep = "\t")
    df_list[[i]] <- df_list[[i]][which(df_list[[i]]$V1 %in% paste0("chr",c(c(1:19),"X","Y"))),]
    df_list[[i]] <- df_list[[i]][which(df_list[[i]]$V5 > depth_threshold),]
    df_list[[i]]$percent <- df_list[[i]]$V4/df_list[[i]]$V5*100
    df_list[[i]]$V3 <- df_list[[i]]$V2
    df_list[[i]] <- df_list[[i]][,c("V1","V2","V3","percent")]
    setDT(df_list[[i]])
    setkey(df_list[[i]], V1, V2, V3) 
    if(search_table$age[i] == "3M"){
      chromHMM <- chromHMM_young
    }else{
      chromHMM <- chromHMM_old
    }
    setDT(chromHMM)  
    setkey(chromHMM, V1, V2, V3) 
    overlaps <- foverlaps(df_list[[i]], chromHMM, type = "any", nomatch = 0L)  
    overlaps <- as.data.frame(overlaps)
    average_percent_by_V4 <- overlaps %>%
      group_by(V4) %>%
      summarise(average_percent = mean(percent, na.rm = TRUE))
    colnames(average_percent_by_V4) <- c("State",search_table$sample_name[i])
    if(nrow(sample_summary)==0){
      sample_summary <- average_percent_by_V4
    }else{
      sample_summary <- merge(sample_summary,average_percent_by_V4,by="State")
    }
  }
  sample_summary$mean <- rowMeans(sample_summary[,-1])
  sample_summary <- sample_summary[,c("State","mean")]
  colnames(sample_summary)[2] <- tissue
  if(nrow(WGBS_tissue_summary)==0){
    WGBS_tissue_summary <- sample_summary
  }else{
    WGBS_tissue_summary <- merge(WGBS_tissue_summary,sample_summary,by="State")
  }
}

### ATAC
ATAC_tissue_summary <- data.frame()
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/ATAC_search_table_batch.csv")
  search_table <- search_table[which(search_table$tissue == tissue),]
  search_table$age <- factor(search_table$age, c("3m","24m"))
  search_table <- search_table[order(search_table$age),]
  file_dir <- paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/split_1k/")  
  files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_young[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
  file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
  chromHMM_young <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
  chromHMM_young <- chromHMM_young[which(chromHMM_young$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  chromHMM_young$V2 <- chromHMM_young$V2+1
  chromHMM_young <- data.table(chromHMM_young)
  
  file_dir <- paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/split_1k/")  
  files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_old[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
  file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
  chromHMM_old <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
  chromHMM_old <- chromHMM_old[which(chromHMM_old$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  chromHMM_old$V2 <- chromHMM_old$V2+1
  chromHMM_old <- data.table(chromHMM_old)
  sample_summary <- data.frame()
  tab <- read.table(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_1kb_bins.counts"),header = T)
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+|DYQ[0-9]+).*"
  colnames(tab)[7:ncol(tab)] <-  gsub(pattern, "\\1",colnames(tab)[7:ncol(tab)])
  rownames(tab) <- tab$Geneid
  counts <- tab[,7:ncol(tab)]
  CPM <- as.data.frame(edgeR::cpm(counts))
  CPM <- merge(tab[,c(2:4)],CPM,by="row.names")
  CPM$Start <- CPM$Start+1
  for(i in c(1:nrow(search_table))){
    df <- CPM[,c("Chr","Start","End",search_table$sample_name[i])]
    colnames(df)[4] <- "CPM"
    setDT(df)
    setkey(df, Chr, Start, End) 
    if(search_table$age[i] == "3m"){
      chromHMM <- chromHMM_young
    }else{
      chromHMM <- chromHMM_old
    }
    setDT(chromHMM)  
    setkey(chromHMM, V1, V2, V3) 
    overlaps <- foverlaps(df, chromHMM, type = "any", nomatch = 0L)  
    overlaps <- as.data.frame(overlaps)
    average_CPM_by_V4 <- overlaps %>%
      group_by(V4) %>%
      summarise(average_CPM = mean(CPM, na.rm = TRUE))
    colnames(average_CPM_by_V4) <- c("State",search_table$sample_name[i])
    if(nrow(sample_summary)==0){
      sample_summary <- average_CPM_by_V4
    }else{
      sample_summary <- merge(sample_summary,average_CPM_by_V4,by="State")
    }
  }
  sample_summary$mean <- rowMeans(sample_summary[,-1])
  sample_summary <- sample_summary[,c("State","mean")]
  colnames(sample_summary)[2] <- tissue
  if(nrow(ATAC_tissue_summary)==0){
    ATAC_tissue_summary <- sample_summary
  }else{
    ATAC_tissue_summary <- merge(ATAC_tissue_summary,sample_summary,by="State")
  }
}
# write.csv(WGBS_tissue_summary,"data/samples/WGBS/all/DNA_methylation_in_15_chromHMM_state.csv")
# write.csv(ATAC_tissue_summary,"data/samples/WGBS/all/ATAC_in_15_chromHMM_state.csv")
WGBS_tissue_summary <- read.csv("data/samples/WGBS/all/DNA_methylation_in_15_chromHMM_state.csv",row.names = 1)
ATAC_tissue_summary <- read.csv("data/samples/WGBS/all/ATAC_in_15_chromHMM_state.csv",row.names = 1)

WGBS_tissue_summary_mean <- WGBS_tissue_summary
WGBS_tissue_summary_mean$mean <- rowMeans(WGBS_tissue_summary_mean[,-1])
WGBS_tissue_summary_mean <- WGBS_tissue_summary_mean[,c("State","mean")]
# WGBS_tissue_summary_mean$delta <- WGBS_tissue_summary_mean$mean - mean(WGBS_tissue_summary_mean$mean)

ATAC_tissue_summary_mean <- ATAC_tissue_summary
ATAC_tissue_summary_mean$mean <- rowMeans(ATAC_tissue_summary_mean[,-1])
ATAC_tissue_summary_mean <- ATAC_tissue_summary_mean[,c("State","mean")]

colnames(WGBS_tissue_summary_mean)[2] <- "DNA methylation"
colnames(ATAC_tissue_summary_mean)[2] <- "Chromosome accessibility"


to_plot <- merge(WGBS_tissue_summary_mean,ATAC_tissue_summary_mean,by="State")
dictionary <- list("E1"=1, "E2"=2, "E3"=3,
                   "E4"=4, "E5"=5, "E6"=7,
                   "E7"=8, "E8"=6, "E9"=9,
                   "E10"=10,"E11"=11,"E12"=15,
                   "E13"=12,"E14"=13,"E15"=14)
keys <- names(dictionary)
values <- unlist(dictionary)
to_plot$State <- values[match(to_plot$State, keys)]
to_plot$State <- paste0("E",to_plot$State)
rownames(to_plot) <- to_plot$State
to_plot$State <- factor(to_plot$State,levels = paste0("E",1:state_num))
to_plot <- to_plot[order(to_plot$State),]
color_palette <- colorRampPalette(c("#a6d0e4","white","#f76b8a"))(100) 
# pheatmap::pheatmap(to_plot[,-1],scale = "column",cluster_cols = F,cluster_rows = F,color = color_palette,filename = "result/figures/chromHMM_15_WGBS_ATAC.pdf",width = 3,height = 6)
# pheatmap::pheatmap(to_plot[,-1],scale = "column",cluster_cols = F,cluster_rows = F,color = color_palette)
breaks <- c(seq(40, 90, length.out = 50))
color_palette <- colorRampPalette(c("white", "#defcf9","#4589C8FF"))(50) 
pheatmap::pheatmap(to_plot[,2,drop=F],cluster_cols = F,cluster_rows = F,breaks = breaks,color = color_palette,filename = "result/figures/chromHMM_15_WGBS.pdf",width = 3,height = 6)


breaks <- c(seq(0.3, 1, length.out = 50))
color_palette <- colorRampPalette(c("white", "#defcf9","#4589C8FF"))(50) 
pheatmap::pheatmap(to_plot[,3,drop=F],cluster_cols = F,cluster_rows = F,breaks = breaks,color = color_palette,filename = "result/figures/chromHMM_15_ATAC.pdf",width = 3,height = 6)


WGBS_tissue_summary_long <- reshape2::melt(WGBS_tissue_summary)
dictionary <- list("E1"=1, "E2"=2, "E3"=3,
                   "E4"=4, "E5"=5, "E6"=7,
                   "E7"=8, "E8"=6, "E9"=9,
                   "E10"=10,"E11"=11,"E12"=15,
                   "E13"=12,"E14"=13,"E15"=14)
keys <- names(dictionary)
values <- unlist(dictionary)
WGBS_tissue_summary_long$State <- values[match(WGBS_tissue_summary_long$State, keys)]
WGBS_tissue_summary_long$State <- paste0("E",WGBS_tissue_summary_long$State)
WGBS_tissue_summary_long$State <- factor(WGBS_tissue_summary_long$State, levels = rev(paste0("E",1:15)))

plot1 <- ggplot(WGBS_tissue_summary_long, aes(x = value, y = State)) +
  geom_violin(fill = "#4589C8FF", color = "black") +
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
  xlim(20,100)+
  labs(x = "DNA Methylation", y= NULL)
ggplot(WGBS_tissue_summary_long[which(WGBS_tissue_summary_long$State %in% c("E10","E11")),], aes(x = State, y = value)) +
  geom_violin(fill = "#4589C8FF", color = "black") +
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
  ylim(0,100)+
  labs(x = "DNA Methylation", y= NULL)
t.test(WGBS_tissue_summary_long$value[which(WGBS_tissue_summary_long$State =="E10")],WGBS_tissue_summary_long$value[which(WGBS_tissue_summary_long$State =="E11")])

ATAC_tissue_summary_long <- reshape2::melt(ATAC_tissue_summary)
ATAC_tissue_summary_long$State <- factor(ATAC_tissue_summary_long$State, levels = rev(paste0("E",1:15)))

plot2 <- ggplot(ATAC_tissue_summary_long,  aes(x = value, y = State)) +
  geom_violin(fill = "#EF7C7AFF", color = "black") +
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
  labs(x = "Chromosome Accessibility", y = NULL)
p <-grid.arrange(plot1, plot2, ncol = 2)
ggsave("result/figures/chromHMM_15_WGBS_ATAC.pdf",p,width = 4,height = 6)
