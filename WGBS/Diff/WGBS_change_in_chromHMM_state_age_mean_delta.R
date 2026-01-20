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
library(dplyr)
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
state_num <- 15
for(tissue in tissues){
  file_dir <- paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/split_1k/")  
  files_to_read <- list.files(path = file_dir, pattern = paste0(tissue,"_young[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
  file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  t_search_table <- search_table[which(search_table$tissue==tissue),]
  chromHMM_young <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
  chromHMM_young <- chromHMM_young[which(chromHMM_young$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  chromHMM_young$V2 <- chromHMM_young$V2+1
  bin <- chromHMM_young
  if(tissue %in% c("ovary","uterus","mammarygland")){
    bin <- bin[which(bin$V1%in% paste0("chr",c(1:19,"X"))),]
    bin$V1 <- factor(bin$V1,levels=paste0("chr",c(1:19,"X")))
  }else{
    bin$V1 <- factor(bin$V1,levels=paste0("chr",c(1:19,"X","Y")))
  }
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
    result$methylation <- result$V4_sum/result$V5_sum*100
    result <- result[,c("V4","methylation")]
    colnames(result) <- c("label",sample)
    if(nrow(summary)==0){
      summary <- result
    }else{
      summary <- merge(summary,result,by="label")
    }
  }
  young_summary <- summary[,c("label",t_search_table$sample_name[which(t_search_table$age=="3M")])]
  old_summary <- summary[,c("label",t_search_table$sample_name[which(t_search_table$age=="24M")])]
  young_summary$young_methylation <- rowMeans(young_summary[,-1])
  old_summary$old_methylation <- rowMeans(old_summary[,-1])
  
  t_tissue_summary <- merge(young_summary,old_summary,by="label")
  t_tissue_summary$delta <- t_tissue_summary$old_methylation - t_tissue_summary$young_methylation
  t_tissue_summary <- t_tissue_summary[,c("label","delta")]
  colnames(t_tissue_summary)[2] <- tissue_label_change(tissue)
  if(nrow(tissue_summary)==0){
    tissue_summary <- t_tissue_summary
  }else{
    tissue_summary <- merge(tissue_summary,t_tissue_summary,by="label",all=T)
  }
}
# write.csv(tissue_summary,"data/samples/WGBS/all_tissues_age_mean_delta_chromHMM_state.csv")
tissue_summary <- read.csv("data/samples/WGBS/all_tissues_age_mean_delta_chromHMM_state.csv",row.names = 1)
dictionary <- list("E1"=1, "E2"=2, "E3"=3,
                   "E4"=4, "E5"=5, "E6"=7,
                   "E7"=8, "E8"=6, "E9"=9,
                   "E10"=10,"E11"=11,"E12"=15,
                   "E13"=12,"E14"=13,"E15"=14)
keys <- names(dictionary)
values <- unlist(dictionary)
tissue_summary$label <- values[match(tissue_summary$label, keys)]
tissue_summary$label <- paste0("state",tissue_summary$label)
tissue_order <- c("Mammary.Gland","Cecum","Thymus","Uterus","iWAT","Stomach","Skin","Spleen",
                  "Muscle","Bone.Marrow","Liver","Ileum","Testis","Cortex","Jejunum","Tongue",
                  "Hippocampus","Colon","Bladder","Aorta","Cerebellum","Lung","Heart","Kidney","BAT","Ovary","Pancreas")
to_plot <- tissue_summary
to_plot <- to_plot[which(to_plot$label=="state11"),]
to_plot <- as.data.frame(t(to_plot))
to_plot$tissue <- rownames(to_plot)
colnames(to_plot)[1] <- "delta"
to_plot <- to_plot[-1,]
to_plot$delta <- as.numeric(to_plot$delta)
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
  ylim(-8,8) +
  ggtitle(NULL) +
  guides(fill = FALSE) +
  geom_hline(yintercept = c(-1, 1), color = "black", linetype = "dashed")

to_plot <- tissue_summary
rownames(to_plot) <- to_plot$label
to_plot <- to_plot[,-1]
to_plot <- to_plot[,tissue_order]
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- c(seq(-5, -2.1, length.out = 40), seq(-2, 2, length.out = 20), seq(2.1, 5, length.out = 40))
pheatmap::pheatmap(to_plot,cluster_rows = T,cluster_cols = F,show_rownames = T,breaks = breaks, color = color_palette,annotation_colors = annotation_color,main = "Whole genome chromHMM state CpG methylation Delta(Old - Young)",filename = "result/Sup_figures/WGBS_change_in_chromHMM_State.pdf",width = 8,height = 6)
