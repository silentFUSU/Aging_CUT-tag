rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(bitmapType="cairo")  
library(grid)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggrepel)
extract_before_bracket <- function(s) {  
  parts <- strsplit(s, "\\(")[[1]]  
  return(parts[1])  
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
motif_count_summary <- data.frame()
motif_summary <- list(increase=data.frame(),decrease=data.frame())
motif_data_frame <- list(increase=data.frame(),decrease=data.frame())
conditions <- c("increase","decrease")
tissues <-  c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
              "thymus","skin","bladder","bonemarrow","Hip","heart",
              "muscle","jejunum","uterus","ovary","liver","tongue",
              "cecum","colon","testis","stomach","pancreas","iWAT","ileum")
data_path <- "data/samples/WGBS/all/snapatac2/mutual_bg/"

for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(conditions))){
    condition <- conditions[j]
    if(file.exists(paste0(data_path,condition,"/enrichment_results_",tissue,"_DMR_",condition,"_delta01.bed.csv"))){
      motif <- read.csv(paste0(data_path,condition,"/enrichment_results_",tissue,"_DMR_",condition,"_delta01.bed.csv"))
      motif$adjusted.p.value[which(motif$log2.fold.change. <= 0 )] <- 1
      motif$adjusted.p.value[which(motif$fg_percent==0)] <- 1
      motif$log2.fold.change.[which(motif$fg_percent>0 & motif$bg_percent ==0)] <- max(motif$log2.fold.change.[is.finite(motif$log2.fold.change.)],na.rm = T)
      motif$id <- ifelse(
        grepl("\\(.*\\)", motif$id),
        sub(".*\\((.*?)\\).*", "\\1", motif$id), 
        sub("^M\\d+_2\\.00\\s*", "", motif$id) 
      )
      
      colnames(motif)[which(colnames(motif)=="adjusted.p.value")] <- tissue_label_change(tissue)
      if(i == 1){
        motif_data_frame[[condition]] <- motif[,c("id",tissue_label_change(tissue))]
      }else{
        motif_data_frame[[condition]] <- merge(motif_data_frame[[condition]],motif[,c("id",tissue_label_change(tissue))],by="id",all=T)
      }
      
      
      motif <- motif[which(motif$log2.fold.change. > 1 & motif[,9] < 0.01),]
      motif_count <- data.frame(count=nrow(motif),
                                tissue=tissue_label_change(tissue),
                                condition=condition)
      motif_count_summary <- rbind(motif_count_summary,motif_count)
      if(nrow(motif) > 0){
        motif <- motif[order(motif[,9],-motif$log2.fold.change.,-motif$fg_percent),]
        motif <- motif[,"id",drop=F]
        motif <- motif[!duplicated(motif$id),,drop=F]
        motif$tissue <- tissue_label_change(tissue)
        motif_summary[[condition]] <- rbind(motif_summary[[condition]],motif) 
      }
    }
  }
}

increase_count <- motif_summary[["increase"]] %>%   
  dplyr::count(id)
increase_tissue <- motif_summary[["increase"]] %>%   
  group_by(id) %>%   
  summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
increase_count <- merge(increase_count,increase_tissue,by="id")

decrease_count <-motif_summary[["decrease"]] %>%   
  dplyr::count(id)
decrease_tissue <- motif_summary[["decrease"]] %>%   
  group_by(id) %>%   
  summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
decrease_count <- merge(decrease_count,decrease_tissue,by="id")

increase_count <- increase_count[order(increase_count$n,decreasing = TRUE),]
decrease_count <- decrease_count[order(decrease_count$n,decreasing = TRUE),]
write.csv(increase_count,"data/samples/WGBS/all/snapatac2/mutual_bg/increase_common_motif_count.csv")
write.csv(decrease_count,"data/samples/WGBS/all/snapatac2/mutual_bg/decrease_common_motif_count.csv")
motif_count_summary$condition <- factor(motif_count_summary$condition,levels=c("increase","decrease"))
motif_count_summary$position <- motif_count_summary$count
motif_count_summary$position[which(motif_count_summary$condition=="decrease")] <- (-motif_count_summary$position[which(motif_count_summary$condition=="decrease")])
summary_by_tissue <- motif_count_summary %>%
  group_by(tissue) %>%
  summarise(count_sum = sum(count))
summary_by_tissue <- summary_by_tissue[order(summary_by_tissue$count_sum),]
motif_count_summary$tissue <- factor(motif_count_summary$tissue,levels=summary_by_tissue$tissue)
tissue_order <-c("Mammary Gland","Cecum","Thymus","Uterus","iWAT","Stomach","Skin","Spleen","Muscle","Bone Marrow","Liver","Ileum","Testis","Cortex","Jejunum","Tongue","Hippocampus","Colon","Bladder",
                 "Aorta","Cerebellum","Lung","Heart","Kidney","BAT","Ovary","Pancreas")
motif_count_summary$tissue <- factor(motif_count_summary$tissue,levels=tissue_order)
ggplot(motif_count_summary, aes(x = tissue, y = ifelse(condition == "increase", count, -count), fill = condition)) +  
  geom_bar(stat = "identity") +  
  # labs(title = paste0("fdr < 0.05 motif count ",label), x = NULL, y = "Count") +
  labs(x = NULL, y = "Count") +
  theme_minimal() +  
  xlab(NULL)+
  scale_y_continuous(labels = abs) +  
  theme(  
    axis.title.x = element_text(size = 14),     
    axis.title.y = element_text(size = 14),    
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),   
    axis.text.y = element_text(size = 12),    
    plot.title = element_text(size = 16, face = "bold")
  ) +  
  geom_text(data =motif_count_summary,   
            aes(label = count, y = position),   
            color = "black", size = 5, vjust = 0.5) + 
  scale_fill_manual(values = c("increase" = "skyblue", "decrease" = "salmon"), name = NULL)  

#### scatter plot
scatter_plot <- merge(increase_count[,c("id","n")],decrease_count[,c("id","n")],by="id",all=T)

colnames(scatter_plot)[2:3] <- c("increase_num","decrease_num")
scatter_plot[is.na(scatter_plot)] <- 0
scatter_plot$total_num <- scatter_plot$increase_num+scatter_plot$decrease_num

scatter_plot$diff <- scatter_plot$increase_num - scatter_plot$decrease_num
write.csv(scatter_plot,"data/samples/WGBS/all/snapatac2/mutual_bg/scatter_plot_DMR_motif_summary.csv")
top_motif <- scatter_plot[which(scatter_plot$id %in% c("Irf2","Irf8","Jund","Nfe2l1","Fos","Stat2","Batf","Dbp","Irx3","Nfic",
                                                       "Fbxl19","Foxn1","Zbtb1","Egr2","Zfp777","Zfp189","E2f4","Lin28a","Tfdp2","Dnmt1")),]
# top_motif <- scatter_plot[which(scatter_plot$id %in% c("Jund","Junb","Fos","E2f4")),]

p <- ggplot(scatter_plot, aes(x = increase_num, y = decrease_num, color=diff)) +
  geom_point(size = 2) +
  labs(x = "# of tissues enriched, Hyper-DMRs", y = "# of tissues enriched, Hypo-DMRs",color = "Total tissue number") +
  theme_minimal()+
  xlim(0,20)+
  ylim(0,20)+
  scale_color_gradient2(low = "#3c5488", high = "#e64b35", midpoint = 0)+
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(hjust = 1, size = 12),  
    axis.text.y = element_text(size = 12),  
    plot.title = element_text(size = 16, face = "bold"),
    panel.background = element_blank(),  
    panel.grid.major = element_line(size = 0.1, linetype = 'solid', color = "grey"),
    panel.grid.minor = element_line(size = 0.1, linetype = 'solid', color = "lightgrey"),
    panel.border = element_rect(color = "black", fill = NA, size = 1) 
  )+
  geom_text_repel(data = top_motif, aes(x = increase_num, y = decrease_num, label = id), max.overlaps=100,
                  size = 5, 
                  nudge_y = 0.2)
ggsave("result/figures/DMR_snapatac2_scatter_plot.pdf",p,width = 7,height = 5.5)
#### heatmap
motifs <- unique(c(increase_count$id[which(increase_count$n >=12)],decrease_count$id[which(decrease_count$n>=12)]))

to_plot <- motif_data_frame[["increase"]]
to_plot <- to_plot[which(to_plot$id %in% motifs),]
rownames(to_plot) <- to_plot$id
to_plot <- to_plot[,-1]
replace_zeros <- function(column) {
  non_zero_min <- min(column[column != 0], na.rm = TRUE)
  if (is.infinite(non_zero_min)) {
    non_zero_max <- 0
  }
  column[column == 0] <- non_zero_min
  return(column)
}
to_plot <- as.data.frame(apply(to_plot, 2, replace_zeros))
to_plot <- as.data.frame(apply(to_plot, 2, function(column) -log10(column)))

rowmeans <- as.data.frame(rowMeans(to_plot))
colnames(rowmeans) <- "means"
rowmeans <- rowmeans[order(rowmeans$means,decreasing = T),,drop=F]
to_plot <- to_plot[rownames(rowmeans),]

colmeans <- as.data.frame(colMeans(to_plot))
colnames(colmeans) <- "means"
colmeans  <- colmeans [order(colmeans$means,decreasing = T),,drop=F]
to_plot <- to_plot[,rownames(colmeans)]
breaks <- c(seq(-15,-(-log10(0.05)+0.0001), length.out = 80),seq(log10(0.05),-log10(0.05), length.out = 40),seq(-log10(0.05)+0.0001, 15, length.out = 80))
# breaks <- c(seq(-8,-0.51, length.out = 80),seq(-0.5,0.5, length.out = 40),seq(0.51, 8, length.out = 80))
color_palette <- c(
  colorRampPalette(c("blue","#defcf9"))(80),
  rep("white", 20), 
  rep("white", 20), 
  colorRampPalette(c("#ffe2e2","red"))(80) 
)
pheatmap::pheatmap(to_plot, cluster_rows =F,cluster_cols = F,color = color_palette,breaks=breaks,show_rownames = T,fontsize = 8)

to_plot <- motif_data_frame[["decrease"]]
rownames(to_plot) <- to_plot$id
to_plot <- to_plot[,-1]
replace_zeros <- function(column) {
  non_zero_min <- min(column[column != 0], na.rm = TRUE)
  if (is.infinite(non_zero_min)) {
    non_zero_max <- 0
  }
  column[column == 0] <- non_zero_min
  return(column)
}
to_plot <- as.data.frame(apply(to_plot, 2, replace_zeros))
to_plot <- as.data.frame(apply(to_plot, 2, function(column) -log10(column)))
to_plot <- to_plot[rownames(rowmeans),]
to_plot <- to_plot[,rownames(colmeans)]
to_plot <- -to_plot
pheatmap::pheatmap(to_plot, cluster_rows =F,cluster_cols = F,color = color_palette,breaks=breaks,show_rownames = T,fontsize = 8)
