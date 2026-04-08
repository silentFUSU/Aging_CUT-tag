rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(bitmapType="cairo")  
library(grid)
library(ggnewscale)
options(scipen=0)
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

top <- 10
conditions <- c("up","down")
tissues <-  sort(c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
              "thymus","skin","bladder","bonemarrow","Hip","heart",
              "muscle","jejunum","uterus","ovary","liver","tongue",
              "cecum","colon","testis","stomach","pancreas","iWAT","ileum"))
data_path <- "data/samples/ATAC/all/ATAC/snapatac2_macs/strict_stable_peaks_summits_spm3/"
p_list <-list()
tissue_p_value <- data.frame()
for(tissue in tissues){
  to_plot <- data.frame()
  max_value <- 0
  for(condition in conditions){
    if(file.exists(paste0(data_path,condition,"/enrichment_results_",tissue,"_summits_spm3.bed.csv"))){
      motif <- read.csv(paste0(data_path,condition,"/enrichment_results_",tissue,"_summits_spm3.bed.csv"))
      motif <- motif[which(motif$log2.fold.change. > 0 & motif$adjusted.p.value < 0.05),]
      
      if(nrow(motif) > 0){
        min_non_zero <- min(motif$adjusted.p.value[motif$adjusted.p.value != 0])
        max_larger_zero_without_Inf <- max(motif$log2.fold.change.[is.finite(motif$log2.fold.change.)])
        motif$id <- sub("^M\\d+_2\\.00\\s*", "", motif$id) 
        motif$id <- sapply(motif$id, function(x) {
          matches <- str_extract_all(x, "\\(([^()]+)\\)")
          matches <- unlist(matches)
          if (length(matches) >= 2) {
            return(paste(matches[1], matches[2], sep = "_"))
          } else {
            return(x)  
          }
        })
        
        motif <- motif[order(motif$adjusted.p.value, -motif$log2.fold.change.),]
        motif <- motif[c(1:min(nrow(motif),top)),]
        motif$adjusted.p.value[motif$adjusted.p.value == 0] <- min_non_zero
        motif$log2.fold.change.[is.infinite(motif$log2.fold.change.)] <- max_larger_zero_without_Inf
        motif$log10fdr <- -log10(motif$adjusted.p.value) 
        if(max(motif$log10fdr) > max_value){
          max_value <- max(motif$log10fdr)
        }
        if(condition == "down"){
          motif$log2.fold.change. <- -motif$log2.fold.change.
        }
        to_plot <- rbind(to_plot,motif)
      }
    }
  }
  # to_plot <- to_plot[,c("id","log2.fold.change.","adjusted.p.value")]
  t_tissue_pvalue <- data.frame(tissue=tissue,pvalue=max_value)
  tissue_p_value <- rbind(tissue_p_value,t_tissue_pvalue)
  if(nrow(to_plot)==0){
    next
  }
  if (nrow(to_plot) < 20) {
    num_missing_rows <- 20 - nrow(to_plot)
    
    empty_rows <- data.frame(
      id = paste0("Empty", 1:num_missing_rows),
      log2.fold.change. = rep(0, num_missing_rows),
      adjusted.p.value = rep(NA, num_missing_rows),
      change_direction = rep(NA, num_missing_rows)
    )
    
    to_plot <- bind_rows(to_plot, empty_rows)
  }
  to_plot$log10fdr_label <- to_plot$log10fdr
  to_plot$log10fdr_label[which(to_plot$log2.fold.change. < 0 )] <- -to_plot$log10fdr_label[which(to_plot$log2.fold.change. < 0 )]
  to_plot <- to_plot[order(!grepl("^Empty", to_plot$id), to_plot$log10fdr_label, to_plot$`log2.fold.change.`, decreasing = TRUE),]
  to_plot$id <- factor(to_plot$id,levels = unique(rev(to_plot$id)))
  to_plot <- to_plot %>%
    mutate(change_direction = ifelse(`log2.fold.change.` > 0, "positive", "negative"))
  to_plot <- to_plot[,c("id","log2.fold.change.","log10fdr","change_direction")]
  p_list[[tissue]] <- ggplot() +
    geom_bar(data = subset(to_plot, change_direction == "negative"), 
             aes(x = `log2.fold.change.`, y = id, 
                 fill = log10fdr), 
             stat = "identity", position = position_dodge2()) 
  
  if (length(unique(to_plot$log10fdr[to_plot$change_direction == "negative"])) == 1) {
    p_list[[tissue]] <- p_list[[tissue]] + scale_fill_gradient(low = "blue", high = "blue")+
      new_scale_fill()
  } else {
    p_list[[tissue]] <- p_list[[tissue]] + scale_fill_gradient(low = "#e2e2ff", high = "blue")+
      new_scale_fill()
  }
  p_list[[tissue]] <-  p_list[[tissue]] + 
    geom_bar(data = subset(to_plot, change_direction == "positive"), 
             aes(x = `log2.fold.change.`, y = id, fill = log10fdr), stat = "identity", position = position_dodge2())
  
  if (length(unique(to_plot$log10fdr[to_plot$change_direction == "positive"])) == 1) {
    p_list[[tissue]] <- p_list[[tissue]] + scale_fill_gradient(low = "red", high = "red")
  } else {
    p_list[[tissue]] <- p_list[[tissue]] + scale_fill_gradient(low = "pink", high = "red")
  }
  p_list[[tissue]] <- p_list[[tissue]]+
    theme_bw() +
    ylab("") +
    xlab("log2(Fold change)")+
    theme(text = element_text(size = 13),
          axis.text.y = element_text(size = 18, face = "bold"),
          plot.title = element_text(size = 18, face = "bold"))+
    ggtitle(tissue_label_change(tissue))+
    scale_y_discrete(labels = function(labels) {
      labels[grepl("^Empty", labels)] <- ""
      labels
    })
}
tissue_p_value <- tissue_p_value[order(tissue_p_value$pvalue,decreasing = T),]
p_list_tissue_order <- p_list[tissue_p_value$tissue]
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
p_list_tissue_order <- Filter(Negate(is.null), p_list_tissue_order)

combined_plot <- plot_a_list(p_list_tissue_order,no_of_rows = 4,no_of_cols = 7)
ggsave("result/all/ATAC/all_tissues_top_motif_macs2_summits_spm3_strict_bg_snapatact2_stable_peaks_background.png",combined_plot,width = 70,height = 35,limitsize = FALSE)

combined_plot <- plot_a_list(p_list_tissue_order[c(1:14)],no_of_rows = 3,no_of_cols = 5)
ggsave("result/all/ATAC/all_tissues_top_motif_macs2_summits_spm3_strict_bg_snapatact2_stable_peaks_background_1_14.png",combined_plot,width = 50,height = 25,limitsize = FALSE)

combined_plot <- plot_a_list(p_list_tissue_order[c(15:length(p_list_tissue_order))],no_of_rows = 3,no_of_cols = 5)
ggsave("result/all/ATAC/all_tissues_top_motif_macs2_summits_spm3_strict_bg_snapatact2_stable_peaks_background_15_27.png",combined_plot,width = 50,height = 25,limitsize = FALSE)
