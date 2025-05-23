rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(data.table)
library(ggplot2)
library(ggalluvial)  

tissue <- "lung"
antibody <- "H3K27me3"

histone_change_in_compartment <- function(tissue, antibody){
  p_list <- list()
  histone <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_10kb_bins_diff_after_remove_batch_effect.csv"))
  compartment <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/compartment_change_50000.csv"))
  compartment$start <- compartment$start+1
  compartment <- as.data.table(compartment)
  setDT(compartment)
  setkey(compartment,chr,start,end)
  
  increase <- histone[which(histone$Significant == "Up"),]
  increase <- increase[,c(2:4)]
  increase$Start <- increase$Start + 1
  increase <- as.data.table(increase)
  setDT(increase)
  setkey(increase,Chr,Start,End)
  overlaps <- foverlaps(increase,compartment, type = "any", nomatch = 0L)  
  to_plot <- overlaps %>%  
    group_by(young,old) %>%  
    summarize(freq = n())  
  
  p_list[["increase"]] <- ggplot(to_plot, aes(axis1 = young, axis2 = old, y = freq)) +  
    geom_alluvium(aes(fill = young)) +  
    geom_stratum() +  
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +  
    theme_minimal() +  
    labs(y = "Count", x = "Compartment Transition", 
         fill = "Comparment")+
    ggtitle(paste0(tissue_label_change(tissue)," ",antibody," increased regions"),"Sankey Plot of State Transitions")
  
  
  
  
  decrease <- histone[which(histone$Significant == "Down"),]
  decrease <- decrease[,c(2:4)]
  decrease$Start <- decrease$Start + 1
  decrease <- as.data.table(decrease)
  setDT(decrease)
  setkey(decrease,Chr,Start,End)
  overlaps <- foverlaps(decrease,compartment, type = "any", nomatch = 0L)  
  to_plot <- overlaps %>%  
    group_by(young,old) %>%  
    summarize(freq = n())  
  
  p_list[["decrease"]] <- ggplot(to_plot, aes(axis1 = young, axis2 = old, y = freq)) +  
    geom_alluvium(aes(fill = young)) +  
    geom_stratum() +  
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +  
    theme_minimal() +  
    labs(y = "Count", x = "Compartment Transition", 
         fill = "Comparment")+
    ggtitle(paste0(tissue_label_change(tissue)," ",antibody," decreased regions"),"Sankey Plot of State Transitions")
  
  return(p_list)
}
tissues <-c("brain","CB", "kidney", "liver", "lung", "bonemarrow", "colon", "heart", "Hip", "mammarygland", "stomach", "thymus")
increase_p_list <- list()
decrease_p_list <- list()
for(tissue in tissues){
  p_list <- histone_change_in_compartment(tissue,"H3K27me3")
  increase_p_list[[tissue]] <- p_list[["increase"]]
  decrease_p_list[[tissue]] <- p_list[["decrease"]]
}
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}
increase_p_combined <- plot_a_list(increase_p_list,no_of_rows = 3,no_of_cols = 4)
ggsave(paste0("result/all/diff/H3K27me3/all_increased_bin_change_in_compartment.png"),increase_p_combined,width = 20,height = 15,type="cairo")
decrease_p_combined <- plot_a_list(decrease_p_list,no_of_rows = 3,no_of_cols = 4)
ggsave(paste0("result/all/diff/H3K27me3/all_decreased_bin_change_in_compartment.png"),decrease_p_combined,width = 20,height = 15,type="cairo")

