rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
options(scipen = 999)  
tissue <- "lung"
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
    }else if(tissue_label == "Mammarygland"){
      tissue_label <- "Mammary Gland"
    }else if(tissue_label == "Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
} 

compartment_change <- function(tissue){
  ages <- c("young","old")
  df_list <- list()
  for(age in ages){
    df <- read.table(paste0("data/samples/HiC/",tissue,"/fanc_compartment/first_eigenvector/fanc_",age,"_1mb.compartment.bed"))
    df <- df[which(df$V1 %in% c(paste0("chr",c(1:19,"X","Y")))),]
    df_new <- data_frame()
    for(j in c(1:nrow(df))){
      if(df[j,3]-df[j,2]+1 > 1000000){
        V1=df[j,1]  
        V4=df[j,4]
        V5=df[j,5]
        V6=df[j,6]
        for(start in seq(df[j,2],df[j,3], by = 1000000)){
          t_df <- data.frame(V1=V1,V2=start,V3=min(start+999999,df[j,3]),V4=V4,V5=V5,V6=V6)
          df_new <- rbind(df_new,t_df)
        }
      }else{
        df_new <- rbind(df_new,df[j,])
      }
    }
    df_list[[age]] <- df_new
    colnames(df_list[[age]])[4] <- age
    df_list[[age]]$label <- paste0(df_list[[age]]$V1,"-",df_list[[age]]$V2,"-",df_list[[age]]$V3)
  }
  if(tissue == "lung"){
    df_list[["old"]][which(df_list[["old"]]$V1=="chr19"),5] <- -df_list[["old"]][which(df_list[["old"]]$V1=="chr19"),5]
    df_list[["old"]] <- df_list[["old"]] %>%  
      mutate(old = ifelse(V1 == "chr19",  
                         case_when(  
                           old == "A" ~ "B",  
                           old == "B" ~ "A",  
                           TRUE ~ old  
                         ),  
                         old)) 
  }
  extracted_columns <- lapply(df_list, function(df) df[, c(4, 7)])  
  merged_data <- Reduce(function(x, y) merge(x, y, by = "label", all = TRUE), extracted_columns)  
  merged_data$condition <- paste0(merged_data$young,"-",merged_data$old)
  to_plot <- as.data.frame(table(merged_data$condition))
  to_plot$Percentage <- to_plot$Freq/sum(to_plot$Freq)*100
  to_plot$Label <- paste0(to_plot$Var1, " (", round(to_plot$Percentage, 1), "%)")
  colors <- read.table("data/samples/7_distinct_color.txt")
  colors <-setNames(colors$V1,to_plot$Label)
  ggplot(to_plot, aes(x = "", y = Freq, fill = Label)) +  
    geom_bar(width = 1, stat = "identity", color = "white") +  
    scale_fill_manual(values = colors)+
    coord_polar("y", start = 0) +  
    theme_void() + 
    labs(fill = NULL) +  
    ggtitle(paste0(tissue_label_change(tissue))) + 
    theme(legend.position = "right",plot.title = element_text(hjust = 0.5),text = element_text(size = 16))
}
