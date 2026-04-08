rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyverse)  

options(scipen = 999)  
tissue <- "pancreas"
resolution <- "50000"
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
check_row <- function(row) {  
  all(row[-1] == row[-1][1])  
}  
compartment_cluster <- function(tissue,resolution){
  df_list <- list(young=list(),old=list())
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  to_plot <- data.frame()
  for(sample in search_table$sample_name){
    if(sample=="WJH-Liver-103"){
      df <- read.table(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_",resolution,".PC1_recorrect.txt"))
    }else{
      df <- read.table(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_",resolution,".PC1.txt"))
    }
    # df <- read.table(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_",resolution,".PC1.txt"))
    # df <- df[which(df$V2 %in% c(paste0("chr",c(1:19,"X","Y")))),]
    df <- df[which(df$V2 %in% c(paste0("chr",c(1:19,"X")))),]
    df <- df[,c(1,6)]
    colnames(df)[2] <- sample
    if(nrow(to_plot)==0){
      to_plot <- df
    }else{
      to_plot <- merge(to_plot,df,by="V1")
    }
  }
  rownames(to_plot) <- to_plot$V1
  to_plot <- to_plot[,-1]
  search_table$age <- factor(search_table$age,levels=c("3M","24M"))
  annotation_col <- search_table[,c("sample_name","age")]
  rownames(annotation_col) <- annotation_col$sample_name
  annotation_col <- annotation_col[,"age",drop=F]
  color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
  breaks <- c(seq(-2, -0.6, length.out = 40), seq(-0.5, 0.5, length.out = 20), seq(0.6, 2, length.out = 40))
  p <- pheatmap::pheatmap(to_plot,cluster_rows = T,show_rownames = F,color = color_palette,breaks = breaks,main = paste(tissue_label_change(tissue),"homer compartment PC1"),annotation_col = annotation_col)
  print(p)
  }

for(tissue in tissues){
  compartment_cluster(tissue,resolution)
}


compartment_change <- function(tissue,resolution){
  ages <- c("young","old")
  df_list <- list(young=list(),old=list())
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  young_samples <- search_table$sample_name[which(search_table$age=="3M")]
  old_samples <- search_table$sample_name[which(search_table$age=="24M")]
  samples_list <- list(young=young_samples,old=old_samples)
  
  for(age in ages){
    for(i in c(1:length(samples_list[[age]]))){
      sample <- samples_list[[age]][i]
      if(sample=="WJH-Liver-103"){
        df <- read.table(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_",resolution,".PC1_recorrect.txt"))
      }else{
        df <- read.table(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_",resolution,".PC1.txt"))
      }
      df <- df[which(df$V2 %in% c(paste0("chr",c(1:19,"X")))),]
      df <- df[,c(1,6)]
      df$compartmeent <- ifelse(df[, 2] > 0, "A", "B")
      colnames(df)[3] <- sample
      df <- df[,c(1,3)]
      df_list[[age]][[i]] <- df
      names(df_list[[age]])[i]<-sample
    }
  }

  young_data <- Reduce(function(x, y) merge(x, y, by = "V1"), df_list[["young"]])  
  old_data <- Reduce(function(x, y) merge(x, y, by = "V1"), df_list[["old"]])  
  
  young_data <- young_data[apply(young_data, 1, check_row), ] 
  young_data <- young_data[,c(1,2)]
  colnames(young_data)[2] <- "young"

  old_data <- old_data[apply(old_data, 1, check_row), ] 
  old_data <- old_data[,c(1,2)]
  colnames(old_data)[2] <- "old"
  
  merged_data <- merge(young_data,old_data,by="V1")
  
  merged_data$condition <- paste0(merged_data$young,"-",merged_data$old)
  to_plot <- as.data.frame(table(merged_data$condition))
  to_plot$Percentage <- to_plot$Freq/sum(to_plot$Freq)*100
  to_plot$Label <- paste0(to_plot$Var1, " (", round(to_plot$Percentage, 1), "%)")
  colors <- read.table("data/samples/7_distinct_color.txt")
  colors <-setNames(colors$V1,to_plot$Label)
  p <- ggplot(to_plot, aes(x = "", y = Freq, fill = Label)) +  
    geom_bar(width = 1, stat = "identity", color = "white") +  
    scale_fill_manual(values = colors)+
    coord_polar("y", start = 0) +  
    theme_void() + 
    labs(fill = NULL) +  
    ggtitle(paste0(tissue_label_change(tissue))) + 
    theme(legend.position = "right",plot.title = element_text(hjust = 0.5),text = element_text(size = 16))
  ggsave(paste0("result/HiC/",tissue,"/compartment/",tissue,"_homer_compartment_change.png"),p,width = 5,height = 4,type="cairo")
  merged_data <- merged_data %>% separate(V1, into = c("chr", "start"), sep = "-")  
  merged_data$start <- as.numeric(merged_data$start)
  merged_data$end <- merged_data$start + as.numeric(resolution)
  merged_data <- merged_data[,c("chr", "start", "end", "young", "old", "condition")]  
  write.csv(merged_data, paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/compartment_change_",resolution,".csv"),row.names = F)
}
for(tissue in tissues){
  compartment_change(tissue,resolution)
}
compartment_change_summary <- data.frame()
tissues <- c("brain","CB", "kidney", "liver", "lung", "bonemarrow", "colon", "heart", "Hip", "mammarygland", "stomach", "thymus","skin","muscle")
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/compartment_change_",resolution,".csv"))  
  df <- as.data.frame(table(df$condition))
  df$Percentage <- df$Freq/sum(df$Freq)*100
  df <- df[which(df$Var1 %in% c("A-B","B-A")), c("Var1","Percentage")]
  df$tissue <- tissue_label_change(tissue)
  compartment_change_summary <- rbind(compartment_change_summary,df)
}
summary_data <- compartment_change_summary %>%  
  group_by(tissue) %>%  
  summarise(Total_Percentage = sum(Percentage)) %>%  
  arrange(desc(Total_Percentage))
compartment_change_summary$Percentage[which(compartment_change_summary$Var1=="A-B")] <- -(compartment_change_summary$Percentage[which(compartment_change_summary$Var1=="A-B")])
compartment_change_summary$tissue <- factor(compartment_change_summary$tissue, levels=summary_data$tissue)
ggplot(compartment_change_summary, aes(x = Percentage, y = tissue, fill = Var1)) +  
  geom_bar(stat = "identity")+
  scale_x_continuous(labels = abs)+
  labs(x = "Percentage(%)", y = NULL, fill = "Comparison") +  
  theme_minimal() + 
  ggtitle("Compartment change proportion") +
  scale_fill_manual(values = c("A-B" = "skyblue", "B-A" =  "salmon")) +  
  theme(  
    axis.title.x = element_text(size = 14),     
    axis.title.y = element_text(size = 14),    
    axis.text.x = element_text(hjust = 1, size = 12),   
    axis.text.y = element_text(size = 12),    
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "bottom"
  ) 

#####generate bed
resolution <- "50000"
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/compartment_change_",resolution,".csv"))
  write.table(df[,c(1:3)],paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/compartment_change_",resolution,".bed"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)
}
