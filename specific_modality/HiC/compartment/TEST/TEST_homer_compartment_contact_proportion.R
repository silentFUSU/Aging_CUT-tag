rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(data.table)
tissue <- "lung"
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

compartment_contact <- function(tissue,resolution){
  compartment_list <- list()
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  for(i in c(1:length(search_table$sample_name))){
    sample <- search_table$sample_name[i]
    df <- read.table(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_",resolution,".PC1.txt"))
    df <- df[which(df$V2 %in% c(paste0("chr",c(1:19,"X","Y")))),]
    df <- df[,c(2:4,6)]
    df$compartment <- ifelse(df[, 4] > 0, "A", "B")
    compartment_list[[i]] <- df
    names(compartment_list)[i]<-sample
  }
  
  contact_list <- list()
  for(i in c(1:length(search_table$sample_name))){
    sample <- search_table$sample_name[i]
    bed <- read.table(paste0("data/samples/HiC/",tissue,"/raw_matrix/",sample,"_",resolution,"_abs.bed"))
    contact <- fread(paste0("data/samples/HiC/",tissue,"/ice_matrix/",sample,"_",resolution,"_iced.matrix"))
    contact <- contact[abs(V1 - V2) > 20] 
    bed$label <- paste0(bed$V1,"-",bed$V2,"-",bed$V3)
    compartment_list[[sample]]$label <- paste0(compartment_list[[sample]]$V2,"-",compartment_list[[sample]]$V3,"-",compartment_list[[sample]]$V4)
    bed <- merge(bed[,c("label","V4")],compartment_list[[sample]][,c("label","compartment")],by="label")
    bed <- as.data.table(bed)
    setkey(bed,V4)
    contact <- contact[bed, on = .(V1 = V4), V1.compartment := i.compartment]
    contact <- contact[bed, on = .(V2 = V4), V2.compartment := i.compartment]
    contact <- contact[!is.na(V1.compartment) & !is.na(V2.compartment)]  
    contact[, compartment_pair := paste0(V1.compartment, "-", V2.compartment)]  
    sum_by_compartment <- contact[, .(V3_sum = sum(V3)), by = compartment_pair]
    sum_by_compartment <- as.data.frame(sum_by_compartment)
    colnames(sum_by_compartment)[2]<- sample
    contact_list[[i]] <- sum_by_compartment
  }
  
  for(i in c(1:length(contact_list))){
    sample <- search_table$sample_name[i]
    sum <- sum(contact_list[[i]][,2])
    contact_list[[i]][,2] <- contact_list[[i]][,2]/sum*100
    contact_list[[i]] <- data.frame(condition=c("A-A","B-B","A-B"),percent=c(contact_list[[i]][which(contact_list[[i]]$compartment_pair=="A-A"),2],
                                                                             contact_list[[i]][which(contact_list[[i]]$compartment_pair=="B-B"),2],
                                                                             contact_list[[i]][which(contact_list[[i]]$compartment_pair=="A-B"),2]+contact_list[[i]][which(contact_list[[i]]$compartment_pair=="B-A"),2]))
    colnames(contact_list[[i]])[2] <- sample
  } 
  
  to_plot <- Reduce(function(x, y) merge(x, y, by = "condition"), contact_list)  
  
  to_plot_melt <- reshape2::melt(to_plot)
  to_plot_melt$condition <- factor(to_plot_melt$condition,c("A-A","B-B","A-B"))
  to_plot_melt$position <- to_plot_melt$value
  to_plot_melt$position[which(to_plot_melt$condition=="A-A")] <- 100
  to_plot_melt$position[which(to_plot_melt$condition=="B-B")] <- 100- to_plot_melt$value[which(to_plot_melt$condition=="A-A")]
  ggplot(to_plot_melt, aes(x = variable, y = value, fill = condition)) +  
    geom_bar(width = 1, stat = "identity", color = "white") +
    theme_minimal()+
    labs(fill = NULL) + 
    xlab(NULL)+
    ylab("Proportion")+
    ggtitle(paste0(tissue_label_change(tissue))) +   
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right",
          plot.title = element_text(hjust = 0.5),
          text = element_text(size = 16))+
    geom_text(data =to_plot_melt,   
              aes(label = paste0(round(to_plot_melt$value,1),"%"), y = position),   
              color = "black", size = 5, vjust = 0.5)
  
  to_plot_heatmap <- to_plot
  to_plot_heatmap$condition <- factor(to_plot_heatmap$condition, levels=c("A-A","B-B","A-B"))
  rownames(to_plot_heatmap) <- to_plot$condition  
  to_plot_heatmap <- to_plot_heatmap[order(to_plot_heatmap$condition),]

  to_plot_heatmap <- to_plot_heatmap[,-1]
  
  annotation <- search_table[,c("sample_name","age")]
  annotation$age <- factor(annotation$age,c("3M","24M"))
  rownames(annotation) <- annotation$sample_name
  annotation <- annotation[,-1,drop=F]
  pheatmap::pheatmap(to_plot_heatmap,main = tissue_label_change(tissue),cluster_cols = F,cluster_rows=F, annotation_col = annotation, scale = "row")

}


compartment_contact_cis <- function(tissue,resolution){
  compartment_list <- list()
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  for(i in c(1:length(search_table$sample_name))){
    sample <- search_table$sample_name[i]
    df <- read.table(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_",resolution,".PC1.txt"))
    df <- df[which(df$V2 %in% c(paste0("chr",c(1:19,"X","Y")))),]
    df <- df[,c(2:4,6)]
    df$compartment <- ifelse(df[, 4] > 0, "A", "B")
    compartment_list[[i]] <- df
    names(compartment_list)[i]<-sample
  }
  
  contact_list <- list()
  for(i in c(1:length(search_table$sample_name))){
    sample <- search_table$sample_name[i]
    bed <- read.table(paste0("data/samples/HiC/",tissue,"/raw_matrix/",sample,"_",resolution,"_abs.bed"))
    contact <- fread(paste0("data/samples/HiC/",tissue,"/ice_matrix/",sample,"_",resolution,"_iced.matrix"))
    contact <- contact[abs(V1 - V2) > 20] 
    bed$label <- paste0(bed$V1,"-",bed$V2,"-",bed$V3)
    compartment_list[[sample]]$label <- paste0(compartment_list[[sample]]$V2,"-",compartment_list[[sample]]$V3,"-",compartment_list[[sample]]$V4)
    bed <- merge(bed[,c("V1","label","V4")],compartment_list[[sample]][,c("label","compartment")],by="label")
    bed <- as.data.table(bed)
    setkey(bed,V4)
    contact <- contact[bed, on = .(V1 = V4), V1.compartment := i.compartment]
    contact <- contact[bed, on = .(V2 = V4), V2.compartment := i.compartment]
    contact <- contact[!is.na(V1.compartment) & !is.na(V2.compartment)]  
    contact <- contact[bed, on = .(V1 = V4), V1.chr := i.V1]
    contact <- contact[bed, on = .(V2 = V4), V2.chr := i.V1]
    contact <- contact[V1.chr==V2.chr] 
    contact[, compartment_pair := paste0(V1.compartment, "-", V2.compartment)]  
    sum_by_compartment <- contact[, .(V3_sum = sum(V3)), by = compartment_pair]
    sum_by_compartment <- as.data.frame(sum_by_compartment)
    colnames(sum_by_compartment)[2]<- sample
    contact_list[[i]] <- sum_by_compartment
  }
  
  for(i in c(1:length(contact_list))){
    sample <- search_table$sample_name[i]
    sum <- sum(contact_list[[i]][,2])
    contact_list[[i]][,2] <- contact_list[[i]][,2]/sum*100
    contact_list[[i]] <- data.frame(condition=c("A-A","B-B","A-B"),percent=c(contact_list[[i]][which(contact_list[[i]]$compartment_pair=="A-A"),2],
                                                                             contact_list[[i]][which(contact_list[[i]]$compartment_pair=="B-B"),2],
                                                                             contact_list[[i]][which(contact_list[[i]]$compartment_pair=="A-B"),2]+contact_list[[i]][which(contact_list[[i]]$compartment_pair=="B-A"),2]))
    colnames(contact_list[[i]])[2] <- sample
  } 
  
  to_plot <- Reduce(function(x, y) merge(x, y, by = "condition"), contact_list)  
  
  to_plot_melt <- reshape2::melt(to_plot)
  to_plot_melt$condition <- factor(to_plot_melt$condition,c("A-A","B-B","A-B"))
  
  to_plot_melt$position <- to_plot_melt$value
  to_plot_melt$position[which(to_plot_melt$condition=="A-A")] <- 100
  to_plot_melt$position[which(to_plot_melt$condition=="B-B")] <- 100- to_plot_melt$value[which(to_plot_melt$condition=="A-A")] 
  ggplot(to_plot_melt, aes(x = variable, y = value, fill = condition)) +  
    geom_bar(width = 1, stat = "identity", color = "white") +
    theme_minimal()+
    labs(fill = NULL) + 
    xlab(NULL)+
    ylab("Proportion")+
    ggtitle(paste0(tissue_label_change(tissue))) +   
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right",
          plot.title = element_text(hjust = 0.5),
          text = element_text(size = 16))+
    geom_text(data =to_plot_melt,   
              aes(label = paste0(round(to_plot_melt$value,1),"%"), y = position),   
              color = "black", size = 5, vjust = 0.5)
  
  to_plot_heatmap <- to_plot
  to_plot_heatmap$condition <- factor(to_plot_heatmap$condition, levels=c("A-A","B-B","A-B"))
  rownames(to_plot_heatmap) <- to_plot$condition  
  to_plot_heatmap <- to_plot_heatmap[order(to_plot_heatmap$condition),]
  
  to_plot_heatmap <- to_plot_heatmap[,-1]
  
  
  annotation <- search_table[,c("sample_name","age")]
  annotation$age <- factor(annotation$age,c("3M","24M"))
  rownames(annotation) <- annotation$sample_name
  annotation <- annotation[,-1,drop=F]
  pheatmap::pheatmap(to_plot_heatmap,main = tissue_label_change(tissue),cluster_cols = F,cluster_rows=F, annotation_col = annotation, scale = "row")
  
}


compartment_contact_trans <- function(tissue,resolution){
  compartment_list <- list()
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  for(i in c(1:length(search_table$sample_name))){
    sample <- search_table$sample_name[i]
    df <- read.table(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_",resolution,".PC1.txt"))
    df <- df[which(df$V2 %in% c(paste0("chr",c(1:19,"X","Y")))),]
    df <- df[,c(2:4,6)]
    df$compartment <- ifelse(df[, 4] > 0, "A", "B")
    compartment_list[[i]] <- df
    names(compartment_list)[i]<-sample
  }
  
  contact_list <- list()
  for(i in c(1:length(search_table$sample_name))){
    sample <- search_table$sample_name[i]
    bed <- read.table(paste0("data/samples/HiC/",tissue,"/raw_matrix/",sample,"_",resolution,"_abs.bed"))
    contact <- fread(paste0("data/samples/HiC/",tissue,"/ice_matrix/",sample,"_",resolution,"_iced.matrix"))
    contact <- contact[abs(V1 - V2) > 1] 
    bed$label <- paste0(bed$V1,"-",bed$V2,"-",bed$V3)
    compartment_list[[sample]]$label <- paste0(compartment_list[[sample]]$V2,"-",compartment_list[[sample]]$V3,"-",compartment_list[[sample]]$V4)
    bed <- merge(bed[,c("V1","label","V4")],compartment_list[[sample]][,c("label","compartment")],by="label")
    bed <- as.data.table(bed)
    setkey(bed,V4)
    contact <- contact[bed, on = .(V1 = V4), V1.compartment := i.compartment]
    contact <- contact[bed, on = .(V2 = V4), V2.compartment := i.compartment]
    contact <- contact[!is.na(V1.compartment) & !is.na(V2.compartment)]  
    contact <- contact[bed, on = .(V1 = V4), V1.chr := i.V1]
    contact <- contact[bed, on = .(V2 = V4), V2.chr := i.V1]
    contact <- contact[V1.chr!=V2.chr] 
    contact[, compartment_pair := paste0(V1.compartment, "-", V2.compartment)]  
    sum_by_compartment <- contact[, .(V3_sum = sum(V3)), by = compartment_pair]
    sum_by_compartment <- as.data.frame(sum_by_compartment)
    colnames(sum_by_compartment)[2]<- sample
    contact_list[[i]] <- sum_by_compartment
  }
  
  for(i in c(1:length(contact_list))){
    sample <- search_table$sample_name[i]
    sum <- sum(contact_list[[i]][,2])
    contact_list[[i]][,2] <- contact_list[[i]][,2]/sum*100
    contact_list[[i]] <- data.frame(condition=c("A-A","B-B","A-B"),percent=c(contact_list[[i]][which(contact_list[[i]]$compartment_pair=="A-A"),2],
                                                                             contact_list[[i]][which(contact_list[[i]]$compartment_pair=="B-B"),2],
                                                                             contact_list[[i]][which(contact_list[[i]]$compartment_pair=="A-B"),2]+contact_list[[i]][which(contact_list[[i]]$compartment_pair=="B-A"),2]))
    colnames(contact_list[[i]])[2] <- sample
  } 
  
  to_plot <- Reduce(function(x, y) merge(x, y, by = "condition"), contact_list)  
  
  to_plot_melt <- reshape2::melt(to_plot)
  to_plot_melt$condition <- factor(to_plot_melt$condition,c("A-A","B-B","A-B"))
  
  to_plot_melt$position <- to_plot_melt$value
  to_plot_melt$position[which(to_plot_melt$condition=="A-A")] <- 100
  to_plot_melt$position[which(to_plot_melt$condition=="B-B")] <- 100- to_plot_melt$value[which(to_plot_melt$condition=="A-A")] 
  ggplot(to_plot_melt, aes(x = variable, y = value, fill = condition)) +  
    geom_bar(width = 1, stat = "identity", color = "white") +
    theme_minimal()+
    labs(fill = NULL) + 
    xlab(NULL)+
    ylab("Proportion")+
    ggtitle(paste0(tissue_label_change(tissue))) +   
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right",
          plot.title = element_text(hjust = 0.5),
          text = element_text(size = 16))+
    geom_text(data =to_plot_melt,   
              aes(label = paste0(round(to_plot_melt$value,1),"%"), y = position),   
              color = "black", size = 5, vjust = 0.5)
  
  to_plot_heatmap <- to_plot
  to_plot_heatmap$condition <- factor(to_plot_heatmap$condition, levels=c("A-A","B-B","A-B"))
  rownames(to_plot_heatmap) <- to_plot$condition  
  to_plot_heatmap <- to_plot_heatmap[order(to_plot_heatmap$condition),]
  
  to_plot_heatmap <- to_plot_heatmap[,-1]
  pheatmap::pheatmap(to_plot_heatmap,main = tissue_label_change(tissue),cluster_cols = F,cluster_rows=F, annotation_col = annotation)
  
}