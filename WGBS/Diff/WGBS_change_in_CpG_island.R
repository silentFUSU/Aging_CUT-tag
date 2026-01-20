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
library(GenomeInfoDb)
library("GenomicRanges")
library(genomation)
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
cpg.file="~/ref_data/for_normal_mapping/mm10/cpgi.mm10.bed.txt"
for(tissue in tissues){
  print(tissue)
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  t_search_table <- search_table[which(search_table$tissue==tissue),]
  cpg.shore.obj=readFeatureFlank(cpg.file,flank = 2000,feature.flank.name=c("CpGi","shores"))
  cpgi <- as.data.frame(cpg.shore.obj@listData$CpGi)
  bin <- cpgi[,1:3]
  colnames(bin) <- c("V1","V2","V3") 
  if(tissue %in% c("ovary","uterus","mammarygland")){
    bin <- bin[which(bin$V1%in% paste0("chr",c(1:19,"X"))),]
    bin$V1 <- factor(bin$V1,levels=paste0("chr",c(1:19,"X")))
  }else{
    bin <- bin[which(bin$V1%in% paste0("chr",c(1:19,"X","Y"))),]
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
    
    result <- overlaps[, .(V4_sum = sum(V4), V5_sum = sum(V5))]
    result <- as.data.frame(result)
    result$methylation <- result$V4_sum/result$V5_sum*100
    result <- result[,c("V4_sum","V5_sum","methylation")]
    result$sample <- sample 
    summary <- rbind(summary,result)
  }
  
  t_tissue_summary <- data.frame(tissue=tissue_label_change(tissue),
                                 young_methylation=mean(summary$methylation[which(summary$sample %in% t_search_table$sample_name[which(t_search_table$age=="3M")])]),
                                 old_methylation=mean(summary$methylation[which(summary$sample %in% t_search_table$sample_name[which(t_search_table$age=="24M")])]),
                                 delta=mean(summary$methylation[which(summary$sample %in% t_search_table$sample_name[which(t_search_table$age=="24M")])])-mean(summary$methylation[which(summary$sample %in% t_search_table$sample_name[which(t_search_table$age=="3M")])]))

  tissue_summary <- rbind(tissue_summary,t_tissue_summary)
}
tissue_order <- c("Mammary Gland","Cecum","Thymus","Uterus","iWAT","Stomach","Skin","Spleen",
                  "Muscle","Bone Marrow","Liver","Ileum","Testis","Cortex","Jejunum","Tongue",
                  "Hippocampus","Colon","Bladder","Aorta","Cerebellum","Lung","Heart","Kidney","BAT","Ovary","Pancreas")
to_plot <- tissue_summary
to_plot$tissue <- factor(to_plot$tissue,levels = tissue_order)
to_plot$condition <- "Up"
to_plot$condition[which(to_plot$delta < 0 )] <- "Down"
color <- setNames(c("#f39b7f","#4dbbd5"),c("Up","Down"))
ggplot(to_plot, aes(x = tissue, y = delta, fill = condition)) +  
  geom_bar(stat = 'identity') +   
  theme_bw() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Delta")+
  ylim(-18,18) +
  ggtitle(NULL) +
  guides(fill = FALSE) +
  geom_hline(yintercept = c(-1, 1), color = "black", linetype = "dashed")
