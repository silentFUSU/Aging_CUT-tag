rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
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
    }
  }
  return(tissue_label)
} 
dup_map <- read.csv("data/samples/all/map_dup_qc.csv")
search_table_cut <- read.csv("data/samples/all/CUTTag_search_table.csv")
search_table_atac <- read.csv("data/samples/all/ATAC_search_table.csv")
search_table <- rbind(search_table_cut,search_table_atac)
colnames(dup_map)[1] <- "sample_name"
search_table$sample_name <- gsub("_", "", search_table$sample_name)  
search_table <- merge(search_table,dup_map[,c(1,3:4)],by="sample_name")

antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me3","H3K4me1","ATAC")
search_table$age <- factor(search_table$age,levels = c("3m","24m"))
for(i in c(1:length(antibodys))){
  antibody <- antibodys[i]
  p1 <- ggplot(search_table[which(search_table$antibody==antibody),],aes(x=tissue,y=map*100,color = age))+    
    geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.7)+
    geom_text(aes(label = mouse_ID), position = position_jitter(width = 0.2), vjust = -1, size = 3) +
    scale_fill_brewer(palette="Set3")+
    geom_hline(yintercept = 80, color = "red", linetype = "dashed") + 
    geom_hline(yintercept = 90, color = "blue", linetype = "dashed") + 
    annotate("text", x = Inf, y = 75, label = "80%", color = "red", vjust = -0.5, hjust = 1.1, size = 5) +   
    annotate("text", x = Inf, y = 85, label = "90%", color = "blue", vjust = -0.5, hjust = 1.1, size = 5) + 
    ggtitle(paste0(antibody," Mapping"))+ylim(0,100)+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")+labs(fill = "", color = "") +ylab("Mapping(%)")
  ggsave(paste0("result/all/QC/map/per_antibody/",antibody,"_map_age_dotplot.png"),p1,width = 15,height = 10,type="cairo")
  p2 <- ggplot(search_table[which(search_table$antibody==antibody),],aes(x=tissue,y=dup*100,color = age))+    
    geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.7)+
    geom_text(aes(label = mouse_ID), position = position_jitter(width = 0.2), vjust = -1, size = 3) +
    scale_fill_brewer(palette="Set3")+
    geom_hline(yintercept = 25, color = "red", linetype = "dashed") + 
    geom_hline(yintercept = 15, color = "blue", linetype = "dashed") + 
    annotate("text", x = Inf, y = 25, label = "25%", color = "red", vjust = -0.5, hjust = 1.1, size = 5) +   
    annotate("text", x = Inf, y = 15, label = "15%", color = "blue", vjust = -0.5, hjust = 1.1, size = 5) + 
    ggtitle(paste0(antibody," Duplication"))+ylim(0,100)+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")+labs(fill = "", color = "") +ylab("Duplication(%)")
  ggsave(paste0("result/all/QC/dup/per_antibody/",antibody,"_dup_age_dotplot.png"),p2,width = 15,height = 10,type="cairo")
}

tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas","Hip","cecum","bonemarrow","ileum","heart","thymus","stomach","skin","aorta","tongue","bladder","CB","jejunum","uterus","ovary")
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  p1 <- ggplot(search_table[which(search_table$tissue==tissue_label_change(tissue)),],aes(x=antibody,y=map*100,color = age))+    
    geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.7)+
    geom_text(aes(label = mouse_ID), position = position_jitter(width = 0.2), vjust = -1, size = 3) +
    scale_fill_brewer(palette="Set3")+
    geom_hline(yintercept = 80, color = "red", linetype = "dashed") + 
    geom_hline(yintercept = 90, color = "blue", linetype = "dashed") + 
    annotate("text", x = Inf, y = 75, label = "80%", color = "red", vjust = -0.5, hjust = 1.1, size = 5) +   
    annotate("text", x = Inf, y = 85, label = "90%", color = "blue", vjust = -0.5, hjust = 1.1, size = 5) + 
    ggtitle(paste0(tissue_label_change(tissue)," Mapping"))+ylim(0,100)+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")+labs(fill = "", color = "") +ylab("Mapping(%)")
  ggsave(paste0("result/all/QC/map/per_tissue/",tissue_label_change(tissue),"_map_age_dotplot.png"),p1,width = 15,height = 10,type="cairo")
  
  p2 <- ggplot(search_table[which(search_table$tissue==tissue_label_change(tissue)),],aes(x=antibody,y=dup*100,color = age))+    
    geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.7)+
    geom_text(aes(label = mouse_ID), position = position_jitter(width = 0.2), vjust = -1, size = 3) +
    scale_fill_brewer(palette="Set3")+
    geom_hline(yintercept = 25, color = "red", linetype = "dashed") + 
    geom_hline(yintercept = 15, color = "blue", linetype = "dashed") + 
    annotate("text", x = Inf, y = 25, label = "25%", color = "red", vjust = -0.5, hjust = 1.1, size = 5) +   
    annotate("text", x = Inf, y = 15, label = "15%", color = "blue", vjust = -0.5, hjust = 1.1, size = 5) + 
    ggtitle(paste0(tissue_label_change(tissue)," Duplication"))+ylim(0,100)+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")+labs(fill = "", color = "") +ylab("Duplication(%)")
  ggsave(paste0("result/all/QC/dup/per_tissue/",tissue_label_change(tissue),"_dup_age_dotplot.png"),p2,width = 15,height = 10,type="cairo")
}

dup_map <- read.csv("data/samples/all/map_dup_qc_RNA.csv")
search_table_RNA <- read.csv("data/samples/all/RNA_search_table.csv")
colnames(dup_map)[1] <- "sample_name"
search_table_RNA <- merge(search_table_RNA,dup_map[,c(1,3:4)],by="sample_name")
search_table_RNA$age <- factor(search_table_RNA$age, levels=c("3m","24m"))
p1 <- ggplot(search_table_RNA,aes(x=tissue,y=map*100,color = age))+    
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.7)+
  geom_text(aes(label = mouse_ID), position = position_jitter(width = 0.2), vjust = -1, size = 3) +
  scale_fill_brewer(palette="Set3")+
  geom_hline(yintercept = 80, color = "red", linetype = "dashed") + 
  geom_hline(yintercept = 90, color = "blue", linetype = "dashed") + 
  annotate("text", x = Inf, y = 75, label = "80%", color = "red", vjust = -0.5, hjust = 1.1, size = 5) +   
  annotate("text", x = Inf, y = 85, label = "90%", color = "blue", vjust = -0.5, hjust = 1.1, size = 5) + 
  ggtitle(paste0("RNA Mapping"))+ylim(0,100)+
  theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")+labs(fill = "", color = "") +ylab("Mapping(%)")
ggsave(paste0("result/all/QC/map/per_antibody/RNA_map_age_dotplot.png"),p1,width = 15,height = 10,type="cairo")
p2 <- ggplot(search_table_RNA,aes(x=tissue,y=dup*100,color = age))+    
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.7)+
  geom_text(aes(label = mouse_ID), position = position_jitter(width = 0.2), vjust = -1, size = 3) +
  scale_fill_brewer(palette="Set3")+
  geom_hline(yintercept = 25, color = "red", linetype = "dashed") + 
  geom_hline(yintercept = 15, color = "blue", linetype = "dashed") + 
  annotate("text", x = Inf, y = 25, label = "25%", color = "red", vjust = -0.5, hjust = 1.1, size = 5) +   
  annotate("text", x = Inf, y = 15, label = "15%", color = "blue", vjust = -0.5, hjust = 1.1, size = 5) + 
  ggtitle(paste0("RNA Duplication"))+ylim(0,100)+
  theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")+labs(fill = "", color = "") +ylab("Duplication(%)")
ggsave(paste0("result/all/QC/dup/per_antibody/RNA_dup_age_dotplot.png"),p2,width = 15,height = 10,type="cairo")
