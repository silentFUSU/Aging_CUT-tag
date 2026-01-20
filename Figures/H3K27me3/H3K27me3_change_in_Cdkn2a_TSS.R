rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(stringr)
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
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
      tissue_label <- "Mammary gland"
    }else if(tissue_label=="Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
}
summary <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_gene_TSS_10kb_diff_after_remove_batch_effect.csv"))
  df <- df[which(df$Geneid=="Cdkn2a"),c("Geneid","LogFC.old.young","Significant")]
  df$tissue <- tissue_label_change(tissue)
  summary <- rbind(summary,df)
  }
summary <- summary[order(summary$LogFC.old.young,decreasing = T),]
summary$tissue <- factor(summary$tissue,levels=summary$tissue)
color <- setNames(c("red","grey","blue"),c("Up","Stable","Down"))
p <- ggplot(summary,mapping = aes(x=LogFC.old.young,y=tissue,fill = Significant))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+ylab("")+
  scale_fill_manual(values=color)+
  theme(text = element_text(size = 13))+xlim(-3,0.5)

tissues_order <- c("Heart","Tongue","Lung","Skin","Uterus","Liver","Ovary","Spleen","Cortex","Bone Marrow","Bladder",
                   "Jejunum","iWAT","Stomach","Muscle","Cecum","Ileum","Thymus","Mammary gland","Testis","Kidney","Hippocampus",
                   "Cerebellum","Aorta","Pancreas","Colon","BAT")

summary$tissue <- factor(summary$tissue,levels=rev(tissues_order))
color <- setNames(c("red","grey","blue"),c("Up","Stable","Down"))
p <- ggplot(summary,mapping = aes(x=LogFC.old.young,y=tissue,fill = Significant))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+ylab("")+
  scale_fill_manual(values=color)+
  theme(text = element_text(size = 13))+xlim(-3,0.5)
p
ggsave("result/figures/Cdkn2a_H3K27me3_diff_in_each_tissue.pdf",p,width = 6,height = 6)
