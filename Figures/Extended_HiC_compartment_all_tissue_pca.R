rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
library(ggplot2)
set.seed(1)
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
    }else if(tissue_label=="Mammarygland"){
      tissue_label <- "Mammary Gland"
    }
  }
  return(tissue_label)
}

tissues <- c("brain","CB","kidney","liver","lung","bonemarrow","colon","heart","Hip","mammarygland","stomach","thymus","skin","muscle","cecum")
df_pca <- data.frame()
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")  
  search_table <- search_table[which(search_table$tissue==tissue),]
  for(i in c(1:length(search_table$sample_name))){
    sample <- search_table$sample_name[i]
    df <- read.table(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_",resolution,".PC1.txt"))
    df <- df[which(df$V2 %in% c(paste0("chr",c(1:19,"X","Y")))),]
    df <- df[,c(1,6)]
    colnames(df)[2] <- sample
    if(nrow(df_pca)==0){
      df_pca <- df
    }else{
      df_pca <- merge(df_pca,df,by="V1")
    }
  }
}
rownames(df_pca) <- df_pca[,1]
df_pca <- df_pca[,-1]
pca <- prcomp(t(df_pca))
to_plot <- data.frame(pca$x)
to_plot$sample_name <- rownames(to_plot)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
search_table <- read.csv("data/samples/all/HiC_search_table.csv")  
to_plot <- merge(to_plot,search_table,by="sample_name")
to_plot$sample_name <- paste0(to_plot$sample_name,"-",to_plot$mouse_ID,"-",to_plot$age)
to_plot$age <- factor(to_plot$age, levels = c("3M","24M"))
to_plot$tissue <- sapply(to_plot$tissue, tissue_label_change)

tissues <- c("mammarygland","lung","liver","kidney","ileum","Hip","skin","bonemarrow","jejunum","colon","ovary","CB","BAT","thymus","testis","heart","stomach","muscle","bladder","aorta","tongue","spleen","pancreas","brain","cecum","uterus","iWAT")
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(sapply(tissues, tissue_label_change)))
color <- color[which(names(color) %in% to_plot$tissue)]
p<- ggplot(to_plot, aes(x=PC1, y=PC2, color=tissue, shape=age)) + 
  geom_point(size=5) +theme_bw()+
  scale_color_manual(values = color) +
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+ 
  ggtitle("HiC")
ggsave("result/Sup_figures/HiC_all_tissues_PCA.pdf",p,width = 8,height = 6)

