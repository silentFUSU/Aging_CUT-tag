rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)


tissue <- "kidney"
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
fdr_distribution <- function(tissue){
  H3K9me3 <- read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_10kb_bins_diff.csv"))
  H3K27me3 <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_10kb_bins_diff.csv"))
  H3K9me3_decrease <- H3K9me3$Geneid[H3K9me3$Significant=="Down"]
  H3K27me3_in_H3K9me3_decrease <- H3K27me3[which(H3K27me3$Geneid %in% H3K9me3_decrease),c("FDR.old.young","LogFC.old.young","Significant")]
  H3K27me3_in_H3K9me3_decrease$tissue <- tissue
  colnames(H3K27me3_in_H3K9me3_decrease)[3] <- "condition"
  return(H3K27me3_in_H3K9me3_decrease)
}

tissues <- c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
             "thymus","skin","bladder","bonemarrow","Hip","heart",
             "muscle","jejunum","uterus","ovary","liver","tongue",
             "cecum","colon","testis","stomach","pancreas")

for (i in c(1:length(1:length(tissues)))){
  if(i == 1){
    df <- fdr_distribution(tissues[i])
  }else{
    df <- rbind(df,fdr_distribution(tissues[i]))
  }
}
df$condition <- factor(df$condition,levels = c("Down","Stable","Up"))
ggplot(df[which(df$condition !="Stable"),], aes(x = tissue, y = -log10(FDR.old.young), fill = condition)) +
  geom_split_violin() + ylim(0,25)

tissue <- "muscle"
fdr_distribution_dot_plot <- function(tissue,data_choose){
  H3K9me3 <- read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_10kb_bins_diff.csv"))
  H3K27me3 <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_10kb_bins_diff.csv"))
  H3K9me3_decrease <- H3K9me3[H3K9me3$Significant=="Down",c("Geneid","FDR.old.young","LogFC.old.young")]
  H3K27me3_in_H3K9me3_decrease <- merge(H3K9me3_decrease,H3K27me3[,c("Geneid","FDR.old.young","LogFC.old.young","Significant")],by="Geneid")
  H3K27me3_in_H3K9me3_decrease$FDR.old.young.x.adjust <- -log10(H3K27me3_in_H3K9me3_decrease$FDR.old.young.x)
  
  H3K27me3_in_H3K9me3_decrease$FDR.old.young.y.adjust <- ifelse(H3K27me3_in_H3K9me3_decrease$LogFC.old.young.y > 0,   
                                                                -log10(H3K27me3_in_H3K9me3_decrease$FDR.old.young.y),   
                                                                log10(H3K27me3_in_H3K9me3_decrease$FDR.old.young.y))  
  H3K27me3_in_H3K9me3_decrease$Significant <- factor(H3K27me3_in_H3K9me3_decrease$Significant,levels = c("Up","Stable","Down"))
  p1 <- ggplot(H3K27me3_in_H3K9me3_decrease,aes(x = FDR.old.young.x.adjust, y = FDR.old.young.y.adjust, color = Significant))+
    geom_point() + theme_bw() + xlab("-log10(fdr) H3K9me3")+ ylab("-log10(fdr) H3K27me3")+
    ggtitle(tissue)+  
    scale_y_continuous(  
      labels = function(y) abs(y)  
    ) +
    theme(text = element_text(size = 18),legend.position = "none")
  
  p2 <- ggplot(H3K27me3_in_H3K9me3_decrease,aes(x = FDR.old.young.x.adjust, y = LogFC.old.young.y, color = Significant))+
    geom_point() + theme_bw() + xlab("-log10(fdr) H3K9me3")+ ylab("log2(Fold Change) H3K27me3")+
    ggtitle(tissue)+
    theme(text = element_text(size = 18),legend.position = "none")
  if(data_choose=="FDR"){
    return(p1)
  }else{
    return(p2)
  }
}
tissues <- c("BAT","mammarygland","CB","lung","aorta","spleen","ovary","Hip",
             "heart","bladder","thymus","bonemarrow","skin","uterus",
             "liver","cecum","brain","tongue","kidney","colon",
             "testis","muscle","jejunum","stomach","pancreas")
p_list <- list()
for (i in c(1:length(tissues))){
  p_list[[i]] <- fdr_distribution_dot_plot(tissues[i],"Foldchange")
}
combined_plot <- plot_a_list(p_list, 4, 7)
ggsave("result/all/H3K27me3_H3K9me3/all_tissue_H3K27me3_foldchange_H3K9me3_fdr.png",combined_plot,width = 35,height = 20,type="cairo")
