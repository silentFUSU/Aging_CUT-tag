rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(MASS)
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
get_density_weight <- function(x, y, density) {  
  ix <- findInterval(x, density$x)  
  iy <- findInterval(y, density$y)  
  return(density$z[cbind(ix, iy)])  
}  
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
antibody <- "H3K27me3"
p_list <- list()
tissues <- sort(c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
                  "thymus","skin","bladder","bonemarrow","Hip","heart",
                  "muscle","jejunum","uterus","ovary","liver","tongue",
                  "cecum","colon","testis","stomach","pancreas","iWAT","ileum"))
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_gene_TSS_diff_after_remove_batch_effect.csv"))
  df <- df[which(df$Significant != "Stable"),]
  RNA <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene.csv"))
  RNA <- RNA[which(RNA$Significant != "Stable"),]
  colnames(RNA)[1] <- "Geneid"
  df <- merge(df[,c("Geneid","LogFC.old.young")],RNA[,c("Geneid","logFC")],by="Geneid")
  if(nrow(df) >=100){
    colnames(df)[2:3] <- c("histone","RNA")
    model <- lm(histone ~ 0+RNA, data = df)
    slope <- coef(model)[1]  
    
    density <- kde2d(df$RNA, df$histone, n = 50)  
    weights <- mapply(get_density_weight, df$RNA, df$histone, MoreArgs = list(density))  
    model_weighted <- lm(histone ~ 0+RNA, data = df, weights = weights)  
    slope_weighted <- coef(model_weighted)[1]  
    fit_data <- data.frame(RNA = range(df$RNA))  
    fit_data$histone <- slope_weighted * fit_data$RNA  
    
    p_list[[tissue]] <- ggplot(df, aes(x = RNA, y = histone)) +  
      geom_point(aes(size = weights), alpha = 0.5) +  
      geom_smooth(method = "lm", col = "red", formula = y ~ 0+x, se = F) +
      geom_line(data = fit_data, aes(x = RNA, y = histone), color = "blue") +
      theme_minimal() +  
      labs(title = paste0(tissue_label_change(tissue)," H3K27me3 relationship with RNA"),  
           x = "RNA log2(Fold Change)",  
           y = "Histone log2(Fold Change)") 
    ggsave(paste0("result/RNA/histone_relationship_with_RNA/H3K27me3/",tissue,"_H3K27me3_relationship_with_RNA_gene_TSS_diff_relationship_with_gene_dot_plot_slope.png"),p_list[[tissue]],width = 6,height = 5,type="cairo")
  }
}
combined_plot <- plot_a_list(p_list,no_of_rows = 4,no_of_cols = 4)
ggsave("result/RNA/histone_relationship_with_RNA/H3K27me3/H3K27me3_relationship_with_RNA_gene_TSS_diff_relationship_with_gene_dot_plot_slope.png",width = 10,height = 10,type="cairo")
