rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(stringr)
library(dplyr)
library(dbplyr)
library(clusterProfiler)
library(GSVA)
library(enrichplot)
library(MASS)  
library(RANSAC)
options(scipen = 0) 
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

tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
                  "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT","ileum"))
# summary_antibody <- read.csv("data/samples/WGBS/all_tissues_delta_in_1000_change_in_high_H3K27me3_level_hmr.csv")
# summary_antibody <- read.csv("data/samples/WGBS/all_tissues_delta_in_1000_change_in_high_EZH2_SUZ12_level_hmr.csv")
summary_antibody <- read.csv("data/samples/WGBS/all_tissues_all_samples_delta_in_top1000_hmr_change_in_high_EZH2_SUZ12_peak_q01_input_level_hmr.csv")
summary_antibody <- summary_antibody[,c("tissue","delta")]

mitotic_nuclear_division <- read.csv("data/public_data/GO_term_summary_0140014.csv")
mitotic_nuclear_division <- unique(mitotic_nuclear_division$Symbol)
target_genes <- mitotic_nuclear_division
GO_id <- "GO:0140014"

description <- "mitotic nuclear division"
search_table <- read.csv("data/samples/all/RNA_search_table.csv")
genelist <- list(score=target_genes)

rpkm <- data.frame()
for(tissue in tissues){
  df <- read.table(paste0("data/samples/RNA/",tissue,"/combined-chrM.counts"),header = T)
  df <- df[,c(1,6,7:ncol(df))]
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|HM[0-9]+).*"
  colnames(df)[-c(1,2)] <- gsub(pattern, "\\1",colnames(df)[-c(1,2)])
  gene_lengths <- df$Length
  total_mapped_reads <- colSums(df[, 3:ncol(df)])
  rpkm_df <- data.frame(Geneid = df$Geneid)
  for (i in 3:ncol(df)) {
    counts <- df[[i]]
    t_rpkm <- (counts / (gene_lengths / 1000)) / (total_mapped_reads[i - 2] / 1e6)
    rpkm_df[[colnames(df)[i]]] <- t_rpkm
  }
  
  if(nrow(rpkm)==0){
    rpkm <- rpkm_df
  }else{
    rpkm <- merge(rpkm,rpkm_df,by="Geneid")
  }
}
rownames(rpkm) <- rpkm$Geneid
rpkm <- rpkm[,-1]
rpkm_matrix <- as.matrix(rpkm)

re <- gsva(rpkm_matrix,genelist , method="ssgsea",ssgsea.norm=TRUE) 
re <- as.data.frame(t(re))

to_plot <- merge(re,search_table,by.x="row.names",by.y="sample_name")
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(unique(as.character(to_plot$tissue))))

average_scores <-aggregate(score ~ tissue, data = to_plot, FUN = mean)
average_scores <- average_scores[order(average_scores$score),]
to_plot$tissue <- factor(to_plot$tissue,levels = average_scores$tissue)
to_plot$age[which(to_plot$age=="3m")] <- "Young"
to_plot$age[which(to_plot$age=="24m")] <- "Old"
to_plot$age <- factor(to_plot$age,levels=c("Young","Old"))
ggplot(to_plot,aes(x=tissue,y=score,color = tissue,shape=age))+    
  geom_jitter(position = position_jitter(width = 0.2), size = 3, alpha = 0.7)+
  scale_color_manual(values = color)+
  ggtitle(paste0(description))+
  theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("Tissues")+labs(fill = "", color = "") 

to_plot$tissue <- as.character(to_plot$tissue)
to_plot$tissue[which(to_plot$tissue=="Mammary gland")] <- "Mammary Gland"
to_plot <- merge(to_plot,summary_antibody,by="tissue")
colnames(to_plot)[ncol(to_plot)] <- "histone"
cor_test <- cor.test(to_plot$score,to_plot$histone,method="spearman")

average_scores$tissue[which(average_scores$tissue=="Mammary gland")] <- "Mammary Gland"
average_scores <- merge(average_scores,summary_antibody,by="tissue")
colnames(average_scores)[ncol(average_scores)] <- "histone"

cor_test <- cor.test(average_scores$score,average_scores$histone,method="spearman")
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(unique(as.character(to_plot$tissue))))
ggplot(to_plot,aes(x=histone,y=score,color = tissue,shape=age))+    
  geom_jitter(size = 3, alpha = 0.7)+
  scale_color_manual(values = color)+
  ggtitle(paste0(description," with DNA methylation"))+
  theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("Delta")+labs(fill = "", color = "")+
  scale_x_reverse() 

p <- ggplot(average_scores[which(!average_scores$tissue%in% c("Mammary Gland","Uterus","Ovary")),],aes(x=histone,y=score))+    
  geom_jitter(size = 3, alpha = 0.7,color="#4DBBD5")+
  geom_smooth(data = average_scores, aes(x = histone, y = score),
              method = "lm", color = "#3C5488", se = TRUE, level = 0.95) +
  # geom_abline(intercept = intercept, slope = slope, color = "#e64b35", size = 1) +
  # geom_ribbon(aes(ymin = lower, ymax = upper))+
  ggtitle(paste0(description))+
  theme_bw()+theme(text = element_text(size = 18))+
  xlab("Delta")+labs(fill = "", color = "")+
  scale_y_continuous(limits = c(2.5, 3.8), breaks = seq(2.5, 3.8, by = 0.3))+
  scale_x_reverse(limits = c(16, -3))
cor_test <- cor.test(average_scores$score[which(!average_scores$tissue%in% c("Mammary Gland","Uterus","Ovary"))],
                     average_scores$histone[which(!average_scores$tissue%in% c("Mammary Gland","Uterus","Ovary"))],method="spearman")
linear_model <- lm(score ~ histone, data = average_scores[which(!average_scores$tissue%in% c("Mammary Gland","Uterus","Ovary")),])
model_summary <- summary(linear_model)
ggsave("result/figures/WGBS_mitotic_nuclear_division_PRC2.pdf",p,width = 6,height = 4)

set.seed(1)
model <- ransac_reg(score ~ histone, data = average_scores, n_min = 10, tol = 0.05,n_iter = 10000,verbose = T)
coefficients <- coef(model)
intercept <- coefficients[1]
slope <- -coefficients[2]

predictions <- average_scores %>%
  mutate(
    fitted = intercept + slope * histone, 
    residuals = score - fitted 
  )

n<- nrow(average_scores)                                
se <- sqrt(sum(predictions$residuals^2) / (n - 2))      
predictions <- predictions %>%
  mutate(
    conf_interval = qt(0.975, df = n - 2) * se
  )
ggplot(average_scores,aes(x=histone,y=score))+    
  geom_jitter(size = 3, alpha = 0.7,color="#f39b7f")+
  geom_abline(intercept = intercept, slope = slope, color = "#e64b35", size = 1) +
  ggtitle(paste0(description))+
  # geom_ribbon(data = predictions, 
  #             aes(x = histone, 
  #                 ymin = fitted - conf_interval, 
  #                 ymax = fitted + conf_interval),
  #             fill = "#e64b35", alpha = 0.2) +
  theme_bw()+theme(text = element_text(size = 18))+
  xlab("Delta")+labs(fill = "", color = "")+
  scale_x_reverse()



