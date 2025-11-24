rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/","/usr/local/lib64/R/library"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(patchwork)
library(stringr)
library(dplyr)
library(tidyr)
library(UpSetR)
library(ChIPseeker)
library(data.table)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
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
      tissue_label <- "IWAT"
    }
  }
  return(tissue_label)
}

tissues <- sort(c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
                  "thymus","skin","bladder","bonemarrow","Hip","heart",
                  "muscle","jejunum","uterus","ovary","liver","tongue",
                  "cecum","colon","testis","stomach","pancreas","iWAT","ileum"))
genes<-data.frame(X = character(),  
                  Significant = character(),  
                  tissue = character(),
                  stringsAsFactors = FALSE)  

for (i in c(1:length(tissues))){
  tissue <- tissues[i]
  df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv")) 
  df <- df[which(df$Significant!="Stable"),c("X","Significant")]
  if(nrow(df) >0){
    df$tissue <- tissue
    genes <- rbind(genes,df)
  }
}
colnames(genes)[1] <-"Geneid"
increase <- genes[which(genes$Significant=="Up"),]
increase_count <- increase %>%   
  count(Geneid)
increase_tissue <- increase %>%   
  group_by(Geneid) %>%   
  summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
increase_count <- merge(increase_count,increase_tissue,by="Geneid")
# increase_count <- merge(increase_count,bin_file,by="Geneid")
write.csv(increase_count[order(increase_count$n,decreasing = T),],"data/samples/all/RNA/gene_expression_increased_common_genes.csv")


colnames(genes)[1] <-"Geneid"
decrease <- genes[which(genes$Significant=="Down"),]
decrease_count <- decrease %>%   
  count(Geneid)
decrease_tissue <- decrease %>%   
  group_by(Geneid) %>%   
  summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
decrease_count <- merge(decrease_count,decrease_tissue,by="Geneid")
write.csv(decrease_count[order(decrease_count$n,decreasing = T),],"data/samples/all/RNA/gene_expression_decreased_common_genes.csv")

PMD_HMD_region <- read.table("data/public_data/PMD_coordinates_mm10.bed")
PMD_HMD_region$V2 <- PMD_HMD_region$V2 + 1
PMD_HMD_region$label <- paste0("bin",c(1:nrow(PMD_HMD_region)))
PMD_HMD_region <- PMD_HMD_region[,c("V1","V2","V3","V5","label")]
PMD_HMD_region$V5[is.na(PMD_HMD_region$V5)] <- "other"
PMD_HMD_region <- as.data.table(PMD_HMD_region)
setDT(PMD_HMD_region)
setkey(PMD_HMD_region,V1,V2,V3)

gene_TSS <- read.table("~/ref_data/for_normal_mapping/mm10/mm10_refseq_TSS.bed")
result <- increase_count
result <- as.data.table(merge(result,gene_TSS[,c(1:3,6)],by.x="Geneid",by.y="V6"))
setDT(result)
setkey(result,V1,V2,V3)
overlaps <- as.data.frame(foverlaps(result, PMD_HMD_region, type = "any", nomatch = 0L))
overlaps$condition <- "larger_10"
overlaps$condition[which(overlaps$n < 10 & overlaps$n >=5)] <- "larger_5"
overlaps$condition[which(overlaps$n < 5)] <- "smaller_5"
to_plot <- overlaps %>%
  group_by(condition, V5) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(condition) %>%
  mutate(total_count = sum(count),
         proportion = count / total_count * 100)
to_plot <- as.data.frame(to_plot)
ggplot(to_plot, aes(x = condition, y = proportion, fill = V5)) +  
  geom_bar(stat = 'identity',color="white") +   
  theme_minimal() +   
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Proportion")+
  ggtitle("common decreased genes")+
  geom_text(aes(label = total_count), 
            y = 100, 
            size = 5, 
            vjust = 0, 
            color = "black")
