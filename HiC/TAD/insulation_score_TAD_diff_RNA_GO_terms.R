rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(stringr)
library(ggrepel)
library(grid)
library(clusterProfiler)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(tidyverse)
library(data.table)
library(GO.db)
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
GO_database <- 'org.Mm.eg.db'
genes <- read.table("~/ref_data/for_normal_mapping/TSS/refBed/mm10_refGene.bed")
genes <- as.data.table(genes[,c(1:3,5)])
setDT(genes)
setkey(genes,V1,V2,V3)

tissue <- "lung"
resolution <- "20000"
TAD_gene_GO_terms <- function(tissue,resolution){
  TAD <- read.csv(paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",tissue,"_redundant_",resolution,"_TAD_diff.csv"))
  colnames(TAD)[1] <- "Geneid"
  TAD <- TAD %>%
    separate(Geneid, into = c("Chr", "Start", "End"), sep = "-")
  TAD$Start <- as.numeric(TAD$Start)
  TAD$End <- as.numeric(TAD$End)
  if(nrow(TAD[which(TAD$Significant=="Up"),])>20){
    increase <- as.data.table(TAD[which(TAD$Significant=="Up"),c("Chr","Start","End")])
    setDT(increase)
    setkey(increase,Chr,Start,End)
    overlaps <- foverlaps(genes, increase, type = "any", nomatch = 0L)  
    genelist_up <- bitr(overlaps$V5,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
    genelist_up_GO <- enrichGO( genelist_up$ENTREZID,
                                OrgDb = GO_database,
                                keyType = "ENTREZID",
                                ont = "BP",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05,
                                readable = T)
    p <- barplot(genelist_up_GO,label_format = 50,title = paste0(tissue_label_change(tissue)," genes within increased TAD GO terms"))
    ggsave(paste0("result/HiC/",tissue,"/differential_analysis/increased_",resolution,"_TAD_genes_GO_terms.png"),p,width=8,height = 4,type="cairo")
    result <- genelist_up_GO@result
    write.csv(result,paste0("result/HiC/",tissue,"/differential_analysis/increased_",resolution,"_TAD_genes_GO_terms.csv"))
  }
  if(nrow(TAD[which(TAD$Significant=="Down"),]) > 20){
    decrease <- as.data.table(TAD[which(TAD$Significant=="Down"),c("Chr","Start","End")])
    setDT(decrease)
    setkey(decrease,Chr,Start,End)
    overlaps <- foverlaps(genes, decrease, type = "any", nomatch = 0L)  
    genelist_down <- bitr(overlaps$V5,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
    genelist_down_GO <- enrichGO( genelist_down$ENTREZID,
                                  OrgDb = GO_database,
                                  keyType = "ENTREZID",
                                  ont = "BP",
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 0.05,
                                  readable = T)
    if(nrow(genelist_down_GO) > 10){
      p <- barplot(genelist_down_GO,label_format = 50,title = paste0(tissue_label_change(tissue)," genes within decreased TAD GO terms"))
      ggsave(paste0("result/HiC/",tissue,"/differential_analysis/decreased_",resolution,"_TAD_genes_GO_terms.png"),p,width=8,height = 4,type="cairo")
      result <- genelist_down_GO@result
      write.csv(result,paste0("result/HiC/",tissue,"/differential_analysis/decreased_",resolution,"_TAD_genes_GO_terms.csv"))
    }
  }
}
tissues <- c("brain","CB","kidney","liver","lung","bonemarrow","colon","heart","Hip","mammarygland","stomach","thymus")
for(tissue in tissues){
  TAD_gene_GO_terms(tissue,resolution)
}

increase_GO_summary <- data.frame()
decrease_GO_summary <- data.frame()
for(tissue in tissues){
  if(file.exists(paste0("result/HiC/",tissue,"/differential_analysis/increased_",resolution,"_TAD_genes_GO_terms.csv"))){
    df <- read.csv(paste0("result/HiC/",tissue,"/differential_analysis/increased_",resolution,"_TAD_genes_GO_terms.csv"))
    df <- df[which(df$p.adjust < 0.05),c(2,3)]
    df$tissue <- tissue_label_change(tissue)
    increase_GO_summary <- rbind(increase_GO_summary,df)
  }
  if(file.exists(paste0("result/HiC/",tissue,"/differential_analysis/decreased_",resolution,"_TAD_genes_GO_terms.csv"))){
    df <- read.csv(paste0("result/HiC/",tissue,"/differential_analysis/decreased_",resolution,"_TAD_genes_GO_terms.csv"))
    df <- df[which(df$p.adjust < 0.05),c(2,3)]
    df$tissue <- tissue_label_change(tissue)
    decrease_GO_summary <- rbind(decrease_GO_summary,df)
  }
}
increase_GO_summary_count <- increase_GO_summary %>%   
  count(ID)
increase_GO_summary_tissue <- increase_GO_summary %>%   
  group_by(ID) %>%   
  summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
increase_GO_summary_count <- merge(increase_GO_summary_count,increase_GO_summary_tissue,by="ID")
increase_GO_summary_count$Description <- sapply(increase_GO_summary_count$ID, function(go_id) {  
  Term(GOTERM[[go_id]])  
})  
to_plot <- as.data.frame(table(increase_GO_summary_count$n))
ggplot(data = to_plot, aes(x = Var1, y = Freq)) +  
  geom_bar(stat = "identity") +  
  labs(  
    title = "Distribution of tissues number in common changed GO terms",  
    x = NULL,  
    y = "Count"  
  ) +  
  geom_text(  
    aes(label = Freq),   
    vjust = 0,   
    size = 3.5  # Adjust text size as needed  
  ) +
  theme_minimal() +  
  xlab(NULL) +  
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
    axis.text.y = element_text(size = 12),  
    plot.title = element_text(size = 16, face = "bold")  
  ) 


decrease_GO_summary_count <- decrease_GO_summary %>%   
  count(ID)
decrease_GO_summary_tissue <- decrease_GO_summary %>%   
  group_by(ID) %>%   
  summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
decrease_GO_summary_count <- merge(decrease_GO_summary_count,decrease_GO_summary_tissue,by="ID")
decrease_GO_summary_count$Description <- sapply(decrease_GO_summary_count$ID, function(go_id) {  
  Term(GOTERM[[go_id]])  
})  
to_plot <- as.data.frame(table(decrease_GO_summary_count$n))
ggplot(data = to_plot, aes(x = Var1, y = Freq)) +  
  geom_bar(stat = "identity") +  
  labs(  
    title = "Distribution of tissues number in common changed GO terms",  
    x = NULL,  
    y = "Count"  
  ) +  
  geom_text(  
    aes(label = Freq),   
    vjust = 0,   
    size = 3.5  # Adjust text size as needed  
  ) +
  theme_minimal() +  
  xlab(NULL) +  
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
    axis.text.y = element_text(size = 12),  
    plot.title = element_text(size = 16, face = "bold")  
  ) 

