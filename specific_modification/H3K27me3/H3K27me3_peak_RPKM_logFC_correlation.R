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
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
                  "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT","ileum")) 
tissue_summary <- data.frame()

for(tissue in tissues){
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  
  tab <- read.table(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_young_old_merge-W5000-G10000-E100.counts"),header = T)
  peak_regions <- as.data.table(tab[,c(1:4)])
  setDT(peak_regions)
  setkey(peak_regions,Chr,Start,End)
  blacklist <- read.table("~/ref_data/mm10-blacklist.v2.bed",sep = "\t")
  blacklist <- as.data.table(blacklist)
  setDT(blacklist)
  setkey(blacklist,V1,V2,V3)
  overlaps <- foverlaps(peak_regions, blacklist, type = "any", nomatch = 0L)
  tab <- tab[which(!tab$Geneid %in% overlaps$Geneid),]
  
  tab_summary <- read.table(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_young_old_merge-W5000-G10000-E100.counts.summary"),header = T)
  
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+|DYQ[0-9]+|HM[0-9]+).*"
  colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
  colnames(tab_summary) <- gsub(pattern,"\\1",colnames(tab_summary))
  
  search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
  counts <- counts[,search_table$sample_name]
  tab_summary <- tab_summary[-2,search_table$sample_name]
  total_reads <- colSums(tab_summary)
  length <- as.numeric(tab$Length)
  rpkm <- sweep(counts,2,total_reads,"/")
  rpkm <- sweep(rpkm,1,length,"/") * 1000000000
  rpkm$mean_all <- rowMeans(rpkm)
  rpkm_young <- rpkm[,search_table$sample_name[which(search_table$age=="3m")]]
  rpkm_old <- rpkm[,search_table$sample_name[which(search_table$age=="24m")]]
  rpkm_young$mean_young <- rowMeans(rpkm_young)
  rpkm_old$mean_old <- rowMeans(rpkm_old)
  rpkm_mean_summary <- merge(rpkm_young[,"mean_young",drop=F],rpkm_old[,"mean_old",drop=F],by="row.names")
  rpkm_mean_summary <- merge(rpkm_mean_summary,rpkm[,"mean_all",drop=F],by.x="Row.names",by.y="row.names")
  
  rpkm_mean_summary$log2FC <- log2(rpkm_mean_summary$mean_old/rpkm_mean_summary$mean_young)
  colnames(rpkm_mean_summary)[1] <- "Geneid"
  rpkm_mean_summary$tissue <- tissue_label_change(tissue)
  
  tissue_summary <- rbind(tissue_summary,rpkm_mean_summary)
  }

p_list <- list()
for(tissue in sort(tissues)){
  to_plot <- tissue_summary[which(tissue_summary$tissue==tissue_label_change(tissue)),]
  to_plot$condition <- NA
  to_plot$condition[which(to_plot$log2FC < 0)] <- "Down"
  to_plot$condition[which(to_plot$log2FC > 0)] <- "Up"
  to_plot$condition <- factor(to_plot$condition,levels=c("Up","Down"))
  x_range <- range(log2(to_plot$mean_young), na.rm = TRUE)  
  y_range <- range(to_plot$log2FC, na.rm = TRUE)  
  x_pos_right <- x_range[2] * 0.9    
  x_pos_left <- x_range[1] * 0.9   
  y_pos_top <- y_range[2] * 0.9    
  y_pos_bottom <- y_range[1] * 0.9 
  ### log2 RPKM log2(O/Y)
  p_list[[tissue]] <- ggplot(to_plot,aes(x=log2(mean_young),y=log2FC,color=condition))+    
    geom_jitter(size = 3, alpha = 0.7)+
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # 添加水平线
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # 添加竖直线
    geom_hline(yintercept = 1, color = "red") +  # 添加水平线
    geom_hline(yintercept = -1, color = "red") +
    theme_bw()+theme(text = element_text(size = 18))+
    xlab("log2(RPKM)")+
    ylab("log2(O/Y)")+
    ggtitle(tissue_label_change(tissue))+
    labs(fill = "", color = "")+
    annotate("text", label = paste0(nrow(to_plot[which(log2(to_plot$mean_young)>0 & to_plot$log2FC>0),])),  
             x = x_pos_right, y = y_pos_top, colour = "#00b8a9", size = 5) +  
    annotate("text", label = paste0(nrow(to_plot[which(log2(to_plot$mean_young)<0 & to_plot$log2FC<0),])),  
             x = x_pos_left, y = y_pos_bottom, colour = "#ff9a00", size = 5) +  
    annotate("text", label = paste0(nrow(to_plot[which(log2(to_plot$mean_young)<0 & to_plot$log2FC>0),])),  
             x = x_pos_left, y = y_pos_top, colour = "#f6416c", size = 5) +  
    annotate("text", label = paste0(nrow(to_plot[which(log2(to_plot$mean_young)>0 & to_plot$log2FC<0),])),  
             x = x_pos_right, y = y_pos_bottom, colour = "#48466d", size = 5) 
  
}
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
combined_plot <- plot_a_list(p_list, 4, 7)
ggsave("tmp.png",combined_plot,width = 35,height = 20,type="cairo")


proportion_summary <- data.frame()
for(tissue in sort(tissues)){
  df <- tissue_summary[which(tissue_summary$tissue==tissue_label_change(tissue)),]  
  df$quadrant <- "other"
  df$quadrant[which(log2(df$mean_young)>0 & df$log2FC>0)] <- "first"
  df <- as.data.frame(table(df$quadrant))
  df$percentage <- df$Freq/sum(df$Freq)*100
  df$tissue <- tissue_label_change(tissue)
  proportion_summary <- rbind(proportion_summary,df)
}
rank <- proportion_summary[which(proportion_summary$Var1=="first"),]
rank <- rank[order(rank$percentage),]
proportion_summary$tissue <- factor(proportion_summary$tissue,levels=rank$tissue) 

ggplot(proportion_summary, aes(x = tissue, y = percentage, fill = Var1)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Tissue", y = "Percentage", title = "Stacked Bar Plot of Tissue Proportions") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

rank <- rank[order(rank$Freq),]
proportion_summary$tissue <- factor(proportion_summary$tissue,levels=rank$tissue) 

ggplot(proportion_summary, aes(x = tissue, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Tissue", y = "# peaks", title = "Stacked Bar Plot of Tissue Proportions") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


to_plot <- tissue_summary[which(tissue_summary$tissue=="Cerebellum"),]
to_plot$condition <- NA
to_plot$condition[which(to_plot$log2FC < 0)] <- "Down"
to_plot$condition[which(to_plot$log2FC > 0)] <- "Up"
to_plot$condition <- factor(to_plot$condition,levels=c("Up","Down"))
ggplot(to_plot, aes(x = condition, y = mean_young,fill=condition)) +
  geom_boxplot(outliers = F) +
  # scale_fill_manual(values = color) +
  labs(x = NULL, y = "RPKM", title = "H3K27me3 signal in diff peak condition") +
  theme_bw()

x_range <- range(log2(to_plot$mean_young), na.rm = TRUE)  
y_range <- range(to_plot$log2FC, na.rm = TRUE)  
x_pos_right <- x_range[2] * 0.9    
x_pos_left <- x_range[1] * 0.9   
y_pos_top <- y_range[2] * 0.9    
y_pos_bottom <- y_range[1] * 0.9 
### log2 RPKM log2(O/Y)
ggplot(to_plot,aes(x=log2(mean_young),y=log2FC,color=condition))+    
  geom_jitter(size = 3, alpha = 0.7)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # 添加水平线
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # 添加竖直线
  geom_hline(yintercept = 1, color = "red") +  # 添加水平线
  geom_hline(yintercept = -1, color = "red") +
  theme_bw()+theme(text = element_text(size = 18))+
  xlab("log2(RPKM)")+
  ylab("log2(O/Y)")+
  ggtitle(tissue_label_change(tissue))+
  labs(fill = "", color = "")+
  annotate("text", label = paste0(nrow(to_plot[which(log2(to_plot$mean_young)>0 & to_plot$log2FC>0),])),  
           x = x_pos_right, y = y_pos_top, colour = "#00b8a9", size = 5) +  
  annotate("text", label = paste0(nrow(to_plot[which(log2(to_plot$mean_young)<0 & to_plot$log2FC<0),])),  
           x = x_pos_left, y = y_pos_bottom, colour = "#ff9a00", size = 5) +  
  annotate("text", label = paste0(nrow(to_plot[which(log2(to_plot$mean_young)<0 & to_plot$log2FC>0),])),  
           x = x_pos_left, y = y_pos_top, colour = "#f6416c", size = 5) +  
  annotate("text", label = paste0(nrow(to_plot[which(log2(to_plot$mean_young)>0 & to_plot$log2FC<0),])),  
           x = x_pos_right, y = y_pos_bottom, colour = "#48466d", size = 5) 

to_plot <- to_plot %>%
  separate(Geneid, into = c("chr", "range"), sep = ":") %>% 
  separate(range, into = c("start", "end"), sep = "-") %>%
  mutate(start = as.numeric(start),
         end = as.numeric(end),
         peak_length = end - start)
to_plot$quadrant <- "first"
to_plot$quadrant[which(log2(to_plot$mean_young)<0 & to_plot$log2FC>0)] <- "second"
to_plot$quadrant[which(log2(to_plot$mean_young)<0 & to_plot$log2FC<0)] <- "third"
to_plot$quadrant[which(log2(to_plot$mean_young)>0 & to_plot$log2FC<0)] <- "fourth"
to_plot$quadrant <- factor(to_plot$quadrant,levels=c("first","second","third","fourth"))
ggplot(to_plot, aes(x = quadrant, y = peak_length,fill=condition)) +
  geom_boxplot(outliers = F) +
  # scale_fill_manual(values = color) +
  labs(x = NULL, y = "length") +
  theme_bw()

df <- to_plot[which(log2(to_plot$mean_young)>0 & to_plot$log2FC>0),]


cortest <- cor.test(to_plot$mean_young,to_plot$log2FC)

ggplot(to_plot,aes(x=log2(mean_young),y=log2FC,color=condition))+    
  geom_jitter(size = 3, alpha = 0.7)+
  theme_bw()+theme(text = element_text(size = 18))+
  xlab("log2(RPKM)")+
  ylab("logFC")+
  labs(fill = "", color = "")

smoothScatter(to_plot$log2FC ~ log2(to_plot$mean_young),
              bandwidth = 0.05,xlab = "log2(RPKM)",ylab="log2(FC)")
abline(h = 0, col = "red", lwd = 2, lty = 2)

color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(unique(to_plot$tissue)))
ggplot(to_plot,aes(x=log2(mean_young),y=log2FC,color=tissue))+    
  geom_jitter(size = 3, alpha = 0.5) +
  scale_color_manual(values = color) +
  theme_bw()+theme(text = element_text(size = 18))+
  xlab("log2(RPKM)")+
  ylab("logFC")+
  labs(fill = "", color = "")

###
tissue_proportion <- tissue_summary
tissue_proportion$quadrant <- "first"
tissue_proportion$quadrant[which(log2(tissue_proportion$mean_young) < 0 & tissue_proportion$log2FC > 0)] <- "second"
tissue_proportion$quadrant[which(log2(tissue_proportion$mean_young) < 0 & tissue_proportion$log2FC < 0)] <- "third"
tissue_proportion$quadrant[which(log2(tissue_proportion$mean_young) > 0 & tissue_proportion$log2FC < 0)] <- "fourth"

tissue_distribution <- tissue_proportion %>%
  group_by(quadrant, tissue) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count))
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(unique(to_plot$tissue)))
ggplot(tissue_distribution, aes(x = quadrant, y = proportion, fill = tissue)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color) +
  labs(x = "Quadrant", y = "Proportion", fill = "Tissue") +
  theme_minimal() +
  ggtitle("Proportion of Tissue Elements by Quadrant")


for(quadrant in c("first","second","third","fourth")){
  t_tissue_distribution <- tissue_distribution[which(tissue_distribution$quadrant==quadrant),]
  print(var(t_tissue_distribution$proportion))
  }

