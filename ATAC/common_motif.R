rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(bitmapType="cairo")  
conditions <- c("up","down")
tissues <-  c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
                        "thymus","skin","bladder","bonemarrow","Hip","heart",
                        "muscle","jejunum","uterus","ovary","liver","tongue",
                        "cecum","colon","testis","stomach","pancreas","iWAT","ileum")
extract_before_bracket <- function(s) {  
  parts <- strsplit(s, "\\(")[[1]]  
  return(parts[1])  
}  
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

motif_count_summary <- data.frame()
motif_summary <- list(up=data.frame(),down=data.frame())
motif_filter_list <- list(up=vector(),down=vector())
top <- 10
motif_data_frame <- list(up=data.frame(),down=data.frame())
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(conditions))){
    condition <- conditions[j]
    if(file.exists(paste0("data/samples/ATAC/ATAC_peak_from_MJ/motif_bg/",condition,"/",tissue,"/knownResults.txt"))){
      motif <- read.delim(paste0("data/samples/ATAC/ATAC_peak_from_MJ/motif_bg/",condition,"/",tissue,"/knownResults.txt"))
      motif$Motif.Name <-  sapply(motif$Motif.Name, extract_before_bracket)  
      motif$Motif.Name <- paste0(motif$Motif.Name,"-",motif$Consensus)
      motif <- motif[!duplicated(motif$Motif.Name),]
      colnames(motif)[which(colnames(motif)=="q.value..Benjamini.")] <- tissue_label_change(tissue)
      if(i == 1){
        motif_data_frame[[condition]] <- motif[,c("Motif.Name",tissue_label_change(tissue))]
      }else{
        motif_data_frame[[condition]] <- merge(motif_data_frame[[condition]],motif[,c("Motif.Name",tissue_label_change(tissue))],by="Motif.Name",all=T)
      }
      
      motif <- motif[which(motif[,5]<0.05),]
      motif_count <- data.frame(count=nrow(motif),
                                tissue=tissue_label_change(tissue),
                                condition=condition)
      motif_count_summary <- rbind(motif_count_summary,motif_count)
      if(nrow(motif) > 0){
        motif <- motif[order(motif[,5]),]
        if(i == 1){
          motif_filter_list[[condition]] <- motif$Motif.Name[1:(max(nrow(motif_filter_list),top))]
        }else{
          motif_filter_list[[condition]] <- union(motif_filter_list[[condition]],motif$Motif.Name[1:(max(nrow(motif_filter_list),top))])
        }
        motif <- motif[,"Motif.Name",drop=F]
        motif <- motif[!duplicated(motif$Motif.Name),,drop=F]
        motif$tissue <- tissue_label_change(tissue)
        motif_summary[[condition]] <- rbind(motif_summary[[condition]],motif)    
      }
    }
  }
}
increase_count <-motif_summary[["up"]] %>%   
  count(Motif.Name)
increase_tissue <- motif_summary[["up"]] %>%   
  group_by(Motif.Name) %>%   
  summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
increase_count <- merge(increase_count,increase_tissue,by="Motif.Name")

decrease_count <-motif_summary[["down"]] %>%   
  count(Motif.Name)
decrease_tissue <- motif_summary[["down"]] %>%   
  group_by(Motif.Name) %>%   
  summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
decrease_count <- merge(decrease_count,decrease_tissue,by="Motif.Name")

increase_count <- increase_count[order(increase_count$n,decreasing = TRUE),]
decrease_count <- decrease_count[order(decrease_count$n,decreasing = TRUE),]
# write.csv(increase_count,"data/samples/all/ATAC/motif_bg/up/all_tissues_ATAC_peaks_increase_motif_count.csv",row.names = F)
# write.csv(decrease_count,"data/samples/all/ATAC/motif_bg/down/all_tissues_ATAC_peaks_decrease_motif_count.csv",row.names = F)

motif_count_summary$condition <- factor(motif_count_summary$condition,levels=c("up","down"))
motif_count_summary$position <- motif_count_summary$count
motif_count_summary$position[which(motif_count_summary$condition=="down")] <- (-motif_count_summary$position[which(motif_count_summary$condition=="down")])
ggplot(motif_count_summary, aes(x = tissue, y = ifelse(condition == "up", count, -count), fill = condition)) +  
  geom_bar(stat = "identity") +  
  labs(title = paste0("fdr < 0.05 motif count (using oppositely changed peaks as background)"), x = NULL, y = "Count") +  
  theme_minimal() +  
  xlab(NULL)+
  scale_y_continuous(labels = abs) +  
  theme(  
    axis.title.x = element_text(size = 14),     
    axis.title.y = element_text(size = 14),    
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),   
    axis.text.y = element_text(size = 12),    
    plot.title = element_text(size = 16, face = "bold")
  ) +  
  geom_text(data =motif_count_summary,   
            aes(label = count, y = position),   
            color = "black", size = 5, vjust = 0.5) + 
  scale_fill_manual(values = c("up" = "skyblue", "down" = "salmon"), name = NULL)  

increase_motif_data_frame <- motif_data_frame[["up"]][which(motif_data_frame[["up"]]$Motif.Name %in% motif_filter_list[["up"]]),]
decrease_motif_data_frame <- motif_data_frame[["down"]][which(motif_data_frame[["down"]]$Motif.Name %in% motif_filter_list[["down"]]),]

rownames_increase_motif_data_frame <- increase_motif_data_frame$Motif.Name
rownames_decrease_motif_data_frame <- decrease_motif_data_frame$Motif.Name

increase_motif_data_frame <- as.data.frame(lapply(increase_motif_data_frame[,-1], function(x) -log10(x)))
rownames(increase_motif_data_frame) <- rownames_increase_motif_data_frame

decrease_motif_data_frame <- as.data.frame(lapply(decrease_motif_data_frame[,-1], function(x) -log10(x)))
rownames(decrease_motif_data_frame) <- rownames_decrease_motif_data_frame
increase_motif_data_frame[] <- lapply(increase_motif_data_frame, function(x) replace(x, is.infinite(x), 4))  
decrease_motif_data_frame[] <- lapply(decrease_motif_data_frame, function(x) replace(x, is.infinite(x), 4))  

breaks <- seq(-log10(0.05), 4, length.out = 101)

pheatmap::pheatmap(increase_motif_data_frame,main = "Motif enriched in increase peaks",filename = "result/all/ATAC/all_tissues_increase_peaks_motif.png",width = 10,height = 20,breaks = breaks,na_col = "grey")
pheatmap::pheatmap(decrease_motif_data_frame,main = "Motif enriched in decrease peaks",filename = "result/all/ATAC/all_tissues_decrease_peaks_motif.png",width = 10,height = 20,breaks = breaks,na_col = "grey")

