rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(bitmapType="cairo")  
library(grid)
library(stringr)
library(dplyr)
library(ggplot2)
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
conditions <- c("up","down")
tissues <-  c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
              "thymus","skin","bladder","bonemarrow","Hip","heart",
              "muscle","jejunum","uterus","ovary","liver","tongue",
              "cecum","colon","testis","stomach","pancreas","iWAT","ileum")
data_path <- "data/samples/ATAC/all/ATAC/snapatac2_macs/strict_stable_peaks_summits_spm3/"

for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(conditions))){
    condition <- conditions[j]
    if(file.exists(paste0(data_path,condition,"/enrichment_results_",tissue,"_summits_spm3.bed.csv"))){
      motif <- read.csv(paste0(data_path,condition,"/enrichment_results_",tissue,"_summits_spm3.bed.csv"))
      motif$adjusted.p.value[which(motif$log2.fold.change. < 0 )] <- 1
      motif$adjusted.p.value[which(motif$fg_percent == 0 )] <- 1
      motif$id <- ifelse(
        grepl("\\(.*\\)", motif$id),
        sub(".*\\((.*?)\\).*", "\\1", motif$id), 
        sub("^M\\d+_2\\.00\\s*", "", motif$id) 
      )
      
      colnames(motif)[which(colnames(motif)=="adjusted.p.value")] <- tissue_label_change(tissue)
      if(i == 1){
        motif_data_frame[[condition]] <- motif[,c("id",tissue_label_change(tissue))]
      }else{
        motif_data_frame[[condition]] <- merge(motif_data_frame[[condition]],motif[,c("id",tissue_label_change(tissue))],by="id",all=T)
      }
      motif$log2.fold.change.[which(motif$fg_percent>0 & motif$bg_percent ==0)] <- max(motif$log2.fold.change.[is.finite(motif$log2.fold.change.)],na.rm = T)
      motif <- motif[which(motif$log2.fold.change. > 0 & motif[,9] < 0.05),]
      motif_count <- data.frame(count=nrow(motif),
                                tissue=tissue_label_change(tissue),
                                condition=condition)
      motif_count_summary <- rbind(motif_count_summary,motif_count)
      if(nrow(motif) > 0){
        motif <- motif[order(motif[,9]),]
        if(i == 1){
          motif_filter_list[[condition]] <- motif$id[1:(max(nrow(motif_filter_list),top))]
        }else{
          motif_filter_list[[condition]] <- union(motif_filter_list[[condition]],motif$id[1:(max(nrow(motif_filter_list),top))])
        }
        motif <- motif[,"id",drop=F]
        motif <- motif[!duplicated(motif$id),,drop=F]
        motif$tissue <- tissue_label_change(tissue)
        motif_summary[[condition]] <- rbind(motif_summary[[condition]],motif) 
      }
    }else{
      motif_count <- data.frame(count=0,
                                tissue=tissue_label_change(tissue),
                                condition=condition)
      motif_count_summary <- rbind(motif_count_summary,motif_count)
    }
  }
}


increase_count <- motif_summary[["up"]] %>%   
  dplyr::count(id)
increase_tissue <- motif_summary[["up"]] %>%   
  group_by(id) %>%   
  summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
increase_count <- merge(increase_count,increase_tissue,by="id")

decrease_count <-motif_summary[["down"]] %>%   
  dplyr::count(id)
decrease_tissue <- motif_summary[["down"]] %>%   
  group_by(id) %>%   
  summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
decrease_count <- merge(decrease_count,decrease_tissue,by="id")

increase_count <- increase_count[order(increase_count$n,decreasing = TRUE),]
decrease_count <- decrease_count[order(decrease_count$n,decreasing = TRUE),]
write.csv(increase_count,paste0(data_path,"increase_common_motif_count.csv"))
write.csv(decrease_count,paste0(data_path,"decrease_common_motif_count.csv"))

motif_count_summary$condition <- factor(motif_count_summary$condition,levels=c("up","down"))
motif_count_summary$position <- motif_count_summary$count
motif_count_summary$position[which(motif_count_summary$condition=="down")] <- (-motif_count_summary$position[which(motif_count_summary$condition=="down")])
summary_by_tissue <- motif_count_summary %>%
  group_by(tissue) %>%
  summarise(count_sum = sum(count))
summary_by_tissue <- summary_by_tissue[order(summary_by_tissue$count_sum),]
motif_count_summary$tissue <- factor(motif_count_summary$tissue,levels=summary_by_tissue$tissue)
ggplot(motif_count_summary, aes(x = tissue, y = ifelse(condition == "up", count, -count), fill = condition)) +  
  geom_bar(stat = "identity") +  
  # labs(title = paste0("fdr < 0.05 motif count ",label), x = NULL, y = "Count") +
  labs(x = NULL, y = "Count") +
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

## top motif heatmap
increase_motif_data_frame <- motif_data_frame[["up"]][which(motif_data_frame[["up"]]$id %in% motif_filter_list[["up"]]),]
decrease_motif_data_frame <- motif_data_frame[["down"]][which(motif_data_frame[["down"]]$id %in% motif_filter_list[["down"]]),]
epsilon <- 1e-20
increase_motif_data_frame[increase_motif_data_frame == 0] <- epsilon
decrease_motif_data_frame[decrease_motif_data_frame == 0] <- epsilon

rownames_increase_motif_data_frame <- increase_motif_data_frame$id
rownames_decrease_motif_data_frame <- decrease_motif_data_frame$id
increase_motif_data_frame <- as.data.frame(lapply(increase_motif_data_frame[,-1], function(x) -log10(x)))
rownames(increase_motif_data_frame) <- rownames_increase_motif_data_frame

decrease_motif_data_frame <- as.data.frame(lapply(decrease_motif_data_frame[,-1], function(x) -log10(x)))
rownames(decrease_motif_data_frame) <- rownames_decrease_motif_data_frame

color_palette <- c(
  rep("white", 20), 
  colorRampPalette(c("#ffe2e2","red"))(80) 
)
breaks <- c(seq(0, -log10(0.05)+0.1, length.out = 20), seq(-log10(0.05)+0.2, 5, length.out = 80))
# p <- pheatmap::pheatmap(increase_motif_data_frame,breaks = breaks, border_color = "grey",  cellwidth = 15,cellheight = 9.5,color =color_palette ,legend_breaks = c(0,1,-log10(0.05),2, 3, 4),legend_labels = c("0", "1","-log10(0.05)", "2", "3", "4"),main = paste0("increased peaks motif ",label))
# p <- pheatmap::pheatmap(increase_motif_data_frame,breaks = breaks, border_color = "grey", color =color_palette ,main = paste0("increased peaks motif ",label),filename = "result/all/ATAC/all_tissues_increase_peaks_motif_snapatac2.png",width = 10,height = 25)
p <- pheatmap::pheatmap(increase_motif_data_frame,breaks = breaks, border_color = "grey", color =color_palette ,main = paste0("increased peaks motif ",label))
color_palette <- c(
  rep("white", 20), 
  colorRampPalette(c("#defcf9","blue"))(80) 
)
# p <- pheatmap::pheatmap(decrease_motif_data_frame,breaks = breaks, border_color = "grey", color =color_palette ,main = paste0("decreased peaks motif ",label),filename = "result/all/ATAC/all_tissues_decrease_peaks_motif_snapatac2.png",width = 10,height = 25)
pheatmap::pheatmap(decrease_motif_data_frame,breaks = breaks, border_color = "grey", color =color_palette ,main = paste0("decreased peaks motif ",label))
## all motif kmeans heatmap
condition <- "down"
motif_summary <- motif_data_frame[[condition]]
rownames(motif_summary) <- motif_summary$id
motif_summary <- motif_summary[,-1]
epsilon <- 1e-20
motif_summary[motif_summary == 0] <- epsilon
motif_summary <- -log10(motif_summary)
motif_summary_na <- motif_summary
if(condition =="up"){
  color <- c("#ffe2e2","red")
}else{
  color <- c("#defcf9","blue")
}
color_palette <- c(
  rep("white", 20), 
  colorRampPalette(color)(80) 
)
breaks <- c(seq(0, -log10(0.05)+0.1, length.out = 20), seq(-log10(0.05)+0.2, 5, length.out = 80))
to_plot <- motif_summary
pheatmap::pheatmap(to_plot, cluster_rows =T,cluster_cols = T,color = color_palette,breaks=breaks,show_rownames = F)

pca_df <- to_plot
pca_df[is.na(pca_df)] <- 0
pca <- prcomp(pca_df)
to_plot_PCA <- data.frame(pca$x)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
ggplot(to_plot_PCA, aes(x=PC1, y=PC2)) + 
  geom_point(size=1) +theme_bw()+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))

wcss <- sapply(1:20, function(k) {
  set.seed(1)
  kmeans_result <- kmeans(pca_df, centers = k, nstart = 25)
  return(kmeans_result$tot.withinss)
})
plot(1:20, wcss, type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters K",
     ylab = "Total within-cluster sum of squares")
silhouette_scores <- sapply(2:20, function(k) {
  cluster_assignment <- kmeans(pca_df, centers = k, nstart = 25)$cluster
  ss <- silhouette(cluster_assignment, dist(pca_df))
  return(mean(ss[, 3]))  
})
plot(2:20, silhouette_scores, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Average silhouette width")
set.seed(1)
umap_result <- umap::umap(as.matrix(pca_df))
umap_df <- data.frame(UMAP1 = umap_result$layout[, 1], UMAP2 = umap_result$layout[, 2], names = rownames(pca_df))
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, label = names)) +
  geom_point() +
  theme_minimal() +
  labs(title = "UMAP Projection")

set.seed(1)
if(condition == "up"){
  k <- 5
}else{
  k <- 4
}

kmeans_df <- to_plot 
kmeans_df[is.na(kmeans_df)] <- 0
kmeans_result <- kmeans(kmeans_df, centers=k)
kmeans_result <- as.data.frame(kmeans_result$cluster)
colnames(kmeans_result) <- "cluster"
umap_df <- merge(umap_df,kmeans_result,by="row.names")
umap_df$cluster <- as.character(umap_df$cluster)
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point() +
  theme_minimal() +
  labs(title = "UMAP Projection")


to_plot <- merge(motif_summary,kmeans_result,by="row.names") 
diff_summary_sorted <- to_plot[order(to_plot$cluster),]
rownames(diff_summary_sorted) <- diff_summary_sorted$Row.names
diff_summary_sorted <- diff_summary_sorted[-which(colnames(diff_summary_sorted) %in% c("Row.names","cluster"))]
data_for_heatmap <- diff_summary_sorted
annotation <- to_plot[,c("Row.names","cluster"),drop=F]
rownames(annotation) <- annotation$Row.names
annotation <- annotation[,-1,drop=F]
annotation$cluster <- as.character(annotation$cluster)
if(condition == "up"){
  annotation_color <- list(cluster=setNames(c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"),c(1:5)))
}else{
  annotation_color <- list(cluster=setNames(c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"),c(1:4)))
}

pheatmap::pheatmap(data_for_heatmap,cluster_rows = F,show_rownames = F,breaks = breaks, annotation_row = annotation,annotation_colors = annotation_color, color = color_palette, clustering_distance_cols="manhattan")

if(condition == "up"){
  to_plot <- as.data.frame(table(increase_count$n))
  label <- "increased"
}else{
  to_plot <- as.data.frame(table(decrease_count$n))
  label <- "decreased"
}

ggplot(data = to_plot, aes(x = Var1, y = Freq)) +  
  geom_bar(stat = "identity") +  
  labs(  
    title = paste0("ATAC ",label," Motif"),  
    x = NULL,  
    y = "Count"  
  ) +  
  theme_minimal() +  
  xlab(NULL) +  
  geom_text(  
    aes(label = Freq),   
    vjust = 0,   
    size = 3.5  # Adjust text size as needed  
  ) +
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
    axis.text.y = element_text(size = 12),  
    plot.title = element_text(size = 16, face = "bold")  
  ) 

