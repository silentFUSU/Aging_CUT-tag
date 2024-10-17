rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(readr)
library(purrr)
library(ggplot2)
library(dplyr)  
te_counts <- read.table("data/public_data/GSE132040/te_counts.txt", header = TRUE, sep = "\t")
colnames(te_counts) <- gsub(".*\\.(SRR\\d+)\\..*", "\\1", colnames(te_counts))
rownames(te_counts) <- te_counts[, 1]
te_counts <- te_counts[, -1]
rows_to_remove <- grep("^ENSMUSG", rownames(te_counts))
te_counts <- te_counts[-rows_to_remove, ]
metadata <- read.delim("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/GSE132040/GSE132040_MACA_Bulk_metadata.csv",sep = ',')
metadata$source.name <- gsub("\\_\\d+", "", metadata$source.name)
colnames(metadata)[3] <- "organ" 

te_counts_Tabula_Muris_Senis <- te_counts[, (colnames(te_counts) %in% metadata$raw.file)]
te_counts_Tabula_Muris_Senis <- te_counts_Tabula_Muris_Senis[, !(colnames(te_counts_Tabula_Muris_Senis) %in% "gene.TE")]

metadata <- metadata[which(metadata$raw.file %in% colnames(te_counts_Tabula_Muris_Senis)),]

metadata <- metadata[order(metadata$raw.file), ]
te_counts_Tabula_Muris_Senis <- te_counts_Tabula_Muris_Senis[, order(colnames(te_counts_Tabula_Muris_Senis))]

norm_factors <- calcNormFactors(te_counts_Tabula_Muris_Senis, method = "TMM")
CPM <- cpm(te_counts_Tabula_Muris_Senis, normalization_factors = norm_factors)

tissues <- unique(metadata$organ)
positive_TE_list  <- list()
negative_TE_list  <- list()
p_value_threshold <- 0.05
for(tissue in tissues){
  files <- metadata$raw.file[which(metadata$organ==tissue & metadata$characteristics..sex=="m")]
  # files <- metadata$raw.file[which(metadata$organ==tissue)]
  organ_cpm <- CPM[, files, drop = FALSE]
  ages <- metadata$characteristics..age[match(colnames(organ_cpm), metadata$raw.file)]
  organ_cpm <- rbind(organ_cpm, ages)
  rownames(organ_cpm)[nrow(organ_cpm)] <- 'age'
  positive_correlation_TE <- character()
  for (i in 1:(nrow(organ_cpm) - 1)) { 
    current_row <- rownames(organ_cpm)[i]
    current_data <- as.numeric(t(organ_cpm[i, -which(rownames(organ_cpm) == "age")]))
    age <- as.numeric(organ_cpm['age', ])
    model_data <- data.frame(age = age, current_row = current_data)
    tryCatch({
      model <- lm(current_row ~ age, data = model_data)
      p_value <- summary(model)$coefficients["age", "Pr(>|t|)"]
      if (coef(model)["age"] > 0 && p_value < p_value_threshold) {
        positive_correlation_TE <- c(positive_correlation_TE, current_row)
      }
    }, error = function(e) {
      cat("ERROR: Unable to fit linear regression for", current_row, ":", conditionMessage(e), "\n")
    })
  }
  positive_TE_list[[tissue]]  <-  positive_correlation_TE
  
  negative_correlation_TE <- character()
  for (i in 1:(nrow(organ_cpm) - 1)) {  
    current_row <- rownames(organ_cpm)[i]
    current_data <- as.numeric(t(organ_cpm[i, -which(rownames(organ_cpm) == "age")]))
    age <- as.numeric(organ_cpm['age', ])
    model_data <- data.frame(age = age, current_row = current_data)
    tryCatch({
      model <- lm(current_row ~ age, data = model_data)
      p_value <- summary(model)$coefficients["age", "Pr(>|t|)"]
      if (coef(model)["age"] < 0 && p_value < p_value_threshold) {
        negative_correlation_TE <- c(negative_correlation_TE, current_row)
      }
    }, error = function(e) {
      cat("ERROR: Unable to fit linear regression for", current_row, ":", conditionMessage(e), "\n")
    })
  }
  negative_TE_list[[tissue]]  <-  negative_correlation_TE
}

count_list <- list()
for (organ in names(positive_TE_list)) {
  count <- length(positive_TE_list[[organ]])
  count_list[[organ]] <- count
}
count_df_positive <- data.frame(organ = names(count_list), count = unlist(count_list))

count_list <- list()
for (organ in names(negative_TE_list)) {
  count <- length(negative_TE_list[[organ]])
  count_list[[organ]] <- count
}
count_df_negative <- data.frame(organ = names(count_list), count = unlist(count_list))

count_df_positive$type <- 'Positive'  
count_df_negative$type <- 'Negative'  
count_df_combined <- rbind(count_df_positive, count_df_negative)  
count_df_combined$type <- factor(count_df_combined$type, levels=c("Positive","Negative"))
count_df_combined <- count_df_combined %>%  
  group_by(organ) %>%  
  mutate(percent = count / sum(count)) %>%  
  ungroup()
count_df_combined$percent <- count_df_combined$percent * 100
tissue_order <- count_df_combined %>%  
  filter(type == "Positive") %>%  
  arrange(percent) %>%  
  distinct(organ) %>%   
  pull(organ)  
count_df_combined$organ <- factor(count_df_combined$organ, levels = tissue_order)
count_df_combined$position <- count_df_combined$count
count_df_combined$position[which(count_df_combined$type=="Negative")] <- (-count_df_combined$position[which(count_df_combined$type=="Negative")])
ggplot(count_df_combined, aes(x = organ, y = ifelse(type == "Positive", count, -count), fill = type)) +  
  geom_bar(stat = "identity") +  
  labs(title = paste0("Count of TEs correlating with aging (p-value threshold = ",p_value_threshold," )"), x = "Organ", y = "Count") +  
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
  geom_text(data =count_df_combined,   
            aes(label = count, y = position),   
            color = "black", size = 5, vjust = 0.5) + 
  scale_fill_manual(values = c("Positive" = "skyblue", "Negative" = "salmon"), name = NULL)  

saveRDS(negative_TE_list, paste0("data/public_data/GSE132040/TE/negative_TE_list_pvalue005.rds"))
saveRDS(positive_TE_list, paste0("data/public_data/GSE132040/TE/positive_TE_list_pvalue005.rds"))
