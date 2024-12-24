rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(ggrepel)
library(reshape2)
args <- commandArgs(trailingOnly = TRUE)  
for (i in seq_along(args)) {  
  if (args[i] == "-t" && (i + 1) <= length(args)) {  
    tissue <- args[i + 1]  
  } else if (args[i] == "-r" && (i + 1) <= length(args)) {  
    resolution <- args[i + 1]  
  }  
}  
cat("Tissue:", tissue, "\n")  
cat("Resolution:", resolution, "\n")  

data <- read.table(paste0("data/samples/HiC/",tissue,"/dynamics/",tissue,"_",resolution,".dynamics"))

# Sort values for plotting
val1 <- sort(abs(data[,4]), decreasing = TRUE)
val2 <- sort(abs(data[,5]), decreasing = TRUE)
val3 <- sort(abs(data[,6]), decreasing = TRUE)

tmp <- 1:length(val1)

data_long <- data.frame(  
  tmp = rep(tmp, 3),  
  values = c(val1, val2, val3),  
  variable = rep(c("val1", "val2", "val3"), each = length(val1))  
)  

p <- ggplot(data, aes(x = tmp)) +  
  geom_line(aes(y = val1, color = 'val1')) +  
  geom_line(aes(y = val2, color = 'val2')) +  
  geom_line(aes(y = val3, color = 'val3')) +  
  scale_color_manual(values = c('val1' = 'black', 'val2' = 'red', 'val3' = 'blue')) +  
  labs(x = 'Index', y = 'Value') + 
  theme_bw() +
  ggtitle(tissue)

ggsave(filename = paste0("data/samples/HiC/",tissue,"/differential_analysis/Wang_diff_analysis_pvalue_check.png"), plot = p, width = 3.5, height = 3,type="cairo",create.dir = T)  

size1 <- nrow(data)  
data_f1 <- data[abs(data[,4]) > abs(data[,5]) & abs(data[,4]) > abs(data[,6]), ]  
data_f2 <- data_f1[abs(data_f1[,4]) >= -1 * log10(0.05), ]  

bg <- data[abs(data[,5]) >= -1 * log10(0.05) | abs(data[,6]) >= -1 * log10(0.05), ]  

re <- data_f2  
re[,2] <- re[,2] + 1 
re[,3] <- re[,3] + 1 

# Calculate specific columns in a loop  
for (i in 1:nrow(data_f2)) {  
  if (floor(i / 2000) * 2000 == i) {  
    print(i)  
  }  
  re[i,5] <- sum(abs(data_f2[i,4]) <= abs(bg[,5])) / size1  
  re[i,6] <- sum(abs(data_f2[i,4]) <= abs(bg[,6])) / size1  
} 

# write.table(re, file = paste0("data/public_data/WANG_cellular_aging_HiC/differential_analysis/G_DS_200000.FDR"), sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)  
# Write to file  
write.table(re, file = paste0("data/samples/HiC/",tissue,"/differential_analysis/",tissue,"_",resolution,".FDR"), sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

