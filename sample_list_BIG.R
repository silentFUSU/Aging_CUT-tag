rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)

search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
search_table <- search_table[which(search_table$tissue=="brain"),]
search_table$tissue <- "Frontal_Cortex"
search_table$age_label <- "Young"
search_table$age_label[which(search_table$age=="24m")] <- "Old"
search_table$rep <- "1"
search_table$rep[which(search_table$batch=="batch2")] <- "2"
search_table$rep[which(search_table$batch=="batch3")] <- "3"

search_table$age_short <- "3"
search_table$age_short[which(search_table$age=="24m")] <- "24"
search_table$sample_label <- paste0(search_table$tissue,"-",search_table$antibody,"-",search_table$age_label,search_table$rep)

data_path <- read.table("data/samples/all/CUTTAG_ATAC_datapath.txt")
search_table <- search_table %>%
  left_join(data_path, by = c("sample_name" = "V1"))

dir.create("data/samples/BIG_list/")
write.csv(search_table,"data/samples/BIG_list/FC_cuttag.csv")
