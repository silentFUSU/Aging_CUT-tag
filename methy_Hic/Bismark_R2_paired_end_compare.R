rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
calculate_m_proportion <- function(cigar_string) {  
  # 提取所有的数字  
  all_numbers <- as.numeric(str_extract_all(cigar_string, "\\d+")[[1]])  
  
  # 提取M之前的数字  
  m_numbers <- as.numeric(str_extract_all(cigar_string, "\\d+(?=M)")[[1]])  
  
  # 计算总和  
  total_sum <- sum(all_numbers)  
  m_sum <- sum(m_numbers)  
  
  # 计算比例  
  proportion <- m_sum / total_sum  
  return(proportion)  
} 
paired <- read.delim("data/raw_data/20240804_WGBS/compare/XX315_R2_paired_end/paired_end_readnames_sort.txt",header = F)
R2 <- read.delim("data/raw_data/20240804_WGBS/compare/XX315_R2_paired_end/R2_readnames_sorted.txt",header=F)

R2_unique <- R2[which(!R2$V1 %in% paired$V1),]
paired_unique <- paired[which(!paired$V1 %in% R2$V1),]
write.table(R2_unique,"data/raw_data/20240804_WGBS/compare/XX315_R2_paired_end/R2_readnames_unique.txt",row.names = F,col.names = F,quote =  F)
write.table(paired_unique,"data/raw_data/20240804_WGBS/compare/XX315_R2_paired_end/paired_end_readnames_unique.txt",row.names = F,col.names = F,quote =  F)
overlay <- R2[which(R2$V1 %in% paired$V1),]
write.table(overlay,"data/raw_data/20240804_WGBS/compare/XX315_R2_paired_end/readnames_overlay.txt",row.names = F,col.names = F,quote =  F)

R2_info <- read.delim("data/raw_data/20240804_WGBS/compare/XX315_R2_paired_end/R2_unique_reads_info.txt")
R2_info$Mpercent <- sapply(R2_info$CIGAR, calculate_m_proportion)
R2_info$Mpercent_value<-"not all map"
R2_info$Mpercent_value[which(R2_info$Mpercent==1)] <- "all map"
write.table(R2_info$ReadName[which(R2_info$Mpercent_value=="all map")],
            "data/raw_data/20240804_WGBS/compare/XX315_R2_paired_end/R2_readnames_unique_all_map.txt",row.names = F,col.names = F,quote =  F)

write.table(R2_info$ReadName[which(R2_info$Mpercent_value=="not all map")],
            "data/raw_data/20240804_WGBS/compare/XX315_R2_paired_end/R2_readnames_unique_not_all_map.txt",row.names = F,col.names = F,quote =  F)


R2_info$mapq_condition <- "larger 40"
R2_info$mapq_condition[which(R2_info$Mapq < 40)] <- "smaller 40"
table(R2_info$mapq_condition)
write.table(R2_info$ReadName[which(R2_info$mapq_condition=="larger 40")],
            "data/raw_data/20240804_WGBS/compare/XX315_R2_paired_end/R2_readnames_unique_mapq_larger40.txt",row.names = F,col.names = F,quote =  F)

write.table(R2_info$ReadName[which(R2_info$mapq_condition=="smaller 40")],
            "data/raw_data/20240804_WGBS/compare/XX315_R2_paired_end/R2_readnames_unique_mapq_smaller40.txt",row.names = F,col.names = F,quote =  F)


