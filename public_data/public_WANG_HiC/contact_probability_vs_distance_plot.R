rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(RColorBrewer)  
library(stringr)

resolution <- 20000
calDistanceProb <- function(tissue,resolution){
  resolution_k <- resolution/1000
  con.list <- list()
  for(sample in c("G1","G2","DS1","DS2")){
    con = read.table(paste0("data/public_data/WANG_cellular_aging_HiC/distance_contact/",sample,".dist.contacts.",resolution_k,"k")) 
    con$sample = sample
    con$prob = con$V2/sum(con$V2)
    con.list[[sample]]= con
  }
  combined = do.call(rbind, con.list)
  colnames(combined)[1:2] = c("distance","counts")
  color <- c("#f38181","#ff2e63","#112d4e","#3f72af")
  color <- setNames(color,c("G1","G2","DS1","DS2"))
  combined$sample <- factor(combined$sample, c("G1","G2","DS1","DS2"))
  ggplot(combined, aes(x=distance,color=sample,y=prob)) + geom_line() + 
    #  geom_smooth()+
    scale_y_log10() + 
    scale_x_log10() +
    scale_color_manual(values = color) +
    theme_bw() + 
    theme(text = element_text(size = 18))+
    xlab("Distance")+
    ylab("Probability")+
    ggtitle(tissue_label_change(tissue),"Frequency distribution of Hi-C contacts")
  
  combined_log2 = combined
  combined_log2$dist_log2 = floor(log2(combined_log2$distance+1))
  combined_log2 = aggregate(prob~dist_log2+sample,combined_log2,sum)
  combined_log2 = combined_log2[which(combined_log2$dist_log2>0),]
  combined_log2$name = combined_log2$sample
  search_table <- data.frame(sample_name=c("G1","G2","DS1","DS2"),age=c("G","G","DS","DS"))
  combined_log2 <- merge(combined_log2,search_table[,c("sample_name","age")],by.x="sample",by.y="sample_name")
  
  
  ave = aggregate(prob~dist_log2+age,combined_log2,mean)
  se = aggregate(prob~dist_log2+age,combined_log2,sd)
  out = merge(ave,se,by=c("dist_log2","age"))
  colnames(out) = c("dist_log2","name","mean","se") 
  breaks = seq(13,27,2)#13:27
  labels = ifelse (2**breaks > 1e6, paste0(round(2**breaks/1e6,1),"M"), 
                   ifelse(2**breaks> 1e3, paste0(round(2**breaks/1e3,1),"K")))
  
  p <- ggplot(combined_log2, aes(dist_log2,color=sample,y=prob,group=sample,shape=age)) + 
    geom_line(size=1.2) +
    geom_point(size=2) +
    scale_color_manual(values = color) +
    ylab("Fraction") +
    scale_x_continuous("Distance(log2)", breaks=breaks,
                       labels=labels ) +
    ylab("Probability")+
    ggtitle("Frequency distribution of Hi-C contacts")+
    theme_bw() +
    theme(text = element_text(size = 18),legend.position = "none")
  print(p)
  ggsave(paste0("data/public_data/WANG_cellular_aging_HiC/distance_contact/contact_probability_vs_distance.png"),p,width = 7,height = 7,type="cairo")
}
for(tissue in c("brain","CB","kidney", "liver", "lung", "bonemarrow", "colon", "heart", "Hip", "mammarygland", "stomach", "thymus")){
  calDistanceProb(tissue)
}
