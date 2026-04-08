#pagerank caculation of each sample of tissues

library(igraph)
library(edgeR)
library(stringr)

load('/path/to/grn/list.rdata')
#This load a list named grn, it's GRN information of different tissues including TF and gene
rna=read.table('/path/to/rna/matrix.txt')
#This file is matrix of raw count RNA expression (gene x sample)
rownames(rna)=gene
colnames(rna)=tissue
rna=cpm(rna,log = F,prior.count = 0)
rna=log1p(rna)
rna=rna[rowMeans(rna)>0,]
rna_scale=t(scale(t(rna)))

regulon=list()
for (i in names(grn)){
  df=grn[[i]]
  pair=unique(paste0(df$TF,':',df$gene))
  edges = data.frame(from=sub(':.*','',pair),to=sub('.*:','',pair))
  #SampleID is samples of processing tissue
  for (j in SampleID){
    rna_sub=rna[,j]
    rna_sub2=rna_scale[,j]
    node_weight = c(rna_sub2[unique(edges$from)],rna_sub2[unique(edges$to)])
    node_weight = exp(node_weight)
    edge_weight = c(rna_sub[edges$from])
    names(edge_weight)=paste0(edges$from,'|',edges$to)
    g <- graph_from_data_frame(edges, directed=TRUE)
    V(g)$weight <- node_weight[V(g)$name]
    E(g)$weight <- sapply(1:ecount(g), function(i) {
      paste0(V(g)[ends(g, i)[1]]$name, "|", V(g)[ends(g, i)[2]]$name)
    }) |> (\(keys) edge_weight[keys])()
    scores <- page_rank(g,
      algo="prpack",
      directed=TRUE,
      damping=0.85,
      personalized=V(g)$weight,
      weights=E(g)$weight)
    score_tf=scores$vector[unique(edges$from)]
    regulon[[j]]=score_tf
  }
}
save(regulon,file = '/path/to/pagerank/result.rdata')

