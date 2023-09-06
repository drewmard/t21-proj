library(data.table)
flst = list.files("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/marker_genes/")
dflst = list()
for (i in 1:length(flst) ) {
  print(i)
  f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/marker_genes/",flst[i])
  dflst[[i]] = fread(f,data.table = F,stringsAsFactors = F)
}

dfall = as.data.frame(do.call(rbind,dflst))
dfall = subset(dfall,pvals_adj < 0.05 & logfoldchanges > 1)
library(dplyr)
dfall = dfall %>% select(type,organ,numerical_labels,cell_type,names,Rank,everything())
fwrite(dfall,"/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/marker_genes/cluster_markers.supp.csv",quote = F,na = "NA",sep = ',',row.names = F,col.names = T)
