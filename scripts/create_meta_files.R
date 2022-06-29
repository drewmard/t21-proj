library(Seurat)
library(data.table)

for (sampletype in c("Femur","Liver")) {
  for (disease_status in c("Healthy","DownSyndrome")) {
    
  disease_status="Healthy"
  f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/10X_",disease_status,"_",sampletype,".umap2d.cells_removed")
  fileName=paste0(f,".rds")
  df <- readRDS(file = fileName)
  colName1 <- paste0("leiden_v",max(as.numeric(substring(colnames(df@meta.data)[grep("leiden_v",colnames(df@meta.data))],nchar("leiden_v")+1)),na.rm=T))
  df@meta.data["leiden_names"] <- df@meta.data[colName1]
  
  f.out=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/meta.10X_",disease_status,"_",sampletype,".umap2d.cells_removed.txt")
  fwrite(df@meta.data,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
  }
}