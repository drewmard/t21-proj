library(Seurat)
library(data.table)

for (disease_status in c("Healthy","DownSyndrome")) {
  for (organ in c("Femur","Liver")) {
    df=readRDS(paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/10X_",disease_status,"_",organ,".umap2d.cells_removed.rds"))
    f.out = (paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/meta.10X_",disease_status,"_",organ,".umap2d.cells_removed.v2.txt"))
    fwrite(df@meta.data,f.out,quote = F,na = "NA",sep = '\t',row.names = T,col.names = T)
  }
}

