library(Seurat)

sampletype="Femur"
disease_status="DownSyndrome"
f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/10X_",disease_status,"_",sampletype,".umap2d.cells_removed")
fileName=paste0(f,".rds")
df <- readRDS(file = fileName)

df@meta.data[c(47839,55579),]
out=(df@assays$RNA[,47839]==df@assays$RNA[,55579])
df@assays$RNA[,c(47839,55579)]
df@meta.data[55579,]

df[,47839]
df[,55579]