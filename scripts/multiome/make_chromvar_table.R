library("Seurat")
library("Signac")
library(data.table)

f = "/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h_v2/Multiome.RNA_ATAC.h.GeneActivity_ChromVAR.rds"
df = readRDS(f)

chromvar_data = df@assays$chromvar@data
chromvar_data = as.data.frame(chromvar_data)
chromvar_data[1:5,1:5]

f.out = "/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h_v2/h.ChromVAR.txt"
fwrite(chromvar_data,f.out,quote = F,na = "NA",sep = '\t',row.names = T,col.names = T)
