library(Seurat)
library(Signac)
library(data.table)

DATASET="DS_Multiome_h"
df <- readRDS(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"/RNA_FindClusters.rds"))
f.out = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"/HVG.txt")
fwrite(data.frame(VariableFeatures(df)),f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
rm(df)

DATASET="DS_Multiome_ds"
df <- readRDS(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"/RNA_FindClusters.rds"))
f.out = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"/HVG.txt")
fwrite(data.frame(VariableFeatures(df)),f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
rm(df)

#############

library(data.table)
DATASET="DS_Multiome_h"
f.out = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"/HVG.txt")
df1=fread(f.out,data.table = F,stringsAsFactors = F)[,1]
DATASET="DS_Multiome_ds"
f.out = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"/HVG.txt")
df2=fread(f.out,data.table = F,stringsAsFactors = F)[,1]
DATASET="DS_Multiome_all"
f.out = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"/HVG.txt")
df3=fread(f.out,data.table = F,stringsAsFactors = F)[,1]

mean(df1 %in% df2)
mean(df1 %in% df3)
mean(df2 %in% df3)
mean(df3 %in% df1)

mean(VariableFeatures(dfrna2) %in% VariableFeatures(dfrna))
