# module load R/4.1.2

library(Seurat)
library(data.table)
library(Matrix.utils)
library( "DESeq2" )
library('variancePartition')
library('edgeR')
library('BiocParallel')

cell_type="HSCs/MPPs"; toggle=FALSE
cell_type="Cycling HSCs/MPPs"; toggle=TRUE
cell_type_filename = gsub("/","_",cell_type)

for (sampletype in c("Liver","Femur")) {
  # disease_status="Healthy"
  for (disease_status in c("Healthy","DownSyndrome")) {
    print(paste(sampletype,disease_status))
    if (toggle) {if (sampletype!="Liver" | disease_status!="DownSyndrome") {next}}
    f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/10X_",disease_status,"_",sampletype,".umap2d.cells_removed")
    fileName=paste0(f,".rds")
    df <- readRDS(file = fileName)
    # df <- NormalizeData(df, normalization.method = "LogNormalize", scale.factor = 10000)
    colName1 <- paste0("leiden_v",max(as.numeric(substring(colnames(df@meta.data)[grep("leiden_v",colnames(df@meta.data))],nchar("leiden_v")+1)),na.rm=T))
    df@meta.data["leiden_names"] <- df@meta.data[colName1]
    df@assays$RNA@key <- "RNA_"
    
    ind.lst <- which(df@meta.data["leiden_names"]==cell_type)
    system("mkdir -p /oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/cell_subset/")
    f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/cell_subset/",disease_status,'_',sampletype,'_',cell_type_filename,".rds")
    saveRDS(df[,ind.lst],file = f.out)
  }
}


# disease_status="DownSyndrome"
# f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/10X_",disease_status,"_",sampletype,".umap2d.cells_removed")
# fileName=paste0(f,".rds")
# df2 <- readRDS(file = fileName)
# # df2 <- NormalizeData(df2, normalization.method = "LogNormalize", scale.factor = 10000)
# colName2 <- paste0("leiden_v",max(as.numeric(substring(colnames(df2@meta.data)[grep("leiden_v",colnames(df2@meta.data))],nchar("leiden_v")+1)),na.rm=T))
# df2@meta.data["leiden_names"] <- df2@meta.data[colName2]
# df2@assays$RNA@key <- "RNA_"
# 
# 
# 
# sampletype="Femur"
# 
# disease_status="Healthy"
# f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/10X_",disease_status,"_",sampletype,".umap2d.cells_removed")
# fileName=paste0(f,".rds")
# df3 <- readRDS(file = fileName)
# # df <- NormalizeData(df, normalization.method = "LogNormalize", scale.factor = 10000)
# colName1 <- paste0("leiden_v",max(as.numeric(substring(colnames(df3@meta.data)[grep("leiden_v",colnames(df3@meta.data))],nchar("leiden_v")+1)),na.rm=T))
# df3@meta.data["leiden_names"] <- df3@meta.data[colName1]
# df3@assays$RNA@key <- "RNA_"
# 
# disease_status="DownSyndrome"
# f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/10X_",disease_status,"_",sampletype,".umap2d.cells_removed")
# fileName=paste0(f,".rds")
# df4 <- readRDS(file = fileName)
# # df2 <- NormalizeData(df2, normalization.method = "LogNormalize", scale.factor = 10000
# colName2 <- paste0("leiden_v",max(as.numeric(substring(colnames(df4@meta.data)[grep("leiden_v",colnames(df4@meta.data))],nchar("leiden_v")+1)),na.rm=T))
# df4@meta.data["leiden_names"] <- df4@meta.data[colName2]
# df4@assays$RNA@key <- "RNA_"
# 
# ind2.lst <- which(df2@meta.data["leiden_names"]=="HSCs/MPPs")
# ind3.lst <- which(df3@meta.data["leiden_names"]=="HSCs/MPPs")
# ind4.lst <- which(df4@meta.data["leiden_names"]=="HSCs/MPPs")
# 
# dfcombined <- merge(df[,ind.lst],
#                     y=df2[,ind2.lst])
# dfcombined <- merge(dfcombined,
#                     df3[,ind3.lst])
# dfcombined <- merge(dfcombined,
#                     df4[,ind4.lst])
# 
# 
#                     
#                     
#                     
#                     