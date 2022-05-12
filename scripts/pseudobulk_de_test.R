library(Seurat)
library(data.table)
library(edgeR)

sampletype="Liver"

disease_status="Healthy"
# f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/10X_",disease_status,"_",sampletype,".umap2d.cells_removed")
f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/subset/data/10X_",disease_status,"_",sampletype,".umap.subset.cells_removed")
fileName=paste0(f,".rds")
df <- readRDS(file = fileName)
df <- NormalizeData(df, normalization.method = "LogNormalize", scale.factor = 10000)
colName1 <- paste0("leiden_v",max(as.numeric(substring(colnames(df@meta.data)[grep("leiden_v",colnames(df@meta.data))],nchar("leiden_v")+1)),na.rm=T))
df@meta.data["leiden_names"] <- df@meta.data[colName1]

disease_status="DownSyndrome"
# f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/10X_",disease_status,"_",sampletype,".umap2d.cells_removed")
f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/10X_",disease_status,"_",sampletype,".umap.subset.cells_removed.rds")
fileName=paste0(f,".rds")
df2 <- readRDS(file = fileName)
df2 <- NormalizeData(df2, normalization.method = "LogNormalize", scale.factor = 10000)
colName2 <- paste0("leiden_v",max(as.numeric(substring(colnames(df2@meta.data)[grep("leiden_v",colnames(df2@meta.data))],nchar("leiden_v")+1)),na.rm=T))
df2@meta.data["leiden_names"] <- df2@meta.data[colName1]

################################################################
# # if simulated data
# i <- as.character(df@meta.data["leiden_names"][,1]) %in% as.character((unique(df@meta.data["leiden_names"])[,1])[1:5])
# df2 <- df[,which(i)]
# df2@meta.data$patient <- sample(c("Pat1","Pat2"),ncol(df2),replace = T)
# df2@meta.data$environment <- "Down Syndrome"
# colName2 <- paste0("leiden_v",max(as.numeric(substring(colnames(df2@meta.data)[grep("leiden_v",colnames(df2@meta.data))],nchar("leiden_v")+1)),na.rm=T))

################################################################

column_to_use="cell_type_groups"
healthy_cells=unique(df@meta.data[column_to_use][,1])
ds_cells=unique(df2@meta.data[column_to_use][,1])
clusters_for_DE <- healthy_cells[healthy_cells %in% ds_cells]
# cell_type=clusters_for_DE[1]
P <- length(clusters_for_DE)

cell_type="Erythroid"
cell_type_filename = gsub("/","_",cell_type)
# 
iter=0
iter = iter + 1
print(paste0(iter,"/",P,": ",cell_type))
dfcombined <- merge(df[,which(df@meta.data[column_to_use][,1]==cell_type)],
                    df2[,which(df2@meta.data[column_to_use][,1]==cell_type)])
# 
# saveRDS(dfcombined,file=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/",sampletype,".sc.",cell_type_filename,".rds"))
# 
df.aggre <- AggregateExpression(dfcombined,assays="RNA",group.by="patient")
df.aggre <- df.aggre[["RNA"]]
# 
fwrite(df.aggre,file=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/",sampletype,".pb.",cell_type_filename,".txt"),quote = F,na = "NA",sep = '\t',row.names = T,col.names = T)

x <- unique(dfcombined@meta.data[,c("patient","environment")]); rownames(x) <- NULL
fwrite(x,file=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/",sampletype,".metadata.txt"),quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

