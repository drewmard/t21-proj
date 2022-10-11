library(ggplot2)
library(Seurat)
library(Signac)
library(data.table)
library(harmony)
library("JASPAR2020")
library("TFBSTools")
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)

####################

# args = commandArgs(trailingOnly=TRUE)
# dir=args[1]
# # start=args[2]
# # end=args[3]
# DATASET=args[4]

dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"

####################################################################

print("Reading ATAC data...")
# f <- paste0(dir,"/output/data/",DATASET,"/round2_FindClusters.rds")
f <- paste0(dir,"/output/data/","DS_Multiome_all","/round2_FindClusters.rds")
ATAC_df <- readRDS(file = f)
y = strsplit(rownames(ATAC_df@meta.data),"_")
ATAC_df@meta.data$cell <- unlist(lapply(y,function(x) x[1]))
ATAC_df@meta.data$dataNum <- unlist(lapply(y,function(x) x[2]))
ATAC_df@meta.data$cell_dataset <- paste(ATAC_df@meta.data$cell,ATAC_df@meta.data$dataset)

DATASET="DS_Multiome_h"
f <- paste0(dir,"/output/data/",DATASET,"_v2","/RNA_FindClusters.rds")
print(paste0("Loading RNA data: ",f,"..."))
RNA_df <- readRDS(file = f)
smallRNA_meta.path=paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"_v2","/meta.091522_subclust.txt")
meta <- read.table(smallRNA_meta.path,stringsAsFactors = F,header=T,sep='\t',row.names=1)
RNA_df@meta.data <- meta
RNA_df1 <- RNA_df
RNA_df1@meta.data$cell_dataset <- paste(RNA_df1@meta.data$cell,RNA_df1@meta.data$dataset)
rm(RNA_df)

DATASET="DS_Multiome_ds"
f <- paste0(dir,"/output/data/",DATASET,"_v2","/RNA_FindClusters.rds")
print(paste0("Loading RNA data: ",f,"..."))
RNA_df <- readRDS(file = f)
smallRNA_meta.path=paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"_v2","/meta.091522_subclust.txt")
meta <- read.table(smallRNA_meta.path,stringsAsFactors = F,header=T,sep='\t',row.names=1)
RNA_df@meta.data <- meta
RNA_df2 <- RNA_df
RNA_df2@meta.data$cell_dataset <- paste(RNA_df2@meta.data$cell,RNA_df2@meta.data$dataset)
rm(RNA_df)


unique(RNA_df1$subclust_v5)[!(unique(RNA_df1$subclust_v5) %in% unique(RNA_df2$subclust_v5))]
unique(RNA_df2$subclust_v5)[!(unique(RNA_df2$subclust_v5) %in% unique(RNA_df1$subclust_v5))]
RNA_df2$subclust_v6 <- RNA_df2$subclust_v5
RNA_df1$subclust_v6 <- RNA_df1$subclust_v5
RNA_df2$subclust_v6[RNA_df2$subclust_v5 %in% c("HSCs1","HSCs2")] <- "HSCs"
RNA_df2$subclust_v6[RNA_df2$subclust_v5 %in% c("No marker")] <- "No markers"

ATAC_df1.sub = subset(ATAC_df,cell_dataset %in% RNA_df1$cell_dataset)
RNA_df1.sub = subset(RNA_df1,cell_dataset %in% ATAC_df$cell_dataset)
ATAC_df2.sub = subset(ATAC_df,cell_dataset %in% RNA_df2$cell_dataset)
RNA_df2.sub = subset(RNA_df2,cell_dataset %in% ATAC_df$cell_dataset)

ATAC_df2.sub = RenameCells(ATAC_df2.sub,old.names=Cells(ATAC_df2.sub),new.names=paste(ATAC_df2.sub@meta.data$cell,as.numeric(ATAC_df2.sub@meta.data$dataNum) - 6,sep="_"))

which(rownames(RNA_df1.sub@meta.data)!=rownames(ATAC_df1.sub@meta.data))
which(rownames(RNA_df2.sub@meta.data)!=rownames(ATAC_df2.sub@meta.data))

dfcombined1 <- ATAC_df1.sub
dfcombined1[["RNA"]] <- RNA_df1.sub[["RNA"]] # CreateAssayObject(data=x)
dfcombined1[["RNA_harmony"]] <- RNA_df1.sub[["harmony"]]
dfcombined1@meta.data <- cbind(dfcombined1@meta.data,
                              RNA_df1.sub@meta.data[,c("nCount_RNA","nFeature_RNA","percent.mt","subclust_v5","subclust_v6")]) # "nCount_RNA","nFeature_RNA" already in ATAC_df

dfcombined2 <- ATAC_df2.sub
dfcombined2[["RNA"]] <- RNA_df2.sub[["RNA"]] # CreateAssayObject(data=x)
dfcombined2[["RNA_harmony"]] <- RNA_df2.sub[["harmony"]]
dfcombined2@meta.data <- cbind(dfcombined2@meta.data,
                              RNA_df2.sub@meta.data[,c("nCount_RNA","nFeature_RNA","percent.mt","subclust_v5","subclust_v6")]) # "nCount_RNA","nFeature_RNA" already in ATAC_df


library(sceasy) #

dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"
DATASET="DS_Multiome_h"
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/Multiome.RNA_ATAC.h.rds")
# saveRDS(dfcombined1,f.out)
dfcombined1 <- readRDS(file = f.out)
# reticulate::use_condaenv("sceasy", required = TRUE)
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/Multiome.RNA_ATAC.h.h5ad")
sceasy::convertFormat(
  dfcombined1, 
  from = "seurat", 
  to = "anndata", 
  outFile = f.out
)

dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"
DATASET="DS_Multiome_ds"
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/Multiome.RNA_ATAC.ds.rds")
# saveRDS(dfcombined2,f.out)
dfcombined2 <- readRDS(file = f.out)
# RETICULATE_PYTHON=/home/amarder/anaconda3/envs/r-reticulate/bin/python
# reticulate::use_condaenv("r-reticulate", required = TRUE)
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/Multiome.RNA_ATAC.ds.h5ad")
sceasy::convertFormat(
  dfcombined2, 
  from = "seurat", 
  to = "anndata", 
  outFile = f.out
)




