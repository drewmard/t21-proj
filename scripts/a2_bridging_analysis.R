# rm(list = ls())
# 
# library(Seurat)
# library(Signac)
# 
# cmd.args <- commandArgs(trailingOnly = TRUE)
# bigRNA.anndata.path <- cmd.args[1]
# # bigRNA.anndata.path <- "../../Multiome_results/test-data/10X_Healthy_Liver_sample.h5ad"
# bigRNA.seurat.path <- cmd.args[2]
# # bigRNA.seurat.path <- "../../Multiome_results/test-data/10X_Healthy_Liver_sample.rds"
# smallRNA.path <- cmd.args[3]
# # smallRNA.path <- "../../Multiome_results/RNA_FindClusters.rds"
# smallATAC.path <- cmd.args[4]
# # smallATAC.path <- "../../Multiome_results/round2_FindClusters-2.rds"
# exclude.list.path <- cmd.args[5]
# # exclude.list.path <- "../../Multiome_results/exclude_cells.txt"
# output.prefix <- cmd.args[6]
# # output.prefix <- "../../Multiome_results/bridged_results"
# output.seurat <- paste0(output.prefix, ".rds")
# output.anndata <- paste0(output.prefix, ".h5ad")
# 
# # Convert scRNA-seq anndata to Seurat
# sceasy::convertFormat(bigRNA.anndata.path, outFile=bigRNA.seurat.path, 
#                       from="anndata", to="seurat")
# 
###################


library(Seurat)
library(Signac)
library(magrittr)
# scRNA-seq data
status="Healthy"
bigRNA.seurat.path="/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/10X_Healthy_Liver.umap2d.cells_removed.rds"
# Multiome data
multiomePath="/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/"
smallRNA.path=paste0(multiomePath,"RNA_FindClusters.rds")
smallRNA_meta.path=paste0(multiomePath,"RNA_meta_v1b.txt")
smallATAC.path=paste0(multiomePath,"round2_FindClusters.GeneActivity_ChromVAR.rds")
output.prefix <- paste0(multiomePath,"multiome_transfer.v2")
output.seurat <- paste0(output.prefix, ".RNA.rds")
output.anndata <- paste0(output.prefix, ".RNA.h5ad")

bigRNA.res <- readRDS(bigRNA.seurat.path)
smallRNA.res <- readRDS(smallRNA.path)
smallATAC.res <- readRDS(smallATAC.path)
meta <- read.table(smallRNA_meta.path,stringsAsFactors = F,header=T,sep='\t',row.names=1)
smallRNA.res@meta.data <- meta
# A list of cells to exclude from multiome data (the weird cluster)
exclude.list.path <- paste0(multiomePath,"/exclude_cells.txt")
cells2exclude <- read.table(exclude.list.path, 
                            sep = "\t", header = TRUE)

# Exclude the weird cluster from multiome scATAC-seq data
print(paste("Initial Multiome scATAC-seq dimension:", 
            paste(dim(smallATAC.res), collapse = ", ")))
smallATAC.res <- smallATAC.res[, !colnames(smallATAC.res) %in% cells2exclude$cell]
print(paste("Multiome scATAC-seq dimension after cell exclusion:", 
            paste(dim(smallATAC.res), collapse = ", ")))

# Filter the multiome scRNA-seq data according to scATAC-seq cells
print(paste("Initial Multiome scRNA-seq dimension:", 
            paste(dim(smallRNA.res), collapse = ", ")))
smallRNA.res <- smallRNA.res[, colnames(smallRNA.res) %in% colnames(smallATAC.res)]
print(paste("Multiome scRNA-seq dimension after aligning with ATAC-seq:", 
            paste(dim(smallRNA.res), collapse = ", ")))

nFeature_RNA.limits=c(500,Inf)
nCount_RNA.limits=c(500,Inf)
percent.mt.limits=c(-Inf,40)
smallRNA.res <- subset(smallRNA.res, subset = (nFeature_RNA > nFeature_RNA.limits[1] & nFeature_RNA < nFeature_RNA.limits[2]) &
                         (nCount_RNA > nCount_RNA.limits[1] & nCount_RNA < nCount_RNA.limits[2]) & 
                         percent.mt > percent.mt.limits[1] & percent.mt < percent.mt.limits[2])
# dfcombined <- subset(dfcombined, subset = (nFeature_RNA > nFeature_RNA.limits[1] & nFeature_RNA < nFeature_RNA.limits[2]) &
#                          (nCount_RNA > nCount_RNA.limits[1] & nCount_RNA < nCount_RNA.limits[2]) & 
#                          percent.mt > percent.mt.limits[1] & percent.mt < percent.mt.limits[2])


# Normalise big scRNA-seq
print("NormalizeData + FindVariableFeatures + ScaleData + RunPCA...")
bigRNA.res@assays$RNA@key <- "rna_"
bigRNA.res <- NormalizeData(bigRNA.res, 
                            normalization.method = "LogNormalize", 
                            scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>% # change to all genes using , features=rownames(pbmc) if using heatmap
  RunPCA()

# Harmonise big & small scRNA-seq data

# dfcombined <- subset(dfcombined, subset = RNA_clusters %in% c(1,2,3,5,6,8,9,11,12,13,14,15,16,17,18,19,20,21,22,23,24,27))
smallRNA.res <- subset(smallRNA.res, subset = 
                         seurat_clusters %in% c(1,2,3,5,6,8,9,11,12,13,14,15,16,17,18,19,20,21,22,23,24,27,4,7,10))
# 0  4  7 10 25 26
leidenNames = colnames(bigRNA.res@meta.data)[max(grep("leiden_v",colnames(bigRNA.res@meta.data)))]
bigRNA.res@meta.data$leidenNames <- bigRNA.res@meta.data[,leidenNames]
bigRNA.res <- subset(bigRNA.res,subset = leidenNames %in% subset(bigRNA.res@meta.data,!(cell_type_groups %in% c("Stroma")))$leidenNames)

dim(bigRNA.res)
dim(smallRNA.res)
table(bigRNA.res@meta.data$leidenNames)

anchors <- FindTransferAnchors(reference = bigRNA.res, query = smallRNA.res,scale = FALSE)
predictions <- TransferData(anchorset = anchors, 
                            refdata = bigRNA.res@meta.data[,leidenNames])
smallRNA.res <- AddMetaData(smallRNA.res, metadata = predictions)

# Save results

saveRDS(smallRNA.res, file = output.seurat)
# library(sceasy)
# sceasy::convertFormat(smallRNA.res, outFile=output.anndata, 
#                       from="seurat", to="anndata")
