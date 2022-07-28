# module load R/4.1.2

library(Seurat)
library(Signac)
library(data.table)

source("~/Documents/Research/Useful_scripts/write_dgCMatrix_csv.R")

headDir="/oak/stanford/groups/smontgom/amarder/t21-proj"
subDirec="data"
suffixDirec="full"

outdir = paste0(headDir,"/out/full/cpdb_Data")

disease_status = "DownSyndrome"
# disease_status = "Healthy"
env = "Liver"

fout=paste0("10X_" , disease_status , "_" , env , ".umap2d.cells_removed.rds")
foutpath=paste0(headDir , "/out/" , suffixDirec , "/" , subDirec , "/" , fout)
print(paste0("\n * Reading Seurat data from file..." , foutpath))
adata=readRDS(foutpath)
colName1 <- paste0("leiden_v",max(as.numeric(substring(colnames(adata@meta.data)[grep("leiden_v",colnames(adata@meta.data))],nchar("leiden_v")+1)),na.rm=T))
adata@meta.data["leiden_names"] <- adata@meta.data[colName1]

print("Done.")

adata@assays$RNA@key <- "rna_"

# we recommend using the normalised non-log transformed data - you can save adata.X before log-transforming in adata.norm for example
adata = NormalizeData(adata, normalization.method="RC",scale.factor=1e6,assay="RNA")
# adata[["RNA"]]["KITLG",1:5]
# adata[["RNA"]][c("KIT","KITLG"),1:5]
x1=(FetchData(subset(adata,leiden_names=="Vascular endothelial cells"),vars="KITLG"))
x2=(FetchData(subset(adata,leiden_names=="HSCs/MPPs"),vars="KIT"))
mean(c(mean(x1[,1]),mean(x2[,1])))


FetchData(adata,vars=c("KIT","KITLG"))
