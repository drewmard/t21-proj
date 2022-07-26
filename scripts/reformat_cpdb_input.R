# module load R/4.1.2
library(data.table)
meta=fread("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cpdb_Data/10X_DownSyndrome_Liver.meta.txt",data.table = F,stringsAsFactors = F)
fields=fread("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cpdb_Data/10X_DownSyndrome_Liver.fields.out",data.table = F,stringsAsFactors = F)
meta.sub <- meta[as.numeric(fields[1,]),]
fwrite(meta.sub,"/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cpdb_Data/10X_DownSyndrome_Liver.meta.sub.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

###########

# module load R/4.1.2
library(data.table)
print("Reading...")
counts=fread("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cpdb_Data/10X_DownSyndrome_Liver.norm_count.txt",data.table = F,stringsAsFactors = F)
meta=fread("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cpdb_Data/10X_DownSyndrome_Liver.meta.txt",data.table = F,stringsAsFactors = F)
# fields=fread("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cpdb_Data/10X_DownSyndrome_Liver.fields.out",data.table = F,stringsAsFactors = F)

print("Modifying...")
counts <- as.matrix(counts,rownames=1)

print("Keep...")
keep = sample(1:nrow(meta),size = 108922,replace = F)
meta.sub <- meta[keep,]
counts.sub <- counts[,keep]

fwrite(counts.sub,"/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cpdb_Data/10X_DownSyndrome_Liver.norm_count.sub.txt",quote = F,na = "NA",sep = '\t',row.names = T,col.names = T)
fwrite(meta.sub,"/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cpdb_Data/10X_DownSyndrome_Liver.meta.sub.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
fwrite(as.data.frame(keep),"/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cpdb_Data/10X_DownSyndrome_Liver.fields.out",quote = F,na = "NA",sep = '\t',row.names = F,col.names = F)


###################

library(Seurat)
library(Signac)
library(data.table)

source("/oak/stanford/groups/smontgom/amarder/etc/write_dgCMatrix_csv.R")

headDir="/oak/stanford/groups/smontgom/amarder/t21-proj"
subDirec="data"
suffixDirec="full"

outdir = paste0(headDir,"/out/full/cpdb_Data")

disease_status = "DownSyndrome"
env = "Liver"

savepath_counts = paste0(outdir , "/10X_", disease_status ,"_" , env,".norm_count.sub.csv")
savepath_meta = paste0(outdir , "/10X_" , disease_status , "_" , env , ".meta.sub.txt")
savepath_keep = paste0(outdir , "/10X_" , disease_status , "_" , env , ".keep.sub.txt")
# data after filtering and normalising
# adata = sc.read(headDir , "/10X_" , disease_status , "_" , env , suffix , ".h5ad")

fout=paste0("10X_" , disease_status , "_" , env , ".umap2d.cells_removed.rds")
foutpath=paste0(headDir , "/out/" , suffixDirec , "/" , subDirec , "/" , fout)
print(paste0("\n * Reading Seurat data from file..." , foutpath))
adata=readRDS(foutpath)
print("Done.")

print("Keep...")
keep = sample(1:nrow(adata@meta.data),size = 108922,replace = F)
keep <- sort(keep)

adata@assays$RNA@key <- "rna_"
adata=adata[,keep]

# we recommend using the normalised non-log transformed data - you can save adata.X before log-transforming in adata.norm for example
adata = NormalizeData(adata, normalization.method="RC",scale.factor=1e6,assay="RNA")

df_meta = data.frame(Cell=rownames(adata@meta.data),cell_type=adata@meta.data[,"leiden_v10"])

# library(data.table)
# fwrite(as.data.frame(adata@assays$RNA@counts),savepath_counts,quote = F,na = "NA",sep = '\t',row.names = T,col.names = T)
# writeMM(adata@assays$RNA@counts,file=savepath_counts)
# fwrite(as.matrix(adata@assays$RNA@counts),paste0(savepath_counts,".v2.txt"),quote = F,na = "NA",sep = '\t',row.names = T,col.names = T)
# write_dgCMatrix_csv(adata@assays$RNA@counts,savepath_counts)
write_dgCMatrix_csv(adata[["RNA"]],savepath_counts)
fwrite(df_meta,savepath_meta,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
fwrite(as.data.frame(keep),savepath_keep,quote = F,na = "NA",sep = '\t',row.names = F,col.names = F)
