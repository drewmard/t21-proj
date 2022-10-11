library(data.table)
library(Seurat)

dir="/oak/stanford/groups/smontgom/amarder/t21-proj"
start=1
end=99
DATASET="DS_Multiome_combined"

f.samp_lst=paste0(dir,"/pipeline/scRNA/2_integrate_and_cluster/input/",DATASET,"/samp_lst")
f.samp_ids=paste0(dir,"/pipeline/scRNA/2_integrate_and_cluster/input/",DATASET,"/samp_ids")
samp_lst <- fread(f.samp_lst,data.table = F,stringsAsFactors = F,header=F,sep='\t')[,1]
samp_ids <- fread(f.samp_ids,data.table = F,stringsAsFactors = F,header=F,sep='\t')[,1]
f.metadata <- paste0(dir,"/pipeline/scRNA/2_integrate_and_cluster/input/",DATASET,"/scRNA.metadata.txt")
metadata <- fread(f.metadata,data.table = F,stringsAsFactors = F,header=T,sep='\t')

print("read_10X_and_saveRNAcounts...")
dfseurat <- list()
i=1
z=0
for (i in 1:length(samp_lst)) {
  samp <- samp_lst[[i]]
  print(paste0("Reading data from ",i,"/",length(samp_lst),": ",samp," ..."))
  f=paste0(samp,"/seurat_obj.rds")
  dfsample <- readRDS(file = f)
  
  z=z + nrow(dfsample@meta.data)
  print(paste0(z," cells..."))
  rm(dfsample)
}

