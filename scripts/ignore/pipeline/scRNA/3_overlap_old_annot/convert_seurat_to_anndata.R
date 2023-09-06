RETICULATE_PYTHON.path="/home/amarder/anaconda3/envs/minimal_env/bin/python"
Sys.setenv(RETICULATE_PYTHON=RETICULATE_PYTHON.path)
library(reticulate)
library(data.table)
library(Seurat)
library(sceasy) #

dir="/oak/stanford/groups/smontgom/amarder/t21-proj"
start=1
end=99
DATASET="DS_Multiome_combined"

f.samp_lst=paste0(dir,"/pipeline/scRNA/2_integrate_and_cluster/input/",DATASET,"/samp_lst")
f.samp_ids=paste0(dir,"/pipeline/scRNA/2_integrate_and_cluster/input/",DATASET,"/samp_ids")
samp_lst <- fread(f.samp_lst,data.table = F,stringsAsFactors = F,header=F,sep='\t')[,1]
samp_ids <- fread(f.samp_ids,data.table = F,stringsAsFactors = F,header=F,sep='\t')[,1]

for (i in 1:length(samp_lst)) {
  samp <- samp_lst[[i]]
  print(paste0("Reading data from ",i,"/",length(samp_lst),": ",samp," ..."))
  f=paste0(samp,"/seurat_obj.rds")
  seurat <- readRDS(file = f)
  
  f.out=paste0(samp,"/seurat_obj.h5ad")
  
  sceasy::convertFormat(
    seurat, 
    from = "seurat", 
    to = "anndata", 
    outFile = f.out
  )
  
}

