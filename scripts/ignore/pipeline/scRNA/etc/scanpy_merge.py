import scanpy as sc
import scarf
import pandas as pd

headDir="/oak/stanford/groups/smontgom/amarder/t21-proj"
start=1
end=99
DATASET="DS_Multiome_combined"

f_samp_lst=headDir+"/pipeline/scRNA/2_integrate_and_cluster/input/"+DATASET+"/samp_lst"
f_samp_ids=headDir+"/pipeline/scRNA/2_integrate_and_cluster/input/"+DATASET+"/samp_ids"
f_metadata=headDir+"/pipeline/scRNA/2_integrate_and_cluster/input/"+DATASET+"/scRNA.metadata.txt"

samp_lst = pd.read_csv(f_samp_lst,header=None)
samp_ids = pd.read_csv(f_samp_ids,header=None)
metadata = pd.read_csv(f_metadata,sep='\t')
num_datasets = len(samp_lst.index)

i=0
f=samp_lst.loc[i,0] + "/seurat_obj.h5ad"
adata=sc.read_h5ad(f)
adata.obs["cell"] = adata.obs.index
adata.obs[""]
metadata_sub = metadata.loc[metadata["ID2"] == samp_ids.loc[i,0],]

for i in range(num_datasets):


dfseurat[[i]]$i <- i
dfseurat[[i]]$cell <- rownames(dfseurat[[i]]@meta.data)
dfseurat[[i]]@meta.data <- cbind(dfseurat[[i]]@meta.data,subset(metadata,ID2==samp_ids[i]))
dfseurat[[i]]$dataset <- dfseurat[[i]]$ID2


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
