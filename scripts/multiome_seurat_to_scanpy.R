library(Seurat)

dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"
DATASET="DS_Multiome_h"
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/Multiome.RNA_ATAC.h.rds")
dfcombined1 <- readRDS(file = f.out)

# only necessary if you want to use a particular python/conda environment:
# reticulate::use_condaenv("sceasy", required = TRUE)

# f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/Multiome.RNA_ATAC.h.h5ad")
# sceasy::convertFormat(
#   dfcombined1, 
#   from = "seurat", 
#   to = "anndata", 
#   outFile = f.out
# )

# to pull out single assay
dfcombined1.rna <- dfcombined1
DefaultAssay(dfcombined1.rna) <- "RNA"
dfcombined1.rna[["ATAC"]] <- NULL
dfcombined1.atac <- dfcombined1
DefaultAssay(dfcombined1.atac) <- "ATAC"
dfcombined1.atac[["RNA"]] <- NULL
dfcombined1.atac = RenameAssays(dfcombined1.atac,ATAC="RNA")

f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/Multiome.RNA_only.h.h5ad")
sceasy::convertFormat(
  dfcombined1.rna, 
  from = "seurat", 
  to = "anndata", 
  outFile = f.out
)
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/Multiome.ATAC_only.h.h5ad")
sceasy::convertFormat(
  dfcombined1.atac, 
  from = "seurat", 
  to = "anndata", 
  outFile = f.out
)

##############################

library(Seurat)

dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"
DATASET="DS_Multiome_ds"
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/Multiome.RNA_ATAC.ds.rds")
dfcombined1 <- readRDS(file = f.out)

# only necessary if you want to use a particular python/conda environment:
# reticulate::use_condaenv("sceasy", required = TRUE)

# f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/Multiome.RNA_ATAC.h.h5ad")
# sceasy::convertFormat(
#   dfcombined1, 
#   from = "seurat", 
#   to = "anndata", 
#   outFile = f.out
# )

# to pull out single assay
dfcombined1.rna <- dfcombined1
DefaultAssay(dfcombined1.rna) <- "RNA"
dfcombined1.rna[["ATAC"]] <- NULL
dfcombined1.atac <- dfcombined1
DefaultAssay(dfcombined1.atac) <- "ATAC"
dfcombined1.atac[["RNA"]] <- NULL
dfcombined1.atac = RenameAssays(dfcombined1.atac,ATAC="RNA")

f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/Multiome.RNA_only.ds.h5ad")
sceasy::convertFormat(
  dfcombined1.rna, 
  from = "seurat", 
  to = "anndata", 
  outFile = f.out
)
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/Multiome.ATAC_only.ds.h5ad")
sceasy::convertFormat(
  dfcombined1.atac, 
  from = "seurat", 
  to = "anndata", 
  outFile = f.out
)

##############################


