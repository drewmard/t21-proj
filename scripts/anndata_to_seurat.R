library(reticulate)
library(Seurat)
library(scater)
library(sceasy)
# ad <- import("anndata", convert = FALSE)

sampletype="Liver"
for (sampletype in c("Liver","Femur")) {
  for (disease_status in c("Healthy","DownSyndrome")) {
    f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/10X_",disease_status,"_",sampletype,".umap2d.cells_removed")
    # f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/subset/data/10X_",disease_status,"_",sampletype,".umap.subset.cells_removed")
    fileName_in=paste0(f,".h5ad")
    fileName_out=paste0(f,".rds")
    sceasy::convertFormat(fileName_in, from="anndata", to="seurat",
                          outFile=fileName_out)
  }
}

