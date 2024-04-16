# conda activate r

library(reticulate)
library(Seurat)
# library(scater)
library(sceasy)
# ad <- import("anndata", convert = FALSE)

use_condaenv("scanpy")

sampletype="Liver"
for (sampletype in c("Femur","Liver")) {
  for (disease_status in c("Healthy","Down Syndrome")) {
    if (sampletype == "Femur" & disease_status =="Healthy") {next}
    f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/cellbender.",disease_status,".",sampletype)
    # f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/10X_",disease_status,"_",sampletype,".umap2d.cells_removed")
    # f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/subset/data/10X_",disease_status,"_",sampletype,".umap.subset.cells_removed")
    fileName_in=paste0(f,".h5ad")
    fileName_out=paste0(f,".rds")
    sceasy::convertFormat(fileName_in, from="anndata", to="seurat",
                          outFile=fileName_out)
  }
}

