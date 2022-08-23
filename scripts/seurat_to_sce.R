library(Seurat)
library(Signac)
f.seurat="/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/round2_FindClusters.rds"
f.sce="/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/round2_FindClusters.sce.rds"

dfatac <- readRDS(f.seurat)
sce <- as.SingleCellExperiment(dfatac)
saveRDS(sce, f.sce)
