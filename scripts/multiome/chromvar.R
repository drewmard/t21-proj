library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(Seurat)
library(Signac)
library(GenomicRanges)
library(future)
library(data.table)
library(harmony)
library("JASPAR2020")
library("TFBSTools")
library(parallel)

dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"
DATASET="DS_Multiome_h"
print("Reading GeneActivity Multiome...")
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/Multiome.RNA_ATAC.h.GeneActivity.rds")
dfcombined1 <- readRDS(file = f.out)

DefaultAssay(dfcombined1) <- 'ATAC'

print("Beginning chromvar analysis...")
print("Get a list of motif position frequency matrices from the JASPAR database...")
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)
print("add motif information...")
dfcombined1 <- AddMotifs(
  object = dfcombined1,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)
print("RunChromVAR...")
dfcombined1 <- RunChromVAR(
  object = dfcombined1,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

print("Saving round 2 results + GeneActivity + ChromVAR...")
f.out <- paste0(dir,"/output/data/",DATASET,"/Multiome.RNA_ATAC.h.GeneActivity_ChromVAR.rds")
saveRDS(dfcombined1,file = f.out)

