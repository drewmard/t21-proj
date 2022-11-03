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

print("Reading GeneActivity Multiome...")

for (i in 2:2) { 
  
  if (i==1) {
    f = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/Multiome.RNA_ATAC.ds.rds")
    f.out = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/Multiome.RNA_ATAC.ds.ChromVAR.rds")
    f.out2 = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/ds.ChromVAR.txt")
  } else {
    f = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h_v2/Multiome.RNA_ATAC.h.rds")
    f.out = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h_v2/Multiome.RNA_ATAC.h.ChromVAR.rds")
    f.out2 = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h_v2/h.ChromVAR.txt")
  }
  
  dfcombined1 <- readRDS(file = f)
  
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
  saveRDS(dfcombined1,file = f.out)
  
  chromvar_data = dfcombined1@assays$chromvar@data
  chromvar_data = as.data.frame(chromvar_data)
  chromvar_data[1:5,1:5]
  
  fwrite(chromvar_data,f.out2,quote = F,na = "NA",sep = '\t',row.names = T,col.names = T)
  
}

