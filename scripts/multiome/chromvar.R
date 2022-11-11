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

for (i in 4:4) { 
  
  if (i==1) {
    f = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/Multiome.RNA_ATAC.ds.rds")
    f.out = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/Multiome.RNA_ATAC.ds.ChromVAR.rds")
    f.out2 = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/ds.ChromVAR.txt")
  } else if (i==2) {
    f = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h_v2/Multiome.RNA_ATAC.h.rds")
    f.out = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h_v2/Multiome.RNA_ATAC.h.ChromVAR.rds")
    f.out2 = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h_v2/h.ChromVAR.txt")
  } else if (i==3) {
    f = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_all/round2_FindClusters.GeneActivity.rds")
    f.out = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_all/round2_FindClusters.GeneActivity.ChromVAR.rds")
    f.out2 = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_all/all.ChromVAR.txt")
  } else if (i==4) {
    f.t21 = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/Multiome.RNA_ATAC.ds.rds")
    f.h = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h_v2/Multiome.RNA_ATAC.h.rds")
    f.out = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/Multiome.RNA_ATAC.h_and_t21.ChromVAR.rds")
    f.out2 = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/h_and_t21.ChromVAR.txt")
  } 
  
  if (i != 4) {
    dfcombined1 <- readRDS(file = f)
  } else {
    dfcombined1 <- readRDS(file = f.t21)
    dfcombined2 <- readRDS(file = f.h)
    dfcombined1 = merge(dfcombined1,dfcombined2)
    
  }
  
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

