# module load R/4.1.2
library(Seurat)
library(Signac)
library(data.table)
library(ggplot2)

dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"
DATASET="DS_Multiome_h"

dfrna <- readRDS(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"/RNA_FindClusters.rds"))
meta2 = fread(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"/RNA_meta_v3.txt"),data.table = F,stringsAsFactors = F)
dfrna@meta.data <- meta2; #rm(meta); rm(meta2)
dfatac <- readRDS(paste0(dir,"/output/data/",DATASET,"/round2_FindClusters.rds"))

rnacluster <- dfrna$Final
dfatac <- AddMetaData(dfatac,metadata = rnacluster,col.name="rnacluster")
# dfatac

Idents(dfrna) <- "Final"

print("cluster_peaks: CallPeaks... (cluster-specific)")
macs2_path <- '/home/amarder/anaconda3/envs/colocalization/bin/macs2'
dfatac <- dfatac[,dfatac$rnacluster!="No markers"]
# table(dfatac$rnacluster)
cluster_peaks <- CallPeaks(
  object = dfatac,
  macs2.path = macs2_path,
  group.by = "rna"
)

saveRDS(cluster_peaks,paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"/RNA_cluster.peaks.rds"))



