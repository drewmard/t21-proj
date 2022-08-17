# srun --account=default --partition=interactive --time=24:00:00 --mem=128G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash

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

DATASET="DS_Multiome_ds"
dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"

clusterData <- fread("/oak/stanford/groups/smontgom/amarder/t21-proj/scripts/subcluster/rna_subclust_v1.csv",data.table = F,stringsAsFactors = F)
genes <- unlist(lapply(strsplit(unique(clusterData$Genes),","),function(x) gsub(" ","",x)))
genes.sub <- genes[genes %in% rownames(dfrna)]

mapping = fread("/oak/stanford/groups/smontgom/amarder/t21-proj/scripts/subcluster/rna_t21_subclust_v1.csv",data.table = F,stringsAsFactors = F)
dfrna <- readRDS(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"/RNA_FindClusters.rds"))
dfatac <- readRDS(paste0(dir,"/output/data/",DATASET,"/round2_FindClusters.GeneActivity.rds"))

cellNames = rownames(dfrna@meta.data)
dfrna@meta.data$i = 1:nrow(dfrna@meta.data)
meta = merge(dfrna@meta.data,mapping,by.x="seurat_clusters",by.y="ClustNum")
meta <- meta[order(meta$i),]
dfrna@meta.data <- meta
dfrna@meta.data[is.na(dfrna$Name),"Name"] <- "No markers"
dfrna@meta.data$Name <- as.factor(dfrna@meta.data$Name)
rownames(dfrna@meta.data) <- cellNames

aggregate(dfrna@meta.data[,c("nCount_RNA","nFeature_RNA","percent.mt")],by=list(Dataset=dfrna@meta.data[,"dataset"]),median)
aggregate(dfrna@meta.data[,c("nCount_RNA","nFeature_RNA","percent.mt")],by=list(Dataset=dfrna@meta.data[,"Name"]),median)

idVar="dataset"
Idents(dfrna) <- idVar
# VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

f.out <- paste0(dir,"/output/data/",DATASET,"/RNA.QC.Violin.",idVar,".pdf")
pdf(f.out,width=17,height=9)
print(VlnPlot(dfrna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) + NoLegend())
dev.off()

idVar="Name"
Idents(dfrna) <- idVar
f.out <- paste0(dir,"/output/data/",DATASET,"/RNA.QC.Violin.",idVar,".pdf")
pdf(f.out,width=24,height=9)
print(VlnPlot(dfrna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) + NoLegend())
dev.off()

rnacluster <- dfrna$Name
dfatac <- AddMetaData(dfatac,metadata = rnacluster,col.name="rnacluster")
dfatac

print("UMAP...")
f.out <- paste0(dir,"/output/data/",DATASET,"/ATAC.RNA_Labels.png")
png(filename=f.out,width = 5000,height=4000,res=500)
p1 <- DimPlot(dfatac, group.by = "rnacluster") & theme(plot.title = element_text(hjust = 0.5)) 
print(LabelClusters(plot = p1, id = "rnacluster"))
dev.off()

print("UMAP...")
f.out <- paste0(dir,"/output/data/",DATASET,"/RNA.RNA_Labels.png")
png(filename=f.out,width = 5000,height=4000,res=500)
p1 <- DimPlot(dfrna, group.by = "Name") & theme(plot.title = element_text(hjust = 0.5)) 
print(LabelClusters(plot = p1, id = "Name"))
dev.off()

system(paste0("rm -r ",dir,"/output/data/",DATASET,"/clust_Name_RNA"))
system(paste0("mkdir -p ",dir,"/output/data/",DATASET,"/clust_Name_RNA"))
rng = seq(1,length(genes.sub));
for (k in rng) {
  # for (k in seq(1,length(genes))) {
  f.out <- paste0(dir,"/output/data/",DATASET,"/clust_Name_RNA/Violin.",genes.sub[k],".pdf")
  print(paste0(k,": ",f.out))
  # png(filename=f.out,width = 5000,height=4000,res=500)
  pdf(f.out,width=17,height=9)
  # print(RidgePlot(dfrna,features=genes[k]) + NoLegend())
  print(VlnPlot(dfrna,features=genes.sub[k],sort=TRUE) + NoLegend())
  dev.off()
}

system(paste0("rm -r ",dir,"/output/data/",DATASET,"/clust_Name_ATAC"))
system(paste0("mkdir -p ",dir,"/output/data/",DATASET,"/clust_Name_ATAC"))
idVar="rnacluster"
Idents(dfatac) <- idVar
# rng = seq(1,length(genes.sub));
rng = seq(58,length(genes.sub));
for (k in rng) {
  # for (k in seq(1,length(genes))) {
  f.out <- paste0(dir,"/output/data/",DATASET,"/clust_Name_ATAC/Violin.",genes.sub[k],".pdf")
  print(paste0(k,": ",f.out))
  # png(filename=f.out,width = 5000,height=4000,res=500)
  pdf(f.out,width=17,height=9)
  # print(RidgePlot(dfrna,features=genes[k]) + NoLegend())
  print(VlnPlot(dfatac,features=genes.sub[k],sort=TRUE) + NoLegend())
  dev.off()
}

subset(dfatac,!(rnacluster %in% c("HSCs2","Early erythroid")

dfatac





DATASET="DS_Multiome_h"
DATASET="DS_Multiome_ds"
dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"

dfrna <- readRDS(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"/RNA_FindClusters.rds"))



idVar="dataset"
Idents(dfrna) <- idVar
# VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

f.out <- paste0(dir,"/output/data/",DATASET,"/RNA.QC.Violin.",idVar,".pdf")
pdf(f.out,width=17,height=9)
print(VlnPlot(dfrna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) + NoLegend())
dev.off()


