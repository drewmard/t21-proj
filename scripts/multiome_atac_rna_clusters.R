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
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/Multiome.RNA_ATAC.h.rds")
dfcombined1 <- readRDS(file = f.out)

DefaultAssay(dfcombined1) <- "ATAC"

print('Calculating GeneActivity...')
dfcombined1[['GeneActivity']] <- CreateAssayObject(counts=GeneActivity(dfcombined1))

print("Normalizing data...")
dfcombined1 <- NormalizeData(
  object = dfcombined1,
  assay = 'GeneActivity',
  normalization.method = 'LogNormalize',
  scale.factor = median(dfcombined1$nCount_GeneActivity)
)

print("Saving round 2 results + GeneActivity...")
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/Multiome.RNA_ATAC.h.GeneActivity.rds")
saveRDS(dfcombined1,file = f.out)

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



df = subset(dfcombined1,subclust_v6=="HSCs")
print("FindTopFeatures...")
df <- FindTopFeatures(df,min.cutoff = 200,assay="ATAC")
# df <- FindTopFeatures(df,min.cutoff = "q75",assay="ATAC")
print("RunTFIDF...")
df <- RunTFIDF(df)
print("RunSVD...")
df <- RunSVD(df)
print("RunHarmony...")
lambda_val=1; df <- RunHarmony(df,"dataset",reduction="lsi",assay.use="ATAC",project.dim=FALSE,lambda=lambda_val)
print("RunUMAP...")
df <- RunUMAP(df, reduction = "harmony", dims = 2:30,reduction.name = "atac_umap")
print("FindNeighbors...")
df <- FindNeighbors(object = df, reduction = 'harmony', dims = 2:30)
print("FindClusters...")
resolution_val=1; algnum=1; df <- FindClusters(object = df, verbose = FALSE, algorithm = algnum,resolution=resolution_val,graph.name="ATAC_snn")
table(df$seurat_clusters)

library(ggplot2)
print("Creating UMAP with seurat clusters (round 2)...")
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/ATAC_HSCs_umap.seurat_clusters.png")
png(filename=f.out,width = 5000,height=4000,res=500)
p1 <- DimPlot(df, group.by = "ATAC_snn_res.1",reduction="atac_umap") & theme(plot.title = element_text(hjust = 0.5)) 
print(LabelClusters(plot = p1, id = "ATAC_snn_res.1"))
dev.off()
print("Creating UMAP with seurat clusters (round 2)...")
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/ATAC_HSCs_umap.GATA1.png")
png(filename=f.out,width = 5000,height=4000,res=500)
p1 <- FeaturePlot(
  object = df,
  features = c('GATA1'),
  reduction="atac_umap",
  pt.size = 0.1,
  max.cutoff = 'q95'
); print(p1)
dev.off()

DefaultAssay(df) <- "GeneActivity"
print("Creating UMAP with seurat clusters (round 2)...")
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/GeneActivity_HSCs_umap.seurat_clusters.png")
png(filename=f.out,width = 5000,height=4000,res=500)
p1 <- DimPlot(df, group.by = "ATAC_snn_res.1",reduction="atac_umap") & theme(plot.title = element_text(hjust = 0.5)) 
print(LabelClusters(plot = p1, id = "ATAC_snn_res.1"))
dev.off()
print("Creating UMAP with seurat clusters (round 2)...")
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/GeneActivity_HSCs_umap.GATA1.png")
png(filename=f.out,width = 5000,height=4000,res=500)
p1 <- FeaturePlot(
  object = df,
  features = c('GATA1'),
  reduction="atac_umap",
  pt.size = 0.1,
  max.cutoff = 'q95'
); print(p1)
dev.off()


######################################3

df = RunPCA(df,assay = "RNA")
print("RunHarmony...")
lambda_val=1; df <- RunHarmony(df,"dataset",reduction = "pca",lambda=lambda_val,reduction.save = "RNA_harmony")

print("RunUMAP...")
df <- RunUMAP(df, reduction = "RNA_harmony", dims = 1:30,reduction.name = "rna_umap")

print("FindNeighbors...")
df <- FindNeighbors(object = df, reduction = 'RNA_harmony', dims = 1:30)

print("FindClusters...")
resolution_val=1; algnum=1; df <- FindClusters(object = df, verbose = FALSE, algorithm = algnum,resolution=resolution_val,graph.name="RNA_snn")

print("Creating UMAP with seurat clusters (round 2)...")
f.out <- paste0(dir,"/output/data/",DATASET,"/RNA_HSCs_umap.seurat_clusters.png")
png(filename=f.out,width = 5000,height=4000,res=500)
p1 <- DimPlot(df, group.by = "RNA_snn_res.1",reduction="rna_umap") & theme(plot.title = element_text(hjust = 0.5)) 
print(LabelClusters(plot = p1, id = "RNA_snn_res.1"))
dev.off()

print("Creating UMAP with seurat clusters (round 2)...")
f.out <- paste0(dir,"/output/data/",DATASET,"/RNA_HSCs_umap.ATAC_clusters.png")
png(filename=f.out,width = 5000,height=4000,res=500)
p1 <- DimPlot(df, group.by = "ATAC_snn_res.1",reduction="rna_umap") & theme(plot.title = element_text(hjust = 0.5)) 
print(LabelClusters(plot = p1, id = "ATAC_snn_res.1"))
dev.off()

print("Creating UMAP with seurat clusters (round 2)...")
f.out <- paste0(dir,"/output/data/",DATASET,"/ATAC_HSCs_umap.RNA_clusters.png")
png(filename=f.out,width = 5000,height=4000,res=500)
p1 <- DimPlot(df, group.by = "RNA_snn_res.1",reduction="atac_umap") & theme(plot.title = element_text(hjust = 0.5)) 
print(LabelClusters(plot = p1, id = "RNA_snn_res.1"))
dev.off()

