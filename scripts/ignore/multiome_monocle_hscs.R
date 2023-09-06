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

f <- paste0(dir,"/output/data/",DATASET,"/Multiome.RNA_ATAC.h.GeneActivity_ChromVAR.rds")
dfcombined1 <- readRDS(file = f)

print("getMatrixByID...")
genelst <- getMatrixByID(JASPAR2020, ID = unique(rownames(dfcombined1@assays$chromvar)))

print("motif_name...")
motif_name <- as.character(
  unlist(
    mclapply(
      genelst,
      function(x) ID(x),
      mc.cores = 4)
  )
)

print("gene_name...")
gene_name <- as.character(
  unlist(
    mclapply(
      genelst,
      function(x) name(x),
      mc.cores = 4)
  )
)

print("motif-to-gene linking...")
tmp = data.frame(motif_name,gene_name)

DefaultAssay(dfcombined1) <- 'ATAC'

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

library(ggplot2)
print("Creating UMAP with seurat clusters (round 2)...")
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/ATAC_HSCs_umap.seurat_clusters.png")
png(filename=f.out,width = 5000,height=4000,res=500)
p1 <- DimPlot(df, group.by = "ATAC_snn_res.1",reduction="atac_umap") & theme(plot.title = element_text(hjust = 0.5)) 
print(LabelClusters(plot = p1, id = "ATAC_snn_res.1"))
dev.off()

print("FindAllMarkers (for ChromVAR)...")
DefaultAssay(df) <- 'chromvar'
res.all <- FindAllMarkers(df,assay = "chromvar",mean.fxn=rowMeans,fc.name = "avg_diff")
print("motif-to-gene linking...")
res.all$i <- 1:nrow(res.all)
res.all.2 <- merge(res.all,
                   data.frame(motif_name,gene_name),
                   by.x='gene',
                   by.y='motif_name')
print("reordering...")
res.all.2 <- res.all.2[order(res.all.2$i),]
res.all.2 = res.all.2[order(res.all.2$p_val),]

f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/FindAllMarkers.HSCs.ChromVAR.txt")
print(paste0("Saving to file: ",f.out))
print("...")
fwrite(res.all.2,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

library(SeuratWrappers)
library(monocle3)

DefaultAssay(dfcombined1) <- "ATAC"
dfcombined1.cds <- as.cell_data_set(dfcombined1)
dfcombined1.cds <- cluster_cells(cds = dfcombined1.cds, reduction_method = "UMAP")
dfcombined1.cds <- learn_graph(dfcombined1.cds, use_partition = TRUE)
hsc = colnames(dfcombined1.cds)[dfcombined1.cds$subclust_v6=="HSCs"]
dfcombined1.cds <- order_cells(dfcombined1.cds, reduction_method = "UMAP", root_cells = hsc)

# plot trajectories colored by pseudotime
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/ATAC.pseudotime.png")
png(filename=f.out,width = 5000,height=4000,res=500)
print(plot_cells(
  cds = dfcombined1.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
))
dev.off()

f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/ATAC.monocle_umap.png")
png(filename=f.out,width = 5000,height=4000,res=500)
print(plot_cells(
  cds = dfcombined1.cds,
  color_cells_by = "subclust_v6",
  show_trajectory_graph = TRUE
))
dev.off()


dfcombined1 <- AddMetaData(
  object = dfcombined1,
  metadata = dfcombined1.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "pseudotime"
)

# plot trajectories colored by pseudotime
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/ATAC.pseudotime2.png")
png(filename=f.out,width = 5000,height=4000,res=500)
print(FeaturePlot(dfcombined1, c("pseudotime"), pt.size = 0.1,reduction = "atac_umap") & scale_color_viridis_c())
dev.off()

# 
library(ggplot2)
print("Creating UMAP with seurat clusters (round 2)...")
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/ATAC_all_umap.seurat_clusters.png")
png(filename=f.out,width = 5000,height=4000,res=500)
p1 <- DimPlot(dfcombined1, group.by = "subclust_v6",reduction="atac_umap") & theme(plot.title = element_text(hjust = 0.5)) 
print(LabelClusters(plot = p1, id = "subclust_v6"))
dev.off()


