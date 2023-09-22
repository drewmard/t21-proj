# srun --account=default --partition=interactive --time=24:00:00 --mem=128G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash

library(ggplot2)
library(Seurat)
library(Signac)

dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"

for (DATASET in c("DS_Multiome_h","DS_Multiome_ds")) { 
  print(DATASET)
  # DATASET="DS_Multiome_h"
  # DATASET="DS_Multiome_ds"
  f <- paste0(dir,"/output/data/",DATASET,"_v2","/WNN_ATAC_RNA.rds")
  df <- readRDS(file = f)
  
  df$subclust_v6 <- df$subclust_v5
  df$subclust_v6 <- df$subclust_v5
  df$subclust_v6[df$subclust_v5 %in% c("HSCs1","HSCs2")] <- "HSCs"
  df$subclust_v6[df$subclust_v5 %in% c("No marker")] <- "No markers"
  
  set.seed(03191995)
  df <- RunUMAP(df, reduction="harmony", reduction.name = "atac.umap", reduction.key = "atacUMAP_",dims=2:30)
  df <- RunUMAP(df, reduction="RNA_harmony", reduction.name = "rna.umap", reduction.key = "rnaUMAP_",dims=1:30)
  
  print("Creating UMAP with seurat clusters...")
  f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/RNA_umap.RNA_clusters.png")
  png(filename=f.out,width = 5000,height=4000,res=500)
  p1 <- DimPlot(df, group.by = "subclust_v6",reduction="rna.umap") & theme(plot.title = element_text(hjust = 0.5)) 
  print(LabelClusters(plot = p1, id = "subclust_v6"))
  dev.off()
  
  print("Creating UMAP with seurat clusters...")
  f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/ATAC_umap.RNA_clusters.png")
  png(filename=f.out,width = 5000,height=4000,res=500)
  p1 <- DimPlot(df, group.by = "subclust_v6",reduction="atac.umap") & theme(plot.title = element_text(hjust = 0.5)) 
  print(LabelClusters(plot = p1, id = "subclust_v6"))
  dev.off()
  
  print("Creating UMAP with seurat clusters...")
  f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/RNA_umap.RNA_clusters.no_label.png")
  png(filename=f.out,width = 5000,height=4000,res=500)
  p1 <- DimPlot(df, group.by = "subclust_v6",reduction="rna.umap") & theme(plot.title = element_text(hjust = 0.5)) 
  # print(LabelClusters(plot = p1, id = "subclust_v6"))
  print(p1)
  dev.off()
  
  print("Creating UMAP with seurat clusters...")
  f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/ATAC_umap.RNA_clusters.no_label.png")
  png(filename=f.out,width = 5000,height=4000,res=500)
  p1 <- DimPlot(df, group.by = "subclust_v6",reduction="atac.umap") & theme(plot.title = element_text(hjust = 0.5)) 
  # print(LabelClusters(plot = p1, id = "subclust_v6"))
  print(p1)
  dev.off()
  
}
