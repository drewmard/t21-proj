library(Seurat)
library(data.table)
# dfcombined <- readRDS("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/RNA_FindClusters.rds")
dfcombined <- readRDS("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/round2_FindClusters.rds")
tmp = dfcombined@meta.data[dfcombined$seurat_clusters %in% c(0,32),]
tmp.out <- data.frame(dataset=tmp$dataset,
                      cell=rownames(tmp),
                      cluster=tmp$seurat_clusters)
fwrite(tmp.out,"/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/exclude_cells.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

####################

library(Seurat)
library(Signac)
library(data.table)
library(ggplot2)
dfrna <- readRDS("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/RNA_FindClusters.rds")
dfatac <- readRDS("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/round2_FindClusters.rds")

rnacluster <- dfrna$seurat_clusters
dfatac <- AddMetaData(dfatac,metadata = rnacluster,col.name="rnacluster")

dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"
DATASET="DS_Multiome_h"
print("Creating ATAC UMAP with RNA clusters (round 2)...")
# p1 <- DimPlot(dfcombined, group.by = "seurat_clusters")
f.out <- paste0(dir,"/output/data/",DATASET,"/round2_umap.RNAclusters.png")
png(filename=f.out,width = 5000,height=4000,res=500)
p1 <- DimPlot(dfatac, group.by = "rnacluster") & theme(plot.title = element_text(hjust = 0.5)) 
print(LabelClusters(plot = p1, id = "rnacluster"))
dev.off()

rna_v_atac <- as.data.frame.matrix(table(dfatac@meta.data[,c("seurat_clusters","rnacluster")]))
f.out <- paste0(dir,"/output/data/",DATASET,"/round2_umap.RNAclusters.csv")
fwrite(rna_v_atac,f.out,quote = F,na = "NA",row.names = T,col.names = T,sep = ',')

table(subset(dfatac@meta.data,seurat_clusters%in%c(0,32))[,c("seurat_clusters","rnacluster")])

table(dfatac@meta.data)[,c("seurat_clusters","rnacluster")])

table(dfatac@meta.data$seurat_clusters)



colnames(dfrna)

aggregate(dfcombined@meta.data[,c("nCount_RNA","nFeature_RNA")],by=list(dataset=dfcombined@meta.data$dataset),mean)
dfcombined

aggregate(dfcombined@meta.data[,c("nCount_RNA","nFeature_RNA")],by=list(dataset=dfcombined@meta.data$dataset),length)


library(Seurat)
dfcombined <- readRDS("~/Documents/Research/neuro-variants/output/data/DS_Multiome_h/round2_FindNeighbors.rds")


library(Seurat)
library(Signac)
dfcombined <- readRDS("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/round2_FindNeighbors.rds")
