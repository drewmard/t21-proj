# module load R/4.1.2
library(Seurat)
library(Signac)
library(data.table)
library(ggplot2)

dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"
DATASET="DS_Multiome_ds"

dfrna <- readRDS(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"/RNA_FindClusters.rds"))
Idents(dfrna) <- "seurat_clusters"

clusterData <- fread("/oak/stanford/groups/smontgom/amarder/t21-proj/scripts/subcluster/rna_subclust_v1.csv",data.table = F,stringsAsFactors = F)
genes <- unlist(lapply(strsplit(unique(clusterData$Genes),","),function(x) gsub(" ","",x)))
genes.sub <- genes[genes %in% rownames(dfrna)]

system(paste0("rm -r ",dir,"/output/data/",DATASET,"/clust"))
system(paste0("mkdir -p ",dir,"/output/data/",DATASET,"/clust"))
rng = seq(1,length(genes.sub));
for (k in rng) {
  # for (k in seq(1,length(genes))) {
  f.out <- paste0(dir,"/output/data/",DATASET,"/clust/Violin.",genes.sub[k],".pdf")
  print(paste0(k,": ",f.out))
  # png(filename=f.out,width = 5000,height=4000,res=500)
  pdf(f.out,width=17,height=9)
  # print(RidgePlot(dfrna,features=genes[k]) + NoLegend())
  print(VlnPlot(dfrna,features=genes.sub[k],sort=TRUE) + NoLegend())
  dev.off()
}


clusterData <- fread("/oak/stanford/groups/smontgom/amarder/t21-proj/scripts/subcluster/rna_subclust_v1.csv",data.table = F,stringsAsFactors = F)
genes <- unlist(lapply(strsplit(unique(clusterData$Genes),","),function(x) gsub(" ","",x)))
genes.sub <- genes[genes %in% rownames(dfrna)]
idVar="Name"
Idents(dfrna) <- idVar

system(paste0("rm -r ",dir,"/output/data/",DATASET,"/clust_v2"))
system(paste0("mkdir -p ",dir,"/output/data/",DATASET,"/clust_v2"))
rng = seq(1,length(genes.sub));
for (k in rng) {
  # for (k in seq(1,length(genes))) {
  f.out <- paste0(dir,"/output/data/",DATASET,"/clust_v2/Violin.",genes.sub[k],".pdf")
  print(paste0(k,": ",f.out))
  # png(filename=f.out,width = 5000,height=4000,res=500)
  pdf(f.out,width=17,height=9)
  # print(RidgePlot(dfrna,features=genes[k]) + NoLegend())
  print(VlnPlot(dfrna,features=genes.sub[k],sort=TRUE) + NoLegend())
  dev.off()
}
