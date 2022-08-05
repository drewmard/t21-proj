# module load R/4.1.2
library(Seurat)
library(Signac)
library(data.table)
library(ggplot2)

dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"
DATASET="DS_Multiome_h"

dfrna <- readRDS("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/RNA_FindClusters.rds")
clusterData <- fread("/oak/stanford/groups/smontgom/amarder/t21-proj/scripts/subcluster/rna_subclust.csv",data.table = F,stringsAsFactors = F)
genes <- unlist(lapply(strsplit(unique(clusterData$Genes),","),function(x) gsub(" ","",x)))
# ind <- dfrna$seurat_clusters %in% c(0)

res.lst <- list(); total_cluster_number=0
set.seed(03191995)
for (g in seq(1,max(clusterData$Group))) {
  print(g)
  tmp <- subset(clusterData,Group==g)
  clusters_to_use <- tmp$Cluster
  Kclust <- unique(tmp$Ncluster)
  
  # kmeans_results <- kmeans(Embeddings(subset(dfrna,seurat_clusters %in% clusters_to_use),reduction="harmony"),Kclust)
  # res.lst[[g]] <- data.frame(kmeans_results$cluster)
  
  # # res.lst[[g]] <- res.lst[[g]] + total_cluster_number
  # # res.lst[[g]] <- res.lst[[g]] - total_cluster_number
  
  Name=tmp$Optional[1]
  Name=ifelse(Name=="",
         paste(tmp$Cluster,collapse=","),
         Name)
         # paste0(Name," - ",paste(tmp$Cluster,collapse=",")))
  celltype = paste0(Name,",",unlist(res.lst[[g]]))
  res.lst[[g]][,1] <- celltype
  
  total_cluster_number = total_cluster_number + Kclust
}

res.df <- as.data.frame(do.call(rbind,res.lst))
colnames(res.df) <- "subclust_v1"

meta = dfrna@meta.data
meta$id  <- 1:nrow(meta)
meta <- transform(merge(meta,res.df,by="row.names"), row.names=Row.names, Row.names=NULL)
meta <- meta[order(meta$id), ]
meta <- meta[,colnames(meta)!="id"]
f.out <- "/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/RNA_meta_v1.txt"
# fwrite(meta,f.out,quote = F,na = "NA",sep = '\t',row.names = T,col.names = T)

# meta <- fread(f.out,data.table = F,stringsAsFactors = F)
# rownames(meta) <- rownames(dfrna@meta.data)

f.out <- "/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/RNA_meta_v1.txt"
meta <- read.table(f.out,stringsAsFactors = F,header=T,sep='\t',row.names=1)
f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/scripts/subcluster/rna_subclust_v1_mapping.csv")
mapping <- fread(f,stringsAsFactors = F,data.table = F)
mapping$ClustName <- paste0(mapping$Number,":",mapping$Name,",",mapping$Mclust)
meta2 <- meta
meta2$id  <- 1:nrow(meta2)
meta2 <- merge(meta2,mapping[,c("Number","ClustName")],by.x="subclust_v1",by.y="Number")
meta2 <- meta2[order(meta2$id), ]
meta2 <- meta2[,colnames(meta2)!="id"]
rownames(meta2) <- rownames(meta)
meta2 <- meta2[,-1]
meta2$subclust_v1 <- meta2$ClustName
meta2 <- meta2[,colnames(meta2)!="ClustName"]
# f.out <- "/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/RNA_meta_v1b.txt"
# fwrite(meta2,f.out,quote = F,na = "NA",sep = '\t',row.names = T,col.names = T)
dfrna@meta.data <- meta2; rm(meta); rm(meta2)

Idents(dfrna) <- "subclust_v1"
print("Creating UMAP with sub clusters...")
f.out <- paste0(dir,"/output/data/",DATASET,"/RNA_umap.subclust_v1b.png")
png(filename=f.out,width = 5000,height=4000,res=500)
p1 <- DimPlot(dfrna, group.by = "subclust_v1",label=FALSE) & theme(plot.title = element_text(hjust = 0.5)) + NoLegend() 
print(LabelClusters(plot = p1, id = "subclust_v1"))
dev.off()

aggregate(FetchData(dfrna,vars="CD163"),by=list(dfrna@meta.data$seurat_clusters),median)

k=15
rng = seq(31,length(genes)); rng <- rng[!rng%in%c(40,65)]
for (k in rng) {
# for (k in seq(1,length(genes))) {
  print(k)
  # system(paste0("mkdir -p ",dir,"/output/data/",DATASET,"/subclust_v1"))
  f.out <- paste0(dir,"/output/data/",DATASET,"/subclust_v1/Violin.",genes[k],".pdf")
  # png(filename=f.out,width = 5000,height=4000,res=500)
  pdf(f.out,width=17,height=9)
  # print(RidgePlot(dfrna,features=genes[k]) + NoLegend())
  print(VlnPlot(dfrna,features=genes[k],sort=TRUE) + NoLegend())
  
  dev.off()
}
genes[c(30,40,65)]
# FetchData(dfrna,vars=genes[k])

print("FindAllMarkers...")
# df <- FindAllMarkers(dfrna,logfc.threshold=0.25,min.pct = 0.1)
df <- FindAllMarkers(dfrna,logfc.threshold=0.5,min.pct = 0.25)
f.out <- paste0(dir,"/output/data/",DATASET,"/FindAllMarkers.RNA.subclust_v1.txt")
print(paste0("Saving to file: ",f.out))
print("...")
fwrite(df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)


# check:
# total_cluster_number
# sum(aggregate(clusterData$Ncluster,by=list(Group=clusterData$Group),mean)$x)

clusters_to_use <- c(0)
Kclust <- 4
