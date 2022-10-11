library(Seurat)
library(data.table)
library(ggplot2)
dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"

DATASET="DS_Multiome_h"
df = readRDS(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"_v2","/RNA_FindClusters.rds"))

# library(data.table)
# f.out = "/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_all/HVG.txt"
# fwrite(data.frame(VariableFeatures(df)),f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
# mean(VariableFeatures(df) %in% VariableFeatures(dfrna2))

DATASET_tmp="DS_Multiome_h"
meta2 = fread(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET_tmp,"/RNA_meta_v3.txt"),data.table = F,stringsAsFactors = F)
rownames(meta2) <- meta2[,1]
meta2 <- meta2[,-1]

DATASET_tmp="DS_Multiome_ds"
meta3 = fread(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET_tmp,"/RNA.meta.clust_v3.txt"),data.table = F,stringsAsFactors = F)
meta3$Name[meta3$Name=="Neutrophils,2"] <- "Granulocyte progenitors"
meta3$Name[meta3$Name=="Neutrophils,1"] <- "Neutrophils"
meta3$Name[meta3$Name=="Pro B cells,1"] <- "pre pro B cells"
meta3$Name[meta3$Name=="Pro B cells,2"] <- "pro B cells"
rownames(meta3) <- meta3[,1]
meta3 <- meta3[,-1]

meta <- df@meta.data
meta$cell <- unlist(lapply(strsplit(rownames(meta),"_"),function(x) x[1]))
meta2$cell <- unlist(lapply(strsplit(rownames(meta2),"_"),function(x) x[1]))
meta3$cell <- unlist(lapply(strsplit(rownames(meta3),"_"),function(x) x[1]))

meta2$Disease <- "H"
meta3$Disease <- "T21"

tmp1 <- meta2[,c("cell","dataset","Disease","Final")]
tmp2 <- meta3[,c("cell","dataset","Disease","Name")]
colnames(tmp2)[4] <- "Final"
meta.split <- as.data.frame(rbind(tmp1,tmp2))

meta$i <- 1:nrow(meta)
meta.mg <- merge(meta,meta.split,by=c("cell","dataset"))
meta.mg <- meta.mg[order(meta.mg$i),]

meta.mg$Final[meta.mg$Final=="Pro B cells"] <- "pro B cells"
meta.mg$Final[meta.mg$Final=="Pre pro B cells"] <- "pre pro B cells"
meta.mg$Final[meta.mg$Final=="Late erythroid cells"] <- "Late erythroid"
meta.mg$Final[meta.mg$Final=="HSC/MPPs"] <- "HSCs"
meta.mg$Final_V2 <- paste0(meta.mg$Disease,":",meta.mg$Final)

rownames(meta.mg) <- rownames(meta)
meta.mg <- meta.mg[,!(colnames(meta.mg) %in% c("i"))]
df@meta.data <- meta.mg

##########
subcluster_function <- function(cellName,clustNum_to_use,kClust) {
  print(cellName)
  dftmp = subset(df,seurat_clusters %in% c(clustNum_to_use))
  kmeans_results <- kmeans(Embeddings(dftmp,reduction="harmony"),kClust)
  kmeans_results1 <- data.frame(kmeans_results$cluster)
  kmeans_results1[,1] <- paste0(cellName,":",paste(clustNum_to_use,collapse=","),".",kmeans_results1[,1])
  return(kmeans_results1)
}
mergecluster_function <- function(cellName,clustNum_to_use){
  print(cellName)
  dftmp = subset(df@meta.data,seurat_clusters %in% c(clustNum_to_use))
  kmeans_results1 <- data.frame(kmeans_results.cluster=rep(paste0(cellName,":",paste(clustNum_to_use,collapse=",")),nrow(dftmp)))
  rownames(kmeans_results1) <- rownames(dftmp)
  return(kmeans_results1)
}

if (DATASET=="DS_Multiome_h") { 
  j=0
  subcluster_results <- list(); 
  j = j + 1; subcluster_results[[j]] = subcluster_function(cellName="HSCs",clustNum_to_use = c(0),kClust = 4)
  j = j + 1; subcluster_results[[j]] = subcluster_function(cellName="MEMPs",clustNum_to_use = c(7,8,11),kClust = 4)
  j = j + 1; subcluster_results[[j]] = subcluster_function(cellName="Early erythroid",clustNum_to_use = c(3),kClust = 4)
  j = j + 1; subcluster_results[[j]] = subcluster_function(cellName="Late erythroid",clustNum_to_use = c(20),kClust = 4)
  j = j + 1; subcluster_results[[j]] = subcluster_function(cellName="Mast cell",clustNum_to_use = c(9),kClust = 4)
  j = j + 1; subcluster_results[[j]] = subcluster_function(cellName="Megakaryocytes",clustNum_to_use = c(1,14),kClust = 4)
  j = j + 1; subcluster_results[[j]] = subcluster_function(cellName="Neutrophils",clustNum_to_use = c(17),kClust = 2)
  j = j + 1; subcluster_results[[j]] = subcluster_function(cellName="Proinflammatory macrophages",clustNum_to_use = c(5),kClust = 4)
  j = j + 1; subcluster_results[[j]] = subcluster_function(cellName="Kupffer cells",clustNum_to_use = c(16),kClust = 2)
  j = j + 1; subcluster_results[[j]] = mergecluster_function(cellName="pDCs",clustNum_to_use = c(19))
  j = j + 1; subcluster_results[[j]] = mergecluster_function(cellName="NK cells",clustNum_to_use = c(10,23))
  j = j + 1; subcluster_results[[j]] = subcluster_function(cellName="?",clustNum_to_use = c(12),kClust = 4)
  j = j + 1; subcluster_results[[j]] = subcluster_function(cellName="Maybe HSCs/Prog",clustNum_to_use = c(2),kClust = 4)
  j = j + 1; subcluster_results[[j]] = subcluster_function(cellName="No markers",clustNum_to_use = c(6),kClust = 4)
  j = j + 1; subcluster_results[[j]] = subcluster_function(cellName="?",clustNum_to_use = c(18),kClust = 4)
  j = j + 1; subcluster_results[[j]] = subcluster_function(cellName="pro B cells",clustNum_to_use = c(4),kClust = 4)
  j = j + 1; subcluster_results[[j]] = subcluster_function(cellName="B cells",clustNum_to_use = c(13,15),kClust = 4)
  j = j + 1; subcluster_results[[j]] = mergecluster_function(cellName="Stroma",clustNum_to_use = c(21,22,24))
} else if (DATASET=="DS_Multiome_ds") {
  j=0
  subcluster_results <- list(); 
  j = j + 1; subcluster_results[[j]] = mergecluster_function(cellName="HSCs",clustNum_to_use = c(4,20))
  j = j + 1; subcluster_results[[j]] = subcluster_function(cellName="MEMPs",clustNum_to_use = c(3,5,7),kClust = 4)
  j = j + 1; subcluster_results[[j]] = subcluster_function(cellName="Early erythroid",clustNum_to_use = c(0,2,10,16,17,19),kClust = 4)
  j = j + 1; subcluster_results[[j]] = mergecluster_function(cellName="Late erythroid",clustNum_to_use = c(13))
  j = j + 1; subcluster_results[[j]] = mergecluster_function(cellName="Cycling erythroid",clustNum_to_use = c(9))
  j = j + 1; subcluster_results[[j]] = mergecluster_function(cellName="Mast cells",clustNum_to_use = c(18))
  j = j + 1; subcluster_results[[j]] = subcluster_function(cellName="Megakaryocytes",clustNum_to_use = c(11,14),kClust = 4)
  j = j + 1; subcluster_results[[j]] = subcluster_function(cellName="Neutrophils",clustNum_to_use = c(21),kClust = 4)
  j = j + 1; subcluster_results[[j]] = subcluster_function(cellName="Proinflammatory macrophages",clustNum_to_use = c(15),kClust = 4)
  j = j + 1; subcluster_results[[j]] = subcluster_function(cellName="Kupffer cells/cDCs",clustNum_to_use = c(6,29),kClust = 4)
  j = j + 1; subcluster_results[[j]] = mergecluster_function(cellName="pDCs",clustNum_to_use = c(23))
  j = j + 1; subcluster_results[[j]] = mergecluster_function(cellName="T cells",clustNum_to_use = c(25))
  j = j + 1; subcluster_results[[j]] = mergecluster_function(cellName="NK cells",clustNum_to_use = c(1,24))
  j = j + 1; subcluster_results[[j]] = subcluster_function(cellName="?",clustNum_to_use = c(8),kClust = 2)
  j = j + 1; subcluster_results[[j]] = subcluster_function(cellName="?",clustNum_to_use = c(26),kClust = 2)
  j = j + 1; subcluster_results[[j]] = subcluster_function(cellName="B lineage",clustNum_to_use = c(22),kClust = 4)
  j = j + 1; subcluster_results[[j]] = subcluster_function(cellName="?",clustNum_to_use = c(12),kClust = 4)
  j = j + 1; subcluster_results[[j]] = mergecluster_function(cellName="?",clustNum_to_use = c(28))
  j = j + 1; subcluster_results[[j]] = mergecluster_function(cellName="Mix-up/endothelial",clustNum_to_use = c(27))
}

df.subclusters = data.frame(do.call(rbind,subcluster_results)); colnames(df.subclusters) <- c("subclust")

meta=df@meta.data
meta[!rownames(meta) %in% rownames(df.subclusters),][1:3,]

meta$i <- 1:nrow(meta)
meta.new <- transform(merge(meta,df.subclusters,by="row.names"), row.names=Row.names, Row.names=NULL)
meta.new = meta.new[order(meta.new$i),]
df@meta.data <- meta.new

f.out = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"_v2","/meta.090822_subclust.txt")
fwrite(df@meta.data,
       f.out,
       quote=F,na = "NA",sep = "\t",row.names = T,col.names = T)

###################################################################

Idents(df) <- "subclust"
system(paste0("rm -r ",dir,"/output/data/",DATASET,"_v2","/090822_subclust"))
system(paste0("mkdir -p ",dir,"/output/data/",DATASET,"_v2","/090822_subclust"))
geneOrder=read.table(text=
                       "CD34, SPINK2, MLLT3
                     GATA1, KLF1, TESPA1
                     AHSP, ALAS2, HBA1, GYPA
                     GATA2, HDC, CPA3
                     ITGA2B, GP9
                     MPO, AZU1, SPI1, LYZ
                     CD14, CD68, S100A9, MNDA, FCN1, VCAN, IL1B
                     CD163, MS4A7, CTSB, STAB1
                     NKG7, PRF1, GZMA
                     IL7R, DHFR, PAX5, MME, IGLL1, IGHM, CD79A, CD19
                     JCHAIN, IRF8, CLEC4C, ILR3A
                     CD1C, CLEC4A, CLEC10A
                     KDR, CDH5, COL1A1, ALB, ACTA2
                     MKI67"
                     ,sep='\n')
geneOrder1=unlist(strsplit(gsub("\\ ","",geneOrder[,1]),","))
s="CD34, PROM1, MLLT3, TESPA1, GATA2, GATA1, BLVRB, ALAS2, HBA1, HDC, ITGA2B, MPO, AZU1, FCN1, VCAN, CTSB, STAB1, IL2RB, DHFR, PAX5, IGHM, IRF8, CLEC4C, IL3RA, ALB, CDH5, MKI67"
geneOrder2=unlist(strsplit(s,", "))
geneOrder = unique(geneOrder1,geneOrder2)
genes.sub <- unique(geneOrder[geneOrder %in% rownames(df)])

system(paste0("rm -r ",dir,"/output/data/",DATASET,"_v2","/090822_clust"))
system(paste0("mkdir -p ",dir,"/output/data/",DATASET,"_v2","/090822_clust"))
rng = seq(1,length(genes.sub));
for (k in rng) {
  f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/090822_clust/Violin.",geneOrder[k],".pdf")
  print(paste0(k,": ",f.out))
  pdf(f.out,width=17,height=9)
  print(VlnPlot(df,features=genes.sub[k],sort=TRUE) + NoLegend() + theme(axis.text.x=element_text(angle=90)))
  dev.off()
}

