library(ggplot2)
library(Seurat)
library(data.table)

DATASET="DS_Multiome_ds"
dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"


# cellOrder=read.table(text=
#                        "HSC/MPP
# Cycling HSC/MPPs
# MEMPs
# Cycling MEMPs
# Early erythroid
# Late erythroid cells
# Mast cells
# Megakaryocytes
# Cycling megakaryocytes
# Granulocyte progenitors
# Monocyte precursors
# Inflammatory macrophages
# Kupffer cells
# NK progenitors
# NK cells
# Pre-pro B cells
# Pro-B cells
# B cells
# pDCs
# Cycling pDCs
# cDC2
# ",sep="\n")[,1]
cellOrder=read.table(text=
                       "HSCs
MEMPs
Early erythroid
Late erythroid cells
Mast cells
Megakaryocytes
Inflammatory macrophages
Kupffer cells
NK cells
Pro B cells
pDCs
cDC2
",sep="\n")[,1]


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
geneOrder=unlist(strsplit(gsub("\\ ","",geneOrder[,1]),","))

mapping = fread("/oak/stanford/groups/smontgom/amarder/t21-proj/scripts/subcluster/rna_t21_subclust_v2.csv",data.table = F,stringsAsFactors = F)
dfrna <- readRDS(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"/RNA_FindClusters.rds"))

cellNames = rownames(dfrna@meta.data)
dfrna@meta.data$i = 1:nrow(dfrna@meta.data)
meta <- dfrna@meta.data
meta = merge(meta,mapping,by.x="seurat_clusters",by.y="ClustNum")
meta <- meta[order(meta$i),]
rownames(meta) <- cellNames
meta[is.na(meta$Name),"Name"] <- "No markers"
# y=as.character(unique(meta$Name))
# meta$Name <- factor(meta$Name,levels=c(cellOrder[cellOrder %in% y],y[!(y %in% cellOrder)]))
dfrna@meta.data <- meta

# subcluster proB cells to 2
kmeans_results <- kmeans(Embeddings(subset(dfrna,Name %in% "Pro B cells"),reduction="harmony"),2)
kmeans_results1 <- data.frame(kmeans_results$cluster)
kmeans_results1[,1] <- paste0("Pro B cells,",kmeans_results1[,1])

kmeans_results <- kmeans(Embeddings(subset(dfrna,Name %in% "Neutrophils"),reduction="harmony"),2)
kmeans_results2 <- data.frame(kmeans_results$cluster)
kmeans_results2[,1] <- paste0("Neutrophils,",kmeans_results2[,1])

df.subclusters = data.frame(rbind(kmeans_results1,kmeans_results2)); colnames(df.subclusters) <- c("subclust")

meta=dfrna@meta.data
i=meta$Name %in% c("Neutrophils","Pro B cells")
meta.sub <- transform(merge(meta[i,],df.subclusters,by="row.names"), row.names=Row.names, Row.names=NULL)
meta.sub <- meta.sub[order(meta.sub$i), ]
meta.sub <- meta.sub[,colnames(meta.sub)!="id"]
meta.sub$Name <- meta.sub$subclust
meta.sub <- meta.sub[,colnames(meta.sub)!="subclust"]
which(rownames(meta[i,])!=rownames(meta.sub))
which(colnames(meta[i,])!=colnames(meta.sub))
meta[i,] <- meta.sub
dfrna@meta.data <- meta

Idents(dfrna) <- "Name"
system(paste0("rm -r ",dir,"/output/data/",DATASET,"/clust_v3"))
system(paste0("mkdir -p ",dir,"/output/data/",DATASET,"/clust_v3"))
rng = seq(1,length(geneOrder));
for (k in rng) {
  f.out <- paste0(dir,"/output/data/",DATASET,"/clust_v3/Violin.",geneOrder[k],".pdf")
  print(paste0(k,": ",f.out))
  pdf(f.out,width=17,height=9)
  print(VlnPlot(dfrna,features=geneOrder[k],sort=TRUE) + NoLegend())
  dev.off()
}

f.out = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"/RNA.meta.clust_v3.txt")
fwrite(dfrna@meta.data,
       f.out,
       quote=F,na = "NA",sep = "\t",row.names = T,col.names = T)

dfrna@meta.data$Name[dfrna@meta.data$Name=="Neutrophils,2"] <- "Granulocyte progenitors"
dfrna@meta.data$Name[dfrna@meta.data$Name=="Neutrophils,1"] <- "Neutrophils"
dfrna@meta.data$Name[dfrna@meta.data$Name=="Pro B cells,1"] <- "pre pro B cells"
dfrna@meta.data$Name[dfrna@meta.data$Name=="Pro B cells,2"] <- "pro B cells"

s="HSCs, MEMPs, Early erythroid, Cycling erythroid cells, Late erythroid cells, Mast cells, Megakaryocytes, Granulocyte progenitors, Neutrophils, Inflammatory macrophages, Kupffer cells, NK cells, T cells, pre pro B cells, pro B cells, pDCs, cDCs"
cellOrder = unlist(strsplit(s,", "))
unique(dfrna@meta.data$Name)[!unique(dfrna@meta.data$Name) %in% cellOrder]

meta <- dfrna@meta.data
y=as.character(unique(meta$Name))
meta$Name <- factor(meta$Name,levels=rev(c(cellOrder[cellOrder %in% y],y[!(y %in% cellOrder)])))
dfrna@meta.data <- meta

Idents(dfrna) <- "Name"
s="CD34, SPINK2, MLLT3, TESPA1, GATA2, GATA1, ALAS2, HBA1, HDC, ITGA2B, MPO, AZU1, SPI1, FCN1, IL1B, VCAN, CTSB, NKG7, PRF1, IL7R, RORC, PAX5, IGHM, IRF8, CLEC4C, CLEC10A, MKI67"
geneOrder = unlist(strsplit(s,", "))

f.out <- paste0(dir,"/output/data/",DATASET,"/RNA_dotplot.final.clust_v3.pdf")
pdf(f.out,width = 17,height=5)
print(
  DotPlot(
    dfrna[,!(dfrna@meta.data$Name %in% c("No markers","Unknown"))],
    features=geneOrder,
    col.min	= 0,
    col.max = 1,
    cluster.idents = FALSE
  ) + 
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))
)
dev.off()


