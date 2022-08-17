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

mapping = fread("/oak/stanford/groups/smontgom/amarder/t21-proj/scripts/subcluster/rna_t21_subclust_v1.csv",data.table = F,stringsAsFactors = F)
dfrna <- readRDS(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"/RNA_FindClusters.rds"))

cellNames = rownames(dfrna@meta.data)
dfrna@meta.data$i = 1:nrow(dfrna@meta.data)
meta <- dfrna@meta.data
meta = merge(meta,mapping,by.x="seurat_clusters",by.y="ClustNum")
meta <- meta[order(meta$i),]
rownames(meta) <- cellNames
meta[is.na(meta$Name),"Name"] <- "No markers"
y=as.character(unique(meta$Name))
meta$Name <- factor(meta$Name,levels=c(cellOrder[cellOrder %in% y],y[!(y %in% cellOrder)]))
dfrna@meta.data <- meta
# dfrna@meta.data$Name <- as.factor(dfrna@meta.data$Name)

genes.sub <- unique(geneOrder[geneOrder %in% rownames(dfrna)])
# levels(dfrna$Name) <- as.character(c(y[y %in% cellOrder],y[!(y %in% cellOrder)]))
Idents(dfrna) <- "Name"
f.out <- paste0(dir,"/output/data/",DATASET,"/RNA_dotplot.final.pdf")
pdf(f.out,width = 17,height=5)
print(
  DotPlot(
    dfrna,
    features=genes.sub,
    col.min	= 0,
    col.max = 1,
    cluster.idents = FALSE
  ) + 
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))
)
dev.off()

