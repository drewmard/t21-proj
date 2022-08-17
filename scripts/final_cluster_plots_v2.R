# module load R/4.1.2
library(Seurat)
# library(Signac)
library(data.table)
library(ggplot2)

dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"
DATASET="DS_Multiome_h"

dfrna <- readRDS("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/RNA_FindClusters.rds")

f <- "/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/RNA_meta_v1.txt"
meta <- read.table(f,stringsAsFactors = F,header=T,sep='\t',row.names=1)
f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/scripts/subcluster/rna_subclust_v1_mapping_v2.csv")
mapping <- fread(f,stringsAsFactors = F,data.table = F)
mapping$ClustName <- paste0(mapping$Number,":",mapping$Name,",",mapping$Mclust)
meta2 <- meta
meta2$id  <- 1:nrow(meta2)
meta2 <- merge(meta2,mapping[,c("Number","ClustName","Final")],by.x="subclust_v1",by.y="Number")
meta2 <- meta2[order(meta2$id), ]
meta2 <- meta2[,colnames(meta2)!="id"]
rownames(meta2) <- rownames(meta)
# meta2 <- meta2[,-1]
meta2$subclust_v1_name <- meta2$ClustName
meta2 <- meta2[,colnames(meta2)!="ClustName"]
meta2$Final[is.na(meta2$Final)] <- "No markers"

cellOrder=read.table(text=
                       "HSC/MPPs
MEMPs
Early erythroid
Late erythroid
Mast cells
Megakaryocytes
Granulocyte progenitors
Monocyte progenitors
Inflammatory macrophages
Kupffer cells
NK cells
Pre pro B cells
Pro B cells
B cells
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
genes.sub <- unique(geneOrder[geneOrder %in% rownames(dfrna)])
y=as.character(unique(meta2$Final))
cellLevels=c(cellOrder[cellOrder %in% y],y[!(y %in% cellOrder)])
meta2$Final <- factor(as.character(meta2$Final),levels=cellLevels)

dfrna@meta.data <- meta2; #rm(meta); rm(meta2)

Idents(dfrna) <- "Final"


print("Creating dotplot with clusters...")
f.out <- paste0(dir,"/output/data/",DATASET,"/RNA_dotplot.final_v2.pdf")
pdf(f.out,width = 17,height=5)
print(
  DotPlot(
    dfrna[,!is.na(dfrna$Final)],
    features=genes.sub,
    col.min	= 0,
    col.max = 1,
    cluster.idents = FALSE
  ) + 
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))
)
dev.off()

