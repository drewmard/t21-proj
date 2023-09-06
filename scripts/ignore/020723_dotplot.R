library(Seurat)
library(Signac)
library(data.table)
library(ggplot2)
dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"
# system(paste0("rm -r ",dir,"/output/data/",DATASET,"_v2","/091522_clust"))

DATASET="DS_Multiome_ds"
df = readRDS("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/Multiome.RNA_ATAC.ds.ChromVAR.rds")

DATASET="DS_Multiome_h"
df = readRDS("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h_v2/Multiome.RNA_ATAC.h.ChromVAR.rds")

s="CD34, SPINK2, PROM1, MLLT3, TESPA1, GATA2, GATA1, ALAS2, HBA1, HDC, ITGA2B, MPO, AZU1, SPI1, FCN1,VCAN,CTSB,NKG7, PRF1, IL2RB, IL7R, PAX5, IGHM, IRF8, CLEC4C, CLEC10A, MKI67"
geneOrder=unlist(strsplit(s,", |,"))
genes.sub <- unique(geneOrder[geneOrder %in% rownames(df@assays$RNA@data)])

DefaultAssay(df) <- "RNA"
# table(df$subclust_v6)
if (DATASET=="DS_Multiome_ds") {
  cellOrder = c("HSCs",
                "MEMPs",
                "Early erythroid",
                "Late erythroid",
                "Cycling erythroid",
                "Mast cells",
                "Megakaryocytes",
                "Granulocyte progenitors",
                "Neutrophils",
                "Pro-inflammatory macrophages",
                "Kupffer cells",
                "NK cells",
                "T cells",
                "B cells",
                "pDCs",
                "cDCs")
  # cellOrder=rev(cellOrder)
}
if (DATASET=="DS_Multiome_h") {
  cellOrder = c("HSCs",
                "MEMPs",
                "Early erythroid",
                "Late erythroid",
                "Mast cells",
                "Megakaryocytes",
                "Granulocyte progenitors",
                "Neutrophils",
                "Pro-inflammatory macrophages",
                "Kupffer cells",
                "NK cells",
                "B cells",
                "pDCs")
}
dim(df)
df = subset(df,subclust_v6 %in% cellOrder)
# df = df[,df@meta.data$subclust_v6 %in% cellOrder]
dim(df)
y=as.character(unique(df@meta.data$subclust_v6))
cellLevels=c(cellOrder[cellOrder %in% y],y[!(y %in% cellOrder)])
cellLevels <- rev(cellLevels)
df@meta.data$subclust_v6b <- factor(as.character(df@meta.data$subclust_v6),levels=cellLevels)
Idents(df) <- "subclust_v6b"

print("Creating dotplot with clusters...")
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/RNA_dotplot.020722.pdf")
pdf(f.out,width = 20,height=10)
print(
  DotPlot(
    df,
    features=genes.sub,
    col.min	= 0,
    col.max = 1,
    cluster.idents = FALSE
  ) + 
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))
)
dev.off()
