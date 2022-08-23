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
f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/scripts/subcluster/rna_h_subclust_v1_mapping_v3.csv")
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

s="HSC/MPPs, Progenitors, MEMPs, Early erythroid, Late erythroid, Mast cells, Megakaryocytes, Granulocyte progenitors, Neutrophils, Monocyte progenitors, Inflammatory macrophages, Kupffer cells, NK cells, Lymphoid progenitors, Pre pro B cells, Pro B cells, B cells, pDCs, cDCs, Hepatocytes, LSECs, No markers"
cellOrder=unlist(strsplit(s,", "))
cellOrder=rev(cellOrder)

s="CD34, PROM1, MLLT3, TESPA1, GATA2, GATA1, BLVRB, ALAS2, HBA1, HDC, ITGA2B, MPO, AZU1, FCN1, VCAN, CTSB, STAB1, IL2RB, DHFR, PAX5, IGHM, IRF8, CLEC4C, IL3RA, ALB, CDH5, MKI67"
geneOrder=unlist(strsplit(s,", "))

genes.sub <- unique(geneOrder[geneOrder %in% rownames(dfrna)])

y=as.character(unique(meta2$Final))
cellLevels=c(cellOrder[cellOrder %in% y],y[!(y %in% cellOrder)])
cellLevels <- rev(cellLevels)
meta2$Final <- factor(as.character(meta2$Final),levels=cellLevels)

dfrna@meta.data <- meta2; #rm(meta); rm(meta2)
fwrite(meta2,"/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/RNA_meta_v3.txt",quote = F,na = "NA",sep = '\t',row.names = T,col.names = T)


Idents(dfrna) <- "Final"
