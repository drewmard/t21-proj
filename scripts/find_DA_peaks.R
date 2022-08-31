library(Seurat)
library(Signac)
library(data.table)
dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"
DATASET="DS_Multiome_h"

f <- paste0(dir,"/output/data/",DATASET,"/RNA_cluster.ATAC.rds")
dfatac <- readRDS(file = f)

Idents(dfatac) <- "rnacluster"

print("Performing QC...")
dfatac <- TSSEnrichment(dfatac)
f.out <- paste0(dir,"/output/data/",DATASET,"/QC_stats_median.txt")
fwrite(aggregate(dfatac@meta.data[,c("nCount_ATAC","nFeature_ATAC","TSS.enrichment")],by=list(dataset=dfatac@meta.data$dataset),median),
       f.out,
       quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

# aggregate(dfatac@meta.data[,c("nCount_ATAC","nFeature_ATAC","TSS.enrichment")],by=list(dataset=dfatac@meta.data$dataset),median)
# aggregate(dfatac@meta.data[,c("nCount_ATAC","nFeature_ATAC","TSS.enrichment")],by=list(dataset=dfatac@meta.data$dataset),mean)
# aggregate(dfatac@meta.data[,c("nCount_ATAC","nFeature_ATAC","TSS.enrichment")],by=list(dataset=dfatac@meta.data$dataset),sd)
dfatac <- subset(x=dfatac,subset = TSS.enrichment > 3 & nCount_ATAC >= (10^3))

print("Normalizing data...")
dfatac <- NormalizeData(
  object = dfatac,
  assay = 'ATAC',
  normalization.method = 'LogNormalize',
  scale.factor = 10000
)

print("FindAllMarkers (for GeneActivity)...")
res.all <- FindAllMarkers(dfatac)
f.out <- paste0(dir,"/output/data/",DATASET,"/FindAllMarkers.RNA_cluster.peaks.txt")
print(paste0("Saving to file: ",f.out))
print("...")
fwrite(res.all,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)



f <- paste0(dir,"/output/data/",DATASET,"/RNA_cluster.peaks.rds")
peaks=readRDS(f)


