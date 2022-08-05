library(ggplot2)
library(Seurat)
library(Signac)
library(data.table)
library(harmony)
library("JASPAR2020")
library("TFBSTools")
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)

####################

# args = commandArgs(trailingOnly=TRUE)
# dir=args[1]
# # start=args[2]
# # end=args[3]
# DATASET=args[4]

dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"
DATASET="DS_Multiome_h"

####################################################################

print("Reading ATAC data...")
f <- paste0(dir,"/output/data/",DATASET,"/round2_FindClusters.GeneActivity_ChromVAR.rds")
ATAC_df <- readRDS(file = f)

f <- paste0(dir,"/output/data/",DATASET,"/RNA_FindClusters.rds")
print(paste0("Loading RNA data: ",f,"..."))
RNA_df <- readRDS(file = f)
smallRNA_meta.path=paste0(dir,"/output/data/",DATASET,"/RNA_meta_v1b.txt")
meta <- read.table(smallRNA_meta.path,stringsAsFactors = F,header=T,sep='\t',row.names=1)
RNA_df@meta.data <- meta

ATAC_df <- ATAC_df[,(rownames(ATAC_df@meta.data) %in% rownames(RNA_df@meta.data))]
RNA_df <- RNA_df[,(rownames(RNA_df@meta.data) %in% rownames(ATAC_df@meta.data))]



# dfcombined <- CreateSeuratObject(GetAssayData(RNA_df))
# Idents(ATAC_df) <- "cluster"
# x <- CreateChromatinAssay(ATAC_df)
# dfcombined[["ATAC"]] <- x

colnames(ATAC_df@meta.data)[colnames(ATAC_df@meta.data)=="seurat_clusters"] <- "ATAC_clusters"
colnames(RNA_df@meta.data)[colnames(RNA_df@meta.data)=="seurat_clusters"] <- "RNA_clusters"

dfcombined <- ATAC_df
dfcombined[["RNA"]] <- RNA_df[["RNA"]] # CreateAssayObject(data=x)
dfcombined[["RNA_harmony"]] <- RNA_df[["harmony"]]
# which(rownames(RNA_df@meta.data)!=rownames(ATAC_df@meta.data))
dfcombined@meta.data <- cbind(dfcombined@meta.data,
             RNA_df@meta.data[,c("percent.mt","RNA_clusters","subclust_v1")]) # "nCount_RNA","nFeature_RNA" already in ATAC_df

f.samp_lst=paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/scripts/single_cell/input/",DATASET,"/samp_lst")
f.samp_ids=paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/scripts/single_cell/input/",DATASET,"/samp_ids")
samp_lst <- fread(f.samp_lst,data.table = F,stringsAsFactors = F,header=F,sep='\t')[,1]
samp_ids <- fread(f.samp_ids,data.table = F,stringsAsFactors = F,header=F,sep='\t')[,1]
macs2_path <- '/home/amarder/anaconda3/envs/colocalization/bin/macs2'

# 5: call peaks within each cluster
print("cluster_peaks: CallPeaks... (cluster-specific)")
cluster_peaks <- CallPeaks(
  object = dfcombined,
  macs2.path = macs2_path,
  group.by = "RNA_clusters",
  combine.peaks	= FALSE
)

f.out = paste0(dir,"/output/data/",DATASET,"/celltype_peak_set.rds")
saveRDS(cluster_peaks,f.out)

library(rtracklayer)
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)

# cluster_peaks = readRDS("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/celltype_peak_set.rds")
dir.create(paste0(dir,"/output/data/",DATASET,"/ATAC_peaks"))
for (ind in 1:length(cluster_peaks)) {
  x <- cluster_peaks[[ind]]
  x <- x[x@seqnames %in% paste0("chr",1:22),]
  bed <- data.frame(chr=x@seqnames,start=x@ranges@start,end=end(x@ranges))
  clustNum <- x$ident[1]
  f.out <- paste0(dir,"/output/data/",DATASET,"/ATAC_peaks/RNA_",clustNum,".ATAC_peaks.hg38.bed")
  fwrite(bed,file=f.out,quote = F,na = "NA",sep = "\t",row.names = F,col.names = F)
  
  seqlevelsStyle(x) = "UCSC"  # necessary
  x19 = liftOver(x, ch)
  x19 = unlist(x19)
  genome(x19) = "hg19"
  bed19 <- data.frame(chr=x19@seqnames,start=x19@ranges@start,end=end(x19@ranges))
  f.out <- paste0(dir,"/output/data/",DATASET,"/ATAC_peaks/RNA_",clustNum,".ATAC_peaks.hg19.bed")
  fwrite(bed19,file=f.out,quote = F,na = "NA",sep = "\t",row.names = F,col.names = F)
}

class(x19)


