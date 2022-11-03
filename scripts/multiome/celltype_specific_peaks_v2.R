library(Signac)
library(Seurat)

f = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h_v2/Multiome.RNA_ATAC.h.rds")

dfcombined1 <- readRDS(file = f)

DefaultAssay(dfcombined1) <- 'ATAC'

# 5: call peaks within each cluster
print("cluster_peaks: CallPeaks... (cluster-specific)")
macs2_path <- '/home/amarder/anaconda3/envs/colocalization/bin/macs2'
cluster_peaks <- CallPeaks(
  object = dfcombined1,
  macs2.path = macs2_path,
  group.by = "subclust_v6",
  combine.peaks = FALSE
)

x <- cluster_peaks[[1]]
x <- x[x@seqnames %in% paste0("chr",1:22),]
bed <- data.frame(chr=x@seqnames,start=x@ranges@start,end=end(x@ranges))
clustNum <- x$ident[1]
f.out <- paste0(dir,"/output/data/",DATASET,"/RNA_",clustNum,".ATAC_peaks.txt")
fwrite(f.out,file=f.out,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)
