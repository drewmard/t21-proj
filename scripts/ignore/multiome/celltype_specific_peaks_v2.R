library(Signac)
library(Seurat)
library(data.table)
for (i in 1:1) { 
  
  if (i==1) {
    f = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/Multiome.RNA_ATAC.ds.rds")
    dir.new = "/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/celltype_peaks"
    dir.create(dir.new)
  } else {
    f = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h_v2/Multiome.RNA_ATAC.h.rds")
    dir.new = "/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h_v2/celltype_peaks"
    dir.create(dir.new)
  }
  
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
  
  for (j in 1:length(cluster_peaks)) {
    
    x <- cluster_peaks[[j]]
    x <- x[as.character(x@seqnames) %in% paste0("chr",1:22),]
    start_end = as.data.frame(x@ranges)
    bed <- data.frame(chr=as.character(x@seqnames),start=start_end$start,end=start_end$end)
    clustID <- x$ident[1]
    f.out = paste0(dir.new,"/",clustID,".peaks.txt")
    fwrite(bed,file=f.out,quote = F,na = "NA",sep = "\t",row.names = F,col.names = F)
    
  }
}