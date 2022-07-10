library(Seurat)
library(data.table)

sampletype="Liver"
disease_status="Healthy"
res.all <- list() 
i = 0
res <- list()
for (sampletype in c("Femur","Liver")) {
  for (disease_status in c("Healthy","DownSyndrome")) {
    f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/meta.10X_",disease_status,"_",sampletype,".umap2d.cells_removed.txt")
    df=fread(f,data.table = F,stringsAsFactors = F)
    cell_type=unique(df$leiden_names)[1]
    for (cell_type in unique(df$leiden_names)) {
      i = i + 1
      tab <- table(subset(df,leiden_names==cell_type)$phase)
      res[[i]] <- as.data.frame.matrix(t(tab/sum(tab)))
      res[[i]]$cell <- cell_type
      res[[i]]$disease_status <- disease_status
      res[[i]]$sampletype <- sampletype
    }
  }
}
res.df <- do.call(dplyr::bind_rows,res)
res.df[is.na(res.df)] <- 0
fwrite(res.df,"/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/CellCycle.stats.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)


