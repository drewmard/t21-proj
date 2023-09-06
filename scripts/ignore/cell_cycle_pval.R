library(Seurat)
library(data.table)

env="Liver"
disease_status="Healthy"
res.all <- list() 
i = 0
df <- list()
for (env in c("Femur","Liver")) {
  for (disease_status in c("Healthy","DownSyndrome")) {
    i = i + 1
    print(i)
    
    f=paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/data/meta.10X_",disease_status,"_",env,".umap2d.cells_removed.txt")
    df[[i]] <- fread(f,data.table = F,stringsAsFactors = F)[,c("organ","environment","patient","sample","sorting","leiden_names","phase")]
  }
}
df.all <- as.data.frame(do.call(rbind,df))
df.all.agg <- aggregate(df.all$phase=="G1",by=list(organ=df.all$organ,
                                                   environment=df.all$environment,
                                                   patient=df.all$patient,
                                                   sample=df.all$sample,
                                                   sorting=df.all$sorting,
                                                   cell=df.all$leiden_names),mean)
df.all.agg <- aggregate(df.all$phase=="G1",by=list(organ=df.all$organ,
                                                   environment=df.all$environment,
                                                   patient=df.all$patient,
                                                   cell=df.all$leiden_names),mean)

dfout.env <- list()
for (env in c("Liver","Femur")) {
  cellCycleStats <- fread("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/data/CellCycle.stats.txt",data.table = F,stringsAsFactors = F)
  df.sub <- subset(cellCycleStats,sampletype==env)
  cell_types_to_use = subset(df.sub,disease_status=="Healthy")$cell[subset(df.sub,disease_status=="Healthy")$cell %in% subset(df.sub,disease_status=="DownSyndrome")$cell]
  
  j=1
  dfout <- list()
  for (j in 1:length(cell_types_to_use)) {
    cell_type = cell_types_to_use[j]
    tmp = subset(df.all.agg,cell == cell_type & organ==env)
    p = wilcox.test(tmp$x[tmp$environment=="Healthy"],
                    tmp$x[tmp$environment=="Down Syndrome"])$p.value
    dfout[[j]] <- data.frame(env,cell_type,p)
  }
  dfout.env[[env]] <- do.call(rbind,dfout)
}
dfout.env <- do.call(rbind,dfout.env)
dfout.env$fdr <- p.adjust(dfout.env$p,method = 'fdr')
dfout.env

f.out = "/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/data/CellCycle.pval.txt"
fwrite(dfout.env,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)







