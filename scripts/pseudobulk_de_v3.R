library(Seurat)
library(data.table)
library(Matrix.utils)
library( "DESeq2" )

sampletype="Liver"
for (sampletype in c("Liver","Femur")) {
  
  disease_status="Healthy"
  f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/10X_",disease_status,"_",sampletype,".umap2d.cells_removed")
  fileName=paste0(f,".rds")
  df <- readRDS(file = fileName)
  df <- NormalizeData(df, normalization.method = "LogNormalize", scale.factor = 10000)
  colName1 <- paste0("leiden_v",max(as.numeric(substring(colnames(df@meta.data)[grep("leiden_v",colnames(df@meta.data))],nchar("leiden_v")+1)),na.rm=T))
  df@meta.data["leiden_names"] <- df@meta.data[colName1]
  
  disease_status="DownSyndrome"
  f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/10X_",disease_status,"_",sampletype,".umap2d.cells_removed")
  fileName=paste0(f,".rds")
  df2 <- readRDS(file = fileName)
  df2 <- NormalizeData(df2, normalization.method = "LogNormalize", scale.factor = 10000)
  colName2 <- paste0("leiden_v",max(as.numeric(substring(colnames(df2@meta.data)[grep("leiden_v",colnames(df2@meta.data))],nchar("leiden_v")+1)),na.rm=T))
  df2@meta.data["leiden_names"] <- df2@meta.data[colName1]
  
  # df <- df[which(rownames(df) %in% rownames(df2)),]
  # df2 <- df2[which(rownames(df2) %in% rownames(df)),]
  
  ################################################################
  # # if simulated data
  # i <- as.character(df@meta.data["leiden_names"][,1]) %in% as.character((unique(df@meta.data["leiden_names"])[,1])[1:5])
  # df2 <- df[,which(i)]
  # df2@meta.data$patient <- sample(c("Pat1","Pat2"),ncol(df2),replace = T)
  # df2@meta.data$environment <- "Down Syndrome"
  # colName2 <- paste0("leiden_v",max(as.numeric(substring(colnames(df2@meta.data)[grep("leiden_v",colnames(df2@meta.data))],nchar("leiden_v")+1)),na.rm=T))
  
  ################################################################
  
  column_to_use="cell_type_groups"
  healthy_cells=unique(df@meta.data[column_to_use][,1])
  ds_cells=unique(df2@meta.data[column_to_use][,1])
  clusters_for_DE <- healthy_cells[healthy_cells %in% ds_cells]
  cell_type=clusters_for_DE[1]
  P <- length(clusters_for_DE)
  
  iter=0; 
  for (cell_type in clusters_for_DE[1]) {
    iter = iter + 1
    print(paste0(iter,"/",P,": ",cell_type))
    cell_type_filename = gsub("/","_",cell_type)
    
    i1=which(as.character(df@meta.data[column_to_use][,1])==cell_type)
    i2=which(as.character(df2@meta.data[column_to_use][,1])==cell_type)
    dfcombined <- merge(df[,i1],
                        y=df2[,i2])

    df.aggre <- aggregate.Matrix(
      t(
        GetAssayData(object = dfcombined, slot = "counts", assay="RNA")
      ),
      groupings=dfcombined$patient,fun="sum")
    df.aggre <- t(df.aggre)
    df.aggre <- as.data.frame(df.aggre)
    system("mkdir -p /oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data_pb_groups")
    fsave=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data_pb_groups/",sampletype,".pb.",cell_type_filename,".txt")
    fwrite(df.aggre,file=fsave,quote = F,na = "NA",sep = '\t',row.names = T,col.names = T)
    
    x <- unique(dfcombined@meta.data[,c("patient","environment")]); 
    rownames(x) <- x$patient
    
    # DESeq2:
    df.aggre <- as.matrix(df.aggre[,match(rownames(x),colnames(df.aggre))])
    dds <- DESeqDataSetFromMatrix(countData=df.aggre, 
                                  colData=x, 
                                  design=~environment)
    dds$environment <- relevel(dds$environment, ref = "Healthy")
    dds <- DESeq(dds)
    res <- results(dds)
    res.df <- as.data.frame(res)
    res.df$names <- rownames(res.df)
    
    system("mkdir -p /oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_cell_type_groups")
    f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_cell_type_groups/",sampletype,".",cell_type_filename,".txt")
    fwrite(res,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
  }
  
  column_to_use="leiden_names"
  healthy_cells=unique(df@meta.data[column_to_use][,1])
  ds_cells=unique(df2@meta.data[column_to_use][,1])
  clusters_for_DE <- healthy_cells[healthy_cells %in% ds_cells]
  # cell_type=clusters_for_DE[1]
  iter=0; for (cell_type in clusters_for_DE[1]) {
    iter = iter + 1
    print(paste0(iter,"/",P,": ",cell_type))
    cell_type_filename = gsub("/","_",cell_type)
    
    dfcombined <- merge(df[,which(as.character(df@meta.data[column_to_use][,1]==cell_type))],
                        df2[,which(as.character(df2@meta.data[column_to_use][,1]==cell_type))])
    
    df.aggre <- aggregate.Matrix(
      t(
        GetAssayData(object = dfcombined, slot = "counts", assay="RNA")
      ),
      groupings=df$patient,fun="sum")
    df.aggre <- t(df.aggre)
    df.aggre <- as.data.frame(df.aggre)
    system("mkdir -p /oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data_pb_leiden")
    fsave=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data_pb_leiden/",sampletype,".pb.",cell_type_filename,".txt")
    fwrite(df.aggre,file=fsave,quote = F,na = "NA",sep = '\t',row.names = T,col.names = T)
    
    x <- unique(dfcombined@meta.data[,c("patient","environment")]); rownames(x) <- NULL
    
    # DESeq2:
    df.aggre <- as.matrix(df.aggre[,match(rownames(x),colnames(df.aggre))])
    dds <- DESeqDataSetFromMatrix(countData=df.aggre, 
                                  colData=x, 
                                  design=~environment)
    dds$environment <- relevel(dds$environment, ref = "Healthy")
    dds <- DESeq(dds)
    res <- results(dds)
    res.df <- as.data.frame(res)
    res.df$names <- rownames(res.df)
    
    system("mkdir -p /oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names")
    f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/",sampletype,".",cell_type_filename,".txt")
    fwrite(res,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
  }
}

