library(Seurat)
library(data.table)
library(edgeR)
library(Matrix.utils)

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
  # cell_type=clusters_for_DE[1]
  P <- length(clusters_for_DE)
  iter=0; for (cell_type in clusters_for_DE) {
    iter = iter + 1
    print(paste0(iter,"/",P,": ",cell_type))
    dfcombined <- merge(df[,which(df@meta.data[column_to_use][,1]==cell_type)],
                        df2[,which(df2@meta.data[column_to_use][,1]==cell_type)])
    
    counts <- 
    pb <- aggregate.Matrix(
      t(
        GetAssayData(object = dfcombined, slot = "counts", assay="RNA")
      ),
      groupings=df$patient,fun="sum")
    pb <- t(pb)
    pb <- as.data.frame(pb)
    
    df.aggre <- AggregateExpression(dfcombined,assays="RNA",group.by="patient")
    df.aggre <- df.aggre[["RNA"]]
    
    x <- unique(dfcombined@meta.data[,c("patient","environment")]); rownames(x) <- NULL
    
    # construct design & contrast matrix
    
    design <- model.matrix(~ 0 + x$environment)
    colnames(design) <- make.names(levels(factor(x$environment)))
    rownames(design) <- x$patient
    
    # A positive FC is increased expression in the DS compared to healthy
    contrast <- makeContrasts(paste0(make.names("Down.Syndrome"),"-Healthy"), levels = design)
    
    y <- DGEList(counts=df.aggre,group=x$environment)
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design)
    fit <- glmFit(y,design)
    lrt <- glmLRT(fit,coef=2)
    res <- topTags(lrt, n=Inf)$table
    res$gene <- rownames(res)
    rownames(res) <- NULL
    
    system("mkdir -p /oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_cell_type_groups")
    cell_type_filename = gsub("/","_",cell_type)
    f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_cell_type_groups/",sampletype,".",cell_type_filename,".txt")
    fwrite(res,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
  }
  
  column_to_use="leiden_names"
  healthy_cells=unique(df@meta.data[column_to_use][,1])
  ds_cells=unique(df2@meta.data[column_to_use][,1])
  clusters_for_DE <- healthy_cells[healthy_cells %in% ds_cells]
  # cell_type=clusters_for_DE[1]
  iter=0; for (cell_type in clusters_for_DE) {
    iter = iter + 1
    print(paste0(iter,"/",P,": ",cell_type))
    
    dfcombined <- merge(df[,which(df@meta.data[column_to_use][,1]==cell_type)],
                        df2[,which(df2@meta.data[column_to_use][,1]==cell_type)])
    
    df.aggre <- AggregateExpression(dfcombined,assays="RNA",group.by="patient")
    df.aggre <- df.aggre[["RNA"]]
    
    x <- unique(dfcombined@meta.data[,c("patient","environment")]); rownames(x) <- NULL
    
    # construct design & contrast matrix
    
    design <- model.matrix(~ 0 + x$environment)
    colnames(design) <- make.names(levels(factor(x$environment)))
    rownames(design) <- x$patient
    
    # A positive FC is increased expression in the DS compared to healthy
    contrast <- makeContrasts(paste0(make.names("Down.Syndrome"),"-Healthy"), levels = design)
    
    y <- DGEList(counts=df.aggre,group=x$environment)
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design)
    fit <- glmFit(y,design)
    lrt <- glmLRT(fit,coef=2)
    res <- topTags(lrt, n=Inf)$table
    res$gene <- rownames(res)
    rownames(res) <- NULL
    
    system("mkdir -p /oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names")
    cell_type_filename = gsub("/","_",cell_type)
    f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/",sampletype,".",cell_type_filename,".txt")
    fwrite(res,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
  }
}

