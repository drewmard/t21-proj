# conda activate seurat

library(data.table)
library(Seurat)
library(Matrix.utils) # for pseudobulking
library(BiocParallel) # for parallel
library(variancePartition) # for limma voom

disease_status = "Healthy"
sampletype = "Femur"
cell_type = "HSCs/MPPs"

# for (sampletype in c("Femur","Liver")) {
for (sampletype in c("Liver")) {
  for (disease_status in c("Healthy","Down Syndrome")) {
    
    print(paste("Reading:",disease_status,sampletype,cell_type))
    f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/cellbender.",disease_status,".",sampletype,".rds")
    dfseurat=readRDS(f)
    dfseurat@assays$RNA@key = "rna_"
    dfseurat@meta.data$i = 1:nrow(dfseurat@meta.data)
    
    f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/meta.10X_",sub(" ","",disease_status),"_",sampletype,".umap2d.cells_removed.v2.txt")
    meta = fread(f,data.table = F,stringsAsFactors = F)
    meta$cell = sub("-.*","",meta$V1)
    
    tmp = dfseurat@meta.data
    tmp$cell = substring(rownames(tmp),nchar(as.character(tmp$Channel))+2)
    tmp$rowname = rownames(tmp)
    # tmp$i = 1:nrow(tmp)
    metainfo = fread("/home/amarder/tmp/metadata.csv",data.table = F,stringsAsFactors = F)
    colnames(metainfo)[1] = "sample"
    tmp = merge(tmp,metainfo[,1:2],by.x="Channel",by.y="CellrangerFile")
    meta$cell_sample = paste(meta$cell,meta$sample)
    meta = meta[!duplicated(meta$cell_sample),]
    tmp = merge(tmp,meta,by=c("cell","sample"))
    sum(duplicated(tmp$rowname))
    tmp = tmp[order(tmp$i,decreasing = F),]
    rownames(tmp) = tmp$rowname
    
    dfseurat = subset(dfseurat,i %in% tmp$i)
    # sum(tmp$i != dfseurat@meta.data$i)
    # sum(rownames(tmp) != rownames(dfseurat@meta.data))
    dfseurat@meta.data = tmp
    
    y=grep("leiden_v",colnames(dfseurat@meta.data))
    y=colnames(dfseurat@meta.data)[y[length(y)]]
    dfseurat@meta.data$leiden_names = dfseurat@meta.data[,y]
    
    # Subset matrix down to cell type of interest
    dfsub = dfseurat[,dfseurat@meta.data[,"leiden_names"]==cell_type]
    print(paste0("Subsetting ",ncol(dfseurat)," cells into ",ncol(dfsub)," cells..."))
    
    dir.create(paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/cellbender_sub"))
    dir.create(paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/cellbender_sub/HSC"))
    saveRDS(dfsub,paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/cellbender_sub/HSC/",disease_status,"_",sampletype,".rds"))
    
    seurat_meta = dfsub@meta.data
    # create pseudobulks:
    df.aggre <- aggregate.Matrix(
      t(
        GetAssayData(object = dfsub, slot = "counts", assay="RNA")
      ),
      groupings=seurat_meta[,"sample"],fun="sum")
    df.aggre <- t(df.aggre)
    df.aggre <- as.data.frame(df.aggre)
    
    # create metadata:
    col_of_interest = c("patient","sample","sorting","environment","organ","age")
    x <- unique(seurat_meta[,col_of_interest])
    rownames(x) <- as.character(x[,"sample"])
    
    # align metadata and pseudobulks:
    metadata_to_use <- x[rownames(x) %in% colnames(df.aggre),]
    df.aggre <- as.data.frame.matrix(df.aggre[,match(rownames(metadata_to_use),colnames(df.aggre))])
    
    # 
    tab <- table(dfsub@meta.data[dfsub@meta.data[,"leiden_names"]==cell_type,"sample"])
    samples_to_keep <- names(tab)[tab >= 10]
    samples_to_keep <- samples_to_keep[samples_to_keep %in% rownames(metadata_to_use)]
    if (length(samples_to_keep)==1) {
      tmprow = rownames(df.aggre)
      df.aggre <- data.frame(sample=df.aggre[,samples_to_keep])
      colnames(df.aggre) = samples_to_keep
      rownames(df.aggre) = tmprow
    } else {
      df.aggre <- df.aggre[,samples_to_keep]
    }
    metadata_to_use <- metadata_to_use[samples_to_keep,]
    
    # save:
    cell_type_filename = gsub("/","_",cell_type)
    
    dir.create("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/cellbender_sub/HSC/pseudobulks")
    f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/cellbender_sub/HSC/pseudobulks/",disease_status,'_',sampletype,".pb.txt")
    fwrite(df.aggre,f.out,na = "NA",sep = "\t",quote = F,row.names = T,col.names = T)
    
    f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/cellbender_sub/HSC/pseudobulks/",disease_status,'_',sampletype,".meta.txt")
    fwrite(metadata_to_use,f.out,na = "NA",sep = "\t",quote = F,row.names = T,col.names = T)
    
  }
}


# for (disease_status in c("Healthy","DownSyndrome")) {
for (disease_status in c("Down Syndrome")) {
  
  # disease_status="Healthy"
  print(paste0("Running analysis: ",disease_status,"..."))
  
  sampletype="Liver"
  f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/cellbender_sub/HSC/pseudobulks/",disease_status,'_',sampletype,".pb.txt")
  df.aggre1<-fread(f,data.table = F,stringsAsFactors = F,header = T)
  f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/cellbender_sub/HSC/pseudobulks/",disease_status,'_',sampletype,".meta.txt")
  meta1<-fread(f,data.table = F,stringsAsFactors = F,header = T)
  rownames(df.aggre1) = df.aggre1[,1]
  df.aggre1 = df.aggre1[,-1]
  rownames(meta1) = meta1[,1]
  meta1 = meta1[,-1]
  
  sampletype="Femur"
  f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/cellbender_sub/HSC/pseudobulks/",disease_status,'_',sampletype,".pb.txt")
  df.aggre2<-fread(f,data.table = F,stringsAsFactors = F,header = T)
  f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/cellbender_sub/HSC/pseudobulks/",disease_status,'_',sampletype,".meta.txt")
  meta2<-fread(f,data.table = F,stringsAsFactors = F,header = T)
  if (nrow(meta2)==1) {
    tmprow = df.aggre2[,1]
    tmpcol = colnames(df.aggre2)[2]
    df.aggre2 <- data.frame(sample=df.aggre2[,2])
    colnames(df.aggre2) = tmpcol
    rownames(df.aggre2) = tmprow
  } else {
    rownames(df.aggre2) = df.aggre2[,1]
    df.aggre2 = df.aggre2[,-1]
  }
  rownames(meta2) = meta2[,1]
  meta2 = meta2[,-1]
  
  # to account for same patients in femur and liver
  rownames(meta1) = paste0(rownames(meta1),".liver")
  rownames(meta2) = paste0(rownames(meta2),".femur")
  colnames(df.aggre1) = paste0(colnames(df.aggre1),".liver")
  colnames(df.aggre2) = paste0(colnames(df.aggre2),".femur")
  
  df.aggre = merge(df.aggre1,df.aggre2,by=0)
  rownames(df.aggre) = df.aggre[,1]
  df.aggre = df.aggre[,-1]
  meta = rbind(meta1,meta2)
  
  # align metadata and pseudobulks:
  metadata_to_use <- meta[rownames(meta) %in% colnames(df.aggre),]
  metadata_to_use$age = as.numeric(substring(metadata_to_use$age,1,2))
  df.aggre <- as.data.frame.matrix(df.aggre[,match(rownames(metadata_to_use),colnames(df.aggre))])
  
  # Specify parallel processing parameters
  # this is used implicitly by dream() to run in parallel
  param = SnowParam(8, "SOCK", progressbar=TRUE)
  
  # The variable to be tested must be a fixed effect
  # for (use_pcw in c(TRUE,FALSE)) {
  form <- ~ organ + sorting
  
  # A positive FC is increased expression in the DS compared to healthy
  metadata_to_use$organ <- factor(metadata_to_use$organ,levels=c("Femur","Liver"))
  if ("sorting" %in% colnames(metadata_to_use)) {
    metadata_to_use$sorting[!(metadata_to_use$sorting %in% c("CD235a-","CD45+"))] <- "Other"
  }
  
  # estimate weights using linear mixed model of dream
  # vobjDream = voomWithDreamWeights( df.aggre, form, metadata_to_use, BPPARAM=param )
  vobjDream = voomWithDreamWeights( df.aggre, form, metadata_to_use, BPPARAM=param ,span="auto")
  
  # Fit the model on each gene
  fitmm = dream( vobjDream, form, metadata_to_use )
  fitmm = eBayes(fitmm)
  
  # reorganize:
  res.df1 <- as.data.frame(topTable( fitmm, number=Inf,coef = "organLiver" ))
  # res.df1 <- as.data.frame(topTable( fitmm, number=Inf,coef = "age" ))
  res.df1$names <- rownames(res.df1)
  res.df1$disease_status <- disease_status
  
  dir.create(paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/","sample","/DE/cellbender"),showWarnings = FALSE)
  # f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/",subset_column,"/DE/",disease_status,'_',cell_type_filename,".de",age_suffix,".txt")
  f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/","sample","/DE/cellbender/",disease_status,'_',cell_type_filename,".de.txt")
  fwrite(res.df1,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
}

### Compare!




