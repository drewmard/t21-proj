library(data.table)
library('variancePartition')
library('edgeR')
library('BiocParallel')

sampletype="Liver"
subset_column="sample"
disease_status="Healthy"
useRandom=FALSE
geneAnnot=fread("/oak/stanford/groups/smontgom/amarder/t21-proj/out/data_small/gene_annot.txt",data.table = F,stringsAsFactors = F)

# for (disease_status in c("DownSyndrome")) {
for (disease_status in c("Healthy","DownSyndrome")) {
  print(disease_status)
  
  print("Reading metadata1...")
  sampletype="Liver"
  f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cellComp/10X_",disease_status,"_",sampletype,".cellComp.csv")
  meta1.full<-fread(f,data.table = F,stringsAsFactors = F)
  cells1 <- unique(meta1.full[,6])
  meta1 <- unique(meta1.full[,c("patient","sample","sorting")])
  meta1$environment <- disease_status
  meta1.full$environment <- disease_status
  meta1$sampletype <- sampletype
  meta1.full$sampletype <- sampletype
  
  print("Reading metadata2...")
  sampletype="Femur"
  f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cellComp/10X_",disease_status,"_",sampletype,".cellComp.csv")
  meta2.full<-fread(f,data.table = F,stringsAsFactors = F)
  cells2 <- unique(meta2.full[,6])
  meta2 <- unique(meta2.full[,c("patient","sample","sorting")])
  meta2$environment <- disease_status
  meta2.full$environment <- disease_status
  meta2$sampletype <- sampletype
  meta2.full$sampletype <- sampletype
  
  colnames(meta1.full)[6] <- "leiden_names"
  colnames(meta2.full)[6] <- "leiden_names"
  meta.all.full <- rbind(meta1.full,meta2.full)
  
  x <- rbind(meta1,meta2)
  rownames(x) <- x[,subset_column]
  x$sorting[!(x$sorting %in% c("CD235a-","CD45+"))] <- "Other"
  
  print("cell types of interest...")
  clusters_for_DE <- cells1[cells1 %in% cells2]
  
  # consider only those that are in all 4 datasets (t21/healthy x femur/liver)

  # cell_type="HSCs/MPPs"
  # cell_type=clusters_for_DE[2]
  for (cell_type in clusters_for_DE) {
  # for (cell_type in clusters_for_DE[15:length(clusters_for_DE)]) {
    print(cell_type)
    # table(meta.all.full[meta.all.full[,"leiden_names"]==cell_type,c('sampletype',"environment")])
    df.aggre.lst <- list()
    metadata_to_use.lst <- list()
    theStopper=FALSE; for (sampletype in c("Liver","Femur")) {
      cell_type_filename = gsub("/","_",cell_type)
      
      f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data_pb_leiden/",sampletype,".pb.",cell_type_filename,".",subset_column,".txt")
      if (!file.exists(f)) {print(paste0("Skipping ",cell_type,"...")); theStopper=TRUE; break}

      df.aggre <- fread(f,data.table = F,stringsAsFactors = F)
      rownames(df.aggre) <- df.aggre[,1]; df.aggre <- df.aggre[,-1]
      
      metadata_to_use <- x[rownames(x) %in% colnames(df.aggre),]
      to_merge <- names(table(metadata_to_use$sorting)[table(metadata_to_use$sorting) < 5])
      df.aggre <- as.matrix(df.aggre[,match(rownames(metadata_to_use),colnames(df.aggre))])
      
      tab <- table(meta.all.full[meta.all.full[,"leiden_names"]==cell_type,'sample'])
      samples_to_keep <- names(tab)[tab >= 10]
      samples_to_keep <- samples_to_keep[samples_to_keep %in% rownames(metadata_to_use)]
      df.aggre <- df.aggre[,samples_to_keep,drop=FALSE]
      metadata_to_use <- metadata_to_use[samples_to_keep,]
      df.aggre.lst[[sampletype]] <- df.aggre
      metadata_to_use.lst[[sampletype]] <- metadata_to_use
    }
    
    if (theStopper) {print("Skipped."); next}
    
    sampletype="Liver"
    i1 <- which(metadata_to_use.lst[[sampletype]]$environment==disease_status)
    metadata_to_use1 <- metadata_to_use.lst[[sampletype]][i1,]
    df.aggre1 <- df.aggre.lst[[sampletype]][,i1,drop=FALSE]
    metadata_to_use1$sampletype <- sampletype
    sampletype="Femur"
    i2 <- which(metadata_to_use.lst[[sampletype]]$environment==disease_status)
    metadata_to_use2 <- metadata_to_use.lst[[sampletype]][i2,]
    df.aggre2 <- df.aggre.lst[[sampletype]][,i2,drop=FALSE]
    metadata_to_use2$sampletype <- sampletype
    
    metadata_to_use <- rbind(metadata_to_use1,metadata_to_use2)
    df.aggre <- cbind(df.aggre1,df.aggre2)
    df.aggre <- df.aggre[rownames(df.aggre) %in% geneAnnot$hgnc_symbol,]
    
    print("Limma/voom...")
    
    # Standard usage of limma/voom
    geneExpr = DGEList( df.aggre )
    keep <- filterByExpr(geneExpr, group=metadata_to_use$environment)
    geneExpr <- geneExpr[keep,]
    
    geneExpr = calcNormFactors( geneExpr )
    
    # Specify parallel processing parameters
    # this is used implicitly by dream() to run in parallel
    param = SnowParam(8, "SOCK", progressbar=TRUE)
    
    # determining the model:
    if (useRandom) {
      form <- ~ sampletype + sorting + (1|patient)
      if (length(unique(metadata_to_use$sorting))==1) {
        form <- ~ sampletype + (1|patient)
      } 
    } else {
      form <- ~ sampletype + sorting #+ patient
    }
    if (max(table(metadata_to_use$patient))<1.5 | nrow(metadata_to_use) < 5) {
      form <- ~ sampletype + sorting
    }
    if ((max(table(metadata_to_use$patient))<1.5 | nrow(metadata_to_use) < 5) & length(unique(metadata_to_use$sorting))==1) {
      form <- ~ sampletype
    }
    if (!useRandom & length(unique(metadata_to_use$sorting))==1) {
      form <- ~ sampletype
    }
    
    metadata_to_use$environment <- factor(metadata_to_use$environment,levels=c("Healthy","DownSyndrome"))
    metadata_to_use$sampletype <- factor(metadata_to_use$sampletype,levels=c("Femur","Liver"))
    metadata_to_use$patient <- factor(metadata_to_use$patient)
    
    # estimate weights using linear mixed model of dream
    vobjDream = voomWithDreamWeights( geneExpr, form, metadata_to_use, BPPARAM=param )
    
    # Fit the dream model on each gene
    # By default, uses the Satterthwaite approximation for the hypothesis test
    fitmm = dream( vobjDream, form, metadata_to_use , BPPARAM = param)
    fitmm = eBayes(fitmm)
    
    # reorganize:
    res.df1 <- as.data.frame(topTable( fitmm, number=Inf,coef = "sampletypeLiver" ))
    res.df1$names <- rownames(res.df1)
    res.df1$disease_status <- disease_status

    # print(head(res.df1))
    # print(subset(res.df1,names%in%c("GATA1","APOC1","MYL4")))
    # mean(res.df1$adj.P.Val < 0.1)
    
    # print(subset(res.df1,names=="GATA1"))
    
    df1 <- merge(res.df1,geneAnnot[,c("hgnc_symbol","chromosome_name","start_position")],by.x="names",by.y="hgnc_symbol")
    # df1$chromosome_name <- factor(df1$chromosome_name,levels=seq(1,22))
    df1$chr21 <- factor(ifelse(df1$chromosome_name=='21','Chr 21','Not Chr 21'),levels=c('Not Chr 21','Chr 21'))
    
    f.out=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_liver_v_femur/",cell_type_filename,'.',disease_status,'.txt')
    if (!useRandom) {f.out=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_liver_v_femur/",cell_type_filename,'.',disease_status,'.FE.txt')}
    fwrite(df1,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
  }
}