library(data.table)
library('variancePartition')
library('edgeR')
library('BiocParallel')

# res.df.lst <- list()
df.aggre.lst <- list()
metadata_to_use.lst <- list()
sampletype="Liver"
subset_column="sample"

for (sampletype in c("Liver","Femur")) { 
  
  print(sampletype)
  
  # 
  print("Reading metadata1...")
  disease_status="Healthy"
  f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cellComp/10X_",disease_status,"_",sampletype,".cellComp.csv")
  meta1.full<-fread(f,data.table = F,stringsAsFactors = F)
  cells1 <- unique(meta1.full[,6])
  meta1 <- unique(meta1.full[,c("patient","sample","sorting")])
  meta1$environment <- disease_status
  
  #
  print("Reading metadata2...")
  disease_status="DownSyndrome"
  f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cellComp/10X_",disease_status,"_",sampletype,".cellComp.csv")
  meta2.full<-fread(f,data.table = F,stringsAsFactors = F)
  cells2 <- unique(meta2.full[,6])
  meta2 <- unique(meta2.full[,c("patient","sample","sorting")])
  meta2$environment <- disease_status
  
  colnames(meta1.full)[6] <- "leiden_names"
  colnames(meta2.full)[6] <- "leiden_names"
  meta.all.full <- rbind(meta1.full,meta2.full)
  
  #
  print("Merge metadata...")
  x <- rbind(meta1,meta2)
  rownames(x) <- x[,subset_column]
  x$sorting[!(x$sorting %in% c("CD235a-","CD45+"))] <- "Other"
  
  
  print("cell types of interest...")
  clusters_for_DE <- cells1[cells1 %in% cells2]
  P <- length(clusters_for_DE)
  
  cell_type="HSCs/MPPs"
  # cell_type="Megakaryocytes"
  # cell_type="Late erythroid cells"
  
  #
  
  cell_type_filename = gsub("/","_",cell_type)
  
  f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data_pb_leiden/",sampletype,".pb.",cell_type_filename,".",subset_column,".txt")
  df.aggre <- fread(f,data.table = F,stringsAsFactors = F)
  rownames(df.aggre) <- df.aggre[,1]; df.aggre <- df.aggre[,-1]
  
  metadata_to_use <- x[rownames(x) %in% colnames(df.aggre),]
  to_merge <- names(table(metadata_to_use$sorting)[table(metadata_to_use$sorting) < 5])
  df.aggre <- as.matrix(df.aggre[,match(rownames(metadata_to_use),colnames(df.aggre))])
  
  tab <- table(meta.all.full[meta.all.full[,"leiden_names"]==cell_type,'sample'])
  samples_to_keep <- names(tab)[tab >= 10]
  samples_to_keep <- samples_to_keep[samples_to_keep %in% rownames(metadata_to_use)]
  df.aggre <- df.aggre[,samples_to_keep]
  metadata_to_use <- metadata_to_use[samples_to_keep,]
  df.aggre.lst[[sampletype]] <- df.aggre
  metadata_to_use.lst[[sampletype]] <- metadata_to_use
}

# disease_status <- "Healthy"
disease_status <- "DownSyndrome"
for (disease_status in c("Healthy","DownSyndrome")) { 
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
  
  # Standard usage of limma/voom
  geneExpr = DGEList( df.aggre )
  keep <- filterByExpr(geneExpr, group=metadata_to_use$environment)
  geneExpr <- geneExpr[keep,]
  
  geneExpr = calcNormFactors( geneExpr )
  
  # Specify parallel processing parameters
  # this is used implicitly by dream() to run in parallel
  param = SnowParam(8, "SOCK", progressbar=TRUE)
  
  form <- ~ sampletype + patient + sorting
  if (length(unique(metadata_to_use$sorting))==1) {
    form <- ~ sampletype + patient
  } 
  if (max(table(metadata_to_use$patient))<1.5) {
    form <- ~ sampletype + sorting
  }
  if (max(table(metadata_to_use$patient))<1.5 & length(unique(metadata_to_use$sorting))==1) {
    form <- ~ sampletype
  }
  
  # form <- ~ environment
  
  # A positive FC is increased expression in the liver compared to femur
  metadata_to_use$sampletype <- factor(metadata_to_use$sampletype,levels=c("Femur","Liver"))
  
  # estimate weights using linear mixed model of dream
  vobjDream = voomWithDreamWeights( geneExpr, form, metadata_to_use, BPPARAM=param )
  
  # Fit the dream model on each gene
  # By default, uses the Satterthwaite approximation for the hypothesis test
  fitmm = dream( vobjDream, form, metadata_to_use )
  fitmm = eBayes(fitmm)
  
  # reorganize:
  # res.df <- as.data.frame(topTable( fitmm, number=Inf))
  res.df <- as.data.frame(topTable( fitmm, number=Inf,coef = "sampletypeLiver" ))
  res.df$names <- rownames(res.df)
  res.df$sampletype <- sampletype
  
  print(head(res.df))
  print(subset(res.df,names=="GATA1"))
  
  res.df.all <- res.df
  annot <- getBM(attributes = c('hgnc_symbol', 'chromosome_name',
                                'start_position', 'end_position'),
                 filters = 'hgnc_symbol', 
                 values = res.df.all$names, 
                 mart = ensembl)
  df1 <- merge(res.df.all,annot[,c("hgnc_symbol","chromosome_name","start_position")],by.x="names",by.y="hgnc_symbol")
  df1 <- df1[df1$chromosome_name %in% c(as.character(seq(1,22)),"X","Y"),]
  # df1$chromosome_name <- as.factor(df1$chromosome_name)
  df1$chr21 <- factor(ifelse(df1$chromosome_name==21,'Chr 21','Not Chr 21'),levels=c('Not Chr 21','Chr 21'))
  res.df.all <- df1
  res.df.all[order(res.df.all$P.Value)[1:5],]
  subset(res.df.all,names=="GATA1")
  # head(res.df.all)
  
  print(aggregate(res.df.all$adj.P.Val < 0.1,by=list(res.df.all$chromosome_name==21),mean))
}


