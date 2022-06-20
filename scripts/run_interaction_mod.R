library(data.table)
library('variancePartition')
library('edgeR')
library('BiocParallel')

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

df.aggre.keep <- list()
metadata_to_use.keep <- list()

sampletype="Femur"
subset_column="sample"

clusters_for_DE <- c("Early erythroid cells","B cells",
                     "NK cells","cDC2"   ,                   
                     "Megakaryocytes","MEMPs"     ,                
                     "pDCs","Pro B cells"    ,           
                     "Pre pro B cells","Mast cells"     ,           
                     "Neutrophils","Cycling MEMPs"     ,        
                     "HSCs/MPPs","Inflammatory macrophages"  ,
                     "Late erythroid cells","Cycling pDCs"      )        

all_results_saved <- list()
do_next=FALSE; for (k in 1:length(clusters_for_DE))    {   
  cell_type <- clusters_for_DE[k]
  print(cell_type)
  # for (celltype in clusters_for_DE[2:length(clusters_for_DE)])    {                 
  for (sampletype in c("Liver","Femur")) {
    print(sampletype)
    
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
    # if (sampletype=="Femur") {clusters_for_DE <- unique(meta.all.full$leiden_names)[unique(meta.all.full$leiden_names) %in% unique(c(meta1.full$leiden_names,meta2.full$leiden_names))]}
    meta.all.full <- rbind(meta1.full,meta2.full)
    
    #
    print("Merge metadata...")
    x <- rbind(meta1,meta2)
    rownames(x) <- x[,subset_column]
    x$sorting[!(x$sorting %in% c("CD235a-","CD45+"))] <- "Other"
    
    print("cell types of interest...")
    cell_type_filename = gsub("/","_",cell_type)
    
    f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data_pb_leiden/",sampletype,".pb.",cell_type_filename,".",subset_column,".txt")
    if (!file.exists(f)) {do_next=TRUE;break}
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
    
    df.aggre.keep[[sampletype]] <- df.aggre
    metadata_to_use.keep[[sampletype]] <- metadata_to_use
    metadata_to_use.keep[[sampletype]]$sampletype <- sampletype
  }
  
  if (do_next) {do_next=FALSE;next}
  
  print(cell_type)
  
  # check
  # which(!(rownames(df.aggre.keep[["Liver"]]) %in% rownames(df.aggre.keep[["Femur"]])))
  # which(!(rownames(df.aggre.keep[["Femur"]]) %in% rownames(df.aggre.keep[["Liver"]])))
  
  df.aggre.keep.all <- do.call(cbind,df.aggre.keep)
  metadata_to_use.keep.all <- do.call(rbind,c(metadata_to_use.keep,make.row.names=FALSE)); rownames(metadata_to_use.keep.all) <- metadata_to_use.keep.all$sample
  
  # to keep consistent
  df.aggre <- df.aggre.keep.all
  metadata_to_use <- metadata_to_use.keep.all
  
  # Standard usage of limma/voom
  geneExpr = DGEList( df.aggre )
  keep <- filterByExpr(geneExpr, group=metadata_to_use$environment) | filterByExpr(geneExpr, group=metadata_to_use$sampletype)
  geneExpr <- geneExpr[keep,]
  
  geneExpr = calcNormFactors( geneExpr )
  
  # Specify parallel processing parameters
  # this is used implicitly by dream() to run in parallel
  param = SnowParam(8, "SOCK", progressbar=TRUE)
  
  form <- ~ sampletype + environment + sampletype*environment + sorting + (1|patient)
  
  metadata_to_use$environment <- factor(metadata_to_use$environment,levels=c("Healthy","DownSyndrome"))
  metadata_to_use$sampletype <- factor(metadata_to_use$sampletype,levels=c("Femur","Liver"))
  metadata_to_use$patient <- factor(metadata_to_use$patient)
  
  # estimate weights using linear mixed model of dream
  vobjDream = voomWithDreamWeights( geneExpr, form, metadata_to_use, BPPARAM=param )
  
  # Fit the dream model on each gene
  # By default, uses the Satterthwaite approximation for the hypothesis test
  fitmm = dream( vobjDream, form, metadata_to_use )
  fitmm = eBayes(fitmm)
  
  # reorganize:
  # res.df <- as.data.frame(topTable( fitmm, number=Inf))
  res.df1 <- as.data.frame(topTable( fitmm, number=Inf,coef = "environmentDownSyndrome" ))
  res.df1$names <- rownames(res.df1)
  # val1=mean(res.df$adj.P.Val < 0.05)
  
  # min(res.df$adj.P.Val)
  # res.df2 <- as.data.frame(topTable( fitmm, number=Inf,coef = "sampletypeLiver" ))
  # res.df2$names <- rownames(res.df2)
  # val2=mean(res.df$adj.P.Val < 0.05)
  
  # min(res.df$adj.P.Val)
  res.df <- as.data.frame(topTable( fitmm, number=Inf,coef = "sampletypeLiver:environmentDownSyndrome" ))
  res.df$names <- rownames(res.df)
  res.df.mg <- merge(res.df,res.df1,by="names")
  
  annot <- getBM(attributes = c('hgnc_symbol', 'chromosome_name',
                                'start_position', 'end_position'),
                 filters = 'hgnc_symbol',
                 values = res.df.mg$names,
                 mart = ensembl)
  res.df.mg <- merge(res.df.mg,annot[,c("hgnc_symbol","chromosome_name","start_position")],by.x="names",by.y="hgnc_symbol")
  
  f.out=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_int/",cell_type_filename,'.txt')
  fwrite(res.df.mg,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
}  


