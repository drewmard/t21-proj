library(data.table)
library('variancePartition')
library('edgeR')
library('BiocParallel')

# sampletype="Liver"
subset_column="sample"

for (sampletype in c("Liver","Femur")) {
  
  print(sampletype)
  
  f = "/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cellComp/10X_Healthy_Liver.cellComp.csv"
  meta1<-fread(f,data.table = F,stringsAsFactors = F)
  healthy_cells <- unique(meta1[,6])
  meta1 <- unique(meta1[,c("patient","sample")])
  meta1$environment <- "Healthy"
  f = "/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cellComp/10X_DownSyndrome_Liver.cellComp.csv"
  meta2<-fread(f,data.table = F,stringsAsFactors = F)
  ds_cells <- unique(meta2[,6])
  meta2 <- unique(meta2[,c("patient","sample")])
  meta2$environment <- "DownSyndrome"
  x <- rbind(meta1,meta2)
  rownames(x) <- x[,subset_column]
  clusters_for_DE <- healthy_cells[healthy_cells %in% ds_cells]
  P <- length(clusters_for_DE)
  
  iter=0; for (cell_type in clusters_for_DE) {
    iter = iter + 1
    print(paste0(iter,"/",P,": ",cell_type))
    cell_type_filename = gsub("/","_",cell_type)
    
    f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data_pb_leiden/",sampletype,".pb.",cell_type_filename,".",subset_column,".txt")
    df.aggre <- fread(f,data.table = F,stringsAsFactors = F)
    rownames(df.aggre) <- df.aggre[,1]; df.aggre <- df.aggre[,-1]
    
    df.aggre <- as.matrix(df.aggre[,match(rownames(x),colnames(df.aggre))])
    
    # Standard usage of limma/voom
    geneExpr = DGEList( df.aggre )
    geneExpr = calcNormFactors( geneExpr )
    
    # Specify parallel processing parameters
    # this is used implicitly by dream() to run in parallel
    param = SnowParam(8, "SOCK", progressbar=TRUE)
    
    # The variable to be tested must be a fixed effect
    form <- ~ environment + (1|patient) 
    
    # estimate weights using linear mixed model of dream
    vobjDream = voomWithDreamWeights( geneExpr, form, x, BPPARAM=param )
    
    # A positive FC is increased expression in the DS compared to healthy
    x$environment <- factor(x$environment,levels=c("Healthy","DownSyndrome"))
    
    # Fit the dream model on each gene
    # By default, uses the Satterthwaite approximation for the hypothesis test
    fitmm = dream( vobjDream, form, x )
    
    res.df <- as.data.frame(topTable( fitmm, number=Inf ))
    res.df$names <- rownames(res.df)
    
    system("mkdir -p /oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names")
    f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/",sampletype,".",cell_type_filename,".",subset_column,".txt")
    fwrite(res.df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
    
  }
}

# cell_type_filename="Early erythroid cells"
# 
# f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data_pb_leiden/",sampletype,".pb.",cell_type_filename,".",subset_column,".txt")
# df.aggre <- fread(f,data.table = F,stringsAsFactors = F)
# rownames(df.aggre) <- df.aggre[,1]; df.aggre <- df.aggre[,-1]
# 
# df.aggre <- as.matrix(df.aggre[,match(rownames(x),colnames(df.aggre))])
# 
# # Standard usage of limma/voom
# geneExpr = DGEList( df.aggre )
# geneExpr = calcNormFactors( geneExpr )
# 
# # Specify parallel processing parameters
# # this is used implicitly by dream() to run in parallel
# param = SnowParam(4, "SOCK", progressbar=TRUE)
# 
# # The variable to be tested must be a fixed effect
# form <- ~ environment + (1|patient) 
# 
# # estimate weights using linear mixed model of dream
# vobjDream = voomWithDreamWeights( geneExpr, form, x, BPPARAM=param )
# 
# # A positive FC is increased expression in the DS compared to healthy
# x$environment <- factor(x$environment,levels=c("Healthy","DownSyndrome"))
# 
# # Fit the dream model on each gene
# # By default, uses the Satterthwaite approximation for the hypothesis test
# fitmm = dream( vobjDream, form, x )
# 
# res.df <- as.data.frame(topTable( fitmm, number=Inf ))
# res.df$names <- rownames(res.df)
# 
# 
# library(biomaRt)
# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
# 
# annot <- getBM(attributes = c('hgnc_symbol', 'chromosome_name',
#                               'start_position', 'end_position'),
#                filters = 'hgnc_symbol', 
#                values = res.df$names, 
#                mart = ensembl)
# df1 <- merge(res.df,annot[,c("hgnc_symbol","chromosome_name","start_position")],by.x="names",by.y="hgnc_symbol")
# df1 <- df1[df1$chromosome_name %in% seq(1,22),]
# df1$chromosome_name <- factor(df1$chromosome_name,levels=seq(1,22))
# df1$chr21 <- factor(ifelse(df1$chromosome_name==21,'Chr 21','Not Chr 21'),levels=c('Not Chr 21','Chr 21'))
# df1.p <- df1
# aggregate(df1[,"logFC"]<0.05,by=list(df1$chr21),mean,na.rm=T)
# aggregate(rank(df1[,"logFC"])/nrow(df1),by=list(df1$chr21),mean,na.rm=T)
# aggregate(df1[,"adj.P.Val"]<0.05,by=list(df1$chr21),mean,na.rm=T)
# 
# sampletype <- "Liver"
# 
# dir <- '/oak/stanford/groups/smontgom/amarder'
# # dir <- '~/Documents/Research'
# 
# # datadir="DE_pb_cell_type_groups"; lfc.col="log2FoldChange"; p.col="padj"
# datadir="DE_pb_leiden_names"; lfc.col="log2FoldChange"; p.col="padj"
# # datadir="DE_cell_type_groups"; lfc.col="logfoldchanges"; p.col="pvals_adj"
# # datadir="DE_leiden_names"; lfc.col="logfoldchanges"; p.col="pvals_adj"
# 
# pathtodir <- paste0(dir,"/t21-proj/out/full/",datadir)
# 
# cell_type=cell_type_filename
# f<- paste0(dir,"/t21-proj/out/full/",datadir,"/",sampletype,".",cell_type,".txt")
# df <- fread(f,data.table = F,stringsAsFactors = F)
# 
# df.mg <- merge(df1,df,by='names')
# 
# cor(df.mg$logFC,df.mg$log2FoldChange,use='complete.obs')
# table(df.mg$adj.P.Val<0.05,df.mg$padj<0.05)
# 
# 
# 
# 
# 
