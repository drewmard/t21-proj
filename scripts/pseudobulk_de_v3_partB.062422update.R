library(data.table)
library('variancePartition')
library('edgeR')
library('BiocParallel')

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
  iter=0; for (cell_type in clusters_for_DE) {
  # for (cell_type in clusters_for_DE[c(19,25:length(clusters_for_DE))]) {
    iter = iter + 1
    print(paste0(iter,"/",P,": ",cell_type))
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
    
    # saveRDS(list(df.aggre,metadata_to_use),"/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data_pb_leiden/Liver.pb.HSCs_MPPs.txt")
    
    
    # need to remodel sorting probably:
    # [1] "19/28: Schwann cells"
    # Fixed effect model, using limma directly...
    # Error in `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]) : 
    #   contrasts can be applied only to factors with 2 or more levels
    # Calls: voomWithDreamWeights ... model.matrix -> model.matrix.default -> contrasts<-
      
    # Standard usage of limma/voom
    geneExpr = DGEList( df.aggre )
    keep <- filterByExpr(geneExpr, group=metadata_to_use$environment)
    geneExpr <- geneExpr[keep,]
    
    geneExpr = calcNormFactors( geneExpr )
    
    # Specify parallel processing parameters
    # this is used implicitly by dream() to run in parallel
    param = SnowParam(8, "SOCK", progressbar=TRUE)
    
    # The variable to be tested must be a fixed effect
    form <- ~ environment + sorting #+ patient
    if (length(unique(metadata_to_use$sorting))==1 | length(unique(metadata_to_use$sorting)) > nrow(metadata_to_use)/2) {
      form <- ~ environment
    }
    
    # form <- ~ environment + sorting + (1|patient)
    # if (length(unique(metadata_to_use$sorting))==1) {
    #   form <- ~ environment + (1|patient)
    # }
    # if (median(table(metadata_to_use$patient))<1.5) {
    #   form <- ~ environment + patient
    # }
    # if (max(table(metadata_to_use$patient))==1) {
    #   form <- ~ environment
    # }
    # 
    # form <- ~ environment + patient + sorting
    # if (length(unique(metadata_to_use$sorting))==1) {
    #   form <- ~ environment + patient
    # }
    # 
    # A positive FC is increased expression in the DS compared to healthy
    metadata_to_use$environment <- factor(metadata_to_use$environment,levels=c("Healthy","DownSyndrome"))
    
    # estimate weights using linear mixed model of dream
    vobjDream = voomWithDreamWeights( geneExpr, form, metadata_to_use, BPPARAM=param )
    
    # Fit the dream model on each gene
    # By default, uses the Satterthwaite approximation for the hypothesis test
    fitmm = dream( vobjDream, form, metadata_to_use )
    fitmm = eBayes(fitmm)
    
    # reorganize:
    # res.df <- as.data.frame(topTable( fitmm, number=Inf))
    res.df <- as.data.frame(topTable( fitmm, number=Inf,coef = "environmentDownSyndrome" ))
    res.df$names <- rownames(res.df)
    
    mean(res.df$adj.P.Val < 0.05)
    # save:
    system("mkdir -p /oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names")
    # f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/",sampletype,".",cell_type_filename,".",subset_column,".random.txt")
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
