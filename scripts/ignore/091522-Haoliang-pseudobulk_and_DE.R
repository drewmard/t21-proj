# here, dfsub is a scRNA-seq Seurat object containing only femur HSCs (both Healthy and T21) 

df.aggre <- aggregate.Matrix(
  t(
    GetAssayData(object = dfsub, slot = "counts", assay="RNA")
  ),
  groupings=dfsub$sample,fun="sum")

df.aggre <- t(df.aggre)
df.aggre <- as.data.frame(df.aggre)

x <- unique(dfsub@meta.data[,c("patient","sample","sorting","environment","organ")])
rownames(x) <- x[,"sample"]

metadata_to_use <- x[rownames(x) %in% colnames(df.aggre),]
df.aggre <- as.matrix(df.aggre[,match(rownames(metadata_to_use),colnames(df.aggre))])

# Standard usage of limma/voom
geneExpr = DGEList( df.aggre[full_diff_exp_results$names,] )
geneExpr = calcNormFactors( geneExpr )

# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param = SnowParam(8, "SOCK", progressbar=TRUE)

# The variable to be tested must be a fixed effect
form <- ~ environment + sorting #+ patient
if (length(unique(metadata_to_use$sorting))==1 | length(unique(metadata_to_use$sorting)) > nrow(metadata_to_use)/2) {
  form <- ~ environment
}

# A positive FC is increased expression in the DS compared to healthy
metadata_to_use$environment <- factor(metadata_to_use$environment,levels=c("Healthy","Down Syndrome"))
metadata_to_use$sorting[!(metadata_to_use$sorting %in% c("CD235a-","CD45+"))] <- "Other"

# estimate weights using linear mixed model of dream
vobjDream = voomWithDreamWeights( geneExpr, form, metadata_to_use, BPPARAM=param )

# Fit the model on each gene
fitmm = dream( vobjDream, form, metadata_to_use )
fitmm = eBayes(fitmm)

# reorganize:
res.df1 <- as.data.frame(topTable( fitmm, number=Inf,coef = "environmentDown Syndrome" ))
res.df1$names <- rownames(res.df1)
res.df1$disease_status <- disease_status

df1 <- merge(res.df1,geneAnnot[,c("hgnc_symbol","chromosome_name","start_position")],by.x="names",by.y="hgnc_symbol")
df1$chr21 <- factor(ifelse(df1$chromosome_name=='21','Chr 21','Not Chr 21'),levels=c('Not Chr 21','Chr 21'))

# save:
subset_column='sample'
system("mkdir -p /oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names_downsamp")
f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names_downsamp/","Liver",".",cell_type_filename,".",subset_column,".perm",permNum,".txt")
fwrite(df1,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
