print("Limma/voom...")

f.out=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_liver_v_femur/",cell_type_filename,'.',"DownSyndrome",'.FE.txt')
full_diff_exp_results <- fread(f.out,data.table = F,stringsAsFactors = F)

# Standard usage of limma/voom
# geneExpr = DGEList( df.aggre )
# keep <- filterByExpr(geneExpr, group=metadata_to_use$environment)
# geneExpr <- geneExpr[keep,]
geneExpr = DGEList( df.aggre[full_diff_exp_results$names,] )

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
colnames(res.df1)[colnames(res.df1)=="ID"] <- "names"
res.df1$disease_status <- disease_status

# print(head(res.df1))
print(subset(res.df1,names%in%c("GATA1","APOC1","MYL4")))
mean(res.df1$adj.P.Val < 0.1)

# print(subset(res.df1,names=="GATA1"))

df1 <- merge(res.df1,geneAnnot[,c("hgnc_symbol","chromosome_name","start_position")],by.x="names",by.y="hgnc_symbol")
df1$chr21 <- factor(ifelse(df1$chromosome_name=='21','Chr 21','Not Chr 21'),levels=c('Not Chr 21','Chr 21'))

f.out=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_liver_v_femur/",cell_type_filename,'.',disease_status,'.FE.keep_t21_genes.txt')
fwrite(df1,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
