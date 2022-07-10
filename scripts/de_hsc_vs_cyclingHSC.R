library(data.table)
library(Matrix.utils)
library(Seurat)
library(DESeq2)
library(variancePartition)
library('edgeR')

print("Reading metadata2...")
disease_status="DownSyndrome"
f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cellComp/10X_","DownSyndrome","_","Liver",".cellComp.csv")
meta2.full<-fread(f,data.table = F,stringsAsFactors = F)
colnames(meta2.full)[6] <- "leiden_names"
cells2 <- unique(meta2.full[,6])
meta2 <- unique(meta2.full[,c("patient","sample","sorting")])
meta2$environment <- disease_status
meta2$sorting[!(meta2$sorting %in% c("CD235a-","CD45+"))] <- "Other"
x <- meta2
rownames(x) <- x[,"sample"]
rm(meta2)

########

iter=5
cell_type="Cycling HSCs/MPPs";
cell_type_filename = gsub("/","_",cell_type)
disease_status="DownSyndrome"
sampletype="Liver"
f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/cell_subset/",disease_status,'_',sampletype,'_',cell_type_filename,".rds")
df <- readRDS(file = f.out)

df.aggre <- aggregate.Matrix(
  t(
    GetAssayData(object = df, slot = "counts", assay="RNA")
  ),
  groupings=df@meta.data[,"sample"],fun="sum")
df.aggre <- t(df.aggre)
df.aggre <- as.data.frame(df.aggre)

metadata_to_use <- x[rownames(x) %in% colnames(df.aggre),]
df.aggre <- as.matrix(df.aggre[,match(rownames(metadata_to_use),colnames(df.aggre))])

tab <- table(meta2.full[meta2.full[,"leiden_names"]==cell_type,'sample'])
samples_to_keep <- names(tab)[tab >= 10]
samples_to_keep <- samples_to_keep[samples_to_keep %in% rownames(metadata_to_use)]
df.aggre <- df.aggre[,samples_to_keep]
metadata_to_use <- metadata_to_use[samples_to_keep,]
df.aggre1=df.aggre
metadata_to_use1 <- metadata_to_use
metadata_to_use1$leiden_names=cell_type
colnames(df.aggre1) <- paste0(colnames(df.aggre1),".cyc")
rownames(metadata_to_use1) <- paste0(rownames(metadata_to_use1),".cyc")

##### 

cell_type="HSCs/MPPs"
cell_type_filename = gsub("/","_",cell_type)
f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data_pb_leiden/","Liver",".pb.",cell_type_filename,".","sample",".txt")

df.aggre <- fread(f,data.table = F,stringsAsFactors = F)
rownames(df.aggre) <- df.aggre[,1]; df.aggre <- df.aggre[,-1]

metadata_to_use <- x[rownames(x) %in% colnames(df.aggre),]
to_merge <- names(table(metadata_to_use$sorting)[table(metadata_to_use$sorting) < 5])
df.aggre <- as.matrix(df.aggre[,match(rownames(metadata_to_use),colnames(df.aggre))])

tab <- table(meta2.full[meta2.full[,"leiden_names"]==cell_type,'sample'])
samples_to_keep <- names(tab)[tab >= 10]
samples_to_keep <- samples_to_keep[samples_to_keep %in% rownames(metadata_to_use)]
df.aggre <- df.aggre[,samples_to_keep]
metadata_to_use <- metadata_to_use[samples_to_keep,]
colnames(df.aggre) <- paste0(colnames(df.aggre),".hsc")
rownames(metadata_to_use) <- paste0(rownames(metadata_to_use),".hsc")
metadata_to_use$leiden_names <- cell_type
df.aggre <- cbind(df.aggre[rownames(df.aggre) %in% rownames(df.aggre1),],df.aggre1)
metadata_to_use <- rbind(metadata_to_use,metadata_to_use1)
##### 

geneExpr = DGEList( df.aggre )
keep <- filterByExpr(geneExpr, group=metadata_to_use$environment)
geneExpr <- geneExpr[keep,]
geneExpr = calcNormFactors( geneExpr )

# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param = SnowParam(8, "SOCK", progressbar=TRUE)

# The variable to be tested must be a fixed effect
form <- ~ leiden_names + sorting #+ patient

# A positive FC is increased expression in the DS compared to healthy
metadata_to_use$leiden_names <- factor(metadata_to_use$leiden_names,levels=c("HSCs/MPPs","Cycling HSCs/MPPs"))

# estimate weights using linear mixed model of dream
vobjDream = voomWithDreamWeights( geneExpr, form, metadata_to_use, BPPARAM=param )

# Fit the dream model on each gene
# By default, uses the Satterthwaite approximation for the hypothesis test
fitmm = dream( vobjDream, form, metadata_to_use )
fitmm = eBayes(fitmm)

# reorganize:
# res.df <- as.data.frame(topTable( fitmm, number=Inf))
res.df <- as.data.frame(topTable( fitmm, number=Inf,coef = "leiden_namesCycling HSCs/MPPs" ))
res.df$names <- rownames(res.df)

# save:
f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/HSC_v_CyclingHSC.txt")
fwrite(res.df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)


