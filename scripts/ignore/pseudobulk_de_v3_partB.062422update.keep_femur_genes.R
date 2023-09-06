library(data.table)
library('variancePartition')
library('edgeR')
library('BiocParallel')

cell_type="HSCs/MPPs"

dir <- '/oak/stanford/groups/smontgom/amarder'
disease_status="DownSyndrome"
sampletype="Liver"
datadir="DE_pb_leiden_names";
cell_type_filename = gsub("/","_",cell_type)
f<- paste0(dir,"/t21-proj/out/full/",datadir,"/",sampletype,".",cell_type_filename,".sample.txt")
# f.out=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_liver_v_femur/",cell_type_filename,'.',disease_status,'.FE.txt')
full_diff_exp_results <- fread(f,data.table = F,stringsAsFactors = F)

sampletype="Femur"
subset_column="sample"

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
iter=0; 
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

# need to remodel sorting probably:
# [1] "19/28: Schwann cells"
# Fixed effect model, using limma directly...
# Error in `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]) : 
#   contrasts can be applied only to factors with 2 or more levels
# Calls: voomWithDreamWeights ... model.matrix -> model.matrix.default -> contrasts<-

# Standard usage of limma/voom
# geneExpr = DGEList( df.aggre )
# keep <- filterByExpr(geneExpr, group=metadata_to_use$environment)
# geneExpr <- geneExpr[keep,]
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
f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/",sampletype,".",cell_type_filename,".",subset_column,".keep_t21_genes.txt")
fwrite(res.df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
