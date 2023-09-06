library(data.table)
# library(Seurat)
library(Matrix.utils) # for pseudobulking
library(BiocParallel) # for parallel
library(variancePartition) # for limma voom

disease_status="Healthy"
sampletype="Liver"
subset_column="sample"
print(paste0("Reading data: ",disease_status," ",subset_column,"..."))
# cell_type = "HSCs/MPPs"
cell_type = "Cycling HSC"
cell_type_filename = gsub("/","_",cell_type)
print("Reading metadata2...")
f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cellComp/10X_",disease_status,"_",sampletype,".cellComp.csv")
meta2.full<-fread(f,data.table = F,stringsAsFactors = F)
colnames(meta2.full)[6] <- "leiden_names"
meta2.full$cellname = unlist(lapply(strsplit(meta2.full$V1,"-"),function(x) paste(x[1:2],collapse = "-")))
integrated_labels = fread("/oak/stanford/groups/smontgom/amarder/t21-proj/out/a3_tab.obs.withannot.csv",data.table = F,stringsAsFactors = F)
integrated_labels$cellname=unlist(lapply(strsplit(integrated_labels$barcodekey,"-"),function(x) paste(x[1:2],collapse = "-")))
integrated_labels$sample=unlist(lapply(strsplit(as.character(integrated_labels$patient_sample)," "),function(x) paste(x[2],collapse = "-")))
integrated_labels = integrated_labels[!duplicated(integrated_labels[,c("cellname","sample","combi_annot")]),] # why??
nrow(meta2.full)
meta2.full = merge(meta2.full,integrated_labels,by=c("sample","cellname"),all.x=TRUE)
nrow(meta2.full)

cluster_label="combi_annot"
cell_type = "Cycling HSC"
suffix = ".integ1"
# cluster_label="leiden_latest"
# cell_type = "Cycling HSCs/MPPs"
# cell_type_filename = gsub("/","_",cell_type)
# suffix = ""
f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/",subset_column,"/pb/",disease_status,'_',sampletype,'_',cell_type_filename,".pb",suffix,".txt")
df.aggre1<-fread(f,data.table = F,stringsAsFactors = F,header = T)
f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/",subset_column,"/meta/",disease_status,'_',sampletype,'_',cell_type_filename,".meta",suffix,".txt")
meta1<-fread(f,data.table = F,stringsAsFactors = F,header = T)
rownames(df.aggre1) = df.aggre1[,1]
df.aggre1 = df.aggre1[,-1]
rownames(meta1) = meta1[,1]
meta1 = meta1[,-1]

tab <- table(meta2.full[meta2.full[,cluster_label]==cell_type,'sample'])
samples_to_keep <- names(tab)[tab >= 10]
df.aggre1 <- df.aggre1[,samples_to_keep]
meta1 <- meta1[samples_to_keep,]

cell_type = "HSCs"
cell_type_filename = "HSCs_MPPs"
cluster_label="combi_annot"
suffix = ".integ1"
# cluster_label="leiden_latest"
# cell_type = "HSCs/MPPs"
# cell_type_filename = gsub("/","_",cell_type)
# suffix = ""
# cell_type="Cycling HSCs/MPPs";
# cell_type_filename = gsub("/","_",cell_type)
# f.out=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data_pb_leiden/Liver.pb.",cell_type_filename,".txt")
# df.aggre2=fread(f.out,data.table = F,stringsAsFactors = F)
# f.out=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data_pb_leiden/Liver.pb.",cell_type_filename,".meta.txt")
# meta2=fread(f.out,data.table = F,stringsAsFactors = F)
f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/",subset_column,"/pb/",disease_status,'_',sampletype,'_',cell_type_filename,".pb",suffix,".txt")
df.aggre2<-fread(f,data.table = F,stringsAsFactors = F,header = T)
f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/",subset_column,"/meta/",disease_status,'_',sampletype,'_',cell_type_filename,".meta",suffix,".txt")
meta2<-fread(f,data.table = F,stringsAsFactors = F,header = T)
rownames(df.aggre2) = df.aggre2[,1]
df.aggre2 = df.aggre2[,-1]
rownames(meta2) = meta2[,1]
meta2 = meta2[,-1]

tab <- table(meta2.full[meta2.full[,cluster_label]==cell_type,'sample'])
samples_to_keep <- names(tab)[tab >= 10]
df.aggre2 <- df.aggre2[,samples_to_keep]
meta2 <- meta2[samples_to_keep,]

# to account for same patients in femur and liver
rownames(meta1) = paste0(rownames(meta1),".cyc")
rownames(meta2) = paste0(rownames(meta2),".hsc")
colnames(df.aggre1) = paste0(colnames(df.aggre1),".cyc")
colnames(df.aggre2) = paste0(colnames(df.aggre2),".hsc")

df.aggre = merge(df.aggre1,df.aggre2,by=0)
rownames(df.aggre) = df.aggre[,1]
df.aggre = df.aggre[,-1]
meta1$celltype = "cyc"
meta2$celltype = "hsc"
meta = rbind(meta1,meta2)

# align metadata and pseudobulks:
metadata_to_use <- meta[rownames(meta) %in% colnames(df.aggre),]
# metadata_to_use$age = as.numeric(substring(metadata_to_use$age,1,2))
df.aggre <- as.data.frame.matrix(df.aggre[,match(rownames(metadata_to_use),colnames(df.aggre))])

# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param = SnowParam(8, "SOCK", progressbar=TRUE)

# The variable to be tested must be a fixed effect
# for (use_pcw in c(TRUE,FALSE)) {
use_pcw = FALSE

form <- ~ celltype + sorting

# A positive FC is increased expression in the DS compared to healthy
metadata_to_use$celltype <- factor(metadata_to_use$celltype,levels=c("hsc","cyc"))
if ("sorting" %in% colnames(metadata_to_use)) {
  metadata_to_use$sorting[!(metadata_to_use$sorting %in% c("CD235a-","CD45+"))] <- "Other"
}

# preprocessing:
library('edgeR')
geneExpr = DGEList( df.aggre )
keep <- filterByExpr(geneExpr, group=metadata_to_use$environment)
geneExpr <- geneExpr[keep,]
geneExpr = calcNormFactors( geneExpr )

# estimate weights using linear mixed model of dream
vobjDream = voomWithDreamWeights( geneExpr, form, metadata_to_use, BPPARAM=param )

# Fit the model on each gene
fitmm = dream( vobjDream, form, metadata_to_use )
fitmm = eBayes(fitmm)

# reorganize:
res.df1 <- as.data.frame(topTable( fitmm, number=Inf,coef = "celltypecyc" ))
# res.df1 <- as.data.frame(topTable( fitmm, number=Inf,coef = "age" ))
res.df1$names <- rownames(res.df1)
res.df1$disease_status <- disease_status

age_suffix = ifelse(use_pcw,".age","")
dir.create(paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/",subset_column,"/DE"),showWarnings = FALSE)
f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/",subset_column,"/DE/",disease_status,'_',"cyc_vs_hsc",".de",suffix,".txt")
fwrite(res.df1,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
# dir.create(paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/test/"),showWarnings = FALSE)
# f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/test/cyc_vs_hsc_newpipe.rds")
# saveRDS(list(df.aggre,metadata_to_use),f.out)



