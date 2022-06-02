# module load R/4.1.2
# module load R/4.0.0

library(Seurat)
library(data.table)
library(Matrix.utils)
library( "DESeq2" )

sampletype="Femur"
subset_column="sample"

print(sampletype)

# 
print("Reading metadata1...")
disease_status="Healthy"
f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cellComp/10X_",disease_status,"_",sampletype,".cellComp.csv")
meta1.full<-fread(f,data.table = F,stringsAsFactors = F)
meta1.full$environment <- disease_status
cells1 <- unique(meta1.full[,6])
meta1 <- unique(meta1.full[,c("patient","sample","sorting","environment")])

# 
print("Reading metadata2...")
disease_status="DownSyndrome"
f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cellComp/10X_",disease_status,"_",sampletype,".cellComp.csv")
meta2.full<-fread(f,data.table = F,stringsAsFactors = F)
meta2.full$environment <- disease_status
cells2 <- unique(meta2.full[,6])
meta2 <- unique(meta2.full[,c("patient","sample","sorting","environment")])

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

# mean(df$nCounts_RNA[df$patient_sample=="15633 F15633C"])
# aggregate(df$nCounts_RNA,by=list(df$patient_sample),mean)
# 
# subset(df@meta.data,sample=="F15633S" & leiden_v9==cell_type)
# subset(df@meta.data,sample=="F15657F" & leiden_v9==cell_type)

tmp <- subset(meta.all.full,leiden_names==cell_type & environment=="Healthy"); 
healthy_summary <- aggregate(tmp$patient,by=list(patient=tmp$patient,sample=tmp$sample),length)
tmp <- subset(meta.all.full,leiden_names==cell_type & environment=="DownSyndrome"); 
t21_summary <- aggregate(tmp$patient,by=list(patient=tmp$patient,sample=tmp$sample),length)

sampletype="Liver"
subset_column="sample"
print(sampletype)
print("Reading metadata1...")
disease_status="Healthy"
f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cellComp/10X_",disease_status,"_",sampletype,".cellComp.csv")
meta1.full<-fread(f,data.table = F,stringsAsFactors = F)
meta1.full$environment <- disease_status
cells1 <- unique(meta1.full[,6])
meta1 <- unique(meta1.full[,c("patient","sample","sorting")])
# 
print("Reading metadata2...")
disease_status="DownSyndrome"
f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cellComp/10X_",disease_status,"_",sampletype,".cellComp.csv")
meta2.full<-fread(f,data.table = F,stringsAsFactors = F)
meta2.full$environment <- disease_status
cells2 <- unique(meta2.full[,6])
meta2 <- unique(meta2.full[,c("patient","sample","sorting")])

colnames(meta1.full)[6] <- "leiden_names"
colnames(meta2.full)[6] <- "leiden_names"
meta.all.full <- rbind(meta1.full,meta2.full)

########################################################################################################

tmp <- subset(meta.all.full,leiden_names==cell_type & environment=="DownSyndrome"); 
new_t21_summary <- aggregate(tmp$patient,by=list(patient=tmp$patient,sample=tmp$sample),length)

# if patient in both femur and liver, and cell count liver > femur, then use X samples from liver
start=TRUE
patient_in_both <- unique(t21_summary$patient[t21_summary$patient %in% new_t21_summary$patient])
patient_of_interest=patient_in_both[1]
for (patient_of_interest in patient_in_both) {
  t21_summary.sub <- subset(t21_summary,patient==patient_of_interest)
  new_t21_summary.sub <- subset(new_t21_summary,patient==patient_of_interest)
  tmp <- new_t21_summary.sub[1:nrow(t21_summary.sub),]
  if (sum(t21_summary.sub$x > tmp$x) > 0) {print("YIkes.");break}
  
  tmp$x <- t21_summary.sub$x
  if (start) {
    t21_overview <- tmp
    start=FALSE
  } else {
    t21_overview <- rbind(t21_overview,tmp)
  }
}

#

tmp <- subset(meta.all.full,leiden_names==cell_type & environment=="Healthy"); 
new_healthy_summary <- aggregate(tmp$patient,by=list(patient=tmp$patient,sample=tmp$sample),length)

# if patient in both femur and liver, and cell count liver > femur, then use X samples from liver
start=TRUE
patient_in_both <- unique(healthy_summary$patient[healthy_summary$patient %in% new_healthy_summary$patient])
patient_of_interest=patient_in_both[1]
for (patient_of_interest in patient_in_both) {
  healthy_summary.sub <- subset(healthy_summary,patient==patient_of_interest)
  new_healthy_summary.sub <- subset(new_healthy_summary,patient==patient_of_interest)
  tmp <- new_healthy_summary.sub[1:nrow(healthy_summary.sub),]
  # if (sum(healthy_summary.sub$x > tmp$x) > 0) {print("YIkes.");break}
  # healthy_summary.sub$x < tmp$x
  tmp$x[healthy_summary.sub$x < tmp$x] <- healthy_summary.sub$x[healthy_summary.sub$x < tmp$x]
  if (start) {
    healthy_overview <- tmp
    start=FALSE
  } else {
    healthy_overview <- rbind(healthy_overview,tmp)
  }
}
# nrow(healthy_overview); nrow(healthy_summary)

# if (patient_in_both)
# subset(new_t21_summary,patient%in% t21_summary$patient)

# if patient not in femur and liver or cell count liver < femur:
# if patient has 2 samples in femur, then need to find a liver patient with 2 samples as well

# if 1 sample, then just need to select liver patient that has not been selected yet

# create running list of available patients

###########

# cell_type_filename = gsub("/","_",cell_type)

sampletype="Liver"

disease_status="Healthy"
f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/10X_",disease_status,"_",sampletype,".umap2d.cells_removed")
fileName=paste0(f,".rds")
df <- readRDS(file = fileName)
# df <- NormalizeData(df, normalization.method = "LogNormalize", scale.factor = 10000)
colName1 <- paste0("leiden_v",max(as.numeric(substring(colnames(df@meta.data)[grep("leiden_v",colnames(df@meta.data))],nchar("leiden_v")+1)),na.rm=T))
df@meta.data["leiden_names"] <- df@meta.data[colName1]
df@assays$RNA@key <- "RNA_"

disease_status="DownSyndrome"
f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/10X_",disease_status,"_",sampletype,".umap2d.cells_removed")
fileName=paste0(f,".rds")
df2 <- readRDS(file = fileName)
# df2 <- NormalizeData(df2, normalization.method = "LogNormalize", scale.factor = 10000)
colName2 <- paste0("leiden_v",max(as.numeric(substring(colnames(df2@meta.data)[grep("leiden_v",colnames(df2@meta.data))],nchar("leiden_v")+1)),na.rm=T))
df2@meta.data["leiden_names"] <- df2@meta.data[colName2]

permNum=1
print(permNum)

ind.lst <- c()
for (i in 1:nrow(healthy_overview)) {
  ind <- unname(which(df$sample==healthy_overview$sample[i] & df$patient==healthy_overview$patient[i] & df$leiden_names==cell_type))
  ind.lst <- c(ind.lst,sample(ind,healthy_overview$x[i]))
}
ind2.lst <- c()
for (i in 1:nrow(t21_overview)) {
  ind <- unname(which(df2$sample==t21_overview$sample[i] & df2$patient==t21_overview$patient[i] & df2$leiden_names==cell_type))
  ind2.lst <- c(ind2.lst,sample(ind,t21_overview$x[i]))
}

dftmp <- df[,ind.lst]
dfcombined <- merge(df[,ind.lst],
                    y=df[,ind.lst])

column_to_use="leiden_names"
i1=which(as.character(df@meta.data[column_to_use][,1])==cell_type)
i2=which(as.character(df2@meta.data[column_to_use][,1])==cell_type)
dfcombined <- merge(df[,i1],
                    y=df2[,i2])

#

for (permNum in 1:10) {
  print(permNum)
  
  ind.lst <- c()
  for (i in 1:nrow(healthy_overview)) {
    ind <- unname(which(df$sample==healthy_overview$sample[i] & df$patient==healthy_overview$patient[i] & df$leiden_names==cell_type))
    ind.lst <- c(ind.lst,sample(ind,healthy_overview$x[i]))
  }
  ind2.lst <- c()
  for (i in 1:nrow(t21_overview)) {
    ind <- unname(which(df2$sample==t21_overview$sample[i] & df2$patient==t21_overview$patient[i] & df2$leiden_names==cell_type))
    ind2.lst <- c(ind2.lst,sample(ind,t21_overview$x[i]))
  }
  
  column_to_use="leiden_names"
  i1=which(as.character(df@meta.data[column_to_use][,1])==cell_type)
  i2=which(as.character(df2@meta.data[column_to_use][,1])==cell_type)
  dfcombined <- merge(df[,i1],
                      y=df2[,i2])
  
  dfcombined <- merge(df[,ind.lst],
                      y=df2[,ind2.lst])
  dfcombined <- merge(df[,ind.lst],
                      y=df[,ind.lst])
  
  dftmp <- df[,ind.lst]
  dfcombined <- merge(dftmp,
                      y=dftmp)
  
  dfcombined <- merge(df2[,ind2.lst],
                      y=df2[,ind2.lst])
  
  
  dfcombined <- df[,ind.lst]
  df.aggre <- aggregate.Matrix(
    t(
      GetAssayData(object = dfcombined, slot = "counts", assay="RNA")
    ),
    groupings=dfcombined$sample,fun="sum")
  
  df.aggre <- t(df.aggre)
  df.aggre <- as.data.frame(df.aggre)
  
  metadata_to_use <- x[rownames(x) %in% colnames(df.aggre),]
  df.aggre <- as.matrix(df.aggre[,match(rownames(metadata_to_use),colnames(df.aggre))])
  tab <- table(meta.all.full[meta.all.full[,"leiden_names"]==cell_type,'sample'])
  samples_to_keep <- healthy_overview$sample[healthy_overview$x >= 10]
  samples_to_keep <- samples_to_keep[samples_to_keep %in% rownames(metadata_to_use)]
  df.aggre <- df.aggre[,samples_to_keep,drop=FALSE]
  metadata_to_use <- metadata_to_use[samples_to_keep,]
  
  # Standard usage of limma/voom
  geneExpr = DGEList( df.aggre )
  keep <- filterByExpr(geneExpr, group=metadata_to_use$environment)
  geneExpr <- geneExpr[keep,]
  
  geneExpr = calcNormFactors( geneExpr )
  
  # Specify parallel processing parameters
  # this is used implicitly by dream() to run in parallel
  param = SnowParam(8, "SOCK", progressbar=TRUE)
  
  # The variable to be tested must be a fixed effect
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
  
  # save:
  system("mkdir -p /oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names_downsamp")
  f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names_downsamp/",sampletype,".",cell_type_filename,".",subset_column,".perm",permNum,".txt")
  fwrite(res.df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
  
}



