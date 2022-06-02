library(data.table)
library('variancePartition')
library('edgeR')
library('BiocParallel')

res.df.lst <- list()
df.aggre.lst <- list()

sampletype="Liver"
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

# which(metadata_to_use[,3] %in% c("L15633R","L15781A"))
# metadata_to_use <- metadata_to_use[c(4,6,8,27,42),]
# df.aggre <- df.aggre[,c(4,6,8,27,42)]

# Standard usage of limma/voom
geneExpr = DGEList( df.aggre )
keep <- filterByExpr(geneExpr, group=metadata_to_use$environment)
geneExpr <- geneExpr[keep,]

geneExpr = calcNormFactors( geneExpr )

# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param = SnowParam(8, "SOCK", progressbar=TRUE)

form <- ~ environment + patient + sorting
if (length(unique(metadata_to_use$sorting))==1) {
  form <- ~ environment + patient
}
# form <- ~ environment

# A positive FC is increased expression in the DS compared to healthy
metadata_to_use$environment <- factor(metadata_to_use$environment,levels=c("Healthy","DownSyndrome"))

#######################################
# 1

design <- model.matrix(form,data=metadata_to_use)
colnames(design)[2] <- 'environment'
# colnames(design) <- make.names(levels(factor(metadata_to_use$environment)))
rownames(design) <- metadata_to_use$patient

# A positive FC is increased expression in the DS compared to healthy
# contrast <- makeContrasts(paste0(make.names("Down.Syndrome"),"-Healthy"), levels = design)

y <- DGEList(counts=df.aggre,group=metadata_to_use$environment)
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
res <- topTags(lrt, n=Inf)$table
res$gene <- rownames(res)
rownames(res) <- NULL
res.df <- res
res.df$sampletype <- sampletype

#######################################
# 1

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
res.df$sampletype <- sampletype

# save:
res.df.lst[[sampletype]] <- res.df

#########################

# res.df.all <- as.data.frame(do.call(rbind,res.df.lst))

# if using limma:
res.df.all <- merge(res.df.lst[["Liver"]][,c("names","logFC","adj.P.Val")],res.df.lst[["Femur"]][,c("names","logFC","adj.P.Val")],by='names',all=T)
# if using edgeR:
res.df.all <- merge(res.df.lst[["Liver"]][,c("gene","logFC","FDR")],res.df.lst[["Femur"]][,c("gene","logFC","FDR")],by='gene',all=T)
colnames(res.df.all)[1] <- "names"

res.df.all.p <- readRDS("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/res.Femur.p.rds")
res.df.all <- merge(res.df.all,res.df.all.p[[1]][,c("names","chromosome_name")],by='names',all=T)

aggregate(res.df.all[,c(3,5)]<0.1,by=list(res.df.all$chromosome_name==21),mean,na.rm=T)
aggregate(res.df.all[,c(2,4)],by=list(chr21=res.df.all$chromosome_name==21),mean,na.rm=T)
aggregate(res.df.all[,c(2,4)],by=list(chr21=res.df.all$chromosome_name==21),median,na.rm=T)

aggregate(res.df.all[,c(2,4)],by=list(chr21=res.df.all$chromosome_name==21,sig=res.df.all[,3]<0.1),median,na.rm=T)
aggregate(res.df.all[,c(2,4)],by=list(chr21=res.df.all$chromosome_name==21,sig=res.df.all[,3]<0.1),mean,na.rm=T)

res.sub <- subset(res.df.all,chromosome_name==21)
cor(res.sub[,c(2:5)],use="na.or.complete")
cor.test(res.sub[,2],res.sub[,4])
t.test(res.sub[,2],res.sub[,4],paired=T)
res.sub <- subset(res.df.all,chromosome_name!=21)
cor(res.sub[,c(2:5)],use="na.or.complete")
t.test(res.sub[,2],res.sub[,4],paired=T)

res.df.all[order(res.df.all$adj.P.Val.x)[1:3],]

(cpm(df.aggre))["SOD1",]


library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/",sampletype,".",cell_type_filename,".",subset_column,".random.txt")
# f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/",sampletype,".",cell_type_filename,".",subset_column,".txt")
fwrite(res.df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

res.df.all[order(res.df.all)[logFC.x]]

annot <- getBM(attributes = c('hgnc_symbol', 'chromosome_name',
                              'start_position', 'end_position'),
               filters = 'hgnc_symbol', 
               values = res.df.all$names, 
               mart = ensembl)
df1 <- merge(res.df.all,annot[,c("hgnc_symbol","chromosome_name","start_position")],by.x="names",by.y="hgnc_symbol")
df1 <- df1[df1$chromosome_name %in% seq(1,22),]
df1$chromosome_name <- factor(df1$chromosome_name,levels=seq(1,22))
df1$chr21 <- factor(ifelse(df1$chromosome_name==21,'Chr 21','Not Chr 21'),levels=c('Not Chr 21','Chr 21'))
res.df.all <- df1


