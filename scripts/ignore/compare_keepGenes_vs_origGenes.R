library(data.table)

sampletype="Liver"
cell_type_filename="HSCs_MPPs"
subset_column="sample"

fileOut <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/data_small/gene_annot.txt")
geneAnnot <- fread(fileOut,data.table = F,stringsAsFactors = F)
geneAnnot <- geneAnnot[!duplicated(geneAnnot$hgnc_symbol),]

# f.out=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_liver_v_femur/",cell_type_filename,'.',disease_status,'.FE.txt')
f.out=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/Femur.HSCs_MPPs.sample.keep_t21_genes.txt")
df <- fread(f.out,data.table=F,stringsAsFactors = F)
df1.full <- df[!duplicated(df$names),]
df1.full <- merge(df1.full,geneAnnot,by.x="names",by.y="hgnc_symbol")
aggregate(df1.full[,"adj.P.Val"]<0.1,by=list(df1.full$chromosome_name==21),mean,na.rm=T)
df1.full.lfc <- df1.full[,c("names","chromosome_name","logFC")]
df1.full.p <- df1.full[,c("names","chromosome_name","P.Value")]
df1.full.fdr <- df1.full[,c("names","chromosome_name","adj.P.Val")]
colnames(df1.full.lfc)[3] <- paste0(colnames(df1.full.lfc)[3],".keep")
colnames(df1.full.p)[3] <- paste0(colnames(df1.full.p)[3],".keep")
colnames(df1.full.fdr)[3] <- paste0(colnames(df1.full.fdr)[3],".keep")

f.out=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/Femur.HSCs_MPPs.sample.txt")
df <- fread(f.out,data.table=F,stringsAsFactors = F)
df1.full.tmp <- df[!duplicated(df$names),]
df1.full.lfc.tmp <- df1.full.tmp[,c("names","logFC")]
df1.full.p.tmp <- df1.full.tmp[,c("names","P.Value")]
df1.full.fdr.tmp <- df1.full.tmp[,c("names","adj.P.Val")]
colnames(df1.full.lfc.tmp)[2] <- paste0(colnames(df1.full.lfc.tmp)[2],".orig")
colnames(df1.full.p.tmp)[2] <- paste0(colnames(df1.full.p.tmp)[2],".orig")
colnames(df1.full.fdr.tmp)[2] <- paste0(colnames(df1.full.fdr.tmp)[2],".orig")
df1.full.lfc = merge(df1.full.lfc,df1.full.lfc.tmp,by='names',all=T)
df1.full.p = merge(df1.full.p,df1.full.p.tmp,by='names',all=T)
df1.full.fdr = merge(df1.full.fdr,df1.full.fdr.tmp,by='names',all=T)

cor(-log10(df1.full.p$P.Value.orig),-log10(df1.full.p$P.Value.keep),use='na.or.complete')
cor(-log10(df1.full.p$P.Value.orig),-log10(df1.full.p$P.Value.keep),use='na.or.complete')
subset(df1.full.p,!is.na(P.Value.keep))$P.Value.orig

which.min(subset(df1.full.p,is.na(P.Value.orig))$P.Value.keep)
subset(df1.full.p,is.na(P.Value.orig))[1513,]
df1.full.fdr$use_orig <- is.na(df1.full.fdr$adj.P.Val.orig)
df1.full.fdr$use_keep <- is.na(df1.full.fdr$adj.P.Val.keep)
mean(df1.full.fdr$use_orig<0.05)
mean(df1.full.fdr$adj.P.Val.orig<0.05,na.rm=T)
aggregate(df1.full.fdr$adj.P.Val.orig<0.05,by=list(keep=df1.full.fdr$use_keep),mean,na.rm=T)
aggregate(df1.full.fdr$adj.P.Val.keep<0.05,by=list(keep=df1.full.fdr$use_orig),sum,na.rm=T)


df1.full.p$use_orig <- is.na(df1.full.p$P.Value.orig)
df1.full.p$use_keep <- is.na(df1.full.p$P.Value.keep)
aggregate(df1.full.p$P.Value.orig<0.05,by=list(keep=df1.full.p$use_keep),mean,na.rm=T)
aggregate(df1.full.p$P.Value.keep<0.05,by=list(keep=df1.full.p$use_orig),mean,na.rm=T)



mean(df1.full.fdr$adj.P.Val.keep<0.05,na.rm=T)

df1.full.tmp <- merge(df1.full.tmp,geneAnnot,by.x="names",by.y="hgnc_symbol")
