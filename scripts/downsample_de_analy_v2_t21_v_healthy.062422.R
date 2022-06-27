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

f.out=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/Femur.HSCs_MPPs.sample.txt")
df.tmp <- fread(f.out,data.table=F,stringsAsFactors = F)
df1.full.tmp <- df.tmp[!duplicated(df.tmp$names),]
df1.full.tmp <- merge(df1.full.tmp,geneAnnot,by.x="names",by.y="hgnc_symbol")
df1.full.tmp <- subset(df1.full.tmp,adj.P.Val < 0.05 & !(names %in% df1.full$names))
df1.full <- as.data.frame(rbind(df1.full,df1.full.tmp))

aggregate(df1.full[,"adj.P.Val"]<0.1,by=list(df1.full$chromosome_name==21),mean,na.rm=T)
df1.full.lfc <- df1.full[,c("names","chromosome_name","logFC")]
df1.full.p <- df1.full[,c("names","chromosome_name","P.Value")]
df1.full.fdr <- df1.full[,c("names","chromosome_name","adj.P.Val")]
colnames(df1.full.lfc)[3] <- paste0(colnames(df1.full.lfc)[3],".fem")
colnames(df1.full.p)[3] <- paste0(colnames(df1.full.p)[3],".fem")
colnames(df1.full.fdr)[3] <- paste0(colnames(df1.full.fdr)[3],".fem")

f.out="/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/Liver.HSCs_MPPs.sample.txt"
df <- fread(f.out,data.table=F,stringsAsFactors = F)
df1.full.tmp <- df[!duplicated(df$names),]
df1.full.lfc.tmp <- df1.full.tmp[,c("names","logFC")]
df1.full.p.tmp <- df1.full.tmp[,c("names","P.Value")]
df1.full.fdr.tmp <- df1.full.tmp[,c("names","adj.P.Val")]
colnames(df1.full.lfc.tmp)[2] <- paste0(colnames(df1.full.lfc.tmp)[2],".liv")
colnames(df1.full.p.tmp)[2] <- paste0(colnames(df1.full.p.tmp)[2],".liv")
colnames(df1.full.fdr.tmp)[2] <- paste0(colnames(df1.full.fdr.tmp)[2],".liv")
df1.full.lfc = merge(df1.full.lfc,df1.full.lfc.tmp,by='names',all.x=TRUE)
df1.full.p = merge(df1.full.p,df1.full.p.tmp,by='names',all.x=TRUE)
df1.full.fdr = merge(df1.full.fdr,df1.full.fdr.tmp,by='names',all.x=TRUE)

# subset(df1.full.fdr,names=="EPX")
# df1.full.tmp <- merge(df1.full.tmp,geneAnnot,by.x="names",by.y="hgnc_symbol")
# aggregate(df1.full.tmp[,"adj.P.Val"]<0.1,by=list(df1.full.tmp$chromosome_name==21),mean,na.rm=T)


#############

permNum=1
aggre.res <- list()
for (permNum in 1:100) {
  print(permNum)
  
  f <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names_downsamp/","Liver",".",cell_type_filename,".",subset_column,".perm",permNum,".txt")
  res.df <- fread(f,data.table=F,stringsAsFactors=F)
  res.df <- res.df[!duplicated(res.df$names),]
  
  df1 <- res.df[!duplicated(res.df$names),]
  
  tmp <- df1[,c("names","logFC")]
  colnames(tmp)[2] <- paste0("sim",permNum)
  df1.full.lfc <- merge(df1.full.lfc,tmp,by="names",all.x=TRUE)
  tmp <- df1[,c("names","P.Value")]
  colnames(tmp)[2] <- paste0("sim",permNum)
  df1.full.p <- merge(df1.full.p,tmp,by="names",all.x=TRUE)
  tmp <- df1[,c("names","adj.P.Val")]
  colnames(tmp)[2] <- paste0("sim",permNum)
  df1.full.fdr <- merge(df1.full.fdr,tmp,by="names",all.x=TRUE)
  
  # aggre.res[[permNum]] <- aggregate(df1[,"logFC"],by=list(df1$chr21),mean,na.rm=T)
  # aggre.res[[permNum]] <- aggregate(df1[,"adj.P.Val"]<0.05,by=list(df1$chr21),mean,na.rm=T)
  # aggre.res[[permNum]]$permNum <- permNum
  
}

df1.full.lfc[is.na(df1.full.lfc)] <- 0
df1.full.p[is.na(df1.full.p)] <- 1
df1.full.fdr[is.na(df1.full.fdr)] <- 1

# df1.full.fdr[order(df1.full.fdr$adj.P.Val.liv),][1,]
# df1.full.p[order(df1.full.p$P.Value.liv),][1,]
# 
# subset(df1.full.lfc,names=="GATA1")
# subset(df1.full.lfc,names=="APOC1")
# subset(df1.full.p,names=="APOC1")
# 
# subset(df1.full.lfc,names=="MYL4")

permNum="ALL"
f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names_downsamp/","Liver",".",cell_type_filename,".",subset_column,".perm",permNum,".lfc.txt")
fwrite(df1.full.lfc,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names_downsamp/","Liver",".",cell_type_filename,".",subset_column,".perm",permNum,".p.txt")
fwrite(df1.full.p,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names_downsamp/","Liver",".",cell_type_filename,".",subset_column,".perm",permNum,".fdr.txt")
fwrite(df1.full.fdr,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)






