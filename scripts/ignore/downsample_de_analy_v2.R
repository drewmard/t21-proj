library(data.table)
sampletype="Liver"
cell_type_filename="HSCs_MPPs"
subset_column="sample"

fileOut <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/data_small/gene_annot.txt")
geneAnnot <- fread(fileOut,data.table = F,stringsAsFactors = F)
geneAnnot <- geneAnnot[!duplicated(geneAnnot$hgnc_symbol),]

f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/",sampletype,".",cell_type_filename,".",subset_column,".txt")
df <- fread(f.out,data.table=F,stringsAsFactors = F)
df1.full <- merge(df,geneAnnot[,c("hgnc_symbol","chromosome_name","start_position")],by.x="names",by.y="hgnc_symbol")
df1.full <- df1.full[!duplicated(df1.full$names),]
aggregate(df1.full[,"adj.P.Val"]<0.1,by=list(df1.full$chromosome_name==21),mean,na.rm=T)

df1.full.lfc <- df1.full[,c("names","chromosome_name","logFC")]
df1.full.p <- df1.full[,c("names","chromosome_name","P.Value")]
df1.full.fdr <- df1.full[,c("names","chromosome_name","adj.P.Val")]

# subset(df1,names=="MRPS5")
# subset(df1.full.p,names=="MRPS5")
# 
df1.full.fdr[duplicated(df1.full.fdr$names),]
res.df[duplicated(res.df$names),]
df1[duplicated(df1$names),]

permNum=1
aggre.res <- list()
for (permNum in 1:50) {
  print(permNum)
  
  f <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names_downsamp/",sampletype,".",cell_type_filename,".",subset_column,".perm",permNum,".txt")
  res.df <- fread(f,data.table=F,stringsAsFactors=F)
  res.df <- res.df[!duplicated(res.df$names),]
  
  df1 <- merge(res.df,geneAnnot[,c("hgnc_symbol","chromosome_name","start_position")],by.x="names",by.y="hgnc_symbol")
  df1 <- df1[!duplicated(df1$names),]
  
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

# df1.full.fdr[order(df1.full.fdr$adj.P.Val),][1,]
# df1.full.p[order(df1.full.p$P.Value),][1,]
# df1.full.lfc[order(df1.full.p$P.Value),][1,]
# 
# cor(df1.full.lfc$logFC[df1.full.p$P.Value < 0.05],df1.full.lfc$sim1[df1.full.p$P.Value < 0.05])
# t.test(abs(df1.full.lfc$logFC[df1.full.p$P.Value < 0.05]),abs(df1.full.lfc$sim1[df1.full.p$P.Value < 0.05]),paired=T)
# t.test(abs(df1.full.lfc$logFC[df1.full.fdr$adj.P.Val < 0.1]),abs(df1.full.lfc$sim1[df1.full.fdr$adj.P.Val < 0.1]),paired=T)

df1.full.fdr$med <- apply(df1.full.fdr[,3:ncol(df1.full.fdr)],1,median)
df1.full.p$med <- apply(df1.full.p[,3:ncol(df1.full.p)],1,median)
df1.full.lfc$med <- apply(df1.full.lfc[,3:ncol(df1.full.lfc)],1,median)

cor(df1.full.lfc$logFC[df1.full.p$P.Value < 0.05],df1.full.lfc$med[df1.full.p$P.Value < 0.05])
t.test(abs(df1.full.lfc$logFC[df1.full.p$P.Value < 0.05]),abs(df1.full.lfc$med[df1.full.p$P.Value < 0.05]),paired=T)
t.test(abs(df1.full.lfc$logFC[df1.full.fdr$adj.P.Val < 0.1]),abs(df1.full.lfc$med[df1.full.fdr$adj.P.Val < 0.1]),paired=T)




aggre.all <- do.call(rbind,aggre.res)
