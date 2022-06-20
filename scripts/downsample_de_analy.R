library(data.table)

subset_column="sample"
sampletype="Liver"
cell_type="HSCs/MPPs"
cell_type_filename = gsub("/","_",cell_type)

df.aggre <- fread("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data_pb_leiden/Liver.pb.HSCs_MPPs.sample.txt",data.table = F,stringsAsFactors = F) 
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
annot <- getBM(attributes = c('hgnc_symbol', 'chromosome_name',
                              'start_position', 'end_position'),
               filters = 'hgnc_symbol',
               values = df.aggre$V1,
               mart = ensembl)


permNum=1
aggre.res <- list()
for (permNum in 1:50) {
  f <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names_downsamp/",sampletype,".",cell_type_filename,".",subset_column,".perm",permNum,".txt")
  res.df <- fread(f,data.table=F,stringsAsFactors=F)
  
  df1 <- merge(res.df,annot[,c("hgnc_symbol","chromosome_name","start_position")],by.x="names",by.y="hgnc_symbol")
  df1 <- df1[df1$chromosome_name %in% seq(1,22),]
  df1$chromosome_name <- factor(df1$chromosome_name,levels=seq(1,22))
  df1$chr21 <- factor(ifelse(df1$chromosome_name==21,'Chr 21','Not Chr 21'),levels=c('Not Chr 21','Chr 21'))
  # aggre.res[[permNum]] <- aggregate(df1[,"logFC"],by=list(df1$chr21),mean,na.rm=T)
  aggre.res[[permNum]] <- aggregate(df1[,"adj.P.Val"]<0.05,by=list(df1$chr21),mean,na.rm=T)
  aggre.res[[permNum]]$permNum <- permNum
  
}

aggre.all <- do.call(rbind,aggre.res)
# subset(aggre.all,Group.1=="Chr 21")
# subset(aggre.all,Group.1!="Chr 21")

mean(subset(aggre.all,Group.1=="Not Chr 21")$x<=0.000635055)
mean(subset(aggre.all,Group.1=="Not Chr 21")$x)
mean(subset(aggre.all,Group.1=="Not Chr 21")$x)/0.000635055
mean(subset(aggre.all,Group.1=="Chr 21")$x<=0)
mean(subset(aggre.all,Group.1=="Chr 21")$x)




