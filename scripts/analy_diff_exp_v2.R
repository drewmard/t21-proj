library(data.table)
library(ggplot2)
library(biomaRt)

# df <- fread("~/Documents/Research/t21-proj/out/full/DE_cell_type_groups/HSC_Progenitors.txt",data.table = F,stringsAsFactors = F)
# df <- fread("~/Documents/Research/t21-proj/out/full/DE_cell_type_groups/Erythroid.txt",data.table = F,stringsAsFactors = F)
# df <- fread("~/Documents/Research/t21-proj/out/full/DE_leiden_names/Early erythroid cells.txt",data.table = F,stringsAsFactors = F)
df <- fread("~/Documents/Research/t21-proj/out/full/DE_cell_type_groups/Liver.Erythroid.txt",data.table = F,stringsAsFactors = F)
head(df)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
annot <- getBM(attributes = c('hgnc_symbol', 'chromosome_name',
                     'start_position', 'end_position'),
      filters = 'hgnc_symbol', 
      values = df$names, 
      mart = ensembl)
df1 <- merge(df,annot,by.x="names",by.y="hgnc_symbol")
df1 <- df1[df1$chromosome_name %in% seq(1,22),]
df1$chromosome_name <- factor(df1$chromosome_name,levels=seq(1,22))
df1$chr21 <- factor(ifelse(df1$chromosome_name==21,'Chr 21','Not Chr 21'),levels=c('Not Chr 21','Chr 21'))

gene_info <- fread("~/Documents/Research/data/grch38/protein_coding_genes.txt",data.table = F,stringsAsFactors = F)
df2 <- merge(df,gene_info[,c("Chromosome/scaffold name","Gene name")],by.x="names",by.y="Gene name")
df2 <- df2[df2$`Chromosome/scaffold name` %in% seq(1,22),]
df2$`Chromosome/scaffold name` <- factor(df2$`Chromosome/scaffold name`,levels=seq(1,22))
df2$chr21 <- factor(ifelse(df2$`Chromosome/scaffold name`==21,'Chr 21','Not Chr 21'),levels=c('Not Chr 21','Chr 21'))

dim(df1)
dim(df2)

subset(df2,!(names %in% df1$names))[1,]
subset(df1,!(names %in% df2$names))[1:23,]

ggplot(df1,aes(x=chr21,y=rank(logfoldchanges)/nrow(df1),fill=chr21)) + 
  geom_boxplot() + 
  theme_bw() +
  labs(x="Chromosome",y='Log Fold Change Percentile (T21 vs Healthy)') +
  scale_fill_brewer(palette = "Set2") + guides(fill=F) + ylim(0,1)

ggplot(df2,aes(x=chr21,y=rank(logfoldchanges)/nrow(df2),fill=chr21)) + 
  geom_boxplot() + 
  theme_bw() +
  labs(x="Chromosome",y='Log Fold Change Percentile (T21 vs Healthy)') +
  scale_fill_brewer(palette = "Set2") + guides(fill=F)





df_pb <- fread("~/Documents/Research/t21-proj/out/full/DE_pb_cell_type_groups/Liver.Erythroid.txt",data.table = F,stringsAsFactors = F)
df.mg <- merge(df,df_pb,by.x="names",by.y="gene")
cor(df.mg$logfoldchanges,df.mg$logFC)
cor(-log10(df.mg$pvals+1e-216),-log10(df.mg$PValue+1e-216),use="complete.obs")

gene_info <- fread("~/Documents/Research/data/grch38/protein_coding_genes.txt",data.table = F,stringsAsFactors = F)

df <- merge(df.mg,gene_info[,c("Chromosome/scaffold name","Gene name")],by.x="names",by.y="Gene name")

df <- df[df$`Chromosome/scaffold name` %in% seq(1,22),]
df$`Chromosome/scaffold name` <- factor(df$`Chromosome/scaffold name`,levels=seq(1,22))
# ggplot(df,aes(x=`Chromosome/scaffold name`,y=logfoldchanges)) + geom_boxplot()

# aggregate(df$logfoldchanges,by=list(df$`Chromosome/scaffold name`),median)

# ggplot(df,aes(x=`Chromosome/scaffold name`==21,y=logfoldchanges)) + geom_boxplot()

source("~/Documents/Research/Useful_scripts/rntransform.R")
aggregate(df$logfoldchanges,by=list(df$`Chromosome/scaffold name`==21),median)

df$chr21 <- factor(ifelse(df$`Chromosome/scaffold name`==21,'Chr 21','Not Chr 21'),levels=c('Not Chr 21','Chr 21'))
ggplot(df,aes(x=chr21,y=rank(logfoldchanges)/nrow(df),fill=chr21)) + 
  geom_boxplot() + 
  theme_bw() +
  labs(x="Chromosome",y='Log Fold Change Percentile (T21 vs Healthy)') +
  scale_fill_brewer(palette = "Set2") + guides(fill=F)
ggplot(df,aes(x=chr21,y=rank(logFC)/nrow(df),fill=chr21)) + 
  geom_boxplot() + 
  theme_bw() +
  labs(x="Chromosome",y='Log Fold Change Percentile (T21 vs Healthy)') +
  scale_fill_brewer(palette = "Set2") + guides(fill=F)

ggplot(df,aes(x=`Chromosome/scaffold name`,y=logfoldchanges)) + geom_boxplot()
ggplot(df,aes(x=`Chromosome/scaffold name`,y=logFC)) + geom_boxplot()

aggregate(df$pvals_adj<0.01,by=list(df$`Chromosome/scaffold name`),mean)

###################

# module load R/4.0
library(data.table)
cell_type_filename="Erythroid"
sampletype="Liver"
f <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data_pb_groups/",sampletype,".pb.",cell_type_filename,".txt")
df.aggre <- fread(f,data.table = F,stringsAsFactors = F,header = T)
rownames(df.aggre) <- df.aggre[,1]; df.aggre <- df.aggre[,-1]
f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/",sampletype,".metadata.txt")
x <- fread(f,data.table = F,stringsAsFactors = F)
x <- data.frame(x,row.names=x$patient)
# x <- x[match(colnames(df.aggre),rownames(x)),]
df.aggre <- as.matrix(df.aggre[,match(rownames(x),colnames(df.aggre))])

# library(seurat)
# cell_type_filename="Erythroid"
# sampletype="Liver"
# f <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/",sampletype,".sc.",cell_type_filename,".rds")
# dfcombined <- readRDS(f)

library( "DESeq2" )
dds <- DESeqDataSetFromMatrix(countData=df.aggre, 
                              colData=x, 
                              design=~environment)
dds$environment <- relevel(dds$environment, ref = "Healthy")
dds <- DESeq(dds)
res <- results(dds)
res.df <- as.data.frame(res)
res.df$names <- rownames(res.df)

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
annot <- getBM(attributes = c('hgnc_symbol', 'chromosome_name',
                              'start_position', 'end_position'),
               filters = 'hgnc_symbol', 
               values = res.df$names, 
               mart = ensembl)
df1 <- merge(res.df,annot,by.x="names",by.y="hgnc_symbol")
df1 <- df1[df1$chromosome_name %in% seq(1,22),]
df1$chromosome_name <- factor(df1$chromosome_name,levels=seq(1,22))
df1$chr21 <- factor(ifelse(df1$chromosome_name==21,'Chr 21','Not Chr 21'),levels=c('Not Chr 21','Chr 21'))
aggregate(df1$log2FoldChange,by=list(df1$chr21),median,na.rm=T)
aggregate(rank(df1$log2FoldChange)/nrow(df1),by=list(df1$chr21),median,na.rm=T)

df <- fread("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_cell_type_groups/Liver.Erythroid.txt",data.table = F,stringsAsFactors = F)
res.df <- merge(res.df,df,by='names')
df1 <- merge(res.df,annot,by.x="names",by.y="hgnc_symbol")
df1 <- df1[df1$chromosome_name %in% seq(1,22),]
df1$chromosome_name <- factor(df1$chromosome_name,levels=seq(1,22))
df1$chr21 <- factor(ifelse(df1$chromosome_name==21,'Chr 21','Not Chr 21'),levels=c('Not Chr 21','Chr 21'))

cor(res.df$logfoldchanges,res.df$log2FoldChange,use='complete.obs')
cor(-log10(res.df$pvals+1e-216),-log10(res.df$pvalue+1e-216),use="complete.obs")

aggregate(df1$padj < 0.1,by=list(df1$chr21),mean,na.rm=T)
aggregate(df1$pvals_adj < 0.1,by=list(df1$chr21),mean,na.rm=T)

aggregate(df1$padj < 0.1 & df1$pvals_adj < 0.1,by=list(df1$chr21),mean,na.rm=T)


aggregate(df1$log2FoldChange,by=list(df1$chr21),median,na.rm=T)
aggregate(rank(df1$log2FoldChange)/nrow(df1),by=list(df1$chr21),median,na.rm=T)
aggregate(df1$logfoldchanges,by=list(df1$chr21),median,na.rm=T)
aggregate(rank(df1$logfoldchanges)/nrow(df1),by=list(df1$chr21),median,na.rm=T)





