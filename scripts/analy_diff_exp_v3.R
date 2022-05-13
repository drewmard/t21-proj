library(data.table)
library(biomaRt)
# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

#########33

dir <- '/oak/stanford/groups/smontgom/amarder'
dir <- '~/Documents/Research'

# datadir="DE_pb_cell_type_groups"; lfc.col="log2FoldChange"; p.col="padj"
# datadir="DE_pb_leiden_names"; lfc.col="log2FoldChange"; p.col="padj"
# datadir="DE_cell_type_groups"; lfc.col="logfoldchanges"; p.col="pvals_adj"
datadir="DE_leiden_names"; lfc.col="logfoldchanges"; p.col="pvals_adj"

pathtodir <- paste0(dir,"/t21-proj/out/full/",datadir)
flist <- list.files(pathtodir)
flist <- flist[grep("Liver",flist)]
cell_type_groups <- substring(flist,nchar("Liver")+2,nchar(flist)-4)
# cell_type_groups <- c("Erythroid","B cells","NK_T cells","Myeloid","Megakaryocytes","HSC_Progenitors","Stroma","Mast cells")

for (i in 1:length(cell_type_groups)) {
  cell_type <- cell_type_groups[i]
  f<- paste0(dir,"/t21-proj/out/full/",datadir,"/Liver.",cell_type,".txt")
  df <- fread(f,data.table = F,stringsAsFactors = F)
  x <- df[,c("names",p.col)]
  colnames(x)[2] <- cell_type
  if (i==1) {
    res.df <- x
  } else {
    res.df <- merge(res.df,x,by='names')
  }
}

annot <- getBM(attributes = c('hgnc_symbol', 'chromosome_name',
                              'start_position', 'end_position'),
               filters = 'hgnc_symbol', 
               values = res.df$names, 
               mart = ensembl)
df1 <- merge(res.df,annot[,c("hgnc_symbol","chromosome_name","start_position")],by.x="names",by.y="hgnc_symbol")
df1 <- df1[df1$chromosome_name %in% seq(1,22),]
df1$chromosome_name <- factor(df1$chromosome_name,levels=seq(1,22))
df1$chr21 <- factor(ifelse(df1$chromosome_name==21,'Chr 21','Not Chr 21'),levels=c('Not Chr 21','Chr 21'))
df1.p <- df1
aggregate(df1[,cell_type_groups]<0.05,by=list(df1$chr21),mean,na.rm=T)

########3

for (i in 1:length(cell_type_groups)) {
  cell_type <- cell_type_groups[i]
  f<- paste0(dir,"/t21-proj/out/full/",datadir,"/Liver.",cell_type,".txt")
  df <- fread(f,data.table = F,stringsAsFactors = F)
  x <- df[,c("names",lfc.col)]
  colnames(x)[2] <- cell_type
  if (i==1) {
    res.df <- x
  } else {
    res.df <- merge(res.df,x,by='names')
  }
}
annot <- getBM(attributes = c('hgnc_symbol', 'chromosome_name',
                              'start_position', 'end_position'),
               filters = 'hgnc_symbol', 
               values = res.df$names, 
               mart = ensembl)
df1 <- merge(res.df,annot[,c("hgnc_symbol","chromosome_name","start_position")],by.x="names",by.y="hgnc_symbol")
df1 <- df1[df1$chromosome_name %in% seq(1,22),]
df1$chromosome_name <- factor(df1$chromosome_name,levels=seq(1,22))
df1$chr21 <- factor(ifelse(df1$chromosome_name==21,'Chr 21','Not Chr 21'),levels=c('Not Chr 21','Chr 21'))
df1.lfc <- df1

aggregate(apply(df1[,cell_type_groups],2,rank)/nrow(df1),by=list(df1$chr21),median,na.rm=T)
aggregate(df1[,cell_type_groups],by=list(df1$chr21),median,na.rm=T)

df1.p[df1.p$names=="U2AF1",]
df1.lfc[df1.lfc$names=="U2AF1",]


aggregate(df1.p[,cell_type_groups]<0.01,by=list(df1.p$chr21),mean,na.rm=T)

df1.p.melt <- melt(aggregate(df1.p[,cell_type_groups]<0.01,by=list(df1.p$chromosome_name),mean,na.rm=T))
library(ggplot2)
ggplot(df1.p.melt,aes(x=Group.1,y=variable,fill=value)) + 
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  theme_minimal()+
  theme(
        plot.title=element_text(hjust=0.5)) + 
  # guides(fill=F) +
  coord_fixed() +
  labs(title='% DEG',x='Chromosome',y='Cell type',fill="% sig")

df1.lfc.perc <- df1.lfc
df1.lfc.perc[,cell_type_groups] <- apply(df1.lfc.perc[,cell_type_groups],2,rank)/nrow(df1.lfc.perc)
df1.lfc.perc.melt <- melt(df1.lfc.perc[,c(cell_type_groups,"chr21")],id.vars=c("chr21"))
ggplot(df1.lfc.perc.melt,aes(x=variable,fill=chr21,y=value)) + 
  geom_boxplot() + 
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=55,hjust = 1)
  ) +
  labs(x="Cell type",y="LFC (Percentile) (T21 vs Healthy)")




