library(data.table)
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

#########33

sampletype <- "Femur"

dir <- '/oak/stanford/groups/smontgom/amarder'
# dir <- '~/Documents/Research'

# datadir="DE_pb_cell_type_groups"; lfc.col="log2FoldChange"; p.col="padj"
# datadir="DE_pb_leiden_names"; lfc.col="log2FoldChange"; p.col="padj"
datadir="DE_pb_leiden_names"; lfc.col="logFC"; p.col="adj.P.Val"

# datadir="DE_cell_type_groups"; lfc.col="logfoldchanges"; p.col="pvals_adj"
# datadir="DE_leiden_names"; lfc.col="logfoldchanges"; p.col="pvals_adj"

pathtodir <- paste0(dir,"/t21-proj/out/full/",datadir)
flist <- list.files(pathtodir)
flist <- flist[grep(sampletype,flist)]
cell_type_groups <- substring(flist,nchar(sampletype)+2,nchar(flist)-4)
cell_type_groups <- cell_type_groups[!grepl("sample",cell_type_groups)]
# cell_type_groups <- c("Erythroid","B cells","NK_T cells","Myeloid","Megakaryocytes","HSC_Progenitors","Stroma","Mast cells")

for (i in 1:length(cell_type_groups)) {
  cell_type <- cell_type_groups[i]
  # f<- paste0(dir,"/t21-proj/out/full/",datadir,"/",sampletype,".",cell_type,".txt")
  f<- paste0(dir,"/t21-proj/out/full/",datadir,"/",sampletype,".",cell_type,".sample.txt")
  df <- fread(f,data.table = F,stringsAsFactors = F)
  x <- df[,c("names",p.col)]
  colnames(x)[2] <- cell_type
  if (i==1) {
    res.df <- x
  } else {
    res.df <- merge(res.df,x,by='names',all = TRUE)
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
  # f<- paste0(dir,"/t21-proj/out/full/",datadir,"/",sampletype,".",cell_type,".txt")
  f<- paste0(dir,"/t21-proj/out/full/",datadir,"/",sampletype,".",cell_type,".sample.txt")
  df <- fread(f,data.table = F,stringsAsFactors = F)
  x <- df[,c("names",lfc.col)]
  colnames(x)[2] <- cell_type
  if (i==1) {
    res.df <- x
  } else {
    res.df <- merge(res.df,x,by='names',all = TRUE)
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
aggregate(apply(df1.lfc[,cell_type_groups],2,rank)/nrow(df1.lfc),by=list(df1.lfc$chr21),median,na.rm=T)
aggregate(df1.lfc[,cell_type_groups],by=list(df1.lfc$chr21),median,na.rm=T)


#######################################################

aggregate(apply(df1.lfc[,cell_type_groups],2,rank)/nrow(df1.lfc),by=list(df1.lfc$chr21),median,na.rm=T)
aggregate(df1.lfc[,cell_type_groups],by=list(df1.lfc$chr21),median,na.rm=T)

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
  labs(x="Cell type",y="LFC (Percentile) (T21 vs Healthy)",fill='') 

########################################################################################
########################################################################################
########################################################################################

head(df1.lfc)
df1.lfc.sort <- df1.lfc[order(df1.lfc$chromosome_name,df1.lfc$start_position),]
df1.lfc.sort$rk <- 1:nrow(df1.lfc.sort)
g <- ggplot(df1.lfc.sort,aes(x=rk,y=Megakaryocytes)) + 
  geom_line() +
  geom_smooth(method='loess',span=0.001) + geom_abline(slope=0,col='red') 

for (chrNum in 1:22) {
  idx <- which(df1.lfc.sort$chromosome_name==chrNum)
  g <- g + 
    geom_vline(xintercept = idx[1],col='red') + 
    geom_vline(xintercept = idx[length(idx)],col='red')
}
g

##############################################

df.obs <- list()
for (chrNum in seq(1:22)) {
  print(chrNum)
  df1.lfc.sub <- subset(df1.lfc,chromosome_name==chrNum)
  df1.p.sub <- subset(df1.p,chromosome_name==chrNum)
  df1.lfc.sort <- df1.lfc.sub[order(df1.lfc.sub$chromosome_name,df1.lfc.sub$start_position),]
  df1.p.sort <- df1.p.sub[order(df1.lfc.sub$chromosome_name,df1.lfc.sub$start_position),]
  # df1.lfc.sort$Megakaryocytes <- rep(c(-1,1),length.out=nrow(df1.lfc.sub))
  flip_ct.lst <- list()
  for (j in 1:length(cell_type_groups)) {
    cell_type = cell_type_groups[j]
    flip_ct = 0; prev <- NULL
    for (i in 1:nrow(df1.lfc.sort)) {
      curr <- df1.lfc.sort[i,cell_type]
      curr.p <- df1.p.sort[i,cell_type]
      if (is.na(curr.p)) {next}
      if (i==1 | is.null(prev)) {
        prev <- curr
      } else {
        if (sign(curr)!=sign(prev)) {
          flip_ct = flip_ct + 1
        }
      }
    }
    flip_ct.lst[[cell_type]] <- flip_ct
  }
  df.obs[[chrNum]] <- as.data.frame(flip_ct.lst)
  df.obs[[chrNum]]$chr <- chrNum
}
df.obs.all <- as.data.frame(do.call(rbind,df.obs))
rownames(df.obs.all) <- NULL

df.perm.chr <- list()
for (chrNum in 1:22) {
  print(paste0("chr: ",chrNum))
  df1.lfc.sub <- subset(df1.lfc,chromosome_name==chrNum)#[1:20,]
  df1.p.sub <- subset(df1.p,chromosome_name==chrNum)#[1:20,]
  df.perm <- list()
  for (k in 1:100) {
    print(paste0("iter: ",k))
    idx <- sample(1:nrow(df1.lfc.sub),size = nrow(df1.lfc.sub),replace = F)
    df1.lfc.sort <- df1.lfc.sub[idx,]
    df1.p.sort <- df1.p.sub[idx,]
    flip_ct.lst <- list()
    # for (j in 8:8) {
    for (j in 1:length(cell_type_groups)) {
      cell_type = cell_type_groups[j]
      flip_ct = 0; prev <- NULL
      for (i in 1:nrow(df1.lfc.sort)) {
        # for (i in 1:10) {
        curr <- df1.lfc.sort[i,cell_type]
        curr.p <- df1.p.sort[i,cell_type]
        if (is.null(curr.p) | is.na(curr.p)) {next}
        if (i==1 | is.null(prev)) {
          prev <- curr
        } else {
          if (sign(curr)!=sign(prev)) {
            flip_ct = flip_ct + 1
            prev <- curr
          }
        }
      }
      flip_ct.lst[[cell_type]] <- flip_ct
    }
    df.perm[[k]] <- as.data.frame(flip_ct.lst)
  }
  df.perm.chr[[chrNum]] <- do.call(rbind,df.perm)
}


p.perm.df <- list()
for (chrNum in 1:22) {
  m = ncol(df.perm.chr[[chrNum]])
  p.perm <- rep(NA,m)
  for (j in 1:m) {
    cell_type = colnames(df.perm.chr[[chrNum]])[j]
    p.perm[j] <- mean(df.obs.all[df.obs.all$chr==chrNum,cell_type] > df.perm.chr[[chrNum]][,cell_type])
  }
  p.perm.df[[chrNum]] <- p.perm
}
p.perm.df.all <- as.data.frame(do.call(rbind,p.perm.df))
colnames(p.perm.df.all) <- colnames(df.perm.chr[[chrNum]])
p.perm.df.all$chr <- 1:nrow(p.perm.df.all)
subset(p.perm.df.all,chr==21)

p.perm.df.all.melt <- melt(p.perm.df.all,id.vars="chr")
library(ggplot2)
ggplot(p.perm.df.all.melt,aes(x=chr,y=variable,fill=value<0.05)) + 
  geom_tile(color = "white")+
  scale_fill_manual(values=c("white","red3")) +
  scale_x_continuous(breaks=seq(1,22)) +
  # scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  theme_minimal()+
  theme(
    plot.title=element_text(hjust=0.5),panel.grid = element_blank()) + 
  # guides(fill=F) +
  coord_fixed() +
  labs(title='',x='Chromosome',y='Cell type',fill="% sig")



#########

g <- ggplot(df1.lfc.sort,aes(x=rk,y=Megakaryocytes)) + 
  # geom_line() +
  geom_abline(slope=0,col='red')
for (chrNum in 1:22) {
  idx <- which(df1.lfc.sort$chromosome_name==chrNum)
  g <- g + 
    geom_vline(xintercept = idx[1],col='red') + 
    geom_vline(xintercept = idx[length(idx)],col='red')
  g <- g + 
    geom_smooth(data=subset(df1.lfc.sort,chromosome_name==chrNum),method='loess',span=0.03) 
}
g + theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5),
        axis.text.x = element_blank()) + 
  labs(x="Genome Position",y="Log Fold Change",title="Megakaryocytes")
g

# g <- ggplot(subset(df1.lfc.sort,chromosome_name==21),aes(x=rk,y=Erythroid)) + 
g <- ggplot(subset(df1.lfc.sort,chromosome_name==21),aes(x=rk,y=`Mast cells`)) + 
  geom_line() +
  geom_abline(slope=0,col='red') +
  geom_smooth(data=subset(df1.lfc.sort,chromosome_name==21),method='loess',span=0.1) 
g + theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5),
        axis.text.x = element_blank()) + 
  labs(x="Genome Position",y="Log Fold Change",title="Megakaryocytes")

########

cell_type <- "Megakaryocytes"
f<- paste0(dir,"/t21-proj/out/full/",datadir,"/",sampletype,".",cell_type,".txt")
df <- fread(f,data.table = F,stringsAsFactors = F)
df2 <- merge(df,annot[,c("hgnc_symbol","chromosome_name","start_position")],by.x="names",by.y="hgnc_symbol")
df2 <- df2[df2$chromosome_name %in% seq(1,22),]
df2$chromosome_name <- factor(df2$chromosome_name,levels=seq(1,22))
df2$chr21 <- factor(ifelse(df2$chromosome_name==21,'Chr 21','Not Chr 21'),levels=c('Not Chr 21','Chr 21'))
df2 <- df2[!is.na(df2$padj),]
# aggregate(df2$log2FoldChange,by=list(df2$chr21),mean)
# aggregate(df2$padj < 0.01,by=list(df2$chr21),mean)

df2 <- df2[order(df2$chromosome_name,df2$start_position),]
df2$rk <- 1:nrow(df2)
# df2$log2FoldChange[(df2$padj > 0.2)] <- NA
g <- ggplot(df2,aes(x=rk,y=log2FoldChange)) +
  # geom_line() +
  geom_abline(slope=0,col='red')
for (chrNum in 1:22) {
  idx <- which(df2$chromosome_name==chrNum)
  g <- g +
    geom_vline(xintercept = idx[1],col='red') +
    geom_vline(xintercept = idx[length(idx)],col='red')
  g <- g +
    # geom_smooth(data=subset(df2,chromosome_name==chrNum),method='loess',span=0.03)
    geom_smooth(data=subset(df2,chromosome_name==chrNum),method='loess',span=0.25)
  
}
g + theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5),
        axis.text.x = element_blank()) +
  labs(x="Genome Position",y="Log Fold Change",title="Megakaryocytes")









df2<-subset(df2,chr21=="Chr 21")

df2[order(df2$log2FoldChange)[1:25],]


