library(data.table)
library(ggplot2)
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

dir <- '/oak/stanford/groups/smontgom/amarder'
# dir <- '~/Documents/Research'

sampletype <- "Liver"
sampletype <- "Femur"
datadir="DE_pb_leiden_names";

#########33

res.df.all <- list()
for (iter in 1:3) {
  print(iter)
  
  pathtodir <- paste0(dir,"/t21-proj/out/full/",datadir)
  flist <- list.files(pathtodir)
  flist <- flist[grep(sampletype,flist)]
  cell_type_groups <- flist[!grepl("res",flist)]
  cell_type_groups <- substring(cell_type_groups,nchar(sampletype)+2,nchar(cell_type_groups)-4)
  cell_type_groups <- cell_type_groups[!grepl("sample",cell_type_groups)]
  
  # for (i in 1:21) {
  for (i in 1:length(cell_type_groups)) {
    print(i)
    cell_type <- cell_type_groups[i]
    if (iter==1) {
      lfc.col="log2FoldChange"; p.col="padj"
      f<- paste0(dir,"/t21-proj/out/full/",datadir,"/",sampletype,".",cell_type,".txt")
    } else if (iter==2) {
      lfc.col="logFC"; p.col="adj.P.Val"
      f<- paste0(dir,"/t21-proj/out/full/",datadir,"/",sampletype,".",cell_type,".sample.txt")
    } else if (iter==3) {
      lfc.col="logFC"; p.col="adj.P.Val"
      f<- paste0(dir,"/t21-proj/out/full/",datadir,"/",sampletype,".",cell_type,".sample.random.txt")
    }
    df <- fread(f,data.table = F,stringsAsFactors = F)
    x <- df[,c("names",p.col)]
    colnames(x)[2] <- cell_type
    if (i==1) {
      res.df <- x
    } else {
      res.df <- merge(res.df,x,by='names',all = TRUE)
    }
  }
  res.df.all[[iter]] <- res.df
}

#############

for (iter in 1:3) {
  res.df <- res.df.all[[iter]]
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
  
  res.df.all[[iter]] <- df1.p
}
res.df.all.p <- res.df.all
f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".p.rds")
saveRDS(res.df.all.p,file=f.out)


###############

res.df.all <- list()
for (iter in 1:3) {
  print(iter)
  
  for (i in 1:length(cell_type_groups)) {
    print(i)
    cell_type <- cell_type_groups[i]
    if (iter==1) {
      lfc.col="log2FoldChange"; p.col="padj"
      f<- paste0(dir,"/t21-proj/out/full/",datadir,"/",sampletype,".",cell_type,".txt")
    } else if (iter==2) {
      lfc.col="logFC"; p.col="adj.P.Val"
      f<- paste0(dir,"/t21-proj/out/full/",datadir,"/",sampletype,".",cell_type,".sample.txt")
    } else if (iter==3) {
      lfc.col="logFC"; p.col="adj.P.Val"
      f<- paste0(dir,"/t21-proj/out/full/",datadir,"/",sampletype,".",cell_type,".sample.random.txt")
    }
    df <- fread(f,data.table = F,stringsAsFactors = F)
    x <- df[,c("names",lfc.col)]
    colnames(x)[2] <- cell_type
    if (i==1) {
      res.df <- x
    } else {
      res.df <- merge(res.df,x,by='names',all = TRUE)
    }
  }
  res.df.all[[iter]] <- res.df
}

#############

for (iter in 1:3) {
  res.df <- res.df.all[[iter]]
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
  
  res.df.all[[iter]] <- df1.p
}

res.df.all.lfc <- res.df.all

f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".lfc.rds")
saveRDS(res.df.all.lfc,file=f.out)

###############

##############

sampletype="Femur"
res.df.all.p <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".p.rds"))
res.df.all.lfc <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".lfc.rds"))
cell_type_groups <- colnames(res.df.all.p[[1]])[(which(colnames(res.df.all.p[[1]])=="names")+1):(which(colnames(res.df.all.p[[1]])=="chromosome_name")-1)]

tab <- lapply(res.df.all.p,function(x) {aggregate(x[,cell_type_groups]<0.1,by=list(x$chr21),mean,na.rm=T)})
tab <- lapply(tab,function(x) {rownames(x) <- x[,1]; x <- x[,-1]; x <- as.data.frame(t(x))})
tab <- do.call(cbind,tab)

# fetus pb more similar to per sample pb w/ mixed model than per sample pb w/ fixed effect model
cor(tab[,c(2,4,6)])
pairs(tab[,c(2,4,6)])
cor(tab[,c(1,3,5)])
pairs(tab[,c(1,3,5)])

# generally: more chr 21 dysregulation is linked to more non chr 21 dysregulation
cor(tab)
plot(tab[,c(2,1)])
plot(tab[,c(4,3)])
plot(tab[,c(6,5)])


# does the random effect model increase power on non chr 21?
t.test(tab[,5] - tab[,1])
# no - significantly decreased power!

# does the random effect model increase power on chr 21?
t.test(tab[,6] - tab[,2])
# again, no! P < 0.05 and mean < 0

###############

tab <- lapply(res.df.all.lfc,function(x) {aggregate(apply(x[,cell_type_groups],2,rank)/nrow(x),by=list(x$chr21),median,na.rm=T)})
tab <- lapply(tab,function(x) {rownames(x) <- x[,1]; x <- x[,-1]; x <- as.data.frame(t(x))})
tab <- do.call(cbind,tab)
cor(tab[,c(2,4,6)])
pairs(tab[,c(2,4,6)])
cor(tab[,c(1,3,5)])
pairs(tab[,c(1,3,5)])

# generally, more upregulation on chr 21 is linked to more downregulation on other chr
cor(tab)
par(mfrow=c(1,3))
plot(tab[,c(2,1)])
plot(tab[,c(4,3)])
plot(tab[,c(6,5)])
par(mfrow=c(1,1))


tab <- lapply(res.df.all.lfc,function(x) {aggregate(x[,cell_type_groups],by=list(x$chr21),median,na.rm=T)})
tab <- lapply(tab,function(x) {rownames(x) <- x[,1]; x <- x[,-1]; x <- as.data.frame(t(x))})
tab <- do.call(cbind,tab)
cor(tab[,c(2,4,6)])
pairs(tab[,c(2,4,6)])
cor(tab[,c(1,3,5)])
pairs(tab[,c(1,3,5)])

# generally, more upregulation on chr 21 is linked to more upregulation on other chr
# neg corr seen in femur is due to a single outlier
cor(tab)
par(mfrow=c(1,3))
plot(tab[,c(2,1)])
plot(tab[,c(4,3)])
plot(tab[,c(6,5)])
par(mfrow=c(1,1))

iter=1
df.chr21 <- subset(res.df.all.lfc[[iter]],chromosome_name==21)
df.not_chr21 <- subset(res.df.all.lfc[[iter]],chromosome_name!=21)
cor.mat <- cor(df.chr21[,cell_type_groups],use='complete.obs')
cor.mat[upper.tri(cor.mat,diag = TRUE)] <- NA
mean(cor.mat,na.rm=TRUE)
t.test(cor.mat,na.rm=TRUE)

cor.mat <- cor(df.not_chr21[,cell_type_groups],use='complete.obs')
cor.mat[upper.tri(cor.mat,diag = TRUE)] <- NA
mean(cor.mat,na.rm=TRUE)
t.test(cor.mat,na.rm=TRUE)

# cluster cell types together?

##########

sampletype="Femur"
res.df.all.p <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".p.rds"))
res.df.all.lfc <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".lfc.rds"))
cell_type_groups <- colnames(res.df.all.p[[1]])[(which(colnames(res.df.all.p[[1]])=="names")+1):(which(colnames(res.df.all.p[[1]])=="chromosome_name")-1)]

iter=3
df1.p <- res.df.all.p[[iter]]
df1.p.melt <- reshape2::melt(aggregate(df1.p[,cell_type_groups]<0.1,by=list(df1.p$chromosome_name),mean,na.rm=T),id.vars="Group.1")
ggplot(df1.p.melt,aes(x=Group.1,y=variable,fill=value)) + 
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  theme_minimal()+
  theme(
    plot.title=element_text(hjust=0.5)) + 
  # guides(fill=F) +
  coord_fixed() +
  labs(title='% DEG',x='Chromosome',y='Cell type',fill="% sig")

df1.lfc.perc <- res.df.all.lfc[[iter]]
df1.lfc.perc[,cell_type_groups] <- apply(df1.lfc.perc[,cell_type_groups],2,rank)/nrow(df1.lfc.perc)
df1.lfc.perc.melt <- reshape2::melt(df1.lfc.perc[,c(cell_type_groups,"chr21")],id.vars=c("chr21"))
ggplot(df1.lfc.perc.melt,aes(x=variable,fill=chr21,y=value)) + 
  geom_boxplot() + 
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=55,hjust = 1)
  ) +
  labs(x="Cell type",y="LFC (Percentile) (T21 vs Healthy)",fill='') 


##################

cell_type <- "Early erythroid cells" # cell_type_groups[1]
cell_type

# dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"
dir="~/Documents/Research/neuro-variants"

outdir = "~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/gsea"
system(paste0("mkdir -p ",outdir))

print("Performing GSEA...")
library(fgsea)
gene_set_list <- c("h.all","c2.all","c5.go","c8.all")
gene_col_name="gene"
gene_set="h.all"
subset_to_use="not_chr21_up"
iter=1

sampletype="Femur"
for (sampletype in c("Femur","Liver")) {
  print(sampletype)
  res.df.all.p <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".p.rds"))
  res.df.all.lfc <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".lfc.rds"))
  cell_type_groups <- colnames(res.df.all.p[[1]])[(which(colnames(res.df.all.p[[1]])=="names")+1):(which(colnames(res.df.all.p[[1]])=="chromosome_name")-1)]
  
  for (iter in 1:3) { 
    print(iter)
    df.lfc=res.df.all.lfc[[iter]]
    df.p=res.df.all.p[[iter]]
    
    for (gene_set in gene_set_list) {
      print(gene_set)
      f.gene_set <- paste0(dir,"/scripts/single_cell/bin/",gene_set,".v7.5.symbols.gmt")
      pathways <- gmtPathways(f.gene_set)
      
      # system(paste0("mkdir -p ",dir,"/output/data/",DATASET,"/gsea_RNA"))
      # system(paste0("mkdir -p ",dir,"/output/data/",DATASET,"/gsea_RNA/",gene_set,".v7.5/"))
      outdir = "~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/gsea"
      system(paste0("mkdir -p ",outdir,"/",gene_set,".v7.5"))
      
      for (i in 1:length(cell_type_groups)) {
        
        for (subset_to_use in c("not_chr21_up","not_chr21_down","chr21_up","chr21_down")) {
          print(subset_to_use)
          cell_type = cell_type_groups[i]
          print(paste0("GSEA - gene set: ",gene_set,' (',cell_type,', ',i,'/',length(cell_type_groups),')'))
          if (subset_to_use=="not_chr21_up") {
            x <- df.lfc$names[df.lfc[,cell_type] > 0 & df.p[,cell_type] < 0.1 & df.lfc[,"chr21"]=="Not Chr 21"]
            f.out <- paste0(outdir,"/",gene_set,".v7.5/gsea.not_chr21.up.iter",iter,".",cell_type,".",sampletype,".txt")
          } else if (subset_to_use=="not_chr21_down") {
            x <- df.lfc$names[df.lfc[,cell_type] < 0 & df.p[,cell_type] < 0.1 & df.lfc[,"chr21"]=="Not Chr 21"]
            f.out <- paste0(outdir,"/",gene_set,".v7.5/gsea.not_chr21.down.iter",iter,".",cell_type,".",sampletype,".txt")
          } else if (subset_to_use=="chr21_up") {
            x <- df.lfc$names[df.lfc[,cell_type] > 0 & df.p[,cell_type] < 0.1 & df.lfc[,"chr21"]=="Chr 21"]
            f.out <- paste0(outdir,"/",gene_set,".v7.5/gsea.chr21.up.iter",iter,".",cell_type,".",sampletype,".txt")
          } else if (subset_to_use=="chr21_down") {
            x <- df.lfc$names[df.lfc[,cell_type] < 0 & df.p[,cell_type] < 0.1 & df.lfc[,"chr21"]=="Chr 21"]
            f.out <- paste0(outdir,"/",gene_set,".v7.5/gsea.chr21.down.iter",iter,".",cell_type,".",sampletype,".txt")
          }
          
          x <- x[!is.na(x)]
          df.sub <- df.lfc[df.lfc$names %in% x,c(cell_type,"names")]
          ranks_stats <- df.sub[,cell_type]
          names(ranks_stats) <- df.sub[,"names"]
          
          gsea_results2 <- tryCatch({
            gsea_results <- fgsea(pathways=pathways, stats=ranks_stats,scoreType="pos")
            gsea_results2 <- gsea_results[order(gsea_results$padj,gsea_results$pval),]; # gsea_results2[1:5,]
          },error=function(cond) {
            gsea_results2 <- data.frame(pathway=NA,pval=NA,padj=NA,log2err=NA,ES=NA,NES=NA,size=NA,leadingEdge=NA)
          })
          
          print(paste0("Writing: ",f.out))
          fwrite(gsea_results2,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
        }
      }
    }
  }
}
  