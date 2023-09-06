library(data.table)
library(ggplot2)
library(cowplot)
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

for (cell_type in cell_type_groups) {
  keep_genes <- res.df.all.p[[2]][!is.na(res.df.all.p[[2]][,cell_type]),'names']
  res.df.all.p[[1]][!(res.df.all.p[[1]][,'names']%in%keep_genes),cell_type] <- NA
  res.df.all.lfc[[1]][!(res.df.all.lfc[[1]][,'names']%in%keep_genes),cell_type] <- NA
}

tab <- lapply(res.df.all.p,function(x) {aggregate(x[,cell_type_groups]<0.1,by=list(x$chr21),mean,na.rm=T)})
tab <- lapply(tab,function(x) {rownames(x) <- x[,1]; x <- x[,-1]; x <- as.data.frame(t(x))})
tab <- do.call(cbind,tab)

apply(tab,2,median)

# cor(tab[,c(2,4,6)])
# pairs(tab[,c(2,4,6)])
# cor(tab[,c(1,3,5)])
# pairs(tab[,c(1,3,5)])

colnames(tab)[1:2] <- paste0(colnames(tab)[1:2]," 1")
colnames(tab)[3:4] <- paste0(colnames(tab)[3:4]," 2")
colnames(tab)[5:6] <- paste0(colnames(tab)[5:6]," 3")

g1 <- ggplot(tab,aes(x=100*`Chr 21 1`,y=100*`Chr 21 2`)) + 
  geom_point() + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.title=element_text(hjust=0.5)) +
  geom_smooth(method='loess',span=3,se=F,col='purple') + 
  geom_abline(slope=1,intercept = 0,col='red',lty='dashed') +
  labs(x="% DEG per cell type (fetus)",y="% DEG per cell type (sample)",title="Chr 21")
g2 <- ggplot(tab,aes(x=100*`Not Chr 21 1`,y=100*`Not Chr 21 2`)) + 
  geom_point() + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.title=element_text(hjust=0.5)) +
  geom_smooth(method='loess',span=3,se=F,col='purple') + 
  geom_abline(slope=1,intercept = 0,col='red',lty='dashed') +
  labs(x="% DEG per cell type (fetus)",y="% DEG per cell type (sample)",title="Not Chr 21")
plot_grid(g1,g2,ncol=2)

# generally: more chr 21 dysregulation is linked to more non chr 21 dysregulation
# cor(tab)
# plot(tab[,c(2,1)])
# plot(tab[,c(4,3)])
# plot(tab[,c(6,5)])

ggplot(tab,aes(x=100*`Chr 21 2`,y=100*`Not Chr 21 2`)) + 
  geom_point() + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.title=element_text(hjust=0.5)) +
  geom_smooth(method='loess',span=3,se=F,col='purple') + 
  labs(x="% DEG per cell type (chr 21)",y="% DEG per cell type (not chr 21)",title="PB per-sample")

tab <- lapply(res.df.all.lfc,function(x) {aggregate(apply(apply(x[,cell_type_groups],2,rank),2,function(x) x/sum(!is.na(x))),by=list(x$chr21),median,na.rm=T)})
tab <- lapply(res.df.all.lfc,function(x) {aggregate(x[,cell_type_groups],by=list(x$chr21),median,na.rm=T)})
tab <- lapply(tab,function(x) {rownames(x) <- x[,1]; x <- x[,-1]; x <- as.data.frame(t(x))})
tab <- do.call(cbind,tab)

apply(tab,2,median)

colnames(tab)[1:2] <- paste0(colnames(tab)[1:2]," 1")
colnames(tab)[3:4] <- paste0(colnames(tab)[3:4]," 2")
colnames(tab)[5:6] <- paste0(colnames(tab)[5:6]," 3")

g1 <- ggplot(tab,aes(x=`Chr 21 1`,y=`Chr 21 2`)) + 
  geom_point() + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.title=element_text(hjust=0.5)) +
  geom_smooth(method='loess',span=3,se=F,col='purple') + 
  geom_abline(slope=1,intercept = 0,col='red',lty='dashed') +
  labs(x="Median LFC per cell type (fetus)",y="Median LFC per cell type (sample)",title="Chr 21")
g2 <- ggplot(tab,aes(x=`Not Chr 21 1`,y=`Not Chr 21 2`)) + 
  geom_point() + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.title=element_text(hjust=0.5)) +
  geom_smooth(method='loess',span=3,se=F,col='purple') + 
  geom_abline(slope=1,intercept = 0,col='red',lty='dashed') +
  labs(x="Median LFC per cell type (fetus)",y="Median LFC per cell type (sample)",title="Not Chr 21")
plot_grid(g1,g2,ncol=2)
g1

########

# sampletype="Femur"
sampletype="Liver"
res.df.all.p <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".p.rds"))
res.df.all.lfc <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".lfc.rds"))
cell_type_groups <- colnames(res.df.all.p[[1]])[(which(colnames(res.df.all.p[[1]])=="names")+1):(which(colnames(res.df.all.p[[1]])=="chromosome_name")-1)]

for (cell_type in cell_type_groups) {
  keep_genes <- res.df.all.p[[2]][!is.na(res.df.all.p[[2]][,cell_type]),'names']
  res.df.all.p[[1]][!(res.df.all.p[[1]][,'names']%in%keep_genes),cell_type] <- NA
  res.df.all.lfc[[1]][!(res.df.all.lfc[[1]][,'names']%in%keep_genes),cell_type] <- NA
}

tab <- lapply(res.df.all.lfc,function(x) {aggregate(x[,cell_type_groups],by=list(x$chr21),median,na.rm=T)})
tab <- lapply(tab,function(x) {rownames(x) <- x[,1]; x <- x[,-1]; x <- as.data.frame(t(x))})
tab <- do.call(cbind,tab)
tab <- tab[,c(2,4)]
colnames(tab) <- c("Fetus","Sample")
tab$cell <- rownames(tab)
tab <- reshape2::melt(tab,id.vars="cell")

ggplot(tab,aes(x=cell,fill=variable,y=value)) + 
  geom_bar(stat = 'identity',col='black',position = position_dodge2(preserve = 'single')) + 
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=55,hjust = 1)
  ) +
  labs(x="Cell type",y="Median LFC (T21 vs Healthy)",fill='Pseudobulk') 



# melt(tab[,1:2])
# ggplot(tab,aes(x=`Chr 21 1`,y=`Chr 21 2`)) + 
#   geom_boxplot() + 
#   theme_bw() + 
#   theme(panel.grid = element_blank(),
#         plot.title=element_text(hjust=0.5)) +
#   geom_smooth(method='loess',span=3,se=F,col='purple') + 
#   labs(x="% DEG per cell type (chr 21)",y="% DEG per cell type (not chr 21)",title="PB per-sample")
# plot_grid(g1,g2,ncol=2)


# does the random effect model increase power on non chr 21?
t.test(tab[,5] - tab[,1])
# no - significantly decreased power!
# does the fixed effect model increase power on non chr 21?
t.test(tab[,3] - tab[,1])
# no sig diff!

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

apply(tab,2,median)

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

apply(tab,2,median)

colnames(tab)[1:2] <- paste0(colnames(tab)[1:2]," 1")
colnames(tab)[3:4] <- paste0(colnames(tab)[3:4]," 2")
colnames(tab)[5:6] <- paste0(colnames(tab)[5:6]," 3")

g1 <- ggplot(tab,aes(x=`Chr 21 1`,y=`Chr 21 2`)) + 
  geom_point() + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.title=element_text(hjust=0.5)) +
  geom_smooth(method='loess',span=3,se=F,col='purple') + 
  geom_abline(slope=1,intercept = 0,col='red',lty='dashed') +
  labs(x="Median LFC per cell type (fetus)",y="Median LFC per cell type (sample)",title="Chr 21")
g2 <- ggplot(tab,aes(x=`Not Chr 21 1`,y=`Not Chr 21 2`)) + 
  geom_point() + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.title=element_text(hjust=0.5)) +
  geom_smooth(method='loess',span=3,se=F,col='purple') + 
  geom_abline(slope=1,intercept = 0,col='red',lty='dashed') +
  labs(x="Median LFC per cell type (fetus)",y="Median LFC per cell type (sample)",title="Not Chr 21")
plot_grid(g1,g2,ncol=2)

# generally, more upregulation on chr 21 is linked to more upregulation on other chr
# neg corr seen in femur is due to a single outlier
cor(tab)
par(mfrow=c(1,3))
plot(tab[,c(2,1)])
plot(tab[,c(4,3)])
plot(tab[,c(6,5)])
par(mfrow=c(1,1))

iter=2
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

iter=1
df1.p <- res.df.all.p[[iter]]
df1.lfc <- res.df.all.lfc[[iter]]

df1.p.melt <- reshape2::melt(aggregate(df1.p[,cell_type_groups]<0.1,by=list(df1.p$chromosome_name),mean,na.rm=T),id.vars="Group.1")
df1.p.melt <- reshape2::melt(aggregate(df1.p[,cell_type_groups]<0.1 & df1.lfc[,cell_type_groups] > 0,by=list(df1.p$chromosome_name),mean,na.rm=T),id.vars="Group.1")
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

sampletype="Femur"
res.df.all.p <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".p.rds"))
res.df.all.lfc <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".lfc.rds"))
cell_type_groups <- colnames(res.df.all.p[[1]])[(which(colnames(res.df.all.p[[1]])=="names")+1):(which(colnames(res.df.all.p[[1]])=="chromosome_name")-1)]

iter=1
df1.p <- res.df.all.p[[iter]]
df1.lfc <- res.df.all.lfc[[iter]]

tab <- lapply(res.df.all.lfc,function(x) {aggregate(x[,cell_type_groups],by=list(x$chr21),median,na.rm=T)})
tab <- lapply(tab,function(x) {rownames(x) <- x[,1]; x <- x[,-1]; x <- as.data.frame(t(x))})
tab <- do.call(cbind,tab)
tab <- tab[,c(2,4)]
colnames(tab) <- c("Fetus","Sample")
tab$cell <- rownames(tab)
tab <- reshape2::melt(tab,id.vars="cell")

ggplot(tab,aes(x=cell,fill=variable,y=value)) + 
  geom_bar(stat = 'identity',col='black',position = position_dodge2(preserve = 'single')) + 
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=55,hjust = 1)
  ) +
  labs(x="Cell type",y="Median LFC (T21 vs Healthy)",fill='Pseudobulk') 

iter=1
tab <- subset(res.df.all.lfc[[iter]],names=="SOD1")
tab <- reshape2::melt(tab[,cell_type_groups])
iter=2
tab2 <- subset(res.df.all.lfc[[iter]],names=="SOD1")
tab2 <- reshape2::melt(tab2[,cell_type_groups])
tab$sample <- "Fetus"
tab2$sample <- "Sample"
tab <- rbind(tab,tab2)

ggplot(tab,aes(x=variable,fill=sample,y=value)) + 
  geom_bar(stat = 'identity',col='black',position = position_dodge2(preserve = 'single')) + 
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=55,hjust = 1),
    plot.title = element_text(hjust=0.5)
  ) +
  labs(x="Cell type",y="LFC (T21 vs Healthy)",fill='Pseudobulk',title="SOD1")

iter=2
subset(res.df.all.lfc[[iter]],names=="SOD1")

#SOD1 downregulated in late erythroid cells in femur in analysis 1... but upregulated in analysis 2 and 3


##################

# dbs <- listEnrichrDbs()

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
iter=2

sampletype="Femur"
for (sampletype in c("Femur","Liver")) {
  print(sampletype)
  res.df.all.p <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".p.rds"))
  res.df.all.lfc <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".lfc.rds"))
  cell_type_groups <- colnames(res.df.all.p[[1]])[(which(colnames(res.df.all.p[[1]])=="names")+1):(which(colnames(res.df.all.p[[1]])=="chromosome_name")-1)]
  
  for (iter in 2:2) { 
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
            gene_lst <- names(ranks_stats)
            dbs <- c("GO_Molecular_Function_2018", "GO_Biological_Process_2018","ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")
            enriched <- enrichr(gene_lst, dbs)
            y <- enriched[[dbs[1]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y <- subset(y,geneCt>3 & Adjusted.P.value<0.1)
            y <- enriched[[dbs[2]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y <- subset(y,geneCt>3 & Adjusted.P.value<0.1)
            y <- enriched[[dbs[3]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y <- subset(y,geneCt>3 & Adjusted.P.value<0.1)
            
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
  