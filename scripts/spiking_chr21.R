library(data.table)
library(Matrix.utils)
# library(Seurat)
library(DESeq2)
library(variancePartition)
library('edgeR')

print("Reading metadata2...")
disease_status="Healthy"
f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cellComp/10X_",disease_status,"_","Liver",".cellComp.csv")
meta2.full<-fread(f,data.table = F,stringsAsFactors = F)
colnames(meta2.full)[6] <- "leiden_names"
cells2 <- unique(meta2.full[,6])
meta2 <- unique(meta2.full[,c("patient","sample","sorting")])
meta2$environment <- disease_status
meta2$sorting[!(meta2$sorting %in% c("CD235a-","CD45+"))] <- "Other"
x <- meta2
rownames(x) <- x[,"sample"]
rm(meta2)

x = subset(x,patient %in% c("15633","15781"))

#######

cell_type="HSCs/MPPs"
cell_type_filename = gsub("/","_",cell_type)
f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data_pb_leiden/","Liver",".pb.",cell_type_filename,".","sample",".txt")

df.aggre <- fread(f,data.table = F,stringsAsFactors = F)
rownames(df.aggre) <- df.aggre[,1]; df.aggre <- df.aggre[,-1]

metadata_to_use <- x[rownames(x) %in% colnames(df.aggre),]
to_merge <- names(table(metadata_to_use$sorting)[table(metadata_to_use$sorting) < 5])
df.aggre <- as.matrix(df.aggre[,match(rownames(metadata_to_use),colnames(df.aggre))])

tab <- table(meta2.full[meta2.full[,"leiden_names"]==cell_type,'sample'])
samples_to_keep <- names(tab)[tab >= 10]
samples_to_keep <- samples_to_keep[samples_to_keep %in% rownames(metadata_to_use)]
df.aggre <- df.aggre[,samples_to_keep]
metadata_to_use <- metadata_to_use[samples_to_keep,]
colnames(df.aggre) <- paste0(colnames(df.aggre),".hsc")
rownames(metadata_to_use) <- paste0(rownames(metadata_to_use),".hsc")
metadata_to_use$leiden_names <- cell_type

###########

metadata_to_use$patient = as.factor(metadata_to_use$patient)
geneExpr = DGEList( df.aggre )
keep <- filterByExpr(geneExpr, group=metadata_to_use$patient)
geneExpr <- geneExpr[keep,]
geneExpr = calcNormFactors( geneExpr )

# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param = SnowParam(8, "SOCK", progressbar=TRUE)

# The variable to be tested must be a fixed effect
form <- ~ patient #+ patient

# A positive FC is increased expression in the DS compared to healthy
# metadata_to_use$leiden_names <- factor(metadata_to_use$leiden_names,levels=c("HSCs/MPPs","Cycling HSCs/MPPs"))

# estimate weights using linear mixed model of dream
vobjDream = voomWithDreamWeights( geneExpr, form, metadata_to_use, BPPARAM=param )

# Fit the dream model on each gene
# By default, uses the Satterthwaite approximation for the hypothesis test
fitmm = dream( vobjDream, form, metadata_to_use )
fitmm = eBayes(fitmm)

# reorganize:
# res.df <- as.data.frame(topTable( fitmm, number=Inf))
res.df <- as.data.frame(topTable( fitmm, number=Inf,coef = "patient15781" ))
res.df$names <- rownames(res.df)

############

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

annot <- getBM(attributes = c('hgnc_symbol', 'chromosome_name',
                              'start_position', 'end_position'),
               filters = 'hgnc_symbol', 
               values = res.df$names, 
               mart = ensembl)
annot <- annot[!duplicated(annot$hgnc_symbol),]
res.df.mg <- merge(res.df,annot[,c("hgnc_symbol","chromosome_name","start_position")],by.x="names",by.y="hgnc_symbol")
aggregate(res.df.mg$logFC,list(chr21=res.df.mg$chromosome_name),median)

library(ggplot2)

res.df.mg = subset(res.df.mg,chromosome_name %in% seq(1,22))
g=ggplot(res.df.mg,aes(x=factor(chromosome_name,levels =seq(1,22)),y=logFC,fill=chromosome_name)) + geom_boxplot() + 
  geom_point(alpha=0.5) + theme_bw() + theme(panel.grid=element_blank()) +
  guides(fill="none") + labs(x="Chromosome",y="LFC")
f.plot="/home/amarder/tmp/chr_real.pdf"
pdf(f.plot,width=6,height=4)
print(g)
dev.off()

tab = aggregate(res.df.mg$logFC,by=list(chr=res.df.mg$chromosome_name),median)
tab$chr21 = ifelse(tab$chr==21,"Chrom 21","Other chrom")
g <- ggplot(tab,aes(x=chr,fill=chr21,y=x)) + 
  geom_bar(stat = 'identity',col='black',position = position_dodge2(preserve = 'single')) + 
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=55,hjust = 1)
  ) +
  labs(x="Chromosome",y="Median LFC",fill='') + ylim(-0.2,0.6)
f.plot="/home/amarder/tmp/chr_real.box.pdf"
pdf(f.plot,width=6,height=4)
print(g)
dev.off()

# save DEG results from real
f.output="/home/amarder/tmp/real.de.txt"
fwrite(res.df.mg,f.output,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)


##############################


ind = res.df.mg$names[res.df.mg$chromosome_name==21]
df.aggre[ind,3:4] = 1.5*df.aggre[ind,3:4] 

geneExpr = DGEList( df.aggre )
keep <- filterByExpr(geneExpr, group=metadata_to_use$patient)
geneExpr <- geneExpr[keep,]
geneExpr = calcNormFactors( geneExpr )

# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param = SnowParam(8, "SOCK", progressbar=TRUE)

# The variable to be tested must be a fixed effect
form <- ~ patient #+ patient

# A positive FC is increased expression in the DS compared to healthy
# metadata_to_use$leiden_names <- factor(metadata_to_use$leiden_names,levels=c("HSCs/MPPs","Cycling HSCs/MPPs"))

# estimate weights using linear mixed model of dream
vobjDream = voomWithDreamWeights( geneExpr, form, metadata_to_use, BPPARAM=param )

# Fit the dream model on each gene
# By default, uses the Satterthwaite approximation for the hypothesis test
fitmm = dream( vobjDream, form, metadata_to_use )
fitmm = eBayes(fitmm)

# reorganize:
# res.df <- as.data.frame(topTable( fitmm, number=Inf))
res.df <- as.data.frame(topTable( fitmm, number=Inf,coef = "patient15781" ))
res.df$names <- rownames(res.df)

############

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

annot <- getBM(attributes = c('hgnc_symbol', 'chromosome_name',
                              'start_position', 'end_position'),
               filters = 'hgnc_symbol', 
               values = res.df$names, 
               mart = ensembl)
annot <- annot[!duplicated(annot$hgnc_symbol),]
res.df.mg <- merge(res.df,annot[,c("hgnc_symbol","chromosome_name","start_position")],by.x="names",by.y="hgnc_symbol")
aggregate(res.df.mg$logFC,list(chr21=res.df.mg$chromosome_name),median)

library(ggplot2)

res.df.mg = subset(res.df.mg,chromosome_name %in% seq(1,22))
g=ggplot(res.df.mg,aes(x=factor(chromosome_name,levels =seq(1,22)),y=logFC,fill=chromosome_name)) + geom_boxplot() + 
  geom_point(alpha=0.5) + theme_bw() + theme(panel.grid=element_blank()) +
  guides(fill="none") + labs(x="Chromosome",y="LFC")
f.plot="/home/amarder/tmp/chr_fake.pdf"
pdf(f.plot,width=6,height=4)
print(g)
dev.off()

aggregate(res.df.mg$logFC,by=list(chr21=res.df.mg$chromosome_name==21),median)
tab = aggregate(res.df.mg$logFC,by=list(chr=res.df.mg$chromosome_name),median)
tab$chr21 = ifelse(tab$chr==21,"Chrom 21","Other chrom")
g <- ggplot(tab,aes(x=chr,fill=chr21,y=x)) + 
  geom_bar(stat = 'identity',col='black',position = position_dodge2(preserve = 'single')) + 
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=55,hjust = 1)
  ) +
  labs(x="Chromosome",y="Median LFC",fill='') + ylim(-0.2,0.6)
f.plot="/home/amarder/tmp/chr_fake.box.pdf"
pdf(f.plot,width=6,height=4)
print(g)
dev.off()

# save DEG results from fake
f.output="/home/amarder/tmp/fake.de.txt"
fwrite(res.df.mg,f.output,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)

##############################################################################
##############################################################################
##############################################################################

f.output="~/Downloads/real.de.txt"
de1 = fread(f.output,data.table = F,stringsAsFactors = F)
f.output="~/Downloads/fake.de.txt"
de2 = fread(f.output,data.table = F,stringsAsFactors = F)

de.mg = merge(de1,de2,by='names')
library(ggplot2)
de.mg$chr21 = ifelse(de.mg$chromosome_name.x==21,"Chrom 21","Other chrom")
g=ggplot(de.mg,aes(x=logFC.x,y=logFC.y,col=chr21)) + geom_point() + geom_abline(slope=1,intercept=0,lty='dashed') + theme_bw() + 
  labs(x='LFC (real data)',y='LFC (spiked chr21)',col='') + theme(panel.grid = element_blank()) + scale_color_manual(values=c("orange","black"))
pdf("~/Downloads/chr21spike_lfc_corr.pdf",width=5.5,height=4)
print(g)
dev.off()
