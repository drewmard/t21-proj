# library(biomaRt)
library(data.table)
# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",mirror = "useast")
# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

data_snv = fread("~/Downloads/data_snv",data.table=F,stringsAsFactors=F)
colnames(data_snv) = c("chrom","start","end","sample","genotype","snv")
data_snv = subset(data_snv,genotype != "T21_MyeloidPreleuk")

###############

disease_status="DownSyndrome"
pseudobulk_covariate="sample"
suffix=""
f = paste0("~/Downloads/",disease_status,'_',"cyc_vs_hsc",".de",suffix,".txt")
atac=fread(f,data.table = F,stringsAsFactors = F)[,1:7]

###########

cre = fread("~/Downloads/GRCh38-cCREs.bed",data.table=F,stringsAsFactors = F)
# cre = fread("/oak/stanford/groups/smontgom/amarder/data/ENCFF394DBM.bed",data.table=F,stringsAsFactors = F)
colnames(cre)[1:3] = c("chrom",'start','end')

############################################

annot <- getBM(attributes = c('hgnc_symbol','chromosome_name',
                              'start_position', 'end_position'),
               filters = 'hgnc_symbol',
               values = atac$names,
               mart = ensembl)
annot = subset(annot,chromosome_name %in% c(1:22,"X","Y"))
annot = annot[!duplicated(annot$hgnc_symbol),]
colnames(annot)[2:4] = c("chrom","start","end")
annot$chrom = paste0("chr",annot$chrom)

###################

df.mg = as.data.frame(valr::bed_intersect(data_snv,annot))
tmp = df.mg[,c("snv.x","hgnc_symbol.y")]
colnames(tmp) = c("snv","hgnc_symbol")
snv.df = merge(data_snv,tmp,by="snv",all.x=TRUE)
snv.df = snv.df[order(is.na(snv.df$hgnc_symbol)),]
snv.df = snv.df[!duplicated(snv.df$snv),]

# snvs in exons
annot <- getBM(attributes = c('hgnc_symbol','ensembl_gene_id'),
               filters = 'hgnc_symbol', 
               values = atac$names, 
               mart = ensembl)
annot2 <- getBM(attributes = c('ensembl_gene_id','chromosome_name',
                               'exon_chrom_start',"exon_chrom_end"),
                filters = 'ensembl_gene_id', 
                values = annot$ensembl_gene_id, 
                mart = ensembl)
annot2 = subset(annot2,chromosome_name %in% c(1:22,"X","Y"))
annot = merge(annot,annot2,by="ensembl_gene_id")
colnames(annot)[3:5] = c("chrom","start","end")
annot$chrom = paste0("chr",annot$chrom)
fwrite(annot,"~/Downloads/hgnc_ensembl_exon.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

snv.df.exon = as.data.frame(valr::bed_intersect(annot,snv.df))
snv.df.exon2 = snv.df.exon[,c("hgnc_symbol.x","chrom","start.y","end.y","sample.y","genotype.y","snv.y")]
colnames(snv.df.exon2) = c("hgnc_symbol","chrom","start","end","sample","genotype","snv")
ind = !duplicated(snv.df.exon2$snv); sum(!ind); snv.df.exon2=snv.df.exon2[ind,]
snv.df$exon = snv.df$snv %in% snv.df.exon2$snv
snv.df$gene = !is.na(snv.df$hgnc_symbol)
snv.df$intron = snv.df$gene & !snv.df$exon

# snv analysis:
snv.df= subset(snv.df,intron)
# snv.df.use= subset(snv.df,intron)

snv.df=snv.df[,c("chrom","start","end")]
# snv.df=df.sub[,1:3]
colnames(snv.df)[1:3] = c("chrom","start","end")
snv.df$posid = paste0("snv",1:nrow(snv.df))
snv.df$genotype = 1
snv.df$strand = "+"
# snv.df$start=snv.df$start-10
# snv.df$end=snv.df$end+10

df.mg = as.data.frame(valr::bed_intersect(cre,snv.df))
print(nrow(df.mg))
print(nrow(snv.df))

snv.df$snv = paste(snv.df$chrom,snv.df$start,snv.df$end,sep = "-")
df.mg$snv = paste(df.mg$chrom,df.mg$start.y,df.mg$end.y,sep = "-")
data_snv$intron = data_snv$snv %in% snv.df$snv
data_snv$intron_cre = data_snv$snv %in% df.mg$snv
tmp = data_snv[data_snv$intron,]

table(subset(data_snv,intron & genotype=="T21")$intron_cre)
table(subset(data_snv,intron & genotype=="D21")$intron_cre)

fisher.test(tmp$genotype,tmp$intron_cre)


k=1
if (k==1) {
  id="D21"
  df.sub=subset(data_snv,genotype==id)
} else if (k==2) {
  id="T21"
  df.sub=subset(data_snv,genotype==id)
}

ind = !duplicated(df.sub$snv); sum(!ind); df.sub=df.sub[ind,]
print(paste0(k," - ",id))

snv.df=df.sub
colnames(snv.df)[1:3] = c("chrom","start","end")

annot <- getBM(attributes = c('hgnc_symbol','chromosome_name',
                              'start_position', 'end_position'),
               filters = 'hgnc_symbol',
               values = atac$names,
               mart = ensembl)
annot = subset(annot,chromosome_name %in% c(1:22,"X","Y"))
annot = annot[!duplicated(annot$hgnc_symbol),]
colnames(annot)[2:4] = c("chrom","start","end")
annot$chrom = paste0("chr",annot$chrom)


