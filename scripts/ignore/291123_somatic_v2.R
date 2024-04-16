# library(biomaRt)
library(data.table)
# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",mirror = "useast")
# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

data_snv = fread("~/Downloads/data_snv",data.table=F,stringsAsFactors=F)
colnames(data_snv) = c("chrom","start","end","sample","genotype","snv")
data_snv = subset(data_snv,genotype != "T21_MyeloidPreleuk")

###############

# atac = fread("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/HSC_v_CyclingHSC.txt",data.table = F,stringsAsFactors = F)

disease_status="DownSyndrome"
pseudobulk_covariate="sample"
suffix=""
f = paste0("~/Downloads/",disease_status,'_',"cyc_vs_hsc",".de",suffix,".txt")
atac=fread(f,data.table = F,stringsAsFactors = F)[,1:7]
atac = subset(atac,P.Value < 0.05)

f="~/Documents/Research/t21-proj/out/full/DEG_list/liver_v_femur.txt"
atac=fread(f,data.table = F,stringsAsFactors = F)
# atac = subset(atac,class=="t21-induced")
atac = subset(atac,P.Value.t21 < 0.05)

###########

cre = fread("~/Downloads/GRCh38-cCREs.bed",data.table=F,stringsAsFactors = F)
# cre = fread("/oak/stanford/groups/smontgom/amarder/data/ENCFF394DBM.bed",data.table=F,stringsAsFactors = F)
colnames(cre)[1:3] = c("chrom",'start','end')

############################################

# map names to ENST
gene_symbols = fread("~/Downloads/gene_symbols.txt",data.table = F,stringsAsFactors = F)


############################################

# intron = fread("~/Downloads/introns",data.table = F,stringsAsFactors = F)
# colnames(intron)[1:3] = c("chrom",'start','end')
# tmp = as.data.frame(valr::bed_intersect(data_snv,intron))
# tmp = tmp[!duplicated(tmp$snv.x),]
# tmp = tmp[,c("chrom","start.x","end.y","sample.x","genotype.x","snv.x","V4.y")]
# colnames(tmp) = c(colnames(data_snv),"intron_info")
# tmp$gene_info = gsub("^([^\\.]+)\\..*$", "\\1",tmp[,7])

# tmp2 = as.data.frame(valr::bed_intersect(data_snv,cre))
# tmp2 = tmp2[!duplicated(tmp2$snv.x),]

intron = fread("~/Downloads/introns",data.table = F,stringsAsFactors = F)
colnames(intron)[1:3] = c("chrom",'start','end')
tmp = as.data.frame(valr::bed_intersect(data_snv,intron))
tmp = tmp[!duplicated(tmp$snv.x),]
tmp = tmp[,c("chrom","start.x","end.y","sample.x","genotype.x","snv.x","V4.y")]
colnames(tmp) = c(colnames(data_snv),"intron_info")
tmp$gene_info = gsub("^([^\\.]+)\\..*$", "\\1",tmp[,7])
tmp = merge(tmp,gene_symbols,by.x="gene_info",by.y="Transcript stable ID")
tmp <- tmp[, c(2:ncol(tmp), 1)]
tmp$intron = TRUE
data_snv = merge(data_snv,tmp[,c("snv","Gene stable ID","Gene name","gene_info","intron")],by="snv",all.x=TRUE)
data_snv$intron[is.na(data_snv$intron)] = FALSE
data_snv <- data_snv[, c(2:ncol(data_snv), 1)]

exon = fread("~/Downloads/exons",data.table = F,stringsAsFactors = F)
colnames(exon)[1:3] = c("chrom",'start','end')
tmp = as.data.frame(valr::bed_intersect(data_snv,exon))
tmp = tmp[!duplicated(tmp$snv.x),]
tmp = tmp[,c("chrom","start.x","end.y","sample.x","genotype.x","snv.x","V4.y")]
colnames(tmp) = c(colnames(data_snv)[1:5],"snv","exon_info")
tmp$gene_info = gsub("^([^\\.]+)\\..*$", "\\1",tmp[,7])
tmp = merge(tmp,gene_symbols,by.x="gene_info",by.y="Transcript stable ID")
tmp <- tmp[, c(2:ncol(tmp), 1)]
tmp$exon = TRUE
data_snv = merge(data_snv,tmp[,c("snv","exon")],by="snv",all.x=TRUE)
data_snv$exon[is.na(data_snv$exon)] = FALSE

tmp2 = as.data.frame(valr::bed_intersect(data_snv,cre))
tmp2 = tmp2[!duplicated(tmp2$snv.x),]
data_snv$cre = data_snv$snv %in% tmp2$snv.x

#####
PROP=0.5
TPM=0.2
f = "~/Downloads/Healthy_Liver_HSCs_MPPs.pb.txt"
x = fread(f,data.table = F,stringsAsFactors = F)
rownames(x)=x[,1];x=x[,-1]
x = cpm(x)
keep = apply(x,1,function(x){mean(x>TPM) > PROP})
genes_to_keep = rownames(x)[keep]
data_snv$hsc_d21 = data_snv$`Gene name` %in% genes_to_keep

f = "~/Downloads/DownSyndrome_Liver_HSCs_MPPs.pb.txt"
x = fread(f,data.table = F,stringsAsFactors = F)
rownames(x)=x[,1];x=x[,-1]
x = cpm(x)
keep = apply(x,1,function(x){mean(x>TPM) > PROP})
genes_to_keep = rownames(x)[keep]
data_snv$hsc_t21 = data_snv$`Gene name` %in% genes_to_keep

f = "~/Downloads/DownSyndrome_Liver_Cycling HSCs_MPPs.pb.txt"
x = fread(f,data.table = F,stringsAsFactors = F)
rownames(x)=x[,1];x=x[,-1]
x = cpm(x)
keep = apply(x,1,function(x){mean(x>TPM) > PROP})
genes_to_keep = rownames(x)[keep]
data_snv$cychsc_t21 = data_snv$`Gene name` %in% genes_to_keep
# fisher.test(data_snv$hsc_d21,data_snv$hsc_t21)
# fisher.test(data_snv$hsc_t21,data_snv$cychsc_t21)
# fisher.test(data_snv$hsc_d21,data_snv$cychsc_t21)

pb_da_peaks = fread("~/Downloads/pb_de_atac.hsc.txt",data.table = F,stringsAsFactors = F)
y = strsplit(pb_da_peaks$names,"-")
pb_da_peaks$chrom = unlist(lapply(y,function(x)x[[1]]))
pb_da_peaks$start = as.numeric(unlist(lapply(y,function(x)x[[2]])))
pb_da_peaks$end = as.numeric(unlist(lapply(y,function(x)x[[3]])))
pb_da_peaks = pb_da_peaks[,c('chrom','start','end','names','logFC','P.Value')]
tmp = as.data.frame(valr::bed_intersect(data_snv,pb_da_peaks))
tmp = tmp[!duplicated(tmp$snv.x),]
tmp = subset(tmp,intron.x & !exon.x & (cychsc_t21.x | hsc_t21.x))
table(tmp$genotype.x,tmp$P.Value.y < 0.001 & tmp$logFC.y > 0.25)
fisher.test(tmp$genotype.x,tmp$cre.x)

data_snv2 = merge(data_snv,tmp[,c("snv","exon")],by="snv",all.x=TRUE)

colnames(intron)[1:3] = c("chrom",'start','end')
tmp = as.data.frame(valr::bed_intersect(data_snv,intron))


tmp = subset(data_snv,!exon & intron)
fisher.test(tmp$cre,tmp$genotype)

# genes$name = gsub("^([^\\.]+)\\..*$", "\\1",genes$V4)
# genes2 = merge(gene_symbols,genes,by.x="Transcript stable ID",by.y="name")
# genes2 = unique(genes2[,c("Gene name","chrom","start","end")])
genes2 = fread("~/Downloads/gene_chr_pos.txt",data.table = F,stringsAsFactors = F)
genes2 = genes2[!duplicated(genes2$`Gene name`),]

atac2 = merge(atac,genes2,by.x="names",by.y="Gene name")

# tmp = subset(data_snv,!exon & intron & `Gene name` %in% atac2$names)
tmp = subset(data_snv,!exon & intron & `Gene name` %in% atac$names)
tab = table(tmp$cre,tmp$genotype)
binom.test(tab[2,2],tab[1,2]+tab[2,2],(tab[2,1])/(tab[2,1]+tab[1,1]))
# fisher.test(tmp$cre,tmp$genotype)

tmp = subset(data_snv,!exon & intron & !(`Gene name` %in% atac$names))
tab = table(tmp$cre,tmp$genotype)
binom.test(tab[2,2],tab[1,2]+tab[2,2],(tab[2,1])/(tab[2,1]+tab[1,1]))

tmp = subset(data_snv,!exon & intron & nchar(`Gene name`) > 0 )
tmp$diffactivity = tmp$`Gene name` %in% atac$names
# fisher.test(tmp$cre,tmp$diffactivity)
tab = table(tmp$diffactivity,tmp$cre)
binom.test(tab[2,2],tab[1,2]+tab[2,2],(tab[2,1])/(tab[2,1]+tab[1,1]))

tmp = subset(data_snv,!exon & intron & nchar(`Gene name`) > 0 )
# tmp$diffactivity = tmp$`Gene name` %in% subset(atac,adj.P.Val < 0.1)$names
tmp$diffactivity = tmp$`Gene name` %in% subset(atac,P.Value < 0.1)$names
# fisher.test(tmp$cre,tmp$diffactivity)
tab = table(tmp$diffactivity,tmp$cre)
binom.test(tab[2,2],tab[1,2]+tab[2,2],(tab[2,1])/(tab[2,1]+tab[1,1]))

tmp1 = subset(tmp,genotype=="D21")
tmp2 = subset(tmp,genotype=="T21")
fisher.test(tmp1$cre,tmp1$diffactivity)
fisher.test(tmp2$cre,tmp2$diffactivity)
table(tmp1$cre,tmp1$diffactivity)
table(tmp2$cre,tmp2$diffactivity)



atac[!(atac$names %in% atac2$names),][1:13,]

head(gene_symbols)
subset(gene_symbols,`Gene name`=="H2AFX")


# # read gene coordinates
# genes = fread("~/Downloads/genes",data.table = F,stringsAsFactors = F)
# colnames(genes)[1:3] = c("chrom",'start','end')

sum(data_snv$intron)
sum(data_snv$cre)

dfsub = subset(data_snv,intron)
fisher.test(dfsub$genotype,dfsub$cre)
aggregate(cre~genotype,dfsub,mean)
table(dfsub$cre,dfsub$genotype)
binom.test(87,616+87,(103/(912+103)))

sum(duplicated(tmp2$snv.x))

######

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


