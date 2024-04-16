library(data.table)
df = fread("~/Downloads/DownSyndrome_HSC_PeakGeneSets/liver_v_femur.txt",data.table = F,stringsAsFactors = F)
cb1 = fread("~/Downloads/cellbender/Down Syndrome_HSCs_MPPs.de.txt",data.table = F,stringsAsFactors = F)
cb2 = fread("~/Downloads/cellbender/Healthy_HSCs_MPPs.de.txt",data.table = F,stringsAsFactors = F)

cb = merge(cb1,cb2,by="names")

df.mg = merge(df,cb,by="names")
df.mg$adj.P.Val.x = p.adjust(df.mg$P.Value.x,method = 'fdr')
df.mg$adj.P.Val.y = p.adjust(df.mg$P.Value.y,method = 'fdr')
tmp = subset(df.mg,P.Value.t21 < 0.05 | P.Value.h < 0.05)
cor(tmp[,c("logFC.x","logFC.y","logFC.t21","logFC.h")])

tmp = subset(df.mg,adj.P.Val.t21 < 0.1 | adj.P.Val.h < 0.1)
cor(tmp[,c("logFC.x","logFC.y","logFC.t21","logFC.h")])

genes = fread("~/Documents/Research/data/grch38/genes.extended.txt",data.table = F,stringsAsFactors = F)
df.mg2 = merge(df.mg,genes[,c("Gene name","Chromosome/scaffold name")],by.x="names",by.y = "Gene name")

df.mg2 = subset(df.mg2,`Chromosome/scaffold name` != 21)[,colnames(df.mg2)!="Chromosome/scaffold name"]
tmp = subset(df.mg2,adj.P.Val.t21 < 0.1 | adj.P.Val.h < 0.1)
cor(tmp[,c("logFC.x","logFC.y","logFC.t21","logFC.h")])

tmp2 = subset(df.mg2,adj.P.Val.t21 < 0.1)
plot(tmp2$logFC.t21,tmp2$logFC.x)


plot(df.mg$logFC.h,df.mg$logFC.y)
plot(df.mg$logFC.t21,df.mg$logFC.x)

table(df.mg$adj.P.Val.t21 < 0.05,df.mg$adj.P.Val.x < 0.05)
fisher.test(df.mg$adj.P.Val.t21 < 0.05,df.mg$adj.P.Val.x < 0.05)

exact2x2::exact2x2(df.mg$P.Value.t21 < 0.05,df.mg$P.Value.x < 0.05)$p.value
exact2x2::exact2x2(df.mg$P.Value.h < 0.05,df.mg$P.Value.y < 0.05)$p.value

gene_lst <- subset(df.mg,class=="t21-induced" & logFC.t21 > 0)$names
library(enrichR)
dbs <- c("GO_Molecular_Function_2021", "GO_Biological_Process_2021","ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")
enriched <- enrichr(gene_lst, dbs)
y <- enriched[[dbs[1]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); y$Set <- dbs[1]; ykeep <- y
y <- enriched[[dbs[2]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); ytmp = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); ytmp$Set <- dbs[2]; ykeep<-rbind(ykeep,ytmp)
y <- enriched[[dbs[3]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); ytmp = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); ytmp$Set <- dbs[3]; ykeep<-rbind(ykeep,ytmp)
ykeep[,1]


# gene_lst <- subset(df.mg,class=="t21-induced" & logFC.t21 > 0 & logFC.x > 0 & adj.P.Val.x < 0.1)$names
gene_lst <- subset(df.mg,class=="t21-induced" & logFC.t21 > 0)$names; length(gene_lst)
gene_lst <- subset(df.mg,class=="t21-induced" & logFC.t21 > 0 & logFC.x > 0 & P.Value.x < 0.05)$names; length(gene_lst)

library(enrichR)
dbs <- c("GO_Molecular_Function_2021", "GO_Biological_Process_2021","ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")
enriched <- enrichr(gene_lst, dbs)
y <- enriched[[dbs[1]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); y$Set <- dbs[1]; ykeep <- y
y <- enriched[[dbs[2]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); ytmp = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); ytmp$Set <- dbs[2]; ykeep<-rbind(ykeep,ytmp)
y <- enriched[[dbs[3]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); ytmp = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); ytmp$Set <- dbs[3]; ykeep<-rbind(ykeep,ytmp)
ykeep[,1]


