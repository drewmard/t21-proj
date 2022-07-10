library(data.table)
df = fread("~/Documents/Research/t21-proj/manuscript/TFs-on-chr21.csv",data.table = F,stringsAsFactors = F,fill = TRUE)
res <- fread("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.Liver.FE.p.txt",data.table = F,stringsAsFactors = F)
res <- fread("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.Femur.FE.fdr.txt",data.table = F,stringsAsFactors = F)

y1=apply(subset(res,names %in% df$gene.name)[,2:22],2,function(x) sum(x<0.05,na.rm=T))
y2=aggregate(res[,2:22],by=list(res$chromosome_name==21),function(x) mean(x<0.05,na.rm=T))
rownames(y2) <- y2[,1]; y2 <- y2[,-1]
tmp <- as.data.frame(merge(t(y2),t(t(y1)),by='row.names',all=TRUE))
colnames(tmp) <- c("cell","not.chr21","chr21","n.tf")
summary(lm(not.chr21~I(chr21^2) + n.tf,tmp))
summary(lm(not.chr21~I(chr21^2),tmp))
summary(lm(not.chr21~I(chr21^2),tmp))

tmp2 <- subset(res,names %in% df$gene.name)
rownames(tmp2) <- tmp2$names
tmp2[is.na(tmp2)] <- 1
m <- rownames(tmp2)
tmp2 <- apply((tmp2 < 0.05)[,2:22],2,as.numeric)
rownames(tmp2) <- m
dataf = merge(tmp,t(tmp2),by.x='cell',by.y='row.names',all=TRUE)

summary(lm(not.chr21 ~I(chr21^2) + BACH1,dataf))
summary(lm(not.chr21 ~I(chr21^2) + ERG,dataf))
summary(lm(not.chr21 ~I(chr21^2) + ETS2,dataf))
summary(lm(not.chr21 ~I(chr21^2) + GABPA,dataf))
summary(lm(not.chr21 ~I(chr21^2) + OLIG1,dataf))
summary(lm(not.chr21 ~I(chr21^2) + PKNOX1,dataf))
summary(lm(not.chr21 ~I(chr21^2) + PRDM15,dataf))
summary(lm(not.chr21 ~I(chr21^2) + RUNX1,dataf))
summary(lm(not.chr21 ~I(chr21^2) + SON,dataf))
summary(lm(not.chr21 ~I(chr21^2) + ZBTB21,dataf))

apply(tmp2,1,function(x) sum(x,na.rm=T))
apply(tmp2,2,function(x) sum(x,na.rm=T))

