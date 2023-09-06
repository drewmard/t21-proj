library(data.table)
sampletype="Femur"
# res <- fread(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".FE.fdr.txt"),data.table = F,stringsAsFactors=F)
res <- fread(paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".FE.fdr.txt"),data.table = F,stringsAsFactors=F)
cell_type_groups <- colnames(res)[(which(colnames(res)=="names")+1):(which(colnames(res)=="chromosome_name")-1)]

res[is.na(res)] = 1
x <- aggregate(res[,cell_type_groups],by=list(chr21=res$chromosome_name==21),function(x) {mean(x<0.05,na.rm=T)})
rownames(x) <- x[,1]; x <- x[,-1]; x <- as.data.frame(t(x))
colnames(x) <- c("Not Chr 21","Chr 21")
x$cell <- rownames(x)
tab <- reshape2::melt(x,id.vars="cell")

library(ggplot2)
ggplot(x,aes(x=100*`Chr 21`,y=100*`Not Chr 21`)) + 
  geom_point() + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.title=element_text(hjust=0.5)) +
  geom_smooth(method='loess',span=3,se=F,col='purple') + 
  labs(x="% DEG per cell type (chr 21)",y="% DEG per cell type (not chr 21)",title=paste0(sampletype))


library(data.table)
# sampletype="Femur"
res.lfc <- fread(paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".FE.fdr.txt"),data.table = F,stringsAsFactors=F)
res.lfc[is.na(res.lfc)] = 1
res.lfc[,cell_type_groups] = as.numeric(res.lfc[,cell_type_groups] < 0.05)
# res.lfc <- fread(paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".FE.expr.txt"),data.table = F,stringsAsFactors=F)
# res.lfc[is.na(res.lfc)] = 0
# res.lfc <- fread(paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".FE.lfc.txt"),data.table = F,stringsAsFactors=F)
# res.lfc[is.na(res.lfc)] = 0
gene_lst=c("ERG","ETS2","BACH1")
gene_lst=subset(res.lfc,chromosome_name==21)$names
res.lfc.sub=subset(res.lfc,names %in% gene_lst)
rownames(res.lfc.sub) = res.lfc.sub$names
res.lfc.sub = res.lfc.sub[,cell_type_groups]

tab.sub = subset(tab,variable=="Chr 21")
tmp = merge(tab.sub,t(res.lfc.sub),by.x="cell",by.y=0)
pval=rep(NA,length(gene_lst))
for (i in 1:length(gene_lst)) {
  if (sum(!is.na(tmp[,gene_lst[i]])) > 4) {
    pval[i] = cor.test(tmp[,gene_lst[i]],tmp$value)$p.value
  }
}
min(pval,na.rm = T)
restmp = data.frame(gene_lst,pval,fdr=p.adjust(pval,method = 'fdr'))
restmp[order(restmp$pval)[1:3],]

plot(tmp[,gene_lst[i]],tmp$value)

cor.test(tmp[,gene_lst[2]],tmp$value)
cor.test(tmp[,gene_lst[3]],tmp$value)

cor.test(res.lfc.sub)
aggregate(res.sub[,cell_type_groups],by=list(chr21=res.sub$chromosome_name==21),function(x) {mean(x<0.05,na.rm=T)})

