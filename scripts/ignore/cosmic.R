library(data.table)

atac = fread("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/differential_analysis/pb_de_atac.hsc.txt",data.table = F,stringsAsFactors = F)

##############


snv.full = fread("/oak/stanford/groups/smontgom/amarder/data/cosmic/CosmicNCV.acute_myeloid_leukaemia.tsv",data.table = F,stringsAsFactors = F)
# snv = subset(snv.full,V9=="acute_myeloid_leukaemia_associated_with_MDS" | )
# snv = subset(snv.full,V9 %in% c("acute_myeloid_leukaemia","acute_myeloid_leukaemia_associated_with_MDS"))
# snv = subset(snv.full,V9=="acute_myeloid_leukaemia_therapy_related")

# snv.full = fread("/oak/stanford/groups/smontgom/amarder/data/cosmic/CosmicNCV.leukaemia.tsv",data.table = F,stringsAsFactors = F)
# snv=snv.full
# snv.full = fread("/oak/stanford/groups/smontgom/amarder/data/cosmic/CosmicNCV.tsv.gz",data.table = F,stringsAsFactors = F)
# snv = subset(snv.full,`Histology subtype 1`=="acute_myeloid_leukaemia")

y=strsplit(snv[,16],":|-")
snv.df = data.frame(chrom=paste0("chr",unlist(lapply(y,function(x)x[[1]]))),
                    start=as.numeric(unlist(lapply(y,function(x)x[[2]])))-1,
                    end=as.numeric(unlist(lapply(y,function(x)x[[2]]))))
snv.df = unique(snv.df)
y=strsplit(atac$names,"-")
atac.df = data.frame(chrom=paste0(unlist(lapply(y,function(x)x[[1]]))),
                     start=as.numeric(unlist(lapply(y,function(x)x[[2]])))-1,
                     end=as.numeric(unlist(lapply(y,function(x)x[[3]]))))
atac.df = cbind(atac.df,atac)

df.mg = as.data.frame(valr::bed_intersect(atac.df,snv.df))
tab = table(df.mg$names.x)
tab = as.data.frame(tab)

atac.df$snv = atac.df$names %in% df.mg$names.x
atac.df = merge(atac.df,tab,by.x="names",by.y="Var1",all.x=TRUE)
atac.df$Freq[is.na(atac.df$Freq)] = 0
atac.df$span = atac.df$end - atac.df$start

# tmp = subset(atac.df,P.Value > 0.05)
# sum(tmp$Freq)/sum(tmp$span)
# tmp = subset(atac.df,P.Value < 0.05 & logFC > 0)

bootstrap_snv_rate = function(j) {
  {i=sample(1:nrow(tmp),nrow(tmp),replace = T); tmp2 = tmp[i,]; val = sum(tmp2$Freq)/sum(tmp2$span); return(val)}
}

iter = 0; res = list()

iter=iter+1
tmp = subset(atac.df,P.Value > 0.05)
val.open = sum(tmp$Freq)/sum(tmp$span)
boot.open = unlist(parallel::mclapply(1:100,bootstrap_snv_rate,mc.cores = 8))
ci = quantile(boot.open,probs=c(0.025,0.975))
res[[iter]] = data.frame(type="Background (P > 0.05)",est=val.open,l=ci[1],h=ci[2])

iter=iter+1
tmp = subset(atac.df,P.Value < 0.05 & logFC > 0)
val.open = sum(tmp$Freq)/sum(tmp$span)
boot.open = unlist(parallel::mclapply(1:100,bootstrap_snv_rate,mc.cores = 8))
ci = quantile(boot.open,probs=c(0.025,0.975))
res[[iter]] = data.frame(type="Ts21 open (P < 0.05; LFC > 0)",est=val.open,l=ci[1],h=ci[2])

iter=iter+1
tmp = subset(atac.df,P.Value < 0.05 & logFC < 0)
val.open = sum(tmp$Freq)/sum(tmp$span)
boot.open = unlist(parallel::mclapply(1:100,bootstrap_snv_rate,mc.cores = 8))
ci = quantile(boot.open,probs=c(0.025,0.975))
res[[iter]] = data.frame(type="Ts21 closed (P < 0.05; LFC < 0)",est=val.open,l=ci[1],h=ci[2])

iter=iter+1
tmp = subset(atac.df,adj.P.Val < 0.05 & logFC > 0)
val.open = sum(tmp$Freq)/sum(tmp$span)
boot.open = unlist(parallel::mclapply(1:100,bootstrap_snv_rate,mc.cores = 8))
ci = quantile(boot.open,probs=c(0.025,0.975))
res[[iter]] = data.frame(type="Ts21 sig open (FDR < 0.05; LFC > 0)",est=val.open,l=ci[1],h=ci[2])

iter=iter+1
tmp = subset(atac.df,adj.P.Val < 0.05 & logFC < 0)
val.open = sum(tmp$Freq)/sum(tmp$span)
boot.open = unlist(parallel::mclapply(1:100,bootstrap_snv_rate,mc.cores = 8))
ci = quantile(boot.open,probs=c(0.025,0.975))
res[[iter]] = data.frame(type="Ts21 sig closed (FDR < 0.05; LFC > 0)",est=val.open,l=ci[1],h=ci[2])

as.data.frame(do.call(rbind,res))


# fwrite(as.data.frame(do.call(rbind,res)),"~/tmp/cosmic_snv.acute_myeloid_leukaemia.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
# fwrite(as.data.frame(do.call(rbind,res)),"~/tmp/cosmic_snv.leukaemia.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
# fwrite(as.data.frame(do.call(rbind,res)),"~/tmp/cosmic_snv.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)



tmp = subset(atac.df,adj.P.Val < 0.05 & logFC < 0)
val.closed = unlist(parallel::mclapply(1:100,bootstrap_snv_rate,mc.cores = 8))
quantile(val.closed,probs=c(0.025,0.975))
tmp = subset(atac.df,adj.P.Val < 0.05 & logFC > 0)
sum(tmp$Freq)/sum(tmp$span)
tmp = subset(atac.df,adj.P.Val < 0.05 & logFC < 0)
sum(tmp$Freq)/sum(tmp$span)

fisher.test(atac.df$snv,atac.df$logFC > 0)
fisher.test(atac.df$snv,atac.df$P.Value < 0.05)

atac.df.sub = subset(atac.df,P.Value < 0.05)
fisher.test(atac.df.sub$snv,atac.df.sub$logFC > 0)

fisher.test(atac.df$snv,atac.df$logFC > 0)
fisher.test(atac.df$snv,atac.df$P.Value < 0.05)

atac.df.sub = subset(atac.df,adj.P.Val < 0.05)
fisher.test(atac.df.sub$snv,atac.df.sub$logFC > 0)

fisher.
subset(df.mg,P.Value.x < 0.05)