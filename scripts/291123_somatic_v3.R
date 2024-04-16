library(edgeR)
library(data.table)
# f = "~/Downloads/Healthy_Liver_HSCs_MPPs.pb.txt"
# f = "~/Downloads/DownSyndrome_Liver_Cycling HSCs_MPPs.pb.txt"
f = "~/Downloads/DownSyndrome_Liver_HSCs_MPPs.pb.txt"
x = fread(f,data.table = F,stringsAsFactors = F)
rownames(x)=x[,1];x=x[,-1]
x = cpm(x)

# data downloaded from Hasaart et al:
data_snv = fread("~/Downloads/data_snv",data.table=F,stringsAsFactors=F)
colnames(data_snv) = c("chrom","start","end","sample","genotype","snv")
data_snv = subset(data_snv,genotype != "T21_MyeloidPreleuk")

# Ts21 intronic somatics more likely than D21 HSC genes to reside in CREs for cycling and Ts21 HSC genes compared to D21 intronic somatics
# but not for D21 HSC genes (similar distribution)

result_generate=function() {
  keep = apply(x,1,function(x){mean(x>TPM) > PROP})
  genes_to_keep = rownames(x)[keep]
  
  tmp = subset(data_snv,!exon & intron & `Gene name` %in% genes_to_keep)
  tab = table(tmp$cre,tmp$genotype)
  k=tab[2,2];n=tab[1,2]+tab[2,2];p=(tab[2,1])/(tab[2,1]+tab[1,1])
  res=binom.test(k,n,p)
  out = data.frame(PROP,TPM,
                   prop_t21=k/n,
                   prop_d21=p,
                   l=res$conf.int[1],h=res$conf.int[2],
                   pval=res$p.value)
  return(out)
}


PROP=0.5
TPM=0.2

# f = "~/Downloads/Healthy_Liver_HSCs_MPPs.pb.txt"
# f = "~/Downloads/DownSyndrome_Liver_Cycling HSCs_MPPs.pb.txt"
f = "~/Downloads/DownSyndrome_Liver_HSCs_MPPs.pb.txt"
i = 0;glst=list()
pvec=list();tmp1=list();tmp2=list()
set.seed(031995)
for (f in c(
  "~/Downloads/Healthy_Liver_HSCs_MPPs.pb.txt",
  "~/Downloads/DownSyndrome_Liver_HSCs_MPPs.pb.txt",
  "~/Downloads/DownSyndrome_Liver_Cycling HSCs_MPPs.pb.txt"
)) { 
  i = i+1
  x = fread(f,data.table = F,stringsAsFactors = F); rownames(x)=x[,1];x=x[,-1]; x = cpm(x)
  print(result_generate())
  
  keep = apply(x,1,function(x){mean(x>TPM) > PROP})
  genes_to_keep = rownames(x)[keep]
  
  tmp = subset(data_snv,!exon & intron & `Gene name` %in% genes_to_keep & genotype=="D21")
  bootstrap_func = function(i) {ind = sample(1:nrow(tmp),nrow(tmp),replace = T); tmp2 = tmp[ind,]; sum(tmp2$cre)/nrow(tmp2)}
  pvec[[i]] = unlist(lapply(1:10000,bootstrap_func))
  tmp1[[i]] = subset(data_snv,!exon & intron & `Gene name` %in% genes_to_keep & genotype=="T21")
  tmp2[[i]] = subset(data_snv,!exon & intron & `Gene name` %in% genes_to_keep & genotype=="D21")
  library(ggplot2)
  cat("P-value = ",2*(1-mean((sum(tmp1[[i]]$cre)/nrow(tmp1[[i]])) > pvec[[i]])))
  g=ggplot(data.frame(pvec[[i]]),aes(x=pvec[[i]])) + geom_histogram(fill='darkblue',col='white',bins=16) + 
    # theme_bw() +
    ggpubr::theme_pubr() + 
    theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)) + 
    geom_vline(aes(xintercept = sum(tmp1[[i]]$cre)/nrow(tmp1[[i]]),color='Down Syndrome'),lty='dashed',lwd=2) +
    geom_vline(aes(xintercept = sum(tmp2[[i]]$cre)/nrow(tmp2[[i]]),color='Disomic'),lty='dashed',lwd=2) +
    labs(x="Empirical bootstrap distribution\nin disomic somatics",
         title="",
         y='Count') +
    scale_color_manual(name = "Observed external somatic set", values = c(`Down Syndrome` = "yellow3", Disomic = "orange")) +
    xlim(0.04,0.2)
  # data.frame(p=sum(tmp2$cre)/nrow(tmp2),l=quantile(y,probs = 0.025),h=quantile(y,probs = 0.975))
  # print(g)
  glst[[i]]=g
}
# library(cowplot)
# plot_grid(glst[[1]],glst[[2]],glst[[3]],ncol=1)

make_plot = function(i) {
  g=ggplot(data.frame(pvec[[i]]),aes(x=pvec[[i]])) + geom_histogram(fill='darkblue',col='white',bins=16) + 
    # theme_bw() +
    ggpubr::theme_pubr() + 
    theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)) + 
    geom_vline(aes(xintercept = sum(tmp1[[i]]$cre)/nrow(tmp1[[i]]),color='Down Syndrome'),lty='dashed',lwd=2) +
    geom_vline(aes(xintercept = sum(tmp2[[i]]$cre)/nrow(tmp2[[i]]),color='Disomic'),lty='dashed',lwd=2) +
    labs(x="Empirical bootstrap distribution\nin disomic somatics",
         y='Count') +
    scale_color_manual(name = "Observed external somatic set", values = c(`Down Syndrome` = "yellow3", Disomic = "orange")) +
    xlim(0.04,0.2)
  return(g)
}
gplot = plot_grid(make_plot(1)+labs(title="Disomy HSC-expressed Genes"),
          make_plot(2)+labs(title="Ts21 HSC-expressed Genes"),
          make_plot(3)+labs(title="Ts21 Cycling HSC-expressed Genes"),
          ncol=1)
pdf("~/Downloads/intron_cre_enrichment.pdf",width = 7,height = 11)
print(gplot)
dev.off()





disease_status="DownSyndrome"
pseudobulk_covariate="sample"
suffix=""
f = paste0("~/Downloads/",disease_status,'_',"cyc_vs_hsc",".de",suffix,".txt")
de=fread(f,data.table = F,stringsAsFactors = F)[,1:7]
# de = subset(de,P.Value < 0.05)
# de = subset(de,P.Value < 0.01)
de = subset(de,adj.P.Val < 0.001)

# Furthermore, Ts21 variants in CREs more likely to be near differentially expressed genes compared to Ds21 variants in CREs
# f="~/Documents/Research/t21-proj/out/full/DEG_list/liver_v_femur.txt"
# de=fread(f,data.table = F,stringsAsFactors = F)
# # # atac = subset(atac,class=="t21-induced")
# # # de = subset(de,adj.P.Val.t21 < 0.001)
# de = subset(de,class=="none")

de = fread("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_leiden_names/Liver.HSCs_MPPs.sample.txt",data.table = F,stringsAsFactors = F)
# de = subset(de,adj.P.Val < 0.25)
# de = subset(de,P.Value < 0.05)
de = subset(de,P.Value < 0.2)


de = fread("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_leiden_names/Liver.HSCs_MPPs.sample.txt",data.table = F,stringsAsFactors = F)
i=0; res.df = list()
rng = seq(-2.6,-0.6,by=0.2)
for (thres in rng) {
  thres = 10^(thres)
# for (thres in c(0.5,0.45,0.4,0.35,0.3,0.25,0.2,0.15,0.1,0.05,0.01,0.005)) {
  i=i+1
  desub = subset(de,P.Value < thres)
  genes_to_keep = desub$names
  tmp = subset(data_snv,!exon & intron & `Gene name` %in% genes_to_keep)
  tab = table(tmp$cre,tmp$genotype)
  k=tab[2,2];n=tab[1,2]+tab[2,2];p=(tab[2,1])/(tab[2,1]+tab[1,1])
  res=binom.test(k,n,p)
  res.df[[i]] = data.frame(thres,
                           prop_t21=k/n,
                           prop_d21=p,
                           l=res$conf.int[1],h=res$conf.int[2],
                           pval=res$p.value)
}
res.df = as.data.frame(do.call(rbind,res.df))

library(ggplot2)
res.df$enrich=res.df$prop_t21/res.df$prop_d21
res.df$enrich.h=res.df$h/res.df$prop_d21
res.df$enrich.l=res.df$l/res.df$prop_d21

# ggplot(res.df,aes(x=log10(thres),col=log10(thres))) +
  # geom_pointrange(aes(y=prop_t21,ymin=l,ymax=h),lwd=rel(1.5)) +
  # geom_line(aes(y=prop_d21)) +
g=ggplot(subset(res.df,thres >= 0.01),aes(x=-log10(thres),col=log10(thres))) +
  geom_pointrange(aes(y=enrich,ymin=enrich.l,ymax=enrich.h),lwd=rel(1.5)) +
  geom_line(aes(y=enrich)) +
  geom_abline(intercept=1,slope=0,col='red',lty='dashed') +
  scale_color_viridis_c() +  
  scale_x_continuous() +
  ggpubr::theme_pubr() +
  labs(x=expression("Differential Expression Threshold ("~-log[10](italic(P))~")"),y="Ts21 Somatic Enrichment\n(relative to disomy somatics)",title = "Ts21 vs Disomy (Liver HSCs)") +
  guides(col="none") +
  theme(plot.title = element_text(hjust=0.5))
pdf("~/Downloads/Ts21_somatic.ts21_vs_disomy_hsc.pdf",width = 7,height = 4)
print(g)
dev.off()


disease_status="DownSyndrome"
pseudobulk_covariate="sample"
suffix=""
f = paste0("~/Downloads/",disease_status,'_',"cyc_vs_hsc",".de",suffix,".txt")
de=fread(f,data.table = F,stringsAsFactors = F)[,1:7]
i=0; res.df = list()
rng = seq(-3,-1,by=0.2)
for (thres in rng) {
  thres = 10^(thres)
  # for (thres in c(0.5,0.45,0.4,0.35,0.3,0.25,0.2,0.15,0.1,0.05,0.01,0.005)) {
  i=i+1
  desub = subset(de,P.Value < thres)
  genes_to_keep = desub$names
  tmp = subset(data_snv,!exon & intron & `Gene name` %in% genes_to_keep)
  tab = table(tmp$cre,tmp$genotype)
  k=tab[2,2];n=tab[1,2]+tab[2,2];p=(tab[2,1])/(tab[2,1]+tab[1,1])
  res=binom.test(k,n,p)
  res.df[[i]] = data.frame(thres,
                           prop_t21=k/n,
                           prop_d21=p,
                           l=res$conf.int[1],h=res$conf.int[2],
                           pval=res$p.value)
}
res.df = as.data.frame(do.call(rbind,res.df))

library(ggplot2)
res.df$enrich=res.df$prop_t21/res.df$prop_d21
res.df$enrich.h=res.df$h/res.df$prop_d21
res.df$enrich.l=res.df$l/res.df$prop_d21

# ggplot(res.df,aes(x=log10(thres),col=log10(thres))) +
# geom_pointrange(aes(y=prop_t21,ymin=l,ymax=h),lwd=rel(1.5)) +
# geom_line(aes(y=prop_d21)) +
g = ggplot(res.df,aes(x=-log10(thres),col=log10(thres))) +
  geom_pointrange(aes(y=enrich,ymin=enrich.l,ymax=enrich.h),lwd=rel(1.5)) +
  geom_line(aes(y=enrich)) +
  geom_abline(intercept=1,slope=0,col='red',lty='dashed') +
  scale_color_viridis_c() +  
  scale_x_continuous() +
  ggpubr::theme_pubr() +
  labs(x=expression("Differential Expression Threshold ("~-log[10](italic(P))~")"),y="Ts21 Somatic Enrichment\n(relative to disomy somatics)",title = "Ts21 cycling vs less-cycling HSCs") +
  guides(col="none") +
  theme(plot.title = element_text(hjust=0.5))
g
pdf("~/Downloads/Ts21_somatic.cyc_vs_less_cyc.pdf",width = 7,height = 4)
print(g)
dev.off()

# de = subset(de,P.Value < 0.001)




genes_to_keep = de$names
tmp = subset(data_snv,!exon & intron & nchar(`Gene name`) > 0 )
tmp$diffactivity = tmp$`Gene name` %in% genes_to_keep
tab = table(tmp$diffactivity,tmp$cre)
k=tab[2,2];n=tab[1,2]+tab[2,2];p=(tab[2,1])/(tab[2,1]+tab[1,1])
res=binom.test(k,n,p)
tmp1 = subset(tmp,genotype=="D21")
tmp2 = subset(tmp,genotype=="T21")
fisher.test(tmp1$cre,tmp1$diffactivity)
fisher.test(tmp2$cre,tmp2$diffactivity)
tab = table(tmp1$diffactivity,tmp1$cre)
k=tab[2,2];n=tab[1,2]+tab[2,2];p=(tab[2,1])/(tab[2,1]+tab[1,1])
res=binom.test(k,n,p);res
tab = table(tmp2$diffactivity,tmp2$cre)
k=tab[2,2];n=tab[1,2]+tab[2,2];p=(tab[2,1])/(tab[2,1]+tab[1,1])
res=binom.test(k,n,p);res

fisher.test(tmp$genotype,tmp$diffactivity)

res.df[[i]] = data.frame(PROP,TPM,
                         prop_t21=k/n,
                         prop_d21=p,
                         l=res$conf.int[1],h=res$conf.int[2],
                         pval=res$p.value)



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

