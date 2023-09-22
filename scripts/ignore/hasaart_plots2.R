library(data.table)
meta=fread("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/data/hasaart_et_al_2020_scirep/metadata.csv",data.table = F,stringsAsFactors = F)

k=0; dfall = list()
for (id in c("T21","D21")) { #,"T21_MyeloidPreleuk")) {
  metasub=subset(meta,genotype==id)
  for (i in 1:nrow(metasub)) {
    k = k + 1
    
    res = fread(paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/data/hasaart_et_al_2020_scirep/out/",metasub$name[i],".snv_enrich.txt"),data.table = F,stringsAsFactors = F)
    dfall[[k]] = data.frame(type=res$type[2:5],log2ratio=log2((res[2:5,2] + 1e-30)/(res[1,2] + 1e-30)),genotype=id,name=metasub$name[i])
  }
}
dfall=do.call(rbind,dfall)
dfall$rank = rank(dfall$log2ratio)

dfall.sub = subset(dfall,type=="Ts21 open (P < 0.05; LFC > 0)")
ggplot(dfall.sub,aes(x=genotype,y=log2ratio)) + geom_boxplot()
t.test(dfall.sub$log2ratio[dfall.sub$genotype=="D21"],dfall.sub$log2ratio[dfall.sub$genotype=="T21"])
dfall.sub = subset(dfall,type=="Ts21 closed (P < 0.05; LFC < 0)")
ggplot(dfall.sub,aes(x=genotype,y=log2ratio)) + geom_boxplot()
dfall.sub = subset(dfall,type=="Ts21 sig open (FDR < 0.05; LFC > 0)")
ggplot(dfall.sub,aes(x=genotype,y=log2ratio)) + geom_boxplot()
dfall.sub = subset(dfall,type=="Ts21 sig closed (FDR < 0.05; LFC > 0)")
wilcox.test(dfall.sub$log2ratio[dfall.sub$genotype=="D21"],dfall.sub$log2ratio[dfall.sub$genotype=="T21"])
ggplot(dfall.sub,aes(x=genotype,y=log2ratio)) + geom_boxplot()


ggplot(dfall.sub,aes(x=genotype,y=log2ratio)) + geom_boxplot()

  
  geom_pointrange() + labs(x="Peak set",y="Somatic mutation rate",title="AML") + 
  ggpubr::theme_pubr() + geom_abline(slope=0,intercept = res$est[1],col='red',lty='dashed') +
  theme(axis.text.x=element_text(hjust=1,angle = 30,vjust=1),plot.title = element_text(hjust=0.5))
