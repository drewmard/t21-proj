library(ggplot2)
library(data.table)

res = fread("~/Downloads/cosmic_snv.acute_myeloid_leukaemia.txt",data.table = F,stringsAsFactors = F)
ggplot(res,aes(x=type,y=est,ymin = l,ymax=h)) + geom_pointrange() + labs(x="Peak set",y="Somatic mutation rate",title="AML") + 
  ggpubr::theme_pubr() + geom_abline(slope=0,intercept = res$est[1],col='red',lty='dashed') +
  theme(axis.text.x=element_text(hjust=1,angle = 30,vjust=1),plot.title = element_text(hjust=0.5))

res = fread("~/Downloads/cosmic_snv.leukaemia.txt",data.table = F,stringsAsFactors = F)
ggplot(res,aes(x=type,y=est,ymin = l,ymax=h)) + geom_pointrange() + labs(x="Peak set",y="Somatic mutation rate",title="Any leukaemia") + 
  ggpubr::theme_pubr() + geom_abline(slope=0,intercept = res$est[1],col='red',lty='dashed') +
  theme(axis.text.x=element_text(hjust=1,angle = 30,vjust=1),plot.title = element_text(hjust=0.5))

res = fread("~/Downloads/cosmic_snv.txt",data.table = F,stringsAsFactors = F)
ggplot(res,aes(x=type,y=est,ymin = l,ymax=h)) + geom_pointrange() + labs(x="Peak set",y="Somatic mutation rate",title="COSMIC (all non-coding)") + 
  ggpubr::theme_pubr() + geom_abline(slope=0,intercept = res$est[1],col='red',lty='dashed') +
  theme(axis.text.x=element_text(hjust=1,angle = 30,vjust=1),plot.title = element_text(hjust=0.5))


res1 = fread("~/Downloads/cosmic_snv.acute_myeloid_leukaemia.txt",data.table = F,stringsAsFactors = F)
res2 = fread("~/Downloads/cosmic_snv.txt",data.table = F,stringsAsFactors = F)

tmp1=log2(res1[2:5,2:4]/res1[1,2])
tmp2=log2(res2[2:5,2:4]/res2[1,2])
tmp1$type = res1$type[2:5]
tmp2$type = res2$type[2:5]
tmp1$class="AML"
tmp2$class="COSMIC"
tmp = rbind(tmp1,tmp2)

ggplot(tmp,aes(x=type,y=est,ymin = l,ymax=h,col=class)) + 
  geom_pointrange(position = position_dodge(width = .2)) + 
  labs(x="Peak set",y="Log2(Peak set/background set)",title="Somatic mutation rate - Enrichment relative to background",col="Set") + 
  ggpubr::theme_pubr() + geom_abline(slope=0,intercept = 1,col='red',lty='dashed') +
  theme(axis.text.x=element_text(hjust=1,angle = 30,vjust=1),plot.title = element_text(hjust=0.5))

  
