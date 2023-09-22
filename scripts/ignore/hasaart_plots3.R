library(ggplot2)
library(data.table)

res = fread("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/data/hasaart_et_al_2020_scirep/out/D21.cycling.snv_enrich.txt",data.table = F,stringsAsFactors = F)
ggplot(res[1:2,],aes(x=type,y=est,ymin = l,ymax=h)) + geom_pointrange() + labs(x="Differentially expressed gene set",y="Somatic mutation rate",title="D21") + 
  ggpubr::theme_pubr() + geom_abline(slope=0,intercept = res$est[3],col='red',lty='dashed') +
  theme(axis.text.x=element_text(hjust=1,angle = 30,vjust=1),plot.title = element_text(hjust=0.5))

res = fread("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/data/hasaart_et_al_2020_scirep/out/T21.cycling.snv_enrich.txt",data.table = F,stringsAsFactors = F)
ggplot(res[1:2,],aes(x=type,y=est,ymin = l,ymax=h)) + geom_pointrange() + labs(x="Differentially expressed gene set",y="Somatic mutation rate",title="T21") + 
  ggpubr::theme_pubr() + geom_abline(slope=0,intercept = res$est[3],col='red',lty='dashed') +
  theme(axis.text.x=element_text(hjust=1,angle = 30,vjust=1),plot.title = element_text(hjust=0.5)) +
  ylim(res$est[3],max(res$h))

res = fread("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/data/hasaart_et_al_2020_scirep/out/T21.snv_enrich.txt",data.table = F,stringsAsFactors = F)
ggplot(res,aes(x=type,y=est,ymin = l,ymax=h)) + geom_pointrange() + labs(x="Peak set",y="Somatic mutation rate",title="AML") + 
  ggpubr::theme_pubr() + geom_abline(slope=0,intercept = res$est[1],col='red',lty='dashed') +
  theme(axis.text.x=element_text(hjust=1,angle = 30,vjust=1),plot.title = element_text(hjust=0.5))

res = fread("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/data/hasaart_et_al_2020_scirep/out/T21_MyeloidPreleuk.snv_enrich.txt",data.table = F,stringsAsFactors = F)
ggplot(res,aes(x=type,y=est,ymin = l,ymax=h)) + geom_pointrange() + labs(x="Peak set",y="Somatic mutation rate",title="AML") + 
  ggpubr::theme_pubr() + geom_abline(slope=0,intercept = res$est[1],col='red',lty='dashed') +
  theme(axis.text.x=element_text(hjust=1,angle = 30,vjust=1),plot.title = element_text(hjust=0.5))

res = fread("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/data/hasaart_et_al_2020_scirep/out/T21_and_T21_MyeloidPreleuk.snv_enrich.txt",data.table = F,stringsAsFactors = F)
ggplot(res,aes(x=type,y=est,ymin = l,ymax=h)) + geom_pointrange() + labs(x="Peak set",y="Somatic mutation rate",title="AML") + 
  ggpubr::theme_pubr() + geom_abline(slope=0,intercept = res$est[1],col='red',lty='dashed') +
  theme(axis.text.x=element_text(hjust=1,angle = 30,vjust=1),plot.title = element_text(hjust=0.5))


res = fread("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/data/hasaart_et_al_2020_scirep/out/E080416.snv_enrich.txt",data.table = F,stringsAsFactors = F)
ggplot(res,aes(x=type,y=est,ymin = l,ymax=h)) + geom_pointrange() + labs(x="Peak set",y="Somatic mutation rate",title="AML") + 
  ggpubr::theme_pubr() + geom_abline(slope=0,intercept = res$est[1],col='red',lty='dashed') +
  theme(axis.text.x=element_text(hjust=1,angle = 30,vjust=1),plot.title = element_text(hjust=0.5))
res = fread("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/data/hasaart_et_al_2020_scirep/out/MH2.snv_enrich.txt",data.table = F,stringsAsFactors = F)
ggplot(res,aes(x=type,y=est,ymin = l,ymax=h)) + geom_pointrange() + labs(x="Peak set",y="Somatic mutation rate",title="AML") + 
  ggpubr::theme_pubr() + geom_abline(slope=0,intercept = res$est[1],col='red',lty='dashed') +
  theme(axis.text.x=element_text(hjust=1,angle = 30,vjust=1),plot.title = element_text(hjust=0.5))

res = fread("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/data/hasaart_et_al_2020_scirep/out/MH3.snv_enrich.txt",data.table = F,stringsAsFactors = F)
ggplot(res,aes(x=type,y=est,ymin = l,ymax=h)) + geom_pointrange() + labs(x="Peak set",y="Somatic mutation rate",title="AML") + 
  ggpubr::theme_pubr() + geom_abline(slope=0,intercept = res$est[1],col='red',lty='dashed') +
  theme(axis.text.x=element_text(hjust=1,angle = 30,vjust=1),plot.title = element_text(hjust=0.5))
res = fread("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/data/hasaart_et_al_2020_scirep/out/F100916W15.snv_enrich.txt",data.table = F,stringsAsFactors = F)
ggplot(res,aes(x=type,y=est,ymin = l,ymax=h)) + geom_pointrange() + labs(x="Peak set",y="Somatic mutation rate",title="AML") + 
  ggpubr::theme_pubr() + geom_abline(slope=0,intercept = res$est[1],col='red',lty='dashed') +
  theme(axis.text.x=element_text(hjust=1,angle = 30,vjust=1),plot.title = element_text(hjust=0.5))

res1 = fread("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/data/hasaart_et_al_2020_scirep/out/D21.snv_enrich.txt",data.table = F,stringsAsFactors = F)
res2 = fread("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/data/hasaart_et_al_2020_scirep/out/T21.snv_enrich.txt",data.table = F,stringsAsFactors = F)

res2$est/res1$est
(res2$est[2:5]/res2$est[1]) / (res1$est[2:5]/res1$est[1])
(res2$l[2:5]/res2$l[1]) / (res1$l[2:5]/res1$l[1])


  # next: add in protein coding genes!
  


