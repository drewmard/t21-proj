f="~/Downloads/Ts21_vs_Disomic_DE.supp.csv"
df=fread(f,data.table = F,stringsAsFactors = F)

oxidative_stress_genes = fread("~/Downloads/GO_term_summary_20240529_175417.txt",data.table = F,stringsAsFactors = F)
oxidative_stress_genes = unique(toupper(oxidative_stress_genes[,2]))
tmp = subset(df,names %in% oxidative_stress_genes)
tmp2 = aggregate(tmp$logFC,by=list(tmp$cell),mean)
tmp3 = aggregate(tmp$logFC,by=list(tmp$cell),function(x) {sd(x)/sqrt(length(x))})
tmp4 = merge(tmp2,tmp3,by='Group.1')
tmp4$scaled = scale(tmp4[,2])
tmp4[order(tmp4$scaled),]

ggplot(tmp4,aes(x=reorder(Group.1,scaled,mean),y=scaled)) + geom_point(aes(col=-1*scaled)) +
  theme_bw() +
  ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5),plot.title = element_text(hjust=.5)) +
  scale_color_continuous(type=c("viridis")) +
  geom_abline(slope=0,intercept=0,col='red',lty='dashed') +
  labs(x="Cell type",y="Mean log-fold change",title = "Response to reactive oxygen species (228 genes)") +
  guides(col='none')
  
  


