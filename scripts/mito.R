library(data.table)
library(glmmTMB)
library(dplyr)

df = fread("~/Downloads/mtg_exp_230719_lin.csv",data.table = F,stringsAsFactors = F)
f.out = "~/Downloads/mitotracker_results.txt"
pdf.out = "~/Downloads/mitotracker_results.v2.pdf"
y_label = "MitoTracker"

# df = fread("~/Downloads/mitosox_exp_230719_lin.csv",data.table = F,stringsAsFactors = F)
# f.out = "~/Downloads/mitosox_results.txt"
# pdf.out = "~/Downloads/mitosox_results.v2.pdf"
# y_label = "MitoSOX"

df$values = df[,1]
df$Sample = as.factor(df$Sample)
df$Status = recode(df$Status,"Healthy" = "Disomic", "DS" = "Ts21")
df = subset(df,!(Cell_type %in% "FLHSC"))
cell.uniq = unique(df$Cell_type)
res.lst=list()

# aggre = aggregate(df$values,df[,c("Sample","Status","Age","Cell_type")],mean)
# library(ggplot2)
# g = ggplot(aggre,aes(x=Status,y=x,fill=Status)) + 
#   geom_violin() + 
#   geom_point() + facet_grid(cols=vars(Cell_type)) + 
#   ggpubr::theme_pubr() + labs(y=y_label,x='Trisomy status') +
#   scale_fill_brewer(palette="Set2") + guides(fill="none")
# pdf(pdf.out,width=9,height=4)
# print(g)
# dev.off()

aggre = aggregate(df$values,df[,c("Sample","Status","Age","Cell_type")],mean)
aggre=as.data.frame(aggre)
aggre2 = merge(aggregate(x~Cell_type+Status,data=aggre,mean),
               aggregate(x~Cell_type+Status,data=aggre,function(x) sd(x)/sqrt(length(x))),
               by=c("Cell_type","Status"))
# aggre2 = merge(aggregate(values~Cell_type+Status,data=df,mean),
#                aggregate(values~Cell_type+Status,data=df,function(x) sd(x)/sqrt(length(x))),
#                by=c("Cell_type","Status"))

colnames(aggre2)[3:4] = c("mu","se")
aggre2$l = aggre2$mu - 1.96*aggre2$se
aggre2$h = aggre2$mu + 1.96*aggre2$se

aggre2$Cell_type[aggre2$Cell_type == "CD38neg"] = "CD38-"
aggre2$Cell_type[aggre2$Cell_type == "CD38pos"] = "CD38+"
aggre2$Cell_type[aggre2$Cell_type == "Linpos"] = "Lin+"
aggre2$Cell_type = factor(aggre2$Cell_type,levels=c("HSC","CD38-","CD38+","Lin+"))
g = ggplot(aggre2,aes(x=Status,y=mu,fill=Status)) + 
  geom_bar(stat='identity',col='black') + 
  geom_pointrange(aes(ymin=l,ymax=h)) + 
  facet_grid(cols=vars(Cell_type)) + 
  ggpubr::theme_pubr() + labs(y=y_label,x='Trisomy status') +
  scale_fill_brewer(palette="Set2") + guides(fill="none")
g
pdf(pdf.out,width=9,height=4)
print(g)
dev.off()


