library(data.table)
library(glmmTMB)

df = fread("~/Downloads/mtg_exp_230719_lin.csv",data.table = F,stringsAsFactors = F)
f.out = "~/Downloads/mitotracker_results.txt"
pdf.out = "~/Downloads/mitotracker_results.pdf"
y_label = "MitoTracker"

# df = fread("~/Downloads/mitosox_exp_230719_lin.csv",data.table = F,stringsAsFactors = F)
# f.out = "~/Downloads/mitosox_results.txt"
# pdf.out = "~/Downloads/mitosox_results.pdf"
# y_label = "MitoSOX"

df$values = df[,1]
df$Sample = as.factor(df$Sample)
df$Status = recode(df$Status,"Healthy" = "Disomic", "DS" = "Ts21")
df = subset(df,!(Cell_type %in% "FLHSC"))
cell.uniq = unique(df$Cell_type)
res.lst=list()

aggre = aggregate(df$values,df[,c("Sample","Status","Age","Cell_type")],mean)
library(ggplot2)
g = ggplot(aggre,aes(x=Status,y=x,fill=Status)) + 
  geom_violin() + 
  geom_point() + facet_grid(cols=vars(Cell_type)) + 
  ggpubr::theme_pubr() + labs(y=y_label,x='Trisomy status') +
  scale_fill_brewer(palette="Set2") + guides(fill="none")
pdf(pdf.out,width=9,height=4)
print(g)
dev.off()

# ggplot(df,aes(x=Status,y=MitoTracker)) + geom_boxplot() + facet_grid(cols=vars(Cell_type)) + theme_bw()

for (i in 1:length(cell.uniq)) {
  print(i)
  cell = cell.uniq[i]
  df.sub = subset(df,Cell_type==cell)
  
  library(lme4)
  df.sub$Status = factor(df.sub$Status,levels = c("Disomic","Ts21"))
  # min_val = min(df.sub$values)
  # if (min_val < 0) {
  #   df.sub$values = df.sub$values + abs(min_val) + 1
  # }
  # df.sub$values = round(df.sub$values)
  source("~/Documents/Research/Useful_scripts/rntransform.R")
  df.sub$values = rntransform(df.sub$values)
  # print(summary(lm(MitoTracker~Status + Age + (Sample),df.sub))$coef["StatusDS",])
  
  # hist(df.sub$values)
  # model = lmer(MitoTracker ~ Status + (1 | Sample), data = df.sub)
  # summary(model)
  # 
  # model <- glmer(MitoTracker ~ Status + (1 | Sample), data = df.sub)#, family = glmmTMB::nbinom1)
  # summary(model)
  # 
  # model <- glmer(MitoTracker ~ Status + (1 | Sample), data = df.sub, family = glmmTMB::nbinom1)
  # summary(model)
  
  # model <- glmmTMB(MitoTracker~Status + scale(Age) + (1 | Sample), df.sub, family = gaussian, REML = FALSE)
  # summary(model)
  
  # model <- glmmTMB(MitoTracker~Status + Age + (1 | Sample), df.sub, family = nbinom1, REML = FALSE)
  model <- glmmTMB(values~Status + Age + (1 | Sample), df.sub, family = gaussian, REML = FALSE)
  tab=coef(summary(model))[[1]]
  res = tab["StatusTs21",]
  # res
  res.lst[[i]] = res
}

res.df = data.frame(do.call(rbind,res.lst))
colnames(res.df) = c('est','se','stat','pval')
res.df$cell = cell.uniq
fwrite(res.df,f.out,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)

f1 = "~/Downloads/mitotracker_results.txt"
df1=fread(f1,data.table = F,stringsAsFactors = F)
df1$analysis = "mitotracker"
f2 = "~/Downloads/mitosox_results.txt"
df2=fread(f2,data.table = F,stringsAsFactors = F)
df2$analysis = "mitosox"
df.mg = rbind(df1,df2)
df.mg$fdr = p.adjust(df.mg$pval,method='fdr')
df.mg = df.mg[,c("analysis","cell","est","se","stat","pval","fdr")]
f.out = "~/Downloads/mito_results_all.txt"
fwrite(df.mg,f.out,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)



